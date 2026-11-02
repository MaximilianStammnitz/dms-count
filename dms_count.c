/*
 * dms_count - single-pass, C/htslib implementation of deepmutscan's count_variants.py
 * (modules/local/variantcounting/templates/count_variants.py, dev @ eb0d3ce).
 *
 * Same inputs, same arithmetic, same output file, byte for byte:
 *   - coverage per reference base inside the ORF (two modes, as in the template: with --min_flank > 0
 *     every aligned base of every mapped read with quality >= --base_qual, excluding the first and last
 *     min_flank query positions; with --min_flank 0 pysam's count_coverage semantics: A/C/G/T only,
 *     quality >= base_qual, reads flagged unmapped/secondary/QC-fail/duplicate skipped),
 *   - one row per distinct nucleotide-variant combination (GATK AnalyzeSaturationMutagenesis notation),
 *     counts, mean coverage over the changed bases, codon/amino-acid annotation, sorted by leading
 *     reference position then by the variant string.
 * Everything happens in one pass over the BAM; the per-variant coverage average is computed at the end
 * from the complete coverage track, which is what the template does too (its workers see the finished
 * coverage dict).
 *
 * build: gcc -O2 -o dms_count dms_count.c -I<htslib> <htslib>/libhts.a -lz -lbz2 -llzma -lpthread -lm
 * usage: dms_count --bam in.bam --fasta ref.fa --reading_frame 352-1383 [--base_qual 40] [--min_flank 2] [--skip_indel_reads]
 *                  [--threads N] [--columns variant_counts_columns.tsv] --id SAMPLE --out out.tsv
 *        dms_count --version
 * The column sheet (--columns, default variant_counts_columns.tsv in the working directory) is the one the
 * Nextflow module publishes next to the count table.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <stdarg.h>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/thread_pool.h"

#define DMS_COUNT_VERSION "0.1.0"
#define die(...) do { fprintf(stderr, "ERROR: " __VA_ARGS__); fprintf(stderr, "\n"); exit(1); } while (0)

/* ---------- genetic code (standard table 1), with Biopython's handling of IUPAC-ambiguous codons:
 * if every expansion of the codon translates to the same residue (stop included) that residue is
 * returned, otherwise 'X'. --------------------------------------------------------------------- */
static const char *AA = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
static int base_idx(char c) { switch (c) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; } return -1; }
static const char *iupac(char c) {
    switch (toupper((unsigned char)c)) {
        case 'A': return "A"; case 'C': return "C"; case 'G': return "G"; case 'T': return "T"; case 'U': return "T";
        case 'R': return "AG"; case 'Y': return "CT"; case 'S': return "CG"; case 'W': return "AT"; case 'K': return "GT";
        case 'M': return "AC"; case 'B': return "CGT"; case 'D': return "AGT"; case 'H': return "ACT"; case 'V': return "ACG";
        case 'N': return "ACGT";
    }
    return NULL;
}
static char translate(const char *codon) {
    const char *e0 = iupac(codon[0]), *e1 = iupac(codon[1]), *e2 = iupac(codon[2]);
    if (!e0 || !e1 || !e2) return 'X';
    char aa = 0;
    for (const char *a = e0; *a; a++) for (const char *b = e1; *b; b++) for (const char *c = e2; *c; c++) {
        char r = AA[base_idx(*a) * 16 + base_idx(*b) * 4 + base_idx(*c)];
        if (!aa) aa = r; else if (aa != r) return 'X';
    }
    return aa;
}

/* ---------- a small string-keyed hash map for the per-variant rows ---------- */
typedef struct {
    char *base_mut, *codon_mut, *aa_mut, *pos_mut;
    int64_t counts; int varying_bases, varying_codons;
    int *pos; int npos;            /* 0-based reference positions of the changed bases (for the coverage mean) */
    int sort_pos;
} row_t;
typedef struct { row_t *rows; size_t n, cap; int64_t *slots; size_t nslots; } table_t;
static uint64_t hash_str(const char *s) { uint64_t h = 1469598103934665603ULL; for (; *s; s++) { h ^= (unsigned char)*s; h *= 1099511628211ULL; } return h; }
static void table_init(table_t *t) { t->n = 0; t->cap = 1024; t->rows = malloc(t->cap * sizeof(row_t)); t->nslots = 4096; t->slots = malloc(t->nslots * sizeof(int64_t)); memset(t->slots, 0xff, t->nslots * sizeof(int64_t)); }
static void table_rehash(table_t *t) {
    size_t ns = t->nslots * 2; int64_t *s = malloc(ns * sizeof(int64_t)); memset(s, 0xff, ns * sizeof(int64_t));
    for (size_t i = 0; i < t->n; i++) { uint64_t h = hash_str(t->rows[i].base_mut) & (ns - 1); while (s[h] >= 0) h = (h + 1) & (ns - 1); s[h] = (int64_t)i; }
    free(t->slots); t->slots = s; t->nslots = ns;
}
static row_t *table_find_or_add(table_t *t, const char *key, int *added) {
    if (t->n * 2 >= t->nslots) table_rehash(t);
    uint64_t h = hash_str(key) & (t->nslots - 1);
    while (t->slots[h] >= 0) { if (strcmp(t->rows[t->slots[h]].base_mut, key) == 0) { *added = 0; return &t->rows[t->slots[h]]; } h = (h + 1) & (t->nslots - 1); }
    if (t->n == t->cap) { t->cap *= 2; t->rows = realloc(t->rows, t->cap * sizeof(row_t)); }
    t->slots[h] = (int64_t)t->n; row_t *r = &t->rows[t->n++]; memset(r, 0, sizeof(*r)); r->base_mut = strdup(key); *added = 1; return r;
}

/* ---------- growable string buffer ---------- */
typedef struct { char *s; size_t n, cap; } sb_t;
static void sb_reset(sb_t *b) { b->n = 0; if (b->s) b->s[0] = 0; }
static void sb_cat(sb_t *b, const char *fmt, ...) {
    va_list ap; for (;;) { va_start(ap, fmt); size_t avail = b->cap - b->n; int w = vsnprintf(b->s ? b->s + b->n : NULL, avail, fmt, ap); va_end(ap);
        if (w < 0) die("format"); if ((size_t)w < avail) { b->n += (size_t)w; return; } b->cap = b->cap ? b->cap * 2 + (size_t)w : (size_t)w + 64; b->s = realloc(b->s, b->cap); }
}

/* ---------- FASTA (single record) ---------- */
static char *read_fasta(const char *path, char **name_out) {
    FILE *f = fopen(path, "r"); if (!f) die("cannot open %s", path);
    char line[1 << 16]; sb_t seq = {0}; char *name = NULL; int nrec = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '>') { nrec++; if (nrec == 1) { line[strcspn(line, " \t\r\n")] = 0; name = strdup(line + 1); } continue; }
        for (char *p = line; *p; p++) if (!isspace((unsigned char)*p)) sb_cat(&seq, "%c", (char)toupper((unsigned char)*p));
    }
    fclose(f);
    if (nrec != 1) { fprintf(stderr, "ERROR: --fasta must contain exactly one sequence (found %d). deepmutscan counts variants against a single ORF.\n", nrec); exit(1); }
    *name_out = name; return seq.s ? seq.s : strdup("");
}

static int cmp_rows(const void *a, const void *b) {
    const row_t *x = a, *y = b; if (x->sort_pos != y->sort_pos) return x->sort_pos < y->sort_pos ? -1 : 1; return strcmp(x->base_mut, y->base_mut);
}

typedef struct { int ref_pos, read_idx, alt_qual, codon_idx, pos_in_codon; char ref_base, alt_base; } var_t;

int main(int argc, char **argv) {
    const char *bam_path = NULL, *fasta = NULL, *rf = NULL, *id = "sample", *out = NULL, *columns = "variant_counts_columns.tsv";
    int base_qual = 40, min_flank = 2, threads = 1, skip_indel = 0;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--skip_indel_reads")) { skip_indel = 1; continue; }
        if (!strcmp(argv[i], "--version")) { printf("dms_count %s (htslib %s)\n", DMS_COUNT_VERSION, hts_version()); return 0; }
        if (!strcmp(argv[i], "--columns")) { columns = argv[++i]; continue; }
        if (!strcmp(argv[i], "--bam")) bam_path = argv[++i]; else if (!strcmp(argv[i], "--fasta")) fasta = argv[++i];
        else if (!strcmp(argv[i], "--reading_frame")) rf = argv[++i]; else if (!strcmp(argv[i], "--base_qual")) base_qual = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--min_flank")) min_flank = atoi(argv[++i]); else if (!strcmp(argv[i], "--threads")) threads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--id")) id = argv[++i]; else if (!strcmp(argv[i], "--out")) out = argv[++i];
        else if (!strcmp(argv[i], "--batch_size")) ++i; else die("unknown option %s", argv[i]);
    }
    if (!bam_path || !fasta || !rf || !out) die("usage: dms_count --bam in.bam --fasta ref.fa --reading_frame A-B [--base_qual Q] [--min_flank F] [--skip_indel_reads] [--threads N] [--columns variant_counts_columns.tsv] --id ID --out out.tsv | --version");
    int o1, o2; if (sscanf(rf, "%d-%d", &o1, &o2) != 2) die("bad --reading_frame %s", rf);
    const int orf_start = o1 - 1, orf_end = o2 - 1;          /* 0-based inclusive, as in the template */
    const int region_len = orf_end - orf_start + 1;

    char *ref_name; char *ref = read_fasta(fasta, &ref_name); int ref_len = (int)strlen(ref);
    if (orf_start < 0 || orf_end >= ref_len) die("--reading_frame %s lies outside the reference (%d nt)", rf, ref_len);
    /* reference codons: only complete codons inside the ORF, 1-based codon index, as in the template */
    int n_full = 0; for (int i = 0; i < region_len - 2; i += 3) if (i + 3 <= region_len) n_full++;

    /* columns sheet, as the template writes it */
    { FILE *c = fopen(columns, "w"); if (c) { fputs("column\tdescription\ncounts\tNumber of reads supporting this variant (observations)\ncov\tMean sequencing coverage across the variant's mutated base positions\nmean_length_variant_reads\tPlaceholder (0); mean read length of supporting reads is not computed\nvarying_bases\tNumber of nucleotide substitutions in the variant\nbase_mut\tNucleotide change(s), 1-based reference position(s), e.g. 352:G>A\nvarying_codons\tNumber of codons altered by the variant\ncodon_mut\tCodon change(s), 1-based codon index, e.g. 1:GCT>ACT\naa_mut\tAmino-acid change per codon; type S=synonymous, M=missense, N=nonsense; stop=X, e.g. M:A>T\npos_mut\tCompact per-codon amino-acid change, e.g. A1T (stop written as X)\n", c); fclose(c); } }

    samFile *fp = sam_open(bam_path, "r"); if (!fp) die("cannot open %s", bam_path);
    htsThreadPool pool = {NULL, 0}; if (threads > 1) { pool.pool = hts_tpool_init(threads); hts_set_opt(fp, HTS_OPT_THREAD_POOL, &pool); }
    sam_hdr_t *hdr = sam_hdr_read(fp); if (!hdr) die("cannot read header");
    int tid_ref = sam_hdr_name2tid(hdr, ref_name);

    int64_t *cov = calloc(region_len, sizeof(int64_t));
    table_t table; table_init(&table);
    bam1_t *b = bam_init1();
    var_t *vars = NULL; int nvars_cap = 0;
    sb_t key = {0}, codon_s = {0}, aa_s = {0}, pos_s = {0};
    int64_t n_reads = 0, n_variant_reads = 0;

    while (sam_read1(fp, hdr, b) >= 0) {
        n_reads++;
        const uint16_t flag = b->core.flag;
        if (flag & BAM_FUNMAP) continue;
        if (b->core.tid != tid_ref) continue;
        const int qlen = b->core.l_qseq; const uint8_t *qual = bam_get_qual(b); const uint8_t *seq = bam_get_seq(b);
        const uint32_t *cigar = bam_get_cigar(b); const int ncig = b->core.n_cigar;
        const int rs = (int)b->core.pos;
        /* --skip_indel_reads: ignore reads whose alignment has an insertion, deletion or reference skip,
           as the pipeline's `awk '$6 !~ /[IDN]/'` filter did on the SAM stream */
        if (skip_indel) {
            int has = 0;
            for (int k = 0; k < ncig; k++) { int op = bam_cigar_op(cigar[k]); if (op == BAM_CINS || op == BAM_CDEL || op == BAM_CREF_SKIP) { has = 1; break; } }
            if (has) continue;
        }

        /* ---- coverage ---- */
        if (min_flank > 0) {
            if (qlen > 0 && qual[0] != 0xff) {
                const int lo = min_flank, hi = qlen - min_flank; int rp = rs, qp = 0;
                for (int k = 0; k < ncig; k++) {
                    int op = bam_cigar_op(cigar[k]), ln = bam_cigar_oplen(cigar[k]);
                    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                        for (int i = 0; i < ln; i++, rp++, qp++) {
                            if (rp < orf_start || rp > orf_end) continue; if (qp < lo || qp >= hi) continue; if (qual[qp] < base_qual) continue;
                            cov[rp - orf_start]++;
                        }
                    } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) qp += ln; else if (op == BAM_CDEL || op == BAM_CREF_SKIP) rp += ln;
                }
            }
        } else {
            if (!(flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))) {
                int rp = rs, qp = 0;
                for (int k = 0; k < ncig; k++) {
                    int op = bam_cigar_op(cigar[k]), ln = bam_cigar_oplen(cigar[k]);
                    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                        for (int i = 0; i < ln; i++, rp++, qp++) {
                            if (rp < orf_start || rp > orf_end) continue; if (qual[0] != 0xff && qual[qp] < base_qual) continue;
                            int nt = bam_seqi(seq, qp); if (nt != 1 && nt != 2 && nt != 4 && nt != 8) continue;   /* A C G T only */
                            cov[rp - orf_start]++;
                        }
                    } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) qp += ln; else if (op == BAM_CDEL || op == BAM_CREF_SKIP) rp += ln;
                }
            }
        }

        /* ---- variants (template: read_bam_in_chunk + parse_read + gatk_formating) ---- */
        if (flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        uint8_t *mdp = bam_aux_get(b, "MD"); if (!mdp) continue;
        const char *md = bam_aux2Z(mdp); if (!md) continue;
        { int alldig = 1; for (const char *p = md; *p; p++) if (!isdigit((unsigned char)*p)) { alldig = 0; break; } if (alldig) continue; }
        /* MD -> list of (type, ref_pos, ref_base), as parse_md */
        int nmd = 0, mdcap = 64; var_t *mdv = malloc(mdcap * sizeof(var_t)); int rp = rs;
        for (const char *p = md; *p;) {
            if (isdigit((unsigned char)*p)) { int n = 0; while (isdigit((unsigned char)*p)) n = n * 10 + (*p++ - '0'); rp += n; }
            else if (*p == '^') { p++; while (*p && isupper((unsigned char)*p)) { if (nmd == mdcap) { mdcap *= 2; mdv = realloc(mdv, mdcap * sizeof(var_t)); } mdv[nmd].ref_pos = rp; mdv[nmd].ref_base = *p; mdv[nmd].alt_base = 'D'; nmd++; rp++; p++; } }
            else if (isupper((unsigned char)*p)) { if (nmd == mdcap) { mdcap *= 2; mdv = realloc(mdv, mdcap * sizeof(var_t)); } mdv[nmd].ref_pos = rp; mdv[nmd].ref_base = *p; mdv[nmd].alt_base = 'X'; nmd++; rp++; p++; }
            else p++;
        }
        if (nmd == 0) { free(mdv); continue; }
        /* CIGAR walk, exactly as parse_read (only op 0 'M' carries mismatches; 'D' entries are kept but never reported) */
        int nv = 0; int brp = rs, bqp = 0, mdi = 0;
        for (int k = 0; k < ncig; k++) {
            int op = bam_cigar_op(cigar[k]), ln = bam_cigar_oplen(cigar[k]);
            if (op == BAM_CMATCH) {
                for (int i = 0; i < ln; i++) {
                    if (mdi < nmd && mdv[mdi].ref_pos == brp) {
                        if (nv == nvars_cap) { nvars_cap = nvars_cap ? nvars_cap * 2 : 64; vars = realloc(vars, nvars_cap * sizeof(var_t)); }
                        var_t *v = &vars[nv++]; v->ref_pos = brp; v->ref_base = mdv[mdi].ref_base; v->alt_base = seq_nt16_str[bam_seqi(seq, bqp)]; v->alt_qual = qual[bqp];
                        v->read_idx = bqp; v->codon_idx = (brp - orf_start >= 0 ? (brp - orf_start) / 3 : -((orf_start - brp + 2) / 3)) + 1; v->pos_in_codon = ((brp - orf_start) % 3 + 3) % 3;
                        if (mdv[mdi].alt_base == 'D') v->alt_base = 0; /* a 'D' MD entry met inside an M block: template records var_type 'D', dropped later */
                        mdi++;
                    }
                    brp++; bqp++;
                }
            } else if (op == BAM_CINS) bqp += ln;
            else if (op == BAM_CDEL) { for (int i = 0; i < ln; i++) { if (mdi < nmd && mdv[mdi].ref_pos == brp) mdi++; brp++; } }   /* template keeps 'D' vars; they never reach base_mut */
            else if (op == BAM_CREF_SKIP) brp += ln; else if (op == BAM_CSOFT_CLIP) bqp += ln;
            /* BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF: no movement, as in the template */
        }
        free(mdv);
        if (nv == 0) continue;
        /* filters: inside the ORF, quality, flank; D entries drop out (they are not 'X') */
        int nf = 0;
        for (int i = 0; i < nv; i++) {
            var_t *v = &vars[i]; if (v->alt_base == 0) continue;
            if (v->ref_pos < orf_start || v->ref_pos > orf_end) continue; if (v->alt_qual < base_qual) continue;
            if (min_flank > 0 && (v->read_idx < min_flank || v->read_idx >= qlen - min_flank)) continue;
            if (v->codon_idx < 1 || v->codon_idx > n_full) die("codon %d of read %s is outside the complete codons of the reading frame (the template raises KeyError here)", v->codon_idx, bam_get_qname(b));
            vars[nf++] = *v;
        }
        if (nf == 0) continue;
        n_variant_reads++;
        /* base_mut, in read order (= increasing reference position) */
        sb_reset(&key);
        for (int i = 0; i < nf; i++) sb_cat(&key, "%s%d:%c>%c", i ? ", " : "", vars[i].ref_pos + 1, vars[i].ref_base, vars[i].alt_base);
        int added; row_t *r = table_find_or_add(&table, key.s, &added);
        r->counts++;
        if (!added) continue;
        /* first occurrence: annotate (the other columns are functions of base_mut) */
        r->varying_bases = nf; r->npos = nf; r->pos = malloc(nf * sizeof(int)); for (int i = 0; i < nf; i++) r->pos[i] = vars[i].ref_pos;
        r->sort_pos = vars[0].ref_pos + 1;
        sb_reset(&codon_s); sb_reset(&aa_s); sb_reset(&pos_s); int ncod = 0;
        /* codons in increasing codon index (variants are already in increasing ref position, so codons are too) */
        for (int i = 0; i < nf;) {
            int ci = vars[i].codon_idx; char refc[4], altc[4]; memcpy(refc, ref + orf_start + 3 * (ci - 1), 3); refc[3] = 0; memcpy(altc, refc, 4);
            int j = i; for (; j < nf && vars[j].codon_idx == ci; j++) altc[vars[j].pos_in_codon] = vars[j].alt_base;
            char ra = translate(refc), aa = translate(altc); char mt = (aa == ra) ? 'S' : (aa == '*') ? 'N' : 'M';
            char ras = ra == '*' ? 'X' : ra, aas = aa == '*' ? 'X' : aa;
            sb_cat(&codon_s, "%s%d:%s>%s", ncod ? ", " : "", ci, refc, altc); sb_cat(&aa_s, "%s%c:%c>%c", ncod ? ", " : "", mt, ras, aas); sb_cat(&pos_s, "%s%c%d%c", ncod ? ";" : "", ras, ci, aas);
            ncod++; i = j;
        }
        r->varying_codons = ncod; r->codon_mut = strdup(codon_s.s); r->aa_mut = strdup(aa_s.s); r->pos_mut = strdup(pos_s.s);
    }
    bam_destroy1(b); sam_hdr_destroy(hdr); sam_close(fp); if (pool.pool) hts_tpool_destroy(pool.pool);

    { int64_t tot = 0; for (int i = 0; i < region_len; i++) tot += cov[i];
      if (tot == 0) { fprintf(stderr, "ERROR: no reads pass the quality filters inside the reading frame %s of %s for sample %s.\n  Check that --reading_frame lies inside the reference, that the reads really come from this amplicon, and that --base_qual (currently %d) is not above the instrument quality of the data.\n", rf, ref_name, id, base_qual); return 1; } }

    qsort(table.rows, table.n, sizeof(row_t), cmp_rows);
    FILE *o = fopen(out, "w"); if (!o) die("cannot write %s", out);
    for (size_t i = 0; i < table.n; i++) {
        row_t *r = &table.rows[i]; int64_t s = 0; for (int j = 0; j < r->npos; j++) s += cov[r->pos[j] - orf_start];
        int64_t cov_avg = r->npos ? (int64_t)((double)s / r->npos) : 0;       /* int(mean(...)): truncation of the float mean */
        fprintf(o, "%lld\t%lld\t0\t%d\t%s\t%d\t%s\t%s\t%s\n", (long long)r->counts, (long long)cov_avg, r->varying_bases, r->base_mut, r->varying_codons, r->codon_mut, r->aa_mut, r->pos_mut);
    }
    fclose(o);
    fprintf(stderr, "dms_count: %lld reads, %lld variant-carrying reads, %zu distinct variants\n", (long long)n_reads, (long long)n_variant_reads, table.n);
    return 0;
}
