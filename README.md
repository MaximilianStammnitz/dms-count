# dms-count

```
  >>>>>>>>>>>>>>>>>>>>>
  >>>>>>X>>>>>>>>>>>>>>          \(o_o)/       S77R   |||
  >>>>>>>>>>>>>X>>>>>>>    ==> dms_count ==>   A298X  ||
  >>>>X>>>>>>>>>>>>>>>>                        E335K  ||
  >>>>>>>>>>>>>>>>>>>>>
```

`dms_count` counts the variants in aligned reads from deep mutational scanning (DMS) experiments. It reads a BAM or SAM file of reads aligned to a single open reading frame (ORF) and writes one table row per distinct variant. Each row holds the variant's read count, its mean coverage, and its annotation at nucleotide, codon and amino-acid level.

It is a single C file on top of [HTSlib](https://www.htslib.org/). It is the variant counter of the nf-core pipeline [nf-core/deepmutscan](https://nf-co.re/deepmutscan) from v1.1.0 on, and its output is byte-identical to the pipeline's earlier Python counter.

| counter (1.14 M reads, one thread) | CPU time | peak memory |
| --- | --- | --- |
| deepmutscan v1.0.0 (Python, pysam) | 39.9 s | 155 MB |
| `dms_count` | 0.7 s | 14 MB |

## Installation

From Bioconda:

```bash
conda install -c bioconda -c conda-forge dms-count
```

From source, with HTSlib 1.17 or later installed:

```bash
make -f build.mk HTSLIB_PREFIX=/usr/local     # prefix of an installed HTSlib (include/, lib/)
make -f build.mk HTSLIB=/path/to/htslib       # or a built HTSlib source tree (libhts.a)
bash test/run_test.sh                         # checks the build against the bundled test data
```

## Usage

```bash
dms_count --bam aligned.bam --fasta orf.fasta --reading_frame 352-1383 \
          --id sample1 --out sample1.variant_counts.tsv
```

| option | default | meaning |
| --- | --- | --- |
| `--bam FILE` | required | BAM or SAM; `-` reads from standard input. Need not be sorted or indexed. |
| `--fasta FILE` | required | the reference, with exactly one sequence |
| `--reading_frame A-B` | required | first and last base of the ORF in the reference (1-based, inclusive) |
| `--id NAME` | `sample` | sample name, used in error messages |
| `--out FILE` | required | the variant table |
| `--base_qual Q` | 40 | minimum base quality, for variant bases and for coverage |
| `--min_flank N` | 2 | ignore variants and coverage in the first and last N bases of each read; with 0, coverage counts A/C/G/T bases only and skips QC-fail and duplicate reads |
| `--skip_indel_reads` | off | ignore reads whose alignment has an insertion, deletion or skip |
| `--threads N` | 1 | HTSlib threads for decompressing the input |
| `--columns FILE` | `variant_counts_columns.tsv` | where to write the column descriptions |
| `--version` | | print the version of `dms_count` and HTSlib |

Notes on the input:

- Variants are read from the `MD` tag, as written by BWA and bwa-mem3. Reads without it are counted for coverage but not for variants; `samtools calmd` adds it.
- Unmapped reads are skipped, and secondary and supplementary alignments add no variants. There is no mapping-quality filter; filter beforehand if you need one (e.g., `samtools view -q 30`).
- Streaming works, e.g. `aligner ... | samtools view -u - | dms_count --bam - ...`.

## Output

A tab-separated table without header, one row per distinct combination of nucleotide changes in a read, sorted by the first changed position. `dms_count` prints the column descriptions to `--columns`:

| column | content |
| --- | --- |
| `counts` | reads carrying this variant |
| `cov` | mean coverage over the variant's changed bases |
| `mean_length_variant_reads` | placeholder, always 0 |
| `varying_bases` | number of changed bases |
| `base_mut` | changes as `position:ref>alt` (reference position, 1-based), e.g. `352:G>A` |
| `varying_codons` | number of changed codons |
| `codon_mut` | codon changes as `codon:ref>alt` (codon 1 = first codon of the reading frame), e.g. `1:GCT>ACT` |
| `aa_mut` | amino-acid change per codon with type: `S` synonymous, `M` missense, `N` nonsense; stop = `X`, e.g. `M:A>T` |
| `pos_mut` | compact amino-acid changes, e.g. `A1T;A2H` |

Codons are translated with the standard genetic code. Ambiguous bases give the amino acid all their codons agree on, otherwise `X`. A summary line goes to standard error:

```
dms_count: 4587 reads, 378 variant-carrying reads, 364 distinct variants
```

## Test data

`test/` holds 4,587 reads of the GID1A amplicon from the nf-core/deepmutscan test data, aligned with bwa-mem3, the reference, and the expected outputs. `bash test/run_test.sh [path/to/dms_count]` runs `dms_count` from the file and from standard input and compares both with the expected table.

## Citation

If you use `dms_count`, please cite nf-core/deepmutscan:

> Wehnert B, et al. bioRxiv preprint (in preparation). https://nf-co.re/deepmutscan

## License

MIT, see [LICENSE](LICENSE).
