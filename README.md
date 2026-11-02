# dms_count: the variant counter in C

`dms_count.c` is a one-pass, single-file C/htslib implementation of `modules/local/variantcounting/templates/count_variants.py` (dev @ eb0d3ce). Same command line as the benchmark CLI of the Python script, same `variant_counts.tsv`, byte for byte (checked on the four test libraries at `--min_flank 2` and `0`, on the 1.14 M-read scaled BAM at 1 and 2 threads, and on the end-to-end-aligned BAMs with their additional multi-base variants):

| counter | CPU s (1.14 M reads, 1 thread) | wall s | peak RSS |
| --- | --- | --- | --- |
| dev (v1.0.0 template) | 39.9 | 37.4 | 155 MB |
| Python, patched (this folder's patch) | 9.1 | 8.1 | 158 MB |
| **dms_count** | **0.7** | **0.7** | **14 MB** |

What it does in that one pass: coverage per ORF position with the template's two semantics (flank + quality filter over all aligned bases for `--min_flank > 0`; pysam `count_coverage` semantics, A/C/G/T only and QC-fail/duplicate reads skipped, for `--min_flank 0`), mismatch extraction from CIGAR + MD exactly as `parse_read`, the GATK-style variant strings, codon/amino-acid annotation with Biopython's rules for ambiguous codons, aggregation, the per-variant coverage mean computed at the end from the complete track, and the template's sort order. It also writes `<out>.columns.tsv` and `<out>.versions.yml`.

Build: `make -f build.mk HTSLIB=/path/to/htslib-1.21` (a built source tree) or `make -f build.mk HTSLIB_PREFIX=$CONDA_PREFIX` (an installed htslib). `dms_count --version` prints the tool and htslib versions; `--columns FILE` sets where the column sheet goes (default `variant_counts_columns.tsv`, as the module expects). Also here: `LICENSE` (MIT), `Dockerfile` (local test image), `bioconda/dms-count/` (recipe) and `pipeline_module/` (the replacement module files); the step-by-step is in `../../README_implementation_steps.md`.

## How it could reach the pipeline

nf-core modules install their tools from Bioconda/Biocontainers, so the binary needs a package: a small Bioconda recipe (`deepmutscan-count`, or a broader `deepmutscan-tools` package that could later also hold the error-correction code), built from a tagged release of this file with `htslib` as its only dependency. The replacement module (`pipeline_module/modules/local/variantcounting/`) is already commit 5 of `patch_performance_boost/`: `dms_count` with the same inputs and outputs as the Python module, label `process_single`; what it still needs is the package, so that `nf-core pipelines lint --fix container_configs` can pin its container. The patched Python template (commit 1) is the output-neutral stepping stone kept in the history. Step-by-step in `../../README_implementation_steps.md`.

The same file is also the natural home for the alignment-free / ungapped counting idea of the plan's §2.3 level 3 and for the v2.0.0 variant model (multi-codon variants, indels): both are easier to add here than to the Python template.
