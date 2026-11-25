# Changelog

## 0.2.0

- Fixed: `cov` is the joint coverage of a variant: the reads whose usable interval spans all its changed bases, capped by their smallest per-base coverage and never below the variant's count. For variants with changes in several codons it falls with the distance between them, as in GATK AnalyzeSaturationMutagenesis. 0.1.0 reported the mean per-base coverage instead, about the same at every distance, which overstated the coverage of every multi-base variant whose changes are not covered by all reads together. Single-base variants are unchanged; two or three changes in one codon can be slightly lower. Every other column, and the rows and their order, are unchanged.
- The joint coverage is identical to the counter nf-core/deepmutscan used for its error-correction development (August 2026): byte-identical tables on four BAMs from 4,535 to 1.12 million reads, at `--min_flank` 2 and 0. It adds about 0.05 s of CPU time and 8 MB of memory per million reads.

## 0.1.0

- First release: a single-pass C/HTSlib port of the nf-core/deepmutscan v1.0.0 variant counter, byte-identical to it.
