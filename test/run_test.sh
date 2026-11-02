#!/usr/bin/env bash
# Runs dms_count on the bundled test BAM and compares the output with the expected tables.
# Usage (from the repository root, after building): bash test/run_test.sh [path/to/dms_count]
set -euo pipefail
cd "$(dirname "$0")"
DMS_COUNT="${1:-../dms_count}"
mkdir -p out
"$DMS_COUNT" --version
"$DMS_COUNT" --bam GID1A_test.bam --fasta GID1A.fasta --reading_frame 352-1383 \
    --base_qual 40 --min_flank 2 --skip_indel_reads \
    --id GID1A_test --columns out/columns.tsv --out out/variant_counts.tsv
# the same alignments read from standard input
"$DMS_COUNT" --bam - --fasta GID1A.fasta --reading_frame 352-1383 \
    --base_qual 40 --min_flank 2 --skip_indel_reads \
    --id GID1A_test --columns out/columns_stdin.tsv --out out/variant_counts_stdin.tsv < GID1A_test.bam
cmp out/variant_counts.tsv expected.variant_counts.tsv
cmp out/variant_counts_stdin.tsv expected.variant_counts.tsv
cmp out/columns.tsv expected.columns.tsv
echo "OK: output identical to the expected tables"
