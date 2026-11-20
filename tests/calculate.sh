#!/bin/bash
# Summarise a provisioned output directory into a form a baseline can be compared against:
# the file inventory and each VCF's record count per FILTER value. No sizes, paths, headers
# or timestamps, so the baseline does not drift.
set -o nounset
set -o errexit
set -o pipefail
export LC_ALL=C

cd "$1" || exit 1

find . -type f | sed 's|^\./||' | sort

for vcf in *.vcf.gz *.vcf.bgz; do
    [ -e "$vcf" ] || continue
    echo "## ${vcf}"
    zcat "$vcf" | grep -v '^#' | awk -F'\t' '{ print $7 }' | sort | uniq -c | sort -k2,2
done
