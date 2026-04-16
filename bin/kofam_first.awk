#!/usr/bin/env -S awk -f
# Keep only the first KO group associated with a protein from a KOFamScan result TSV
# Usage:
# awk -f kofam_first.awk KOFAM_TSV
BEGIN {
    FS="\t"
    OFS="\t"
}

{
    # Check if we are on the first KO result
    if ($1 == "*") {
        protein = $2
        ko = $3
        print protein, ko
    }
}
