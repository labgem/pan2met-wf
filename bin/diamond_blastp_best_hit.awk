#!/usr/bin/env -S awk -f

# AWK script to retrieve the best hit based on identity value from diamond blastp TSV output
# Assumes input TSV has columns: query, subject, identity, length, mismatches, gapopens, qstart, qend, sstart, send, evalue, bitscore
# and that the TSV is sorted by query id.
# Usage:
# awk -f diamond_blastp_best_hit.awk DIAMOND_BLASTP_TSV > DIAMOND_BLASTP_TSV_BEST_HITS

BEGIN {
    # Initialize variables
    best_identity = 0
    current_query = ""
    best_hit = ""
}

{
    # Skip header lines
    if ($1 ~ /^#/) next

    query = $1
    subject = $2
    identity = $3

    # If we're moving to a new query
    if (query != current_query) {
        # Print previous best hit if exists
        if (current_query != "") {
            print best_hit
        }

        # Reset for new query
        current_query = query
        best_identity = identity
        best_hit = $0
    } else {
        # For same query, check if this hit has higher identity
        if (identity > best_identity) {
            best_identity = identity
            best_hit = $0
        }
    }
}

# Print the last query's best hit
END {
    if (current_query != "") {
        print best_hit
    }
}
