#!/usr/bin/env -S awk -f
#
# Input: a two column file:
# - ID
# - VALUE
# where ID can appear several times
# Output:
# A file with two columns:
# - ID
# - VALUES
# where VALUES is the concatenation (separatad by ',')
# of all values for ID.

BEGIN {
    FS="\t"
    OFS="\t"
}

{
    id = $1
    value = $2
    if (id in table) {
        table[id] = table[id] "," value
    } else {
        table[id] = value
    }
}

END {
    for (id in table) {
        print id, table[id]
    }
}
