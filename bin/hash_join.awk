#!/usr/bin/env -S awk -f
# Usage:
# awk -f hash_join.awk -v key1=FILE1_KEY_COLUMN -v key2=FILE2_KEY_COLUMN -v value1=FILE1_VALUE_COLUMN -v value2=FILE2_VALUE_COLUMN FILE1 FILE2
# File 1
# Record a hash table with every cross reference
FNR == NR {
    key = $key1
    value = $value1
    values[key] = value
    # print key, value
    next
}

# File 2
{
    key = $key2
    if (key in values) {
        value = values[key]
    } else {
        value = ""
    }
    print $value2, value
}
