#!/usr/bin/env -S awk -f
# Merge two TSV files.
# The first column of both file is the protein ID.
# The second column of both file is the column to merge.
# If the second columns contains "NA" or "" in the first file,
# the column will be populated with the content of the second column
# of the second file. Otherwise, the content will be merged.
# The first column act as a key in the hash table.

# Usage:
#
# awk -f merge_two_file_priority.awk PRIORITY_FILE ALTERNATIVE_FILE

function join(array, sep, start, end)
{
    if (start == "") {
        start = 1
    }
    if (end == "") {
        end = length(array)
    }
    if (sep == "")
       sep = " "
    else if (sep == SUBSEP) # magic value
       sep = ""
    result = array[start]
    for (_i_join = start + 1; _i_join <= end; _i_join++)
        result = result sep array[_i_join]
    return result
}

function find(item, set) {
    start = 1
    end = length(set)
    for (_i_find = start; _i_find <= end; _i_find++) {
        if (set[_i_find] == item)
            return 1
    }
    return 0
}

# Merge two comma-seperated lists,
# to keep each entry once
function union_reactions(set1, set2) {
    split(set1, set1_array, ",")
    split(set2, set2_array, ",")
    for (_i_union_reactions in set2_array) {
        if (!(find(set2_array[_i_union_reactions], set1_array))) {
            set1_array[length(set1_array) + 1] = set2_array[_i_union_reactions]
        }
    }
    return join(set1_array, ",")
}

BEGIN {
    FS="\t"
    OFS="\t"
}

function add(key, reactions) {
    if (reactions == "NA")
        return
    if (key in reaction_array) {
        reactions = union_reactions(reactions, reaction_array[key])
    }
    reaction_array[key] = reactions
}

# Read the first file
NR == FNR {
    add($1, $2)
    priorit_file[$1] = 1
}

# Read the second file
{
    if (!($1 in priority_file)) {
        add($1, $2)
    }
}

END {
    for (key in reaction_array) {
        print key, reaction_array[key]
    }
}
