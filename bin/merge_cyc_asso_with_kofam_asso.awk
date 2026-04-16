#!/usr/bin/env -S awk -f
# Merge two files, with first file having priority
# The first file contains two columns:
# - protein identifier
# - metacyc reactions list
# The second file contains three columns:
# - protein identifier
# - metacyc reactions list
# - ec number list
# lists are comma-separated

# Read the first file and store all Protein to Reaction associations.
# Ensure that all protein identifier is unique in the first column.
#
# Output: a three column file:
# - protein ID
# - MetaCyc reactions
# - EC numbers
# If MetaCyc reaction is set, then EC number is not set.
BEGIN {
    FS="\t"
    OFS="\t"
}

function die(message) {
    print "Error:" message > "/dev/stderr"
    exit(1)
}

# First file
NR == FNR {
    protein = $1
    reactions = $2
    protein_to_reactions[protein] = reactions
    next
}

# Second file
{
    protein = $1
    reactions = $2
    ec_numbers = $3
    if (!(protein in protein_to_reactions)) {
        if (reactions) {
            protein_to_reactions[protein] = reactions
        } else if (ec_numbers) {
            protein_to_ec_numbers[protein] = ec_numbers
        }
    }
}

END {
    for (protein in protein_to_reactions) {
        print protein, protein_to_reactions[protein], ""
    }
    for (protein in protein_to_ec_numbers) {
        print protein, "", protein_to_ec_numbers[protein]
    }
}
