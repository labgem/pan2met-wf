#!/usr/bin/env -S awk -f
# Create PathoLogic .pf file

BEGIN {
    FS="\t"
    OFS="\t"
}

{
    gene=$1
    rcyc=$2
    ec=$3

    if (gene!=bufgene) {
        if (bufgene) print "//"
        print "ID", gene
        print "NAME", gene
        print "STARTBASE", 1
        print "ENDBASE", 9 # Multiple of 3
        print "PRODUCT-TYPE", "P"
    }

    if (rcyc) {
        already_has_rcyc[gene] = 1 # Set a flag, to avoid adding EC numbers,
                                   # when a gene is already linked to a reaction
        split(rcyc, rcyc_array, ",")
        for (i in rcyc_array) {
            print "METACYC", rcyc_array[i]
        }
    }
    if (ec && !(gene in already_has_rcyc)) {
        split(ec, ec_array, ",")
        for (i in ec_array) {
            print "EC", ec_array[i]
        }
    }

    bufgene=gene
}

END {
    print "//"
}
