/** Extract Gene - Protein - Reaction relations
/*
/**/

process MERGE_ASSOCIATION {

    input:
    path '1.tsv'
    path '2.tsv'

    output:
    path 'merged.tsv' , emit: asso

    shell:
    """
    awk -F "\\t" -v OFS="\\t" -f "${baseDir}/bin/merge_two_files_priority.awk" "1.tsv" "2.tsv" > "merged.tsv"
    """
}

process MERGE_ASSOCIATION_WITH_EC {


    input:
    path cyc_association_prot_metacyc_reactions
    path kofam_association_prot_metacyc_reactions_and_ec

    output:
    path 'merged_cyc_and_kofam.tsv', emit: asso

    shell:
    """
    awk -F "\\t" -v OFS="\t" -f "${baseDir}/bin/merge_cyc_asso_with_kofam_asso.awk" "${cyc_association_prot_metacyc_reactions}" "${kofam_association_prot_metacyc_reactions_and_ec}" > 'merged_cyc_and_kofam.tsv'
    """
}
