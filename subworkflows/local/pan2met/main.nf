/** Extract Gene - Protein - Reaction relations
/*
/**/

/* Merge associations from 1.tsv with associations from 2.tsv
/* keeping association from 2.tsv only when there is no association for column 1 key in 1.tsv
/* */
process MERGE_ASSOCIATION {
    input:
    path "1.tsv"
    path "2.tsv"

    output:
    path "*.asso", emit: asso

    script:
    def prefix = task.ext.prefix ?: "merged"
    """
    awk -F "\\t" -v OFS="\\t" -f "${workflow.projectDir}/bin/merge_two_files_priority.awk" "1.tsv" "2.tsv" > "${prefix}.asso"
    """
}

process MERGE_ASSOCIATION_WITH_EC {
    input:
    path cyc_association_prot_metacyc_reactions
    path kofam_association_prot_metacyc_reactions_and_ec

    output:
    path '*.asso', emit: asso

    script:
    def prefix = task.ext.prefix ?: "merged_with_ec"
    """
    awk -F "\\t" -v OFS="\t" -f "${workflow.projectDir}/bin/merge_cyc_asso_with_kofam_asso.awk" \
        "${cyc_association_prot_metacyc_reactions}" "${kofam_association_prot_metacyc_reactions_and_ec}" \
        > "${prefix}.asso"
    """
}
