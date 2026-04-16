
include { DIAMOND_ALIGN_REFERENCE } from "../../../modules/local/diamond"

process ASSOCIATE_PROTEIN_TO_REACTION {

    input:
    path 'monomer_to_reactions.tsv'
    path 'blastp.tsv'

    output:
    path 'protein_reaction.asso', emit: asso

    shell:
    """
    awk -F "\\t" -v OFS="\\t" \
        -f "$baseDir/bin/hash_join.awk" \
        -v key1=1 -v value1=2 \
        -v key2=2 -v value2=1 \
            "monomer_to_reactions.tsv" \
        <(grep -v "^#" 'blastp.tsv')  > "protein_reaction.asso"
    """
}

process REMOVE_METACYC_META_PREFIX {
    input:
    path metacyc_monomers_fasta
    output:
    path "formatted_metacyc_monomers.fasta", emit: fasta
    script:
    """
    sed 's/>gnl|META|/>/g' "${metacyc_monomers_fasta}"  > "formatted_metacyc_monomers.fasta"
    """
}

process REMOVE_ECOCYC_ECO_PREFIX {
    input:
    path ecocyc_monomers_fasta
    output:
    path "formatted_ecocyc_monomers.fasta", emit: fasta
    script:
    """
    sed 's/>gnl|ECO|/>/g' "${ecocyc_monomers_fasta}"  > "formatted_ecocyc_monomers.fasta"
    """
}


workflow CYC_BASED_ASSOCIATION {
 
    take:
    family_proteins
    reference_proteins
    monomer_to_reactions
    coverage_threshold
    identity_threshold

    main:

    DIAMOND_ALIGN_REFERENCE(reference_proteins, family_proteins, coverage_threshold, identity_threshold)
    ASSOCIATE_PROTEIN_TO_REACTION(monomer_to_reactions, DIAMOND_ALIGN_REFERENCE.out.tsv)

    emit:
    asso = ASSOCIATE_PROTEIN_TO_REACTION.out.asso
    versions = DIAMOND_ALIGN_REFERENCE.out.versions
}
