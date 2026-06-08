include { SEQKIT_SPLIT2 as SEQKIT_SPLIT_FASTA } from '../../../modules/nf-core/seqkit/split2/main'
include { KOFAMSCAN } from '../../../modules/nf-core/kofamscan/main'

process KOFAMSCAN_ASSOCIATION {

    input:
    path kofamscan_tsv

    output:
    path "kofamscan_metacyc_ec.asso", emit: asso

    script:
    """
    awk -f "${workflow.projectDir}/bin/kofam_first.awk" "${kofamscan_tsv}" > "kofamscan_first.tsv"
    awk -f "${workflow.projectDir}/bin/merge_cell.awk" "kofamscan_first.tsv" > "kofamscan_uniq_id.tsv"
    sort --field-separator=\$'\\t' --key=2 "kofamscan_uniq_id.tsv" > "kofamscan_uniq_id_sorted_by_ko.tsv"
    sort --field-separator=\$'\\t' --key=1 "${params.kegg_kos_to_metacyc_reactions}" | cut -d\$'\\t' -f1,3 > "sorted_kegg_kos_to_metacyc_reactions.tsv"
    sort --field-separator=\$'\\t' --key=1 "${params.kegg_kos_to_ec_numbers}" > "sorted_kegg_kos_to_ec_numbers.tsv"
    # Left join to add MetaCyc reactions
    join -t \$'\\t' -1 2 -2 1 -a 1 "kofamscan_uniq_id_sorted_by_ko.tsv" "sorted_kegg_kos_to_metacyc_reactions.tsv" > "kofamscan_metacyc.tsv"
    # Left join to add EC-numbers
    join -t \$'\\t' -1 1 -2 1 -a 1 "kofamscan_metacyc.tsv" "sorted_kegg_kos_to_ec_numbers.tsv" | cut -d\$'\\t' -f2,3,4 > "kofamscan_metacyc_ec.asso"
    """
}

process CONCAT {
    input:
    path '*.tsv'

    output:
    path 'concatenation.tsv'

    script:
    """
    cat *.tsv > concatenation.tsv
    """
}

process EXTRACT_FILE {
    input:
    tuple val(meta), path(file)
    output:
    path "$file", emit: file
    script:
    """
    # Nothing to do.
    """
}



workflow KOFAMSCAN_BASED_ASSOCIATION {
    take:
    proteins

    main:
    ch_versions = channel.empty()
    ch_proteins = proteins.map{
        prot -> [[ id: prot.baseName, single_end: true ], prot ]
    }
    SEQKIT_SPLIT_FASTA(ch_proteins)
    // ch_versions = ch_versions.mix(SEQKIT_SPLIT_FASTA.out.versions_seqkit) // TODO handle the special case of nf-core module seqkit version channel
    ch_fasta = SEQKIT_SPLIT_FASTA.out.reads.transpose()
    KOFAMSCAN(ch_fasta, params.kofam_profiles, params.kofam_ko_list)
    ch_versions = ch_versions.mix(KOFAMSCAN.out.versions)
    EXTRACT_FILE(KOFAMSCAN.out.tsv) | collect | CONCAT | KOFAMSCAN_ASSOCIATION

    emit:
    asso = KOFAMSCAN_ASSOCIATION.out.asso
    versions = ch_versions
}
