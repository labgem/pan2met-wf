
process KOFAMSCAN_ANNOTATION {

    label 'process_high'

    conda 'bioconda::kofamscan==1.3.0'

    /* container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2' :
        'biocontainers/kofamscan:1.3.0--hdfd78af_2' }"
     */ // Container raises an issue: how to mount the shared bank ?


    input:
    path fasta_chunk

    output:
    path "${fasta_chunk.baseName}.kofamscan.tsv", emit: tsv

    shell:
    """
    exec_annotation \
        --profile "${params.kofam_db}/profiles" \
        -k "${params.kofam_db}/ko_list" \
        --tmp-dir "./tmp" \
        -o "${fasta_chunk.baseName}.kofamscan.tsv" \
        --format=detail-tsv \
        --create-alignment \
        --cpu="${task.cpus}" \
        "${fasta_chunk}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kofamscan: \$(exec_annotation --version | sed 's/exec_annotation //g')
    END_VERSIONS
    """
}

process KOFAMSCAN_ASSOCIATION {

    input:
    path kofamscan_tsv

    output:
    path "kofamscan_metacyc_ec.tsv", emit: asso

    shell:
    """
    awk -f "${baseDir}/bin/kofam_first.awk" "${kofamscan_tsv}" > "kofamscan_first.tsv"
    awk -f "${baseDir}/bin/merge_cell.awk" "kofamscan_first.tsv" > "kofamscan_uniq_id.tsv"
    sort --field-separator=\$'\\t' --key=2 "kofamscan_uniq_id.tsv" > "kofamscan_uniq_id_sorted_by_ko.tsv"
    sort --field-separator=\$'\\t' --key=1 "${params.kegg_kos_to_metacyc_reactions}" | cut -d\$'\\t' -f1,3 > "sorted_kegg_kos_to_metacyc_reactions.tsv"
    sort --field-separator=\$'\\t' --key=1 "${params.kegg_kos_to_ec_numbers}" > "sorted_kegg_kos_to_ec_numbers.tsv"
    # Left join to add MetaCyc reactions
    join -t \$'\\t' -1 2 -2 1 -a 1 "kofamscan_uniq_id_sorted_by_ko.tsv" "sorted_kegg_kos_to_metacyc_reactions.tsv" > "kofamscan_metacyc.tsv"
    # Left join to add EC-numbers
    join -t \$'\\t' -1 1 -2 1 -a 1 "kofamscan_metacyc.tsv" "sorted_kegg_kos_to_ec_numbers.tsv" | cut -d\$'\\t' -f2,3,4 > "kofamscan_metacyc_ec.tsv"
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

workflow KOFAMSCAN_BASED_ASSOCIATION {
    take:
    proteins

    main:
    proteins.splitFasta(by: params.chunkSize, file: true).set{ ch_fasta }
    KOFAMSCAN_ANNOTATION(ch_fasta) | collect | CONCAT | KOFAMSCAN_ASSOCIATION


    emit:
    asso = KOFAMSCAN_ASSOCIATION.out.asso
    versions = KOFAMSCAN_ANNOTATION.out.versions
}
