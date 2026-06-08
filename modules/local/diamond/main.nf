#!/usr/bin/env nextflow

process DIAMOND_MAKEDB {

    conda "bioconda::diamond==2.1.13"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.10--h5ca1c30_3' :
        'biocontainers/diamond:2.1.13--h13889ed_0' }"


    input:
    path 'reference.faa'

    output:
    path 'diamond.dmnd', emit: diamond_db
    path 'versions.yml', emit: versions
    
    script:
    """
    diamond makedb --threads ${task.cpus} --db "diamond" --in "reference.faa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version | sed 's/diamond version //g')
    END_VERSIONS
    """
}

process DIAMOND_BLASTP {

    conda "bioconda::diamond==2.1.13"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.10--h5ca1c30_3' :
        'biocontainers/diamond:2.1.13--h13889ed_0' }"


    input:
    val coverage_threshold
    val identity_threshold
    path input_fasta
    path diamond_db

    output:
    path "${input_fasta.baseName}.diamond.blastp.tsv", emit: tsv
    path 'versions.yml', emit: versions

    script:
    """
    diamond blastp --threads ${task.cpus} --header \
        --query "${input_fasta}" --db "${diamond_db}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --out "${input_fasta.baseName}.diamond.blastp.tsv" \
        --ultra-sensitive --query-cover "${coverage_threshold}" --subject-cover "${identity_threshold}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version | sed 's/diamond version //g')
    END_VERSIONS
    """
}

workflow DIAMOND_ALIGN_REFERENCE {
    take:
        reference_proteins
        proteins
        coverage_threshold
        identity_threshold
    main:
    ch_versions = channel.empty()
    DIAMOND_MAKEDB(reference_proteins)
    ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)
    DIAMOND_BLASTP(coverage_threshold, identity_threshold, proteins, DIAMOND_MAKEDB.out.diamond_db)
    ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)
    emit:
    tsv = DIAMOND_BLASTP.out.tsv
    versions = ch_versions
}
