
process PPANGGOLIN_WORKFLOW {
    label 'process_large'

    conda "bioconda::ppanggolin=2.2.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.4--h0fa9677_0' :
        'biocontainers/ppanggolin:2.2.4--h0fa9677_0' }"

    input:
    path genomes_list

    output:
    path "pangenome.h5", emit: pangenome
    path "versions.yml", emit: versions

    script:
    """
    ppanggolin workflow --fasta "${genomes_list}" -o "${params.pgdb}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """
}
