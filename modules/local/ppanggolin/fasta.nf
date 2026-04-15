
process PPANGGOLIN_FASTA_FAMILY_PROTEINS {

    conda 'bioconda::ppanggolin==2.2.4'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ppanggolin:2.2.4--h0fa9677_0' :
        'biocontainers/ppanggolin:2.2.4--h0fa9677_0' }"

    input:
    path 'pangenome.h5'

    output:
    path 'proteins/all_protein_families.faa', emit: pangenome_family_proteins
    path 'versions.yml', emit: versions

    """
    ppanggolin fasta -p "pangenome.h5" --output ./proteins --prot_families all -f

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ppanggolin: \$(ppanggolin --version | sed 's/ppanggolin //g')
    END_VERSIONS
    """
}

workflow PPANGGOLIN_EXTRACT_SEQUENCES {

    take:
    pangenome

    main:
    PPANGGOLIN_FASTA_FAMILY_PROTEINS(pangenome)

    emit:
    pangenome_family_proteins = PPANGGOLIN_FASTA_FAMILY_PROTEINS.out.pangenome_family_proteins
    versions = PPANGGOLIN_FASTA_FAMILY_PROTEINS.out.versions
}
