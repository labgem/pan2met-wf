
process PROTEIN_COMPLEX_RULE {

    // conda "pan2met"
    // TODO: add dependency to the library.

    input:
    path 'monomers.list'
    path 'monomer_complex_rules.lp'

    output:
    path 'complexes.list', emit: complex
    path 'versions.yml', emit: versions

    shell:
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pan2met: \$(pan2met --version | sed 's/pan2met //g')
    END_VERSIONS
    
    """
    
}
