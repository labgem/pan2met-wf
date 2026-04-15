process PREPARE_PATHOLOGIC_INPUT {

    input:
    path 'pangenome.gpr'

    output:
    path "${params.pgdb}"

    shell:
    """
    mkdir -p "${params.pgdb}/"
    # Prepare fake fna file
    cat > "${params.pgdb}/${params.pgdb}.fna" <<EOF
    >${params.pgdb}
    AAAAAAAAAA
    TTTTTTTTTT
    EOF
    # Prepare genetic-elements.dat
    cat > "${params.pgdb}/genetic-elements.dat" <<EOF
    ID	${params.pgdb}
    NAME	${params.pgdb}
    TYPE	:CHRSM
    CIRCULAR?	Y
    ANNOT-FILE	${params.pgdb}.pf
    SEQ-FILE	${params.pgdb}.fna
    //
    EOF
    # Prepare organism-params.dat
    cat > "${params.pgdb}/organism-params.dat" <<EOF
    ID	${params.pgdb}
    STORAGE	FILE
    NAME	${params.species}
    ABBREV-NAME ${params.species}
    STRAIN	${params.strain}
    CREATE?	T
    DOMAIN	${params.domain}
    RANK	${params.rank}
    AUTHOR	${params.author}
    EOF
    # Create .pf file
    awk -f "${baseDir}/bin/create_pf.awk" "pangenome.gpr" > "${params.pgdb}/${params.pgdb}.pf"
    """
}

process PREPARE_MODIFIED_PATHWAY_TOOLS {
    input:

    output:
    path "modified-pathway-tools"

    script:
    """
    # Modify pathway-tools to allow a specific PTOOLS_LOCAL_PATH to be used.
    cp \$(which pathway-tools) ./modified-pathway-tools

    # Delete the line that create PTOOLS_LOCAL_PATH since we export it ourself and a check
    sed -i '/PTOOLS_LOCAL_PATH=/d' ./modified-pathway-tools
    sed -i '20i if [[ -z "\${PTOOLS_LOCAL_PATH}" ]]; then echo "Please define PTOOLS_LOCAL_PATH (directory containing the ptools-local directory)"; exit 1; else echo "Using PTOOLS_LOCAL_PATH=\${PTOOLS_LOCAL_PATH}"; fi' ./modified-pathway-tools

    # Force PT_PATCH_FLAG to "-no-patch-download"
    sed -i 's/PT_PATCH_FLAG="  "/PT_PATCH_FLAG="-no-patch-download"/g' ./modified-pathway-tools
    """
}


process PATHOLOGIC {

    input:
    path modified_pathway_tools
    path pgdb_input
    path ptools_local_template

    output:
    path "${pgdb}", emit: ptools_local

    shell:
    """
    # Create a temporary ptools-local path, from pathway-tools ptools-local
    cp -rp "}"

    # Launch Xvfb in background
    DISPLAY=":\$\$"
    export DISPLAY

    Xvfb -ac "\${DISPLAY}" 1>/dev/null 2>&1 &
    PIDXVFB=\$(echo \$!)

    # Run PathoLogic in batch mode


    # Run PathoLogic
    set +e
    ${modified_pathway_tools} -patho "${pgdb_input}" -no-web-cel-overview
    set -e

    # Kill Xxfb
    kill \$PIDXVFB
    """
}

workflow INFER_METABOLIC_NETWORK {

    take:
    gene_protein_reaction

    main:
    PREPARE_PATHOLOGIC_INPUT(gene_protein_reaction)

    // | PATHOLOGIC

}
