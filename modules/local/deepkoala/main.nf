process DEEPKOALA_ASSOCIATION {

    input:
    path(deepkoala_output_csv)

    output:
    path("${deepkoala_output_csv.baseName}_metacyc_ec.tsv"), emit: asso

    shell:
    """
    # Convert the deepkoala output CSV to a TSV
    sed 's/,/\t/g' "${deepkoala_output_csv}" > "deepkoala_output.tsv"
    # Sort
    sort --field-separator=\$'\\t' --key=2 "deepkoala_output.tsv" > "deepkoala_sorted_by_ko.tsv" # Sort on KO
    sort --field-separator=\$'\\t' --key=1 "${params.kegg_kos_to_metacyc_reactions}" | cut -d\$'\\t' -f1,3 > "sorted_kegg_kos_to_metacyc_reactions.tsv"
    sort --field-separator=\$'\\t' --key=1 "${params.kegg_kos_to_ec_numbers}" > "sorted_kegg_kos_to_ec_numbers.tsv"
    # Left join to add MetaCyc reactions
    join -t \$'\\t' -1 2 -2 1 -a 1 "deepkoala_sorted_by_ko.tsv" "sorted_kegg_kos_to_metacyc_reactions.tsv" > "deepkoala_metacyc.tsv"
    # Left join to add EC-numbers
    join -t \$'\\t' -1 1 -2 1 -a 1 "deepkoala_metacyc.tsv" "sorted_kegg_kos_to_ec_numbers.tsv" | cut -d\$'\\t' -f2,3,4 > "${deepkoala_output_csv.baseName}_metacyc_ec.tsv"
    """
}

process DEEPKOALA {

    // conda "environment.yml"

    input:
    path(input_fasta)

    output:
    path("${output_csv}"), emit: csv
    // path 'versions.yml', emit: versions // deepkoala is in beta and does not seem to provide a --version or similar option.
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    output_csv = "${input_fasta.baseName}.deepkoala.csv"
    """
    
    python3 -m deepkoala.cli \\
        -i "${input_fasta}" \\
        -o "${output_csv}" \\
        --model "full" \\
        $args \\
        --date "latest" \\
        --batch_size 8 \\
        --num_workers 4 \\
        --topk 1 
    """
}

workflow DEEPKOALA_BASED_ASSOCIATION {
    take:
    proteins
    
    main:
    DEEPKOALA(proteins)
    DEEPKOALA_ASSOCIATION(DEEPKOALA.out.csv)
    emit:
    asso = DEEPKOALA_ASSOCIATION.out.asso
}


