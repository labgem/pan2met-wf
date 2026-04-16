
include { HMMER_HMMSCAN } from "../../../modules/local/hmmer/hmmscan/main"
include { GUNZIP } from "../../../modules/nf-core/gunzip/main"

process JOIN_EC {

    input:
    path ncbifam_domtbl_gz
    path hmm_pgap_tsv

    output:
    path "${ncbifam_domtbl_gz.baseName}_with_ec.tsv", emit: asso

    shell:
    """
    cut -d\$'\\t' -f1,14 "${hmm_pgap_tsv}" > "ncbifam_to_ec.tsv"
    sort --field-separator=\$'\\t' --key=1 "ncbifam_to_ec.tsv" > "ncbifam_to_ec_sorted.tsv"
    # Convert the domtbl to a proper TSV containing only NCBIFam accession and query name
    awk -v OFS="\t" '/^[^#]/ { print \$2, \$4 }' <(zcat "${ncbifam_domtbl_gz}") > "ncbifam_protein_hit.tsv"
    sort --field-separator=\$'\\t' --key=1 "ncbifam_protein_hit.tsv" > "ncbifam_protein_hit_sorted.tsv"
    join -t \$'\\t' -1 1 -2 1 -a 1 "ncbifam_protein_hit_sorted.tsv" "ncbifam_to_ec_sorted.tsv" > "${ncbifam_domtbl_gz.baseName}_with_ec.tsv"
    """
}

process EXTRACT_FILE {
    input:
    tuple val(meta), path(file)
    output:
    path("${file}"), emit: file
    shell:
    """
    """
}

workflow NCBIFAM_BASED_ASSOCIATION {
    take:
    proteins
    
    main:
    meta = [ id: null ]
    write_align = false
    write_target = true
    write_domain = true
    write_pfam = true
    hmmdb = params.ncbifam_db
    HMMER_HMMSCAN([meta, hmmdb, proteins, write_align, write_target, write_domain, write_pfam])
    EXTRACT_FILE(HMMER_HMMSCAN.out.domain_summary)
    JOIN_EC(EXTRACT_FILE.out.file, params.ncbifam_table)
    
    emit:
    asso = JOIN_EC.out.asso
    versions = HMMER_HMMSCAN.out.versions
}


