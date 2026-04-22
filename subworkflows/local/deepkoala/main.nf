
include { DEEPKOALA } from '../../../modules/local/deepkoala'

process JOIN_KO_METACYC {

    input:
    path deepkoala_csv
    path kegg_kos_to_metacyc_reactions

    script:
    """
    TODO
    """
}

workflow DEEPKOALA_BASED_ASSOCIATION {

    take:
    proteins
    kegg_kos_to_metacyc_reactions

    main:

    ch_versions = Channel.empty()

    DEEPKOALA(proteins)

    

    emit:
    versions = ch_versions
    
}
