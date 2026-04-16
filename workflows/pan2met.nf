/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pan2met-wf_pipeline'

include { PPANGGOLIN_FASTA_FAMILY_PROTEINS } from "../modules/local/ppanggolin/fasta.nf"
include { PPANGGOLIN_WORKFLOW } from "../modules/local/ppanggolin/workflow.nf"
include { REMOVE_METACYC_META_PREFIX } from "../subworkflows/local/pan2met/cyc.nf"
include { REMOVE_ECOCYC_ECO_PREFIX } from "../subworkflows/local/pan2met/cyc.nf"
include { CYC_BASED_ASSOCIATION as METACYC_BASED_ASSOCIATION } from "../subworkflows/local/pan2met/cyc.nf"
include { CYC_BASED_ASSOCIATION as ECOCYC_BASED_ASSOCIATION } from "../subworkflows/local/pan2met/cyc.nf"
include { KOFAMSCAN_BASED_ASSOCIATION } from "../modules/local/kofamscan"
include { NCBIFAM_BASED_ASSOCIATION } from "../modules/local/ncbifam"
include { MERGE_ASSOCIATION as MERGE_ECOCYC_METACYC } from "../subworkflows/local/pan2met"
include { MERGE_ASSOCIATION_WITH_EC as MERGE_CYC_KOFAMSCAN } from "../subworkflows/local/pan2met"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// params.genomes = null;
// params.pangenome = null;
// params.proteome = null;

workflow PAN2MET {

    // take:
    // ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()

    // Extract the reference proteins from a pangenome
    // Or directly from a protein fasta file
    // Or first build the pangenome with PPanGGOLiN, then extract the reference proteins

    
    // First ensure that if params.genomes is set, neither params.pangenome, nor params.proteome is set
    if ((params.genomes && params.pangenome) || (params.pangenome && params.proteome) || (params.genomes && params.proteome)) {
        error "--genomes, --pangenome and --proteome parameters are mutually exclusive."
    }
    
    if (params.genomes) {
        // 
        // Prepare the pangenome h5 running the PPanGGOLiN workflow
        //
        genomes = file(params.genomes)
        pangenome = PPANGGOLIN_WORKFLOW(genomes).out.pangenome
        ch_versions = ch_versions.mix(PPANGGOLIN_WORKFLOW.out.versions)
        PPANGGOLIN_FASTA_FAMILY_PROTEINS(pangenome)
        ch_proteins = PPANGGOLIN_FASTA_FAMILY_PROTEINS.out.family_proteins
        ch_versions = ch_versions.mix(PPANGGOLIN_FASTA_FAMILY_PROTEINS.out.versions)
    } else if (params.pangenome) {
        //
        // Extract pangenome proteins amino-acid sequences and CDS
        //
        pangenome = file(params.pangenome)
        PPANGGOLIN_FASTA_FAMILY_PROTEINS(pangenome)
        ch_proteins = PPANGGOLIN_FASTA_FAMILY_PROTEINS.out.family_proteins
        ch_versions = ch_versions.mix(PPANGGOLIN_FASTA_FAMILY_PROTEINS.out.versions)
    } else if (params.proteome) {
        // By pass the pangenome reference protein extraction step
        ch_proteins = Channel.fromPath(params.proteome)
    } else {
        error "No input provided. Please set one of --genomes, --pangenome or --proteome input parameters."
    }

    // Prepare the association file using only the selected annotation sources.
    annotation_sources = params.annotations.split(",")
    available_sources = ["ecocyc", "metacyc", "kofamscan", "ncbifam"]
    for (method in annotation_sources) {
        if (!available_sources.contains(method)) {
            error "GPR annotation method not handled: " + method
        }
    }

    print ch_proteins

    
    // Duplicate proteins channel into 4 channels, one for each available annotation source
    
    if (annotation_sources.contains("metacyc")) { 
        metacyc_reference_proteins = file(params.metacyc_proteins)
        metacyc_reference_proteins_formatted = REMOVE_METACYC_META_PREFIX(metacyc_reference_proteins)
        metacyc_monomer_reactions = file(params.metacyc_monomers_to_reactions)
        METACYC_BASED_ASSOCIATION(
            ch_proteins,
            metacyc_reference_proteins_formatted,
            metacyc_monomer_reactions,
            params.metacyc_coverage_threshold,
            params.metacyc_identity_threshold)
        ch_versions = ch_versions.mix(METACYC_BASED_ASSOCIATION.out.versions)
    }

    if (annotation_sources.contains("ecocyc")) {    
        ecocyc_reference_proteins = file(params.ecocyc_proteins)
        ecocyc_reference_proteins_formatted = REMOVE_ECOCYC_ECO_PREFIX(ecocyc_reference_proteins)
        ecocyc_monomer_reactions = file(params.ecocyc_monomers_to_metacyc_reactions)
        ECOCYC_BASED_ASSOCIATION(
            ch_proteins,
            ecocyc_reference_proteins_formatted,
            ecocyc_monomer_reactions,
            params.ecocyc_coverage_threshold,
            params.ecocyc_identity_threshold)

        ch_versions = ch_versions.mix(ECOCYC_BASED_ASSOCIATION.out.versions)
    }
    if (annotation_sources.contains("kofamscan")) {
        KOFAMSCAN_BASED_ASSOCIATION(ch_proteins)
        ch_versions = ch_versions.mix(KOFAMSCAN_BASED_ASSOCIATION.out.versions)
    }

    if (annotation_sources.contains("ncbifam")) {
        NCBIFAM_BASED_ASSOCIATION(ch_proteins)
        ch_versions = ch_versions.mix(NCBIFAM_BASED_ASSOCIATION.out.versions)
    }

    // TODO: make merges depend on selected annotation source 
    MERGE_ECOCYC_METACYC(ECOCYC_BASED_ASSOCIATION.out.asso, METACYC_BASED_ASSOCIATION.out.asso)
    MERGE_CYC_KOFAMSCAN(MERGE_ECOCYC_METACYC.out.asso, KOFAMSCAN_BASED_ASSOCIATION.out.asso)

    //
    // Run metabolism prediction
    //
    // INFER_METABOLIC_NETWORK(MERGE_CYC_KOFAMSCAN.out.asso)
    // ch_versions = ch_versions.mix(INFER_METABOLIC_NETWORK.out.versions)
 
    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'pan2met_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
