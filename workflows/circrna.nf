/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowCircrna.initialise(params, log)

// Check modules paramater
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

def defineModuleList() {
    return [
    'circrna_discovery',
    'mirna_prediction',
    'differential_expression'
    ]
}

moduleList = defineModuleList()
module = params.module ? params.module.split(',').collect{it.trim().toLowerCase()} : []
if(!checkParameterList(module, moduleList)) {
    log.error "error: Unknown module selected, please choose one of the following:\n\n  circrna_discovery\n\n  mirna_prediction\n\n  differential_expression\n\nPlease refer to the help documentation for a description of each module."
    System.exit(1)
}

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

ch_phenotype   = params.phenotype && params.module.contains('differential_expression') ? file(params.phenotype, checkIfExists:true) : Channel.empty()
ch_fasta       = params.fasta ? file(params.fasta) : 'null'
ch_gtf         = params.gtf ? file(params.gtf) : 'null'
ch_mature      = params.mature && params.module.contains('mirna_prediction') ? file(params.mature) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SUBWORKFLOWS:
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { PREPARE_GENOME    } from '../subworkflows/local/prepare_genome'
include { CIRCRNA_DISCOVERY } from '../subworkflows/local/circrna_discovery'
include { CIRCRNA_DISCOVERY_CIRIQUANT } from '../subworkflows/local/circrna_discovery_ciriquant'
include { CIRCRNA_DISCOVERY_CIRCRNA_FINDER } from '../subworkflows/local/circrna_discovery_circrna_finder'
include { CIRCRNA_DISCOVERY_STAR_ALIGN } from '../subworkflows/local/circrna_discovery_star_align'
include { CIRCRNA_DISCOVERY_SEGEMEHL } from '../subworkflows/local/circrna_discovery_segemehl'
include { MIRNA_PREDICTION  } from '../subworkflows/local/mirna_prediction'
include { DIFFERENTIAL_EXPRESSION } from '../subworkflows/local/differential_expression'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// MODULES:
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { ANNOTATION                  } from '../modules/local/annotation/full_annotation/main'
include { FASTA                       } from '../modules/local/fasta/main'
// SUBWORKFLOWS:
include { FASTQC_TRIMGALORE } from '../subworkflows/nf-core/fastqc_trimgalore'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CIRCRNA {

    ch_reports  = Channel.empty()
    ch_versions = Channel.empty()

    //
    // 1. Pre-processing
    //

    // SUBWORKFLOW:
    // Validate input samplesheet & phenotype file
    INPUT_CHECK (
        ch_input,
        ch_phenotype
    )
    .reads
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // MODULE:
    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // SUBORKFLOW:
    // Prepare index files &/or use iGenomes if chosen.
    // NOTE packaging an arg here could enable the correct genome behavior
    PREPARE_GENOME (
        ch_fasta,
        ch_gtf
    )

    // Stage the indices via newly created indices, iGenomes or empty list if tool not selected.
    bowtie_index   = params.fasta ? params.bowtie ? Channel.fromPath(params.bowtie) : PREPARE_GENOME.out.bowtie : []
    bowtie2_index  = params.fasta ? params.bowtie2 ? Channel.fromPath(params.bowtie2).map{ it -> [[id:'bowtie2'], it] } : PREPARE_GENOME.out.bowtie2 : []
    bwa_index      = params.fasta ? params.bwa ? Channel.fromPath(params.bwa).map{ it -> [[id:'bwa'], it] } : PREPARE_GENOME.out.bwa : []
    chromosomes    = params.fasta && ( params.tool.contains('mapsplice') || params.tool.contains('find_circ') ) ? PREPARE_GENOME.out.chromosomes : []
    // NOTE
    // temporary fix
    // hisat2_index   = params.fasta ? params.hisat2 && ( params.tool.contains('ciriquant') || params.module.contains('differential_expression') ) ? Channel.fromPath(params.hisat2).map{ [[id:"hisat2"], it]} : PREPARE_GENOME.out.hisat2 : []
    hisat2_index   = Channel.fromPath(params.hisat2).map{ [[id:"hisat2"], it]}
    star_index     = params.fasta ? params.star ? Channel.fromPath(params.star).map{[[id:'star'], it]}: PREPARE_GENOME.out.star : []
    segemehl_index = params.fasta ? params.segemehl ? Channel.fromPath(params.segemehl) : PREPARE_GENOME.out.segemehl : []
    ch_versions    = ch_versions.mix(PREPARE_GENOME.out.versions)

    // MODULE: Run FastQC, trimgalore!
    ch_filtered_reads = Channel.empty()
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )
    ch_filtered_reads = FASTQC_TRIMGALORE.out.reads
    ch_filtered_reads.view()
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_reports  = ch_reports.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    //
    // 2. circRNA Discovery
    //

    // CIRCRNA_DISCOVERY(
        // FASTQC_TRIMGALORE.out.reads,
        // ch_fasta,
        // ch_gtf,
        // bowtie_index,
        // bowtie2_index,
        // bwa_index,
        // chromosomes,
        // hisat2_index,
        // segemehl_index,
        // star_index,
        // params.bsj_reads,
        // params.tool_filter,
        // params.duplicates_fun,
        // params.exon_boundary
    // )

    ch_star_align_sam = Channel.empty()
    ch_star_align_junction = Channel.empty()
    ch_star_align_tab = Channel.empty()
    if ( ( params.tool.contains('circexplorer2') || params.tool.contains('circrna_finder') ) ) {


        CIRCRNA_DISCOVERY_STAR_ALIGN(
            ch_filtered_reads,
            ch_fasta,
            ch_gtf,
            star_index,
            params.bsj_reads
        )

        ch_star_align_sam = CIRCRNA_DISCOVERY_STAR_ALIGN.out.sam
        ch_star_align_junction = CIRCRNA_DISCOVERY_STAR_ALIGN.out.junction
        ch_star_align_tab = CIRCRNA_DISCOVERY_STAR_ALIGN.out.tab
    }


    ch_segemehl_results = Channel.empty()
    if ( params.tool.contains('segemehl') ) {


        CIRCRNA_DISCOVERY_SEGEMEHL(
            ch_filtered_reads,
            ch_fasta,
            segemehl_index,
            params.bsj_reads
        )

        CIRCRNA_DISCOVERY_CIRCRNA_FINDER.out.circrna_finder_results.view()

        ch_circrna_finder_results = CIRCRNA_DISCOVERY_CIRCRNA_FINDER.out.circrna_finder_results
    }

    ch_circrna_finder_results = Channel.empty()
    if ( params.tool.contains('circrna_finder') ) {


        CIRCRNA_DISCOVERY_CIRCRNA_FINDER(
            ch_star_align_sam,
            ch_star_align_junction,
            ch_star_align_tab,
            ch_fasta,
            params.bsj_reads
        )

        CIRCRNA_DISCOVERY_CIRCRNA_FINDER.out.circrna_finder_results.view()

        ch_circrna_finder_results = CIRCRNA_DISCOVERY_CIRCRNA_FINDER.out.circrna_finder_results
    }

    ch_ciriquant_results = Channel.empty()
    if ( params.tool.contains('ciriquant') ) {

        CIRCRNA_DISCOVERY_CIRIQUANT(
            ch_filtered_reads,
            ch_fasta,
            ch_gtf,
            bwa_index,
            hisat2_index,
            params.bsj_reads
        )

        CIRCRNA_DISCOVERY_CIRIQUANT.out.ciriquant_results.view()

        ch_ciriquant_results = CIRCRNA_DISCOVERY_CIRIQUANT.out.ciriquant_results
    }

    ch_biotypes = Channel.fromPath("${projectDir}/bin/unwanted_biotypes.txt")
    ch_annotation = Channel.empty()

    ch_forAnnotation = ch_ciriquant_results.mix(ch_circrna_finder_results).groupTuple( by:0 )

    ANNOTATION( ch_forAnnotation , ch_gtf, ch_biotypes.collect(), params.exon_boundary )

    ch_annotation = ANNOTATION.out.bed

    FASTA( ch_annotation, ch_fasta )

    //
    // 3. miRNA prediction
    //

    // MIRNA_PREDICTION(
        // CIRCRNA_DISCOVERY.out.fasta,
        // CIRCRNA_DISCOVERY.out.circrna_bed12,
        // ch_mature
    // )

    // ch_versions = ch_versions.mix(MIRNA_PREDICTION.out.versions)

    //
    // 4. Differential expression tests
    //

    // ch_ensembl_database_map = params.module.contains('differential_expression') ? Channel.fromPath("${projectDir}/bin/ensembl_database_map.txt") : Channel.empty()

    // DIFFERENTIAL_EXPRESSION(
        // ch_filtered_reads,
        // ch_gtf,
        // ch_fasta,
        // hisat2_index,
        // PREPARE_GENOME.out.splice_sites,
        // ch_phenotype,
        // CIRCRNA_DISCOVERY.out.dea_matrix,
        // CIRCRNA_DISCOVERY.out.clr_matrix,
        // params.species,
        // ch_ensembl_database_map,
        // params.exon_boundary
    // )

    // ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)
    // ch_reports  = ch_reports.mix(DIFFERENTIAL_EXPRESSION.out.reports)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MODULE: MultiQC
    workflow_summary    = WorkflowCircrna.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowCircrna.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
