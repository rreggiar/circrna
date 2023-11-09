include { CIRIQUANT_YML                    } from '../../modules/local/ciriquant/yml/main'
include { CIRIQUANT                        } from '../../modules/local/ciriquant/ciriquant/main'
include { CIRIQUANT_FILTER                 } from '../../modules/local/ciriquant/filter/main'

workflow CIRCRNA_DISCOVERY_CIRIQUANT {

    take:
    reads
    fasta
    gtf
    bwa_index
    hisat2_index
    bsj_reads
    tool_filter
    duplicates_fun
    exon_boundary

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)
    ch_gtf   = Channel.fromPath(gtf)

    ch_fasta.map{ it ->
        meta = [:]
        meta.id = it.simpleName
        return [ meta, [it] ]
    }.set{ fasta_tuple }

    gtf_tuple = ch_gtf.map{ it ->
        meta = [:]
        meta.id = it.simpleName
        return [ meta, [it] ]
    }.collect()

    //
    // CIRIQUANT WORKFLOW:
    //

    // only need path to bwa, only need path to hisat2.
    // do not want to upset the collect declr for all indices just for this.
    CIRIQUANT_YML( gtf, fasta, bwa_index.map{ meta, index -> return index }, hisat2_index.map{ meta, index -> return index } )
    CIRIQUANT( reads, CIRIQUANT_YML.out.yml.collect() )
    CIRIQUANT_FILTER( CIRIQUANT.out.gtf.map{ meta, gtf -> meta.tool = "ciriquant"; return [ meta, gtf ] }, bsj_reads )

    ch_versions = ch_versions.mix(CIRIQUANT.out.versions)
    ch_versions = ch_versions.mix(CIRIQUANT_FILTER.out.versions)

    emit:
    ciriquant_out = CIRIQUANT_FILTER.out
    // circrna_bed12 = ANNOTATION.out.bed
    // fasta = FASTA.out.analysis_fasta
    versions = ch_versions
    // dea_matrix
    // clr_matrix
}
