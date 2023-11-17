include { SEGEMEHL_ALIGN                   } from '../../modules/nf-core/segemehl/align/main'
include { SEGEMEHL_FILTER                  } from '../../modules/local/segemehl/filter/main'

workflow CIRCRNA_DISCOVERY_SEGEMEHL {

    take:
    reads
    fasta
    segemehl_index
    bsj_reads

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)

    //
    // SEGEMEHL WORKFLOW:
    //

    SEGEMEHL_ALIGN( reads, fasta, segemehl_index.collect() )
    segemehl_filter = SEGEMEHL_ALIGN.out.results.map{ meta, results ->  def name = meta.clone(); name.tool = "segemehl"; return [ name, results ] }

    SEGEMEHL_FILTER( segemehl_filter, bsj_reads )

    ch_versions = ch_versions.mix(SEGEMEHL_ALIGN.out.versions)
    ch_versions = ch_versions.mix(SEGEMEHL_FILTER.out.versions)

    emit:
    segemehl_results = SEGEMEHL_FILTER.out.results
    segemehl_matrix = SEGEMEHL_FILTER.out.matrix
    versions = ch_versions
}
