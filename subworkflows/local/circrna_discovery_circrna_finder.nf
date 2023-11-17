include { CIRCRNA_FINDER_FILTER            } from '../../modules/local/circrna_finder/filter/main'

workflow CIRCRNA_DISCOVERY_CIRCRNA_FINDER {

    take:
    sam
    junction
    tab
    fasta
    bsj_reads

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)

    //
    // CIRCRNA_FINDER WORKFLOW:
    //

    circrna_finder_stage = sam.join( junction ).join( tab )
    // circrna_finder_filter = circrna_finder_stage.map{ meta, sam, junction, tab -> meta + [tool: "circrna_finder"] ; return [ meta, sam, junction, tab ] }.groupTuple( by:0 )
    circrna_finder_filter = circrna_finder_stage.map{ meta, sam, junction, tab -> def name = meta.clone(); name.tool = "circrna_finder"; return [ name, sam, junction, tab ] }
    CIRCRNA_FINDER_FILTER( circrna_finder_filter, fasta, bsj_reads )

    ch_versions = ch_versions.mix(CIRCRNA_FINDER_FILTER.out.versions)

    emit:
    circrna_finder_results = CIRCRNA_FINDER_FILTER.out.results
    circrna_finder_matrix = CIRCRNA_FINDER_FILTER.out.matrix
    versions = ch_versions
}
