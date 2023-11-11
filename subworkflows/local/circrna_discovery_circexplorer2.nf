include { CIRCEXPLORER2_REFERENCE as CIRCEXPLORER2_REF } from '../../modules/local/circexplorer2/reference/main'
include { CIRCEXPLORER2_PARSE as CIRCEXPLORER2_PAR } from '../../modules/nf-core/circexplorer2/parse/main'
include { CIRCEXPLORER2_ANNOTATE as CIRCEXPLORER2_ANN } from '../../modules/nf-core/circexplorer2/annotate/main'
include { CIRCEXPLORER2_FILTER as CIRCEXPLORER2_FLT } from '../../modules/local/circexplorer2/filter/main'

workflow CIRCRNA_DISCOVERY_CIRCEXPLORER2 {

    take:
    fasta
    junction
    gtf
    bsj_reads

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)
    ch_gtf   = Channel.fromPath(gtf)

    //
    // CIRCEXPLORER2 WORKFLOW:
    //

    CIRCEXPLORER2_REF( gtf )
    CIRCEXPLORER2_PAR( STAR_2ND_PASS.out.junction )
    CIRCEXPLORER2_ANN( CIRCEXPLORER2_PAR.out.junction, fasta, CIRCEXPLORER2_REF.out.txt )
    circexplorer2_filter = CIRCEXPLORER2_ANN.out.txt.map{ meta, txt -> def name = meta.clone(); name.tool = "circexplorer2"; return [ name, txt ] }
    CIRCEXPLORER2_FLT( circexplorer2_filter, bsj_reads )

    ch_versions = ch_versions.mix(CIRCEXPLORER2_REF.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_PAR.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_ANN.out.versions)
    ch_versions = ch_versions.mix(CIRCEXPLORER2_FLT.out.versions)

    emit:
    circexplorer2_results = CIRCEXPLORER2_FLT.out.results
    circexplorer2_matrix = CIRCEXPLORER2_FLT.out.matrix
    versions = ch_versions
}
