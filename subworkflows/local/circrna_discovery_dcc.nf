include { DCC                              } from '../../modules/local/dcc/dcc/main'
include { STAR_ALIGN as DCC_1ST_PASS       } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_2ND_PASS       } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_SJDB                 } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE1_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE1_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE1_SJDB           } from '../../modules/local/star/sjdb/main'
include { STAR_ALIGN as DCC_MATE2_1ST_PASS } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as DCC_MATE2_2ND_PASS } from '../../modules/nf-core/star/align/main'
include { SJDB as DCC_MATE2_SJDB           } from '../../modules/local/star/sjdb/main'
include { DCC_FILTER                       } from '../../modules/local/dcc/filter/main'

workflow CIRCRNA_DISCOVERY_DCC {

    take:
    reads
    fasta
    gtf
    star_index
    bsj_reads

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)
    ch_gtf   = Channel.fromPath(gtf)

    gtf_tuple = ch_gtf.map{ it ->
        meta = [:]
        meta.id = it.simpleName
        return [ meta, [it] ]
    }.collect()

    // Define variables here, star_ignore_sjdbgtf not supposed to be toggled by user.
    star_ignore_sjdbgtf = true
    seq_center     = params.seq_center ?: ''
    seq_platform   = ''

    //
    // DCC WORKFLOW
    //

    DCC_1ST_PASS( reads, star_index.collect(), gtf_tuple, star_ignore_sjdbgtf, seq_platform, seq_center )
    DCC_SJDB( DCC_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "dcc_sjdb"], it]}, bsj_reads )
    DCC_2ND_PASS( reads, star_index.collect(), DCC_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )


    // this should all actually be a .branch
    // maybe the whole workflow actually ....
    // single_end = reads.map{ meta, reads -> return meta.single_end }

    // if ( single_end ) {

        // dcc_stage = DCC_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, meta, junction]}
            // .map{ id, meta, junction -> return [ meta, junction ]}
            // .groupTuple()

    // } else {

    mate1 = reads.map{ meta, reads -> return [ [id: meta.id, single_end: true], reads[0] ] }
    DCC_MATE1_1ST_PASS( mate1, star_index.collect(), gtf_tuple, star_ignore_sjdbgtf, seq_platform, seq_center )
    DCC_MATE1_SJDB( DCC_MATE1_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "mate1_sjdb"], it]}, bsj_reads )
    DCC_MATE1_2ND_PASS( mate1, star_index.collect(), DCC_MATE1_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )

    mate2 = reads.map{ meta, reads -> return [ [id: meta.id, single_end: true], reads[1] ] }
    DCC_MATE2_1ST_PASS( mate2, star_index.collect(), gtf_tuple, star_ignore_sjdbgtf, seq_platform, seq_center )
    DCC_MATE2_SJDB( DCC_MATE2_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "mate2_sjdb"], it]}, bsj_reads )
    DCC_MATE2_2ND_PASS( mate2, star_index.collect(), DCC_MATE2_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )

    dcc_stage = DCC_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, meta, junction]}
        .join(
            DCC_MATE1_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, junction] }
        )
        .join(
            DCC_MATE2_2ND_PASS.out.junction.map{ meta, junction -> return [ meta.id, junction] }
        )
        .map{ id, meta, junction, mate1, mate2 -> return [ meta, junction, mate1, mate2 ]}
        .groupTuple()

    ch_versions = ch_versions.mix(DCC_MATE1_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE1_SJDB.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE1_2ND_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE2_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE2_SJDB.out.versions)
    ch_versions = ch_versions.mix(DCC_MATE2_2ND_PASS.out.versions)
    // }

    dcc = dcc_stage.map{ it -> def meta = it[0]; if( meta.single_end ){ return [ it[0], it[1], [], [] ] } else { return it } }
    DCC( dcc, fasta, gtf )
    DCC_FILTER( DCC.out.txt.map{ meta, txt -> def name = meta.clone(); name.tool = "dcc"; return [ name, txt ] }, bsj_reads )

    ch_versions = ch_versions.mix(DCC.out.versions)
    ch_versions = ch_versions.mix(DCC_FILTER.out.versions)

    emit:
    dcc_results = DCC_FILTER.out.results
    dcc_matrix = DCC_FILTER.out.matrix
    versions = ch_versions
}
