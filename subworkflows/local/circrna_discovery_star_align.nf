include { STAR_ALIGN as STAR_1ST_PASS      } from '../../modules/nf-core/star/align/main'
include { STAR_ALIGN as STAR_2ND_PASS      } from '../../modules/nf-core/star/align/main'
include { SJDB as STAR_SJDB                } from '../../modules/local/star/sjdb/main'

workflow CIRCRNA_DISCOVERY_STAR_ALIGN {

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

    //
    // STAR WORFKLOW:
    //

    // Define variables here, star_ignore_sjdbgtf not supposed to be toggled by user.
    star_ignore_sjdbgtf = true
    seq_center     = params.seq_center ?: ''
    seq_platform   = ''

    STAR_1ST_PASS( reads, star_index.collect(), gtf_tuple, star_ignore_sjdbgtf, seq_platform, seq_center)
    sjdb = STAR_1ST_PASS.out.tab.map{ meta, tab -> return tab }.collect().map{[[id: "star_sjdb"], it]}
    STAR_SJDB( sjdb, bsj_reads )
    STAR_2ND_PASS( reads, star_index.collect(), STAR_SJDB.out.sjtab, star_ignore_sjdbgtf, seq_platform, seq_center )

    ch_versions = ch_versions.mix(STAR_1ST_PASS.out.versions)
    ch_versions = ch_versions.mix(STAR_2ND_PASS.out.versions)


    emit:
    sam = STAR_2ND_PASS.out.sam
    junction = STAR_2ND_PASS.out.junction
    tab = STAR_2ND_PASS.out.tab
    versions = ch_versions
}
