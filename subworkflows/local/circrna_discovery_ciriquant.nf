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

    main:
    ch_versions = Channel.empty()
    ch_fasta = Channel.fromPath(fasta)
    ch_gtf   = Channel.fromPath(gtf)

    //
    // CIRIQUANT WORKFLOW:
    //

    // only need path to bwa, only need path to hisat2.
    // do not want to upset the collect declr for all indices just for this.
    CIRIQUANT_YML( gtf, fasta, bwa_index.map{ meta, index -> return index }, hisat2_index.map{ meta, index -> return index } )
    CIRIQUANT( reads, CIRIQUANT_YML.out.yml.collect() )
    CIRIQUANT_FILTER( CIRIQUANT.out.gtf.map{ meta, gtf -> meta + [tool: "ciriquant"]; return [ meta, gtf ] }, bsj_reads )

    ch_versions = ch_versions.mix(CIRIQUANT.out.versions)
    ch_versions = ch_versions.mix(CIRIQUANT_FILTER.out.versions)

    emit:
    ciriquant_results = CIRIQUANT_FILTER.out.results
    versions = ch_versions
}
