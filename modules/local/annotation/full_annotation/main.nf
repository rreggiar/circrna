process ANNOTATION {
    tag "${meta.id}:${meta.tool}"
    label 'process_high'

    conda "bioconda::ucsc-gtftogenepred=377 bioconda::ucsc-genepredtobed=377 bioconda::bedtools=2.27.0"
    // RER
    // breaking change!!!
    // this only works on sherlock
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/home/users/rreggiar/rreggiar_dev/singularity/sra/getreads_mamba_AMD.sif' :
        'quay.io/biocontainers/mulled-v2-d7ee3552d06d8acebbc660507b48487c7369e221:07daadbfe8182aa3c974c7b78924d5c8730b922d-0' }"

    input:
    tuple val(meta), path(bed)
    path gtf
    path biotypes
    val exon_boundary


    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    path("*.log")                         , emit: log
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '377'
    """
    # RER
    # breaking change!!!
    # this only works on sherlock
    source /etc/conda/etc/profile.d/conda.sh
    conda activate rreggiar

    grep -vf $biotypes $gtf > filt.gtf
    mv $bed circs.bed
    # parallel; split bed into file per-line
    split_input_bed.sh circs.bed \$(( ${task.cpus} * 1 ))

    # parallel; process each bed file in parallel up to task.cpus (or more?)
    ls -d splitBed/* | parallel -u -j \$(( ${task.cpus} * 1 )) annotate_outputs_parallel.sh ${exon_boundary} {} &> ${prefix}.log

    # clean and combine output
    aggregate_outputs.sh

    mv master_bed12.bed ${prefix}.bed.tmp

    awk -v FS="\t" '{print \$11}' ${prefix}.bed.tmp > mature_len.tmp
    awk -v FS="," '{for(i=t=0;i<NF;) t+=\$++i; \$0=t}1' mature_len.tmp > mature_length

    paste ${prefix}.bed.tmp mature_length > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version | head -n1 | cut -d' ' -f3 | sed 's/,//g' )
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        ucsc: $VERSION
    END_VERSIONS
    """
}
