process SJDB {
    label 'process_single'

    input:
    path(sjdb)
    val(bsj_reads)

    output:
    path("dataset.SJ.out.tab"), emit: sjtab

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat *.tab | awk -v BSJ=${bsj_reads} '(\$7 >= BSJ && \$6==0)' | cut -f1-6 | sort | uniq > dataset.SJ.out.tab
    """
}
