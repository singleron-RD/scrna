process CELL_CALLING {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    tuple val(meta), path(raw_matrix)
    val soloCellFilter

    output:
    tuple val(meta), path("${meta.id}.matrix/")       , emit: matrix
    tuple val(meta), path("${meta.id}.matrix/filtered/barcodes.tsv.gz") , emit: barcodes
    path  "versions.yml"                      , emit: versions

    script:
    def prefix = "${meta.id}"
    def matrix_dir = "${prefix}.matrix"
    def filtered_matrix = "${matrix_dir}/filtered"
    """
    STAR \\
        --runMode soloCellFiltering \\
        ${raw_matrix} \\
        ${filtered_matrix}/ \\
        --soloCellFilter ${soloCellFilter} 

    cp -r -L ${raw_matrix} ${matrix_dir}/raw
    gzip ${matrix_dir}/*/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}