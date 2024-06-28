process STARSOLO_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    input:
    tuple val(meta), path(read_stats), path(summary), path(filtered_matrix)

    output:
    tuple val(meta), path("*.json"), emit: json

    script:

    """
    starsolo_summary.py \\
        --read_stats ${read_stats} \\
        --filtered_matrix ${filtered_matrix} \\
        --summary ${summary} \\
        --sample ${meta.id}
    """
}