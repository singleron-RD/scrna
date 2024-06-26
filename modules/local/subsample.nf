process SUBSAMPLE {
    tag "$meta.id"
    cpus '1'
    memory '32 GB'

    conda 'bioconda::pysam==0.22.1'
    container "biocontainers/pysam:0.22.1--py38h15b938a_0"

    input:
    tuple val(meta), path(bam), path(barcodes)

    output:
    tuple val(meta), path("*.json"), emit: json

    script:

    """
    subsample.py \\
        -b ${bam} \\
        -c ${barcodes} \\
        -s ${meta.id}
    """
}