process PROTOCOL_CMD {
    tag "$meta.id"
    label 'process_single'

    conda 'bioconda::pyfastx=2.1.0'
    container "biocontainers/pyfastx:2.1.0--py39h3d4b85c_0"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    tuple val(meta), path(reads)
    path index
    path assets_dir
    val protocol

    output:
    tuple val(meta), path("${meta.id}.starsolo_cmd.txt"), emit: starsolo_cmd
    tuple val(meta), path('*.json'),  emit: json
    path  "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.id}"

    // separate forward from reverse pairs
    def (forward, reverse) = reads.collate(2).transpose()
    """
    protocol_cmd.py \\
        --sample ${prefix} \\
        --genomeDir ${index} \\
        --fq1 ${forward.join( "," )} \\
        --fq2 ${reverse.join( "," )} \\
        --assets_dir ${assets_dir} \\
        --protocol ${protocol} \\
        --pattern ${params.pattern} \\
        --whitelist \"${params.whitelist}\" \\
        --ext_args \"${args}\" \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyfastx: \$(pyfastx --version | sed -e "s/pyfastx version //g")
    END_VERSIONS
    """
}