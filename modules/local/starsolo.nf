process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    tuple val(meta), path(reads, stageAs: "?/*"), val(starsolo_cmd)
    path index
    path assets_dir

    output:
    tuple val(meta), path("${meta.id}.matrix/")       , emit: matrix
    tuple val(meta), path('*d.out.bam')               , emit: bam
    tuple val(meta), path('*.Solo.out')               , emit: solo_out
    path('*Log.final.out')                            , emit: log_final
    path "*.Solo.out/GeneFull_Ex50pAS/Summary.csv"    , emit: summary
    tuple val(meta), path("*.Solo.out/GeneFull_Ex50pAS/CellReads.stats")    , emit: read_stats
    path "${meta.id}.matrix/filtered/barcodes.tsv.gz" , emit: barcodes
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*out.mate')               , optional:true, emit: unmap
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    def (forward, reverse) = reads.collate(2).transpose()
    def args = task.ext.args ?: ''
    // do not indent
    """
    STAR \\
        ${starsolo_cmd} \\
        --readFilesIn ${reverse.join( "," )} ${forward.join( "," )} \\
        --genomeDir $index \\
        --outFileNamePrefix $prefix. \\
        --runThreadN ${task.cpus} \\
        $args

    if [ -d ${prefix}.Solo.out ]; then
        # Backslashes still need to be escaped (https://github.com/nextflow-io/nextflow/issues/67)
        find ${prefix}.Solo.out \\( -name "*.tsv" -o -name "*.mtx" \\) -exec gzip -f {} \\;
    fi

    mkdir ${prefix}.matrix
    mv ${prefix}.Solo.out/GeneFull_Ex50pAS/{raw,filtered} ./${prefix}.matrix/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}