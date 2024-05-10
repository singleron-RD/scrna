/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                } from '../modules/local/multiqc_sgr/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastqc=0.12.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}.${entry}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads $task.cpus \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

process FILTER_GTF {
    tag "$gtf"
    label 'process_single'

    conda 'conda-forge::python==3.12'
    container "biocontainers/python:3.12"

    input:
    //
    // Input reads are expected to come as: [ meta, [ pair1_read1, pair1_read2, pair2_read1, pair2_read2 ] ]
    // Input array for a sample is created in the same order reads appear in samplesheet as pairs from replicates are appended to array.
    //
    path gtf
    val attributes

    output:
    path "*.filtered.gtf", emit: filtered_gtf
    path "gtf_filter.log", emit: log_file

    script:
    def args = task.ext.args ?: ''

    """
    filter_gtf.py ${gtf} \"${attributes}\"
    """
}

process STAR_GENOME {
    tag "$genome_name"
    label 'process_medium'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    path fasta
    path gtf
    val genome_name

    output:
    path "$genome_name"            , emit: index
    path "versions.yml"            , emit: versions

    script:
    def args        = task.ext.args ?: ''
    def memory      = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 100000000}" : ''
    def include_gtf = gtf ? "--sjdbGTFfile $gtf" : ''
    def fasta_sa = ( Math.log(fasta.size()) / Math.log(2) ) / 2 - 1
    def sa = Math.floor( Math.min(14, fasta_sa) )
    """
    mkdir ${genome_name}
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir ${genome_name}/ \\
        --genomeFastaFiles $fasta \\
        $include_gtf \\
        --runThreadN $task.cpus \\
        --genomeSAindexNbases ${sa} \\
        $memory \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        star: \$(STAR --version | sed -e "s/STAR_//g")
    END_VERSIONS
    """
}

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
    tuple val(meta), path(reads), emit: reads
    path "${meta.id}.starsolo_cmd.txt", emit: starsolo_cmd
    path "${meta.id}.protocol.txt", emit: parsed_protocol
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

process STARSOLO {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::star==2.7.11b"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/star:2.7.11b--h43eeafb_0' :
        'biocontainers/star:2.7.11b--h43eeafb_0' }"

    input:
    tuple val(meta), path(reads)
    path index
    path assets_dir
    val starsolo_cmd

    output:
    tuple val(meta), path("${meta.id}.matrix/")       , emit: matrix
    tuple val(meta), path('*d.out.bam')               , emit: bam
    tuple val(meta), path('*.Solo.out')               , emit: solo_out
    path('*Log.final.out')                            , emit: log_final
    path "*.Solo.out/GeneFull_Ex50pAS/Summary.csv"    , emit: summary
    tuple val(meta), path("*.Solo.out/GeneFull_Ex50pAS/CellReads.stats")    , emit: read_stats
    path "${meta.id}.matrix/filtered/barcodes.tsv.gz" , emit: barcodes
    path  "versions.yml"                      , emit: versions

    tuple val(meta), path('*sortedByCoord.out.bam')  , optional:true, emit: bam_sorted
    tuple val(meta), path('*toTranscriptome.out.bam'), optional:true, emit: bam_transcript
    tuple val(meta), path('*Aligned.unsort.out.bam') , optional:true, emit: bam_unsorted
    tuple val(meta), path('*out.mate')               , optional:true, emit: unmap
    tuple val(meta), path('*.tab')                   , optional:true, emit: tab

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"

    // do not indent
"""
${starsolo_cmd} --runThreadN ${task.cpus}

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


process STARSOLO_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda 'conda-forge::pandas==1.5.2'
    container "biocontainers/pandas:1.5.2"

    input:
    tuple val(meta), path(read_stats)
    path summary
    path barcodes

    output:
    path "${meta.id}.read_stats.json", emit: read_stats_json
    path "${meta.id}.summary.json", emit: summary_json
    path "${meta.id}.umi_count.json", emit: umi_count_json

    script:

    """
    starsolo_summary.py \\
        --read_stats ${read_stats} \\
        --barcodes ${barcodes} \\
        --summary ${summary} \\
        --sample ${meta.id}
    """
}

workflow SCRNA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // STAR genome
    def star_genome = null
    if (params.star_genome) {
        star_genome = params.star_genome
    } else {
        FILTER_GTF(
            params.gtf,
            params.keep_attributes,
        )
        ch_gtf = FILTER_GTF.out.filtered_gtf
        STAR_GENOME(
            params.fasta,
            ch_gtf,
            params.genome_name
        )
        ch_versions = ch_versions.mix(STAR_GENOME.out.versions.first())
        star_genome = STAR_GENOME.out.index
    }

    // create cmd
    PROTOCOL_CMD (
        ch_samplesheet,
        star_genome,
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_versions = ch_versions.mix(PROTOCOL_CMD.out.versions.first())
    ch_reads = PROTOCOL_CMD.out.reads
    ch_starsolo_cmd = PROTOCOL_CMD.out.starsolo_cmd.map { it.text }

    // starsolo
    STARSOLO (
        ch_reads,
        star_genome,
        "${projectDir}/assets/",
        ch_starsolo_cmd,
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO.out.log_final.collect())

    // statsolo summary
    STARSOLO_SUMMARY (
        STARSOLO.out.read_stats,
        STARSOLO.out.summary,
        STARSOLO.out.barcodes,
    )
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO_SUMMARY.out.read_stats_json.collect()).mix(STARSOLO_SUMMARY.out.summary_json.collect()).mix(STARSOLO_SUMMARY.out.umi_count_json.collect())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
