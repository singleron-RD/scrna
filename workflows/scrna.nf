/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FILTER_GTF             } from '../modules/local/filter_gtf'
include { STAR_GENOME            } from '../modules/local/star_genome'
include { PROTOCOL_CMD           } from '../modules/local/protocol_cmd'
include { STARSOLO               } from '../modules/local/starsolo'
include { CELL_CALLING           } from '../modules/local/cell_calling'
include { STARSOLO_SUMMARY       } from '../modules/local/starsolo_summary'
include { SUBSAMPLE              } from '../modules/local/subsample'
include { MULTIQC                } from '../modules/local/multiqc_sgr'

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scrna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow SCRNA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // fastqc
    if (params.run_fastqc) {
        FASTQC (
            ch_samplesheet
        )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

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
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_versions = ch_versions.mix(PROTOCOL_CMD.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(PROTOCOL_CMD.out.json.collect{it[1]})

    // starsolo
    ch_merge = ch_samplesheet.join(PROTOCOL_CMD.out.starsolo_cmd.map{ [it[0], it[1].text] })
    STARSOLO (
        ch_merge,
        star_genome,
        "${projectDir}/assets/",
    )
    ch_versions = ch_versions.mix(STARSOLO.out.versions.first())

    // cell-calling
    CELL_CALLING (
        STARSOLO.out.raw_matrix,
        params.soloCellFilter,
    )
    ch_versions = ch_versions.mix(CELL_CALLING.out.versions.first())

    // statsolo summary
    ch_merge = STARSOLO.out.read_stats.join(STARSOLO.out.summary).join(CELL_CALLING.out.filtered_matrix)           
    STARSOLO_SUMMARY ( ch_merge )
    ch_multiqc_files = ch_multiqc_files.mix(STARSOLO_SUMMARY.out.json.collect{it[1]})

    // subsample
    if (params.run_subsample) {
        ch_merge = STARSOLO.out.bam_sorted.join(CELL_CALLING.out.barcodes)                
        SUBSAMPLE ( ch_merge )
        ch_multiqc_files = ch_multiqc_files.mix(SUBSAMPLE.out.json.collect{it[1]})
    }

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
        "${projectDir}/multiqc_sgr/singleron_logo.png",
        "${projectDir}/multiqc_sgr/",
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
