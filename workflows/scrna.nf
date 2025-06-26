/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
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


process MKREF {
    cpus params.thread
    memory params.limitBAMsortRAM

    container "quay.io/singleron-rd/celescope:v2.6.1"

    input:
    path fasta
    path gtf
    val genome_name

    output:
    path "${genome_name}_filtered", emit: genomeDir


    script:
    """
    celescope utils mkgtf ${gtf} ${gtf}.filtered
    celescope rna mkref \
    --genome_name ${genome_name}_filtered \
    --fasta ${fasta} \
    --gtf ${gtf}.filtered \
    --mt_gene_list mt_gene_list.txt \
    --thread ${params.thread}

    mkdir ${genome_name}_filtered && find . -maxdepth 1 -type f -exec mv {} ${genome_name}_filtered/ \\;
    """

}


process CELESCOPE {
    tag "$meta.id"
    cpus params.thread
    memory params.limitBAMsortRAM

    container "quay.io/singleron-rd/celescope:v2.6.1"

    input:
    tuple val(meta), path(reads, stageAs: "?/*")
    path genomeDir

    output:
    tuple val(meta), path("${meta.id}/${meta.id}_report.html"), emit: report
    tuple val(meta), path("${meta.id}/outs/filtered"), emit: filtered_matrix

    script:
    def (r1, r2) = reads.collate(2).transpose()
    r1 = r1.join(",")
    r2 = r2.join(",")
    
    """
    celescope rna sample --outdir ./${meta.id}/00.sample --sample ${meta.id} --chemistry ${params.chemistry}  --fq1 ${r1}
    celescope rna starsolo --outdir .//${meta.id}/01.starsolo --sample ${meta.id} --thread ${params.thread} --chemistry ${params.chemistry} --adapter_3p AAAAAAAAAAAA --genomeDir ${genomeDir} --outFilterMatchNmin ${params.outFilterMatchNmin} --soloCellFilter "${params.soloCellFilter}" --limitBAMsortRAM ${params.limitBAMsortRAM} --soloFeatures "GeneFull_Ex50pAS Gene" --soloCBmatchWLtype 1MM --report_soloFeature GeneFull_Ex50pAS  --fq1 ${r1} --fq2 ${r2}
    celescope rna analysis --outdir .//${meta.id}/02.analysis --sample ${meta.id} --thread ${params.thread} --genomeDir ${genomeDir} --matrix_file .//${meta.id}/outs/filtered 
    """
}


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
    def genomeDir = null
    if (params.genomeDir) {
        genomeDir = params.genomeDir
    }  else {
        MKREF (
            file(params.fasta, checkIfExists: true),
            file(params.gtf, checkIfExists: true),
            params.genome_name,
        )
        genomeDir = MKREF.out.genomeDir
    }

    // celescope
    CELESCOPE(
        ch_samplesheet,
        genomeDir,
    )


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
