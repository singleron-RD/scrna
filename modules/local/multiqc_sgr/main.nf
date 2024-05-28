process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc==1.21"
    container "biocontainers/multiqc:1.21--pyhdfd78af_0"
    containerOptions '--env HOME=/tmp'

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
    path(multiqc_plugin)

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''
    def is_docker = workflow.profile.tokenize(',').intersect(['docker', 'singularity']).size() >= 1
    def pipuser = is_docker ? '--user' : ''
    def pythonuser = is_docker ? 'export PYTHONNOUSERSITE=0' : 'export PYTHONNOUSERSITE=1'
    """
    cp -r -L ./${multiqc_plugin} ./multiqc.plugin
    if [ -d ./multiqc.plugin/build ]; then rm -r ./multiqc.plugin/build; fi
    pip install ./multiqc.plugin --no-cache-dir --no-deps --force-reinstall $pipuser
    $pythonuser
    python -m site
    multiqc \\
        --force \\
        $args \\
        $config \\
        $extra_config \\
        $logo \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
