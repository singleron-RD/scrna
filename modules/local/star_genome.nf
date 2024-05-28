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