#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.fasta = null
params.fastq = null
params.index = null
params.outdir = 'results'
params.threads = 4
params.primary_only = true
params.config = null

process PARALLAX_ALIGN {
    tag "${fastq.baseName}"
    cpus params.threads
    
    input:
    path fasta
    path fastq
    path index
    path config
    
    output:
    tuple val(fastq.baseName), path("${fastq.baseName}.sorted.bam"), emit: bam
    
    script:
    def primary_flag = params.primary_only ? '-p' : ''
    def config_flag = config.name != 'NO_CONFIG' ? "-c ${config}" : ''
    def index_flag = index.name != 'NO_INDEX' ? "-x ${index}" : ''
    """
    parallax align \\
        ${fasta} \\
        ${fastq} \\
        ${primary_flag} \\
        ${index_flag} \\
        ${config_flag} \\
        -t ${task.cpus} \\
        | samtools sort -@ ${task.cpus} -o ${fastq.baseName}.sorted.bam -
    """
}

process SAMTOOLS_INDEX {
    tag "${sample_id}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    
    output:
    tuple val(sample_id), path(bam), path("${bam}.bai"), emit: indexed_bam
    
    script:
    """
    samtools index ${bam}
    """
}

workflow {
    // Validate required parameters
    if (!params.fasta) {
        error "Please specify --fasta"
    }
    if (!params.fastq) {
        error "Please specify --fastq"
    }
    
    // Create channels
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true)
    fastq_ch = channel.fromPath(params.fastq, checkIfExists: true)
    
    // Handle optional index directory
    index_ch = params.index 
        ? channel.fromPath(params.index, type: 'dir', checkIfExists: true)
        : channel.fromPath('NO_INDEX')
    
    // Handle optional config file
    config_ch = params.config
        ? channel.fromPath(params.config, checkIfExists: true)
        : channel.fromPath('NO_CONFIG')
    
    // Run pipeline
    PARALLAX_ALIGN(fasta_ch, fastq_ch, index_ch, config_ch)
    SAMTOOLS_INDEX(PARALLAX_ALIGN.out.bam)
}
