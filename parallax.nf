#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Input parameters
params.fasta = null
params.fastq = null
params.samplesheet = null  // CSV with: sample,fastq,platform,library,lane
params.index = null
params.outdir = 'results'
params.threads = 4
params.primary_only = true
params.config = null

// Read group parameters (for single-sample mode)
params.sample = null      // --rg-sm
params.platform = null    // --rg-pl (ILLUMINA, PACBIO, ONT)
params.library = null     // --rg-lb
params.platform_unit = null  // --rg-pu
params.center = null      // --rg-cn

// Helper function to build read group arguments
def buildRgArgs(meta) {
    def args = []
    if (meta.id)       args << "--rg-id ${meta.id}"
    if (meta.sample)   args << "--rg-sm ${meta.sample}"
    if (meta.platform) args << "--rg-pl ${meta.platform}"
    if (meta.library)  args << "--rg-lb ${meta.library}"
    if (meta.pu)       args << "--rg-pu ${meta.pu}"
    if (meta.center)   args << "--rg-cn ${meta.center}"
    return args.join(' ')
}

process PARALLAX_ALIGN {
    tag "${meta.id}"
    cpus params.threads
    
    input:
    tuple val(meta), path(fastq)
    path fasta
    path index
    path config
    
    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    
    script:
    def primary_flag = params.primary_only ? '-p' : ''
    def config_flag = config.name != 'NO_CONFIG' ? "-c ${config}" : ''
    def index_flag = index.name != 'NO_INDEX' ? "-x ${index}" : ''
    def rg_args = buildRgArgs(meta)
    """
    parallax align \\
        ${fasta} \\
        ${fastq} \\
        ${primary_flag} \\
        ${index_flag} \\
        ${config_flag} \\
        ${rg_args} \\
        -t ${task.cpus} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.sorted.bam -
    """
}

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: indexed_bam
    
    script:
    """
    samtools index ${bam}
    """
}

// Parse samplesheet CSV into channel of [meta, fastq] tuples
def parseSamplesheet(samplesheet) {
    return channel
        .fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.sample + (row.lane ? "_${row.lane}" : ""),
                sample: row.sample,
                platform: row.platform ?: 'UNKNOWN',
                library: row.library ?: row.sample,
                pu: row.lane ? "${row.sample}.${row.lane}" : null,
                center: row.center ?: null
            ]
            def fastq = file(row.fastq, checkIfExists: true)
            return [ meta, fastq ]
        }
}

// Create meta map from params (single-sample mode)
def metaFromParams(fastq_path) {
    def fastq = file(fastq_path)
    def sample_name = params.sample ?: fastq.baseName
    return [
        id: sample_name,
        sample: sample_name,
        platform: params.platform,
        library: params.library ?: sample_name,
        pu: params.platform_unit,
        center: params.center
    ]
}

workflow {
    // Validate parameters
    if (!params.fasta) {
        error "Please specify --fasta"
    }
    if (!params.fastq && !params.samplesheet) {
        error "Please specify --fastq or --samplesheet"
    }
    if (params.fastq && params.samplesheet) {
        error "Please specify either --fastq or --samplesheet, not both"
    }
    
    // Create input channel based on mode
    if (params.samplesheet) {
        // Batch mode: read from samplesheet
        reads_ch = parseSamplesheet(params.samplesheet)
    } else {
        // Single-sample mode: use params
        def meta = metaFromParams(params.fastq)
        reads_ch = channel.of([ meta, file(params.fastq, checkIfExists: true) ])
    }
    
    // Reference and optional files
    fasta_ch = channel.fromPath(params.fasta, checkIfExists: true).first()
    
    index_ch = params.index 
        ? channel.fromPath(params.index, type: 'dir', checkIfExists: true).first()
        : channel.fromPath('NO_INDEX').first()
    
    config_ch = params.config
        ? channel.fromPath(params.config, checkIfExists: true).first()
        : channel.fromPath('NO_CONFIG').first()
    
    // Run pipeline
    PARALLAX_ALIGN(reads_ch, fasta_ch, index_ch, config_ch)
    SAMTOOLS_INDEX(PARALLAX_ALIGN.out.bam)
}
