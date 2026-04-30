#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────────────────
params.reference        = null
params.fastq            = null        // Path to fastq (or ubam) file(s), glob ok
params.samplesheet      = null        // CSV: sample,file  (alternative to --fastq)
params.index            = null        // Pre-built parallax index directory
params.outdir           = 'compare_results'
params.threads          = 4
params.parallax_config  = null
params.parallax         = "${projectDir}/target/release/parallax"

// ─── Index helpers ────────────────────────────────────────────────────────
def indexExists(index_path) {
    if (!index_path) return false
    def idx_dir = file(index_path)
    if (idx_dir.isDirectory() && file("${index_path}/chrom_info.json").exists()) return true
    if (idx_dir.isDirectory() && file("${index_path}/index/chrom_info.json").exists()) return true
    return false
}

def resolveIndexDir(index_path) {
    if (file("${index_path}/chrom_info.json").exists()) return index_path
    if (file("${index_path}/index/chrom_info.json").exists()) return "${index_path}/index"
    return index_path
}

// ─── Processes ────────────────────────────────────────────────────────────

process BUILD_INDEX {
    tag "index"
    cpus params.threads
    memory '30 GB'

    input:
    path reference

    output:
    path "index", type: 'dir', emit: index

    script:
    """
    ${params.parallax} index \\
        ${reference} \\
        -o index \\
        -p \\
        -t ${task.cpus}
    """
}

process ALIGN_PARALLAX {
    tag "plx_${meta.id}"
    cpus params.threads
    memory '30 GB'

    input:
    tuple val(meta), path(reads)
    path reference
    path index_dir

    output:
    tuple val(meta), path("${meta.id}.plx.bam"), emit: bam

    script:
    def config_flag = params.parallax_config ? "-c ${params.parallax_config}" : ''
    """
    ${params.parallax} align \\
        ${reference} \\
        ${reads} \\
        -x ${index_dir} \\
        -p \\
        ${config_flag} \\
        -t ${task.cpus} \\
        | samtools sort -n -@ ${task.cpus} -o ${meta.id}.plx.bam -
    """
}

process ALIGN_MINIMAP2 {
    tag "mm2_${meta.id}"
    cpus params.threads
    memory '30 GB'

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}.mm2.bam"), emit: bam

    script:
    // Detect ubam vs fastq by extension; minimap2 handles both natively
    def fmt = (reads.name ==~ /.*\.(bam|ubam)/) ? '-a -x map-hifi --secondary=yes' : '-a -x map-hifi --secondary=yes'
    """
    minimap2 ${fmt} -t ${task.cpus} \\
        ${reference} ${reads} \\
        | samtools sort -n -@ ${task.cpus} -o ${meta.id}.mm2.bam -
    """
}

process COMPARE_ALIGNMENTS {
    tag "cmp_${meta.id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(plx_bam), path(mm2_bam)
    path compare_script

    output:
    tuple val(meta), path("${meta.id}.compare.tsv"), emit: tsv

    script:
    """
    python3 ${compare_script} ${plx_bam} ${mm2_bam} -o ${meta.id}.compare.tsv
    """
}

// ─── Input helpers ────────────────────────────────────────────────────────

def parseSamplesheet(csv) {
    return channel
        .fromPath(csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample ]
            def f = file(row.file, checkIfExists: true)
            return [ meta, f ]
        }
}

// ─── Workflow ─────────────────────────────────────────────────────────────
workflow {
    if (!params.reference) error "Please specify --reference"
    if (!params.fastq && !params.samplesheet) error "Please specify --fastq or --samplesheet"
    if (params.fastq && params.samplesheet) error "Specify --fastq or --samplesheet, not both"

    reference_ch = channel.fromPath(params.reference, checkIfExists: true).first()

    // Build reads channel
    if (params.samplesheet) {
        reads_ch = parseSamplesheet(params.samplesheet)
    } else {
        reads_ch = channel
            .fromPath(params.fastq, checkIfExists: true)
            .map { f -> [ [id: f.simpleName], f ] }
    }

    // Resolve or build the parallax index
    def index_param = params.index
    if (indexExists(index_param)) {
        def idx_dir = resolveIndexDir(index_param)
        log.info "Using existing index at ${idx_dir}"
        index_ch = channel.fromPath(idx_dir, type: 'dir').first()
    } else {
        log.info "Building parallax index"
        BUILD_INDEX(reference_ch)
        index_ch = BUILD_INDEX.out.index
    }

    // Align with both tools (name-sorted output)
    ALIGN_PARALLAX(reads_ch, reference_ch, index_ch)
    ALIGN_MINIMAP2(reads_ch, reference_ch)

    // Join by sample id and compare read by read
    compare_ch = ALIGN_PARALLAX.out.bam
        .join(ALIGN_MINIMAP2.out.bam)

    compare_script = file("${projectDir}/scripts/compare_alignments_per_read.py")
    COMPARE_ALIGNMENTS(compare_ch, compare_script)
}
