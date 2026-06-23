#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────────────────
params.reference        = null
params.fastq            = null        // Path to fastq file(s), glob ok
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
    storeDir "${params.outdir}/index"

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

// parallax emits SAM; convert to unsorted BAM directly
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
        | samtools view -b -o ${meta.id}.plx.bam -
    """
}

// minimap2 emits SAM; convert to unsorted BAM
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
    """
    minimap2 -a -x map-hifi --secondary=yes -t ${task.cpus} \\
        ${reference} ${reads} \\
        | samtools view -b -o ${meta.id}.mm2.bam -
    """
}

// Name-sort
process SORT_BY_NAME {
    tag "nsort_${meta.id}_${bam.baseName}"
    cpus params.threads
    memory '16 GB'
    storeDir "${params.outdir}/bam"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.baseName}.nsorted.bam"), emit: bam

    script:
    """
    samtools sort -n -@ ${task.cpus} -o ${bam.baseName}.nsorted.bam ${bam}
    """
}

// Coordinate-sort and index
process SORT_AND_INDEX {
    tag "csort_${meta.id}_${bam.baseName}"
    cpus params.threads
    memory '16 GB'
    storeDir "${params.outdir}/bam"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.baseName}.sorted.bam"),
                     path("${bam.baseName}.sorted.bam.bai"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${bam.baseName}.sorted.bam ${bam}
    samtools index ${bam.baseName}.sorted.bam
    """
}

process COMPARE_ALIGNMENTS {
    tag "cmp_${meta.id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(meta), path(plx_bam), path(mm2_bam)

    output:
    tuple val(meta), path("${meta.id}.compare.txt"), emit: txt

    script:
    """
    python3 ${projectDir}/scripts/compare_alignments_2.py ${plx_bam} ${mm2_bam} > ${meta.id}.compare.txt
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

    if (params.samplesheet) {
        reads_ch = parseSamplesheet(params.samplesheet)
    } else {
        reads_ch = channel
            .fromPath(params.fastq, checkIfExists: true)
            .map { f -> [ [id: f.simpleName], f ] }
    }

    if (indexExists(params.index)) {
        def idx_dir = resolveIndexDir(params.index)
        log.info "Using existing index at ${idx_dir}"
        index_ch = channel.fromPath(idx_dir, type: 'dir').first()
    } else {
        log.info "Building parallax index"
        BUILD_INDEX(reference_ch)
        index_ch = BUILD_INDEX.out.index
    }

    // Align -> unsorted BAM
    ALIGN_PARALLAX(reads_ch, reference_ch, index_ch)
    ALIGN_MINIMAP2(reads_ch, reference_ch)

    // Name-sort for comparison
    SORT_BY_NAME(ALIGN_PARALLAX.out.bam.mix(ALIGN_MINIMAP2.out.bam))

    // Split name-sorted channel back into plx / mm2 by filename suffix
    plx_nsorted = SORT_BY_NAME.out.bam.filter { meta, bam -> bam.name.contains('.plx.') }
    mm2_nsorted  = SORT_BY_NAME.out.bam.filter { meta, bam -> bam.name.contains('.mm2.') }

    COMPARE_ALIGNMENTS(plx_nsorted.join(mm2_nsorted))

    // Coordinate-sort + index for IGV
    SORT_AND_INDEX(ALIGN_PARALLAX.out.bam.mix(ALIGN_MINIMAP2.out.bam))
}
