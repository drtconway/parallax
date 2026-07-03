#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────────────────
params.reference        = null
params.fastq            = null        // Path to fastq file(s), glob ok
params.samplesheet      = null        // CSV: sample,file  (alternative to --fastq)
params.index            = null        // Pre-built parallax index directory
params.outdir           = 'curation'
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
    publishDir "${params.outdir}/bam", mode: 'copy'

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
    publishDir "${params.outdir}/bam", mode: 'move'

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

// Extract reads that disagree between parallax and minimap2 for curation
process EXTRACT_CURATION_READS {
    tag "extract_${meta.id}"
    cpus 1
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(meta), path(compare_txt), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.curation.fq.gz"), emit: fastq

    script:
    """
    grep -v '^read_id' ${compare_txt} | awk '\$3 < 1.0 { print \$1 }' > disagree.txt
    python3 ${projectDir}/scripts/extract-reads.py disagree.txt < <(gzip -dc ${fastq}) \
        | gzip > ${meta.id}.curation.fq.gz
    """
}

// Align curation-subset reads with parallax, emitting both alignments and seeds
process ALIGN_PARALLAX_CURATION {
    tag "plx_curation_${meta.id}"
    cpus params.threads
    memory '30 GB'
    publishDir "${params.outdir}/bam", mode: 'move'

    input:
    tuple val(meta), path(reads)
    path reference
    path index_dir

    output:
    tuple val(meta), path("${meta.id}.curation.plx.sorted.bam"),
                     path("${meta.id}.curation.plx.sorted.bam.bai"),  emit: bam
    tuple val(meta), path("${meta.id}.curation.seeds.sorted.bam"),
                     path("${meta.id}.curation.seeds.sorted.bam.bai"), emit: seeds

    script:
    def config_flag = params.parallax_config ? "-c ${params.parallax_config}" : ''
    """
    # Write a TOML config that enables seed dumping to a known path
    cat > curation.toml <<'TOML'
[seeding]
debug_seeds_sam = "seeds.sam"
TOML

    # Merge with any user-supplied config by passing both -c flags (last wins
    # for duplicate keys; seeds path is only in curation.toml)
    ${params.parallax} align \\
        ${reference} \\
        ${reads} \\
        -x ${index_dir} \\
        -p \\
        ${config_flag} \\
        -c curation.toml \\
        -t ${task.cpus} \\
        | samtools sort -@ ${task.cpus} -o ${meta.id}.curation.plx.sorted.bam -
    samtools index ${meta.id}.curation.plx.sorted.bam

    # seeds.sam may not exist if all reads were unmapped; create an empty one
    [ -f seeds.sam ] || samtools view -H ${meta.id}.curation.plx.sorted.bam > seeds.sam
    samtools sort -@ ${task.cpus} -o ${meta.id}.curation.seeds.sorted.bam seeds.sam
    samtools index ${meta.id}.curation.seeds.sorted.bam
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

    // Extract disagreeing reads and run parallax on them with seed dumping
    curation_input = COMPARE_ALIGNMENTS.out.txt.join(reads_ch)
    EXTRACT_CURATION_READS(curation_input)
    ALIGN_PARALLAX_CURATION(EXTRACT_CURATION_READS.out.fastq, reference_ch, index_ch)
}
