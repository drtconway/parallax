#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────────────────
params.reference   = '/Users/tom.conway/data/hg38/hg38_primary.fasta'
params.index       = null                     // Pre-built parallax index directory (default: projectDir/hg38_idx)
params.outdir      = 'bench_results'
params.seed        = 42
params.num_reads   = 10000
params.mean_length = 15000
params.std_dev     = 3000.0
params.error_rate  = 0.0
params.threads     = 4
params.tolerance   = 50
params.parallax_config = null          // Optional TOML config for parallax

// Path to the parallax binary (release build)
params.parallax    = "${projectDir}/target/release/parallax"

process SIMULATE_READS {
    tag "seed_${seed}"
    publishDir "${params.outdir}", mode: 'copy', pattern: '*.fq.gz'

    input:
    val seed
    val project_dir

    output:
    tuple val(seed), path("reads_s${seed}.fq.gz"), emit: reads

    script:
    """
    cargo run --manifest-path ${project_dir}/Cargo.toml \\
        --release --example simulate_reads -- \\
        --reference ${params.reference} \\
        --output reads_s${seed}.fq \\
        --num-reads ${params.num_reads} \\
        --mean-length ${params.mean_length} \\
        --std-dev ${params.std_dev} \\
        --error-rate ${params.error_rate} \\
        --seed ${seed} \\
        --primary-only

    gzip reads_s${seed}.fq
    """
}

process ALIGN_MINIMAP2 {
    tag "mm2_s${seed}"
    cpus params.threads
    memory '30 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(seed), path(fastq)

    output:
    tuple val(seed), path("mm2_s${seed}.sam"), emit: sam

    script:
    """
    minimap2 -a -x map-hifi -t ${task.cpus} \\
        ${params.reference} ${fastq} \\
        > mm2_s${seed}.sam
    """
}

process ALIGN_PARALLAX {
    tag "plx_s${seed}"
    cpus params.threads
    memory '30 GB'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(seed), path(fastq)
    val project_dir

    output:
    tuple val(seed), path("plx_s${seed}.sam"), emit: sam

    script:
    def config_flag = params.parallax_config ? "-c ${params.parallax_config}" : ''
    def index_dir = params.index ?: "${project_dir}/hg38_idx"
    """
    ${params.parallax} align \\
        ${params.reference} \\
        ${fastq} \\
        -x ${index_dir} \\
        -p \\
        ${config_flag} \\
        -t ${task.cpus} \\
        > plx_s${seed}.sam
    """
}

process COMPARE {
    tag "cmp_s${seed}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(seed), path(plx_sam), path(mm2_sam)
    path compare_script

    output:
    tuple val(seed), path("compare_s${seed}.tsv"),   emit: tsv
    tuple val(seed), path("compare_s${seed}.log"),   emit: log

    script:
    """
    python3 ${compare_script} \\
        ${plx_sam} ${mm2_sam} \\
        --name-a parallax --name-b minimap2 \\
        -t ${params.tolerance} \\
        -v \\
        > compare_s${seed}.tsv \\
        2> compare_s${seed}.log
    """
}

// ─── Workflow ──────────────────────────────────────────────────────────────
workflow {
    seed_ch = channel.of(params.seed)
    compare_script = file("${projectDir}/scripts/compare_simulated_alignments.py")

    SIMULATE_READS(seed_ch, projectDir)
    ALIGN_MINIMAP2(SIMULATE_READS.out.reads)
    ALIGN_PARALLAX(SIMULATE_READS.out.reads, projectDir)

    // Join parallax and minimap2 SAMs by seed
    compare_ch = ALIGN_PARALLAX.out.sam
        .join(ALIGN_MINIMAP2.out.sam)

    COMPARE(compare_ch, compare_script)
}
