#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────────────────
params.reference   = '/Users/tom.conway/data/hg38/hg38_primary.fasta'
params.index       = "hg38_idx"                     // Pre-built parallax index directory (default: projectDir/hg38_idx)
params.outdir      = 'bench_results'
params.seed        = 42
params.num_reads   = 10000
params.mean_length = 15000
params.std_dev     = 3000.0
params.error_rate  = 0.0
params.threads     = 4
params.tolerance   = 50
params.parallax_config = null          // Optional TOML config for parallax
params.reads_cache = "${projectDir}/cached_reads"  // Persistent store for simulated reads

// Path to the parallax binary (release build)
params.parallax    = "${projectDir}/target/release/parallax"

// ─── Index helpers (mirrors parallax.nf) ──────────────────────────────────
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

process SETUP_PYTHON {
    input:
    val project_dir

    output:
    env 'PYTHON3', emit: python3

    script:
    """
    VENV="${project_dir}/.venv"
    # Resolve the correct Python interpreter, honouring .python-version via pyenv
    if command -v pyenv >/dev/null 2>&1 && [ -f "${project_dir}/.python-version" ]; then
        WANT_VER=\$(cat "${project_dir}/.python-version")
        PYTHON_BIN=\$(PYENV_VERSION=\$WANT_VER pyenv which python3 2>/dev/null || which python3)
    else
        PYTHON_BIN=\$(which python3)
    fi
    # Recreate the venv if it doesn't exist or was built with a different Python
    GOT_VER=\$("\$VENV/bin/python3" --version 2>/dev/null | cut -d' ' -f2 || echo none)
    EXP_VER=\$("\$PYTHON_BIN" --version | cut -d' ' -f2)
    if [ "\$GOT_VER" != "\$EXP_VER" ]; then
        rm -rf "\$VENV"
        "\$PYTHON_BIN" -m venv "\$VENV"
    fi
    PYTHON3="\$VENV/bin/python3"
    \$PYTHON3 -c "import pysam" 2>/dev/null || \$PYTHON3 -m pip install --quiet pysam
    """
}

process BUILD_INDEX {
    tag "index"
    cpus params.threads
    memory '30 GB'
    publishDir "${params.outdir}/${params.index}", mode: 'move'

    input:
    val project_dir

    output:
    path "index", type: 'dir', emit: index

    script:
    """
    ${params.parallax} index \\
        ${params.reference} \\
        -o index \\
        -p \\
        -t ${task.cpus}
    """
}

process SIMULATE_READS {
    tag "seed_${seed}"

    // Cache reads by the simulation parameters that affect output.
    // As long as these are unchanged, the process is skipped on subsequent runs.
    storeDir "${params.reads_cache}/n${params.num_reads}_l${params.mean_length}_s${params.std_dev}_e${params.error_rate}_seed${seed}"

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
    path index_dir

    output:
    tuple val(seed), path("plx_s${seed}.sam"), emit: sam

    script:
    def config_flag = params.parallax_config ? "-c ${params.parallax_config}" : ''
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
    val python3

    output:
    tuple val(seed), path("compare_s${seed}.tsv"),   emit: tsv
    tuple val(seed), path("compare_s${seed}.log"),   emit: log

    script:
    """
    ${python3} ${compare_script} \\
        ${plx_sam} ${mm2_sam} \\
        --name-a parallax --name-b minimap2 \\
        -t ${params.tolerance} \\
        --fasta ${params.reference} \\
        -v \\
        > compare_s${seed}.tsv \\
        2> compare_s${seed}.log
    """
}

// ─── Workflow ──────────────────────────────────────────────────────────────
workflow {
    seed_ch = channel.of(params.seed)
    compare_script = file("${projectDir}/scripts/compare_simulated_alignments.py")

    // Resolve or build the parallax index
    def index_param = params.index ?: "${projectDir}/hg38_idx"
    if (indexExists(index_param)) {
        def idx_dir = resolveIndexDir(index_param)
        log.info "Using existing index at ${idx_dir}"
        index_ch = channel.fromPath(idx_dir, type: 'dir').first()
    } else {
        log.info "Building index (no existing index found at ${index_param})"
        BUILD_INDEX(projectDir)
        index_ch = BUILD_INDEX.out.index
    }

    SETUP_PYTHON(projectDir)
    SIMULATE_READS(seed_ch, projectDir)
    ALIGN_MINIMAP2(SIMULATE_READS.out.reads)
    ALIGN_PARALLAX(SIMULATE_READS.out.reads, index_ch)

    // Join parallax and minimap2 SAMs by seed
    compare_ch = ALIGN_PARALLAX.out.sam
        .join(ALIGN_MINIMAP2.out.sam)

    COMPARE(compare_ch, compare_script, SETUP_PYTHON.out.python3)
}
