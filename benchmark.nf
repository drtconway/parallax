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
params.vcf             = null          // Optional VCF of structural variants to apply
params.global_sampling = false         // With --vcf, use unbiased global sampling instead of variant-biased

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
    storeDir "${params.reads_cache}/n${params.num_reads}_l${params.mean_length}_s${params.std_dev}_e${params.error_rate}_seed${seed}${params.vcf ? '_vcf' + file(params.vcf).baseName : ''}${params.global_sampling ? '_global' : ''}"

    input:
    val seed
    val project_dir

    output:
    tuple val(seed), path("reads_s${seed}.fq.gz"), emit: reads

    script:
    def vcf_flag = params.vcf ? "--vcf ${file(params.vcf)}" : ''
    def global_flag = params.global_sampling ? '--global-sampling' : ''
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
        --primary-only \\
        ${vcf_flag} ${global_flag}

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
    tuple val(seed), path("parallax-stats.tsv"), emit: stats

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
        plx_s${seed}.sam
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
    tuple val(seed), path("compare_s${seed}.md"),    emit: md

    script:
    """
    ${python3} ${compare_script} \\
        ${plx_sam} ${mm2_sam} \\
        --name-a parallax --name-b minimap2 \\
        -t ${params.tolerance} \\
        --fasta ${params.reference} \\
        -v \\
        > compare_s${seed}.tsv \\
        2> compare_s${seed}.md
    """
}

process READS_TO_BED {
    tag "bed_s${seed}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(seed), path(fastq)

    output:
    tuple val(seed), path("truth_s${seed}.bed"), emit: bed

    script:
    """
    gunzip -c ${fastq} \
      | awk 'NR % 4 == 1 { split(\$1, a, ":"); split(a[2], segs, ","); for (i in segs) { n=split(segs[i], f, "_"); strand=f[n]; end=f[n-1]; start=f[n-2]; chrom=""; for(j=1;j<=n-3;j++){chrom=chrom (j>1?"_":"") f[j]}; print chrom "\t" start "\t" end "\t" substr(a[1],2) "\t0\t" strand } }' \
      | sort -k1,1 -k2,2n \
      > truth_s${seed}.bed
    """
}

process SAM_TO_BAM {
    tag "${name}_s${seed}"
    cpus 2
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(seed), path(sam), val(name)

    output:
    tuple val(seed), path("${name}_s${seed}.bam"), path("${name}_s${seed}.bam.bai"), val(name), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${name}_s${seed}.bam ${sam}
    samtools index ${name}_s${seed}.bam
    """
}

process IGV_SESSION {
    tag "igv_s${seed}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(seed), path(plx_bam), path(plx_bai), path(mm2_bam), path(mm2_bai), path(truth_bed)

    output:
    tuple val(seed), path("igv_session_s${seed}.xml"), emit: session

    script:
    def outdir = file(params.outdir).toAbsolutePath()
    def ref    = file(params.reference).toAbsolutePath()
    def plx_path = "${outdir}/${plx_bam}"
    def mm2_path = "${outdir}/${mm2_bam}"
    def bed_path = "${outdir}/${truth_bed}"
    """
    cat > igv_session_s${seed}.xml <<'EOF'
<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<Session genome="${ref}" version="8">
    <Resources>
        <Resource path="${plx_path}" type="bam"/>
        <Resource path="${mm2_path}" type="bam"/>
        <Resource path="${bed_path}" type="bed"/>
    </Resources>
    <Panel name="PanelPlx" width="1775">
        <Track attributeKey="${plx_bam} Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="10" id="${plx_path}_coverage" name="${plx_bam} Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track attributeKey="${plx_bam} Junctions" autoScale="false" clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" groupByStrand="false" height="60" id="${plx_path}_junctions" maxdepth="50" name="${plx_bam} Junctions" visible="false"/>
        <Track attributeKey="${plx_bam}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="EXPANDED" experimentType="THIRD_GEN" fontSize="10" id="${plx_path}" name="${plx_bam}" visible="true">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel name="PanelMm2" width="1775">
        <Track attributeKey="${mm2_bam} Coverage" autoScale="true" clazz="org.broad.igv.sam.CoverageTrack" fontSize="10" id="${mm2_path}_coverage" name="${mm2_bam} Coverage" snpThreshold="0.2" visible="true">
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="10.0" minimum="0.0" type="LINEAR"/>
        </Track>
        <Track attributeKey="${mm2_bam} Junctions" autoScale="false" clazz="org.broad.igv.sam.SpliceJunctionTrack" fontSize="10" groupByStrand="false" height="60" id="${mm2_path}_junctions" maxdepth="50" name="${mm2_bam} Junctions" visible="false"/>
        <Track attributeKey="${mm2_bam}" clazz="org.broad.igv.sam.AlignmentTrack" color="185,185,185" displayMode="EXPANDED" experimentType="THIRD_GEN" fontSize="10" id="${mm2_path}" name="${mm2_bam}" visible="true">
            <RenderOptions/>
        </Track>
    </Panel>
    <Panel height="70" name="FeaturePanel" width="1775">
        <Track attributeKey="${truth_bed}" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" displayMode="COLLAPSED" fontSize="10" id="${bed_path}" name="Truth regions" visible="true"/>
        <Track attributeKey="Reference sequence" clazz="org.broad.igv.track.SequenceTrack" fontSize="10" id="Reference sequence" name="Reference sequence" sequenceTranslationStrandValue="+" shouldShowTranslation="false" visible="true"/>
    </Panel>
    <PanelLayout dividerFractions="0.006880733944954129,0.533256880733945,0.9139908256880734"/>
    <HiddenAttributes>
        <Attribute name="DATA FILE"/>
        <Attribute name="DATA TYPE"/>
        <Attribute name="NAME"/>
    </HiddenAttributes>
</Session>
EOF
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

    // Extract truth regions from simulated read names
    READS_TO_BED(SIMULATE_READS.out.reads)

    // Convert SAMs to sorted BAMs for IGV
    plx_bam_ch = ALIGN_PARALLAX.out.sam.map { seed, sam -> [seed, sam, 'plx'] }
    mm2_bam_ch = ALIGN_MINIMAP2.out.sam.map { seed, sam -> [seed, sam, 'mm2'] }
    SAM_TO_BAM(plx_bam_ch.mix(mm2_bam_ch))

    // Join BAMs and truth BED by seed, then build IGV session
    igv_ch = SAM_TO_BAM.out.bam
        .branch {
            plx: it[3] == 'plx'
            mm2: it[3] == 'mm2'
        }
    igv_input = igv_ch.plx
        .map { seed, bam, bai, _name -> [seed, bam, bai] }
        .join(igv_ch.mm2.map { seed, bam, bai, _name -> [seed, bam, bai] })
        .join(READS_TO_BED.out.bed)
    // igv_input: [seed, plx_bam, plx_bai, mm2_bam, mm2_bai, truth_bed]
    IGV_SESSION(igv_input)
}
