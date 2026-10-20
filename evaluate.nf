#!/usr/bin/env nextflow

// Truth-free comparison of parallax and minimap2 alignments.
//
// Each read set is aligned by each aligner to each reference (hg38, CHM13),
// along with its reverse complement and fixed-length tiles of it.  Metrics:
//   - per-alignment / per-read statistics (scripts/eval/alignment_stats.py)
//   - invariance under reverse complement and fragmentation (scripts/eval/invariance.py)
//   - pileup allele fractions, depth and candidate false hets (parallax eval-pileup)
// all gathered into one report (scripts/eval/report.py).
//
//   nextflow run evaluate.nf -profile slurm --samplesheet reads.csv \
//       --hg38 hg38.fa --chm13 chm13v2.fa [--hg38_index dir] [--chm13_index dir]
//
// Samplesheet columns: sample,fastq,platform[,parallax_config]
// platform is HIFI, ONT or ONT_HQ (selects the minimap2 preset).

nextflow.enable.dsl = 2

// ─── Parameters ────────────────────────────────────────────────────────────
params.samplesheet      = null
params.hg38             = null        // uncompressed FASTA (.fai made if missing)
params.hg38_index       = null        // prebuilt parallax index (built if absent)
params.strat_hg38       = null        // TSV of name<TAB>bed (absolute paths)
params.chm13            = null
params.chm13_index      = null
params.strat_chm13      = null
params.index_store      = "${projectDir}/eval_indexes"  // where built indexes are kept
params.subsample        = 0           // use only the first N reads of each set (0 = all)
params.chunk_size       = 50000       // reads per alignment job
params.fragment_lengths = '1000,3000,10000'  // '' to skip fragment invariance
params.rc               = true        // reverse-complement invariance
params.parallax         = "${projectDir}/target/release/parallax"
params.parallax_config  = null
params.parallax_args    = ''
params.minimap2_args    = '--eqx -Y --secondary=yes'
params.pileup_chroms    = '^(chr)?([0-9]+|X|Y|M)$'  // chromosomes to pile up
params.min_mapq         = 20
params.publish_bams     = true

// minimap2 preset for each samplesheet platform
def mm2Presets() {
    [HIFI: 'map-hifi', ONT: 'map-ont', ONT_HQ: 'lr:hq']
}

// ─── Index helpers (as in compare.nf) ─────────────────────────────────────
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

def stem(meta) {
    "${meta.id}.${meta.ref}.${meta.aligner}"
}

def chunkKey(meta) {
    [ meta.id, meta.ref, meta.aligner, meta.chunk ]
}

// ─── Processes ────────────────────────────────────────────────────────────

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
    \$PYTHON3 -c "import pysam, pandas, matplotlib" 2>/dev/null \\
        || \$PYTHON3 -m pip install --quiet pysam pandas matplotlib
    """

    stub:
    """
    PYTHON3=python3
    """
}

process FAIDX {
    tag "${ref_id}"
    label 'small'

    input:
    tuple val(ref_id), path(fasta)

    output:
    tuple val(ref_id), path(fasta), path("${fasta}.fai"), emit: ref

    script:
    """
    samtools faidx ${fasta}
    """

    stub:
    """
    printf 'chr1\\t1000\\t6\\t60\\t61\\n' > ${fasta}.fai
    """
}

process BUILD_INDEX {
    tag "${ref_id}"
    label 'index'
    storeDir "${params.index_store}/${ref_id}"

    input:
    tuple val(ref_id), path(fasta)

    output:
    tuple val(ref_id), path("index", type: 'dir'), emit: index

    script:
    """
    ${params.parallax} index ${fasta} -o index -p -t ${task.cpus}
    """

    stub:
    """
    mkdir index
    """
}

process MM2_INDEX {
    tag "${ref_id}_${preset}"
    label 'index'
    storeDir "${params.index_store}/${ref_id}"

    input:
    tuple val(ref_id), val(preset), path(fasta)

    output:
    tuple val(ref_id), val(preset), path("${ref_id}.${preset.replace(':', '_')}.mmi"), emit: mmi

    script:
    """
    minimap2 -x ${preset} -t ${task.cpus} -d ${ref_id}.${preset.replace(':', '_')}.mmi ${fasta}
    """

    stub:
    """
    touch ${ref_id}.${preset.replace(':', '_')}.mmi
    """
}

// Optionally subsample, then split into chunks of params.chunk_size reads.
process PREP_READS {
    tag "${meta.id}"
    label 'small'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.chunk*.fq.gz"), emit: chunks

    script:
    """
    gzip -cdf ${fastq} | awk -v n=${params.chunk_size} -v max=${params.subsample} -v p="${meta.id}" '
        BEGIN { cur = -1 }
        NR % 4 == 1 {
            r++
            if (max > 0 && r > max) exit
            c = int((r - 1) / n)
            if (c != cur) {
                if (out != "") close(out)
                cur = c
                out = sprintf("gzip -1 > %s.chunk%04d.fq.gz", p, c)
            }
        }
        { print | out }'
    """

    stub:
    """
    touch ${meta.id}.chunk0000.fq.gz ${meta.id}.chunk0001.fq.gz
    """
}

process MAKE_VARIANTS {
    tag "${meta.id}_${meta.chunk}"
    label 'small'

    input:
    tuple val(meta), path(reads)
    val python3
    path scripts, stageAs: 'eval_scripts/*'

    output:
    tuple val(meta), path("variants/*.fq.gz"), emit: variants

    script:
    def rc_flag = params.rc ? "--rc variants/rc.fq.gz" : ''
    def frag_flags = fragmentLengths().collect { "--frag ${it}:variants/frag${it}.fq.gz" }.join(' ')
    """
    mkdir variants
    ${python3} eval_scripts/make_variants.py ${reads} ${rc_flag} ${frag_flags}
    """

    stub:
    """
    mkdir variants
    ${params.rc ? 'touch variants/rc.fq.gz' : ''}
    ${fragmentLengths().collect { "touch variants/frag${it}.fq.gz" }.join('\n')}
    """
}

process ALIGN_PARALLAX {
    tag "${stem(meta)}.${meta.chunk}.${meta.variant}"
    label 'align'

    input:
    tuple val(meta), path(reads), path(fasta), path(fai), path(index_dir)

    output:
    tuple val(meta), path("${stem(meta)}.${meta.chunk}.${meta.variant}.bam"), emit: bam

    script:
    def cfg = meta.parallax_config ?: params.parallax_config
    def config_flag = cfg ? "-c ${cfg}" : ''
    """
    ${params.parallax} align2 \\
        ${fasta} \\
        ${reads} \\
        -x ${index_dir} \\
        -p \\
        ${config_flag} \\
        -t ${task.cpus} \\
        ${params.parallax_args} \\
        | samtools view -b -o ${stem(meta)}.${meta.chunk}.${meta.variant}.bam -
    """

    stub:
    """
    touch ${stem(meta)}.${meta.chunk}.${meta.variant}.bam
    """
}

process ALIGN_MINIMAP2 {
    tag "${stem(meta)}.${meta.chunk}.${meta.variant}"
    label 'align'

    input:
    tuple val(meta), path(reads), path(mmi)

    output:
    tuple val(meta), path("${stem(meta)}.${meta.chunk}.${meta.variant}.bam"), emit: bam

    script:
    """
    minimap2 -a -x ${meta.preset} ${params.minimap2_args} -t ${task.cpus} ${mmi} ${reads} \\
        | samtools view -b -o ${stem(meta)}.${meta.chunk}.${meta.variant}.bam -
    """

    stub:
    """
    touch ${stem(meta)}.${meta.chunk}.${meta.variant}.bam
    """
}

process ALIGNMENT_STATS {
    tag "${stem(meta)}.${meta.chunk}"
    label 'metrics'
    publishDir "${params.outdir}/metrics/alignments", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(fasta), path(fai), val(strata)
    val python3
    path scripts, stageAs: 'eval_scripts/*'

    output:
    path "${stem(meta)}.${meta.chunk}.*", emit: files

    script:
    def p = "${stem(meta)}.${meta.chunk}"
    def strata_flag = strata ? "--strata ${strata}" : ''
    """
    samtools sort -n -@ ${task.cpus} -o ns.bam ${bam}
    ${python3} eval_scripts/alignment_stats.py ns.bam \\
        --prefix ${p} --reference ${fasta} ${strata_flag}
    rm ns.bam
    """

    stub:
    def p = "${stem(meta)}.${meta.chunk}"
    """
    touch ${p}.alignments.tsv.gz ${p}.reads.tsv.gz ${p}.indels.tsv ${p}.errpos.tsv ${p}.outliers.tsv.gz
    """
}

process INVARIANCE {
    tag "${stem(meta)}.${meta.chunk}.${meta.variant}"
    label 'metrics'
    publishDir "${params.outdir}/metrics/invariance", mode: 'copy'

    input:
    tuple val(meta), path(orig_bam), path(variant_bam), val(strata)
    val python3
    path scripts, stageAs: 'eval_scripts/*'

    output:
    path "${stem(meta)}.${meta.chunk}.${meta.variant}.invariance.tsv.gz", emit: files

    script:
    def kind = meta.variant == 'rc' ? 'rc' : 'frag'
    def strata_flag = strata ? "--strata ${strata}" : ''
    """
    samtools sort -n -@ ${task.cpus} -o orig.ns.bam ${orig_bam}
    samtools sort -n -@ ${task.cpus} -o variant.ns.bam ${variant_bam}
    ${python3} eval_scripts/invariance.py orig.ns.bam variant.ns.bam \\
        --kind ${kind} --min-mapq ${params.min_mapq} ${strata_flag} \\
        --out ${stem(meta)}.${meta.chunk}.${meta.variant}.invariance.tsv.gz
    rm orig.ns.bam variant.ns.bam
    """

    stub:
    """
    touch ${stem(meta)}.${meta.chunk}.${meta.variant}.invariance.tsv.gz
    """
}

// Merge the original-read chunks of one sample/ref/aligner; coordinate sort.
process MERGE_SORT {
    tag "${stem(meta)}"
    label 'sort'
    publishDir "${params.outdir}/bam", mode: 'copy', enabled: params.publish_bams

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${stem(meta)}.bam"), path("${stem(meta)}.bam.bai"), emit: bam

    script:
    """
    samtools cat -o - ${bams} \\
        | samtools sort -@ ${task.cpus} -m 2G -o ${stem(meta)}.bam -
    samtools index ${stem(meta)}.bam
    """

    stub:
    """
    touch ${stem(meta)}.bam ${stem(meta)}.bam.bai
    """
}

process PILEUP {
    tag "${stem(meta)}.${chrom}"
    label 'pileup'
    publishDir "${params.outdir}/metrics/pileup", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), val(chrom), path(fasta), path(fai), val(strata)

    output:
    path "${stem(meta)}.${chrom}.*", emit: files

    script:
    def strata_flag = strata ? "--strata ${strata}" : ''
    """
    ${params.parallax} eval-pileup ${bam} -r ${fasta} --chrom ${chrom} \\
        --min-mapq ${params.min_mapq} ${strata_flag} -o ${stem(meta)}.${chrom}
    """

    stub:
    def p = "${stem(meta)}.${chrom}"
    """
    touch ${p}.hist.tsv ${p}.totals.tsv ${p}.windows.tsv ${p}.sites.bed
    """
}

process REPORT {
    label 'metrics'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path inputs, stageAs: 'inputs/*'
    val python3
    path scripts, stageAs: 'eval_scripts/*'

    output:
    path "report", emit: report

    script:
    """
    ${python3} eval_scripts/report.py inputs --out report
    """

    stub:
    """
    mkdir report
    touch report/summary.tsv report/report.html
    """
}

// ─── Input helpers ────────────────────────────────────────────────────────

def fragmentLengths() {
    def spec = params.fragment_lengths?.toString() ?: ''
    return spec.tokenize(',').collect { it.trim() as int }
}

def parseSamplesheet(csv) {
    return channel
        .fromPath(csv, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || row.sample.contains('.')) error "Sample name must be non-empty with no '.': '${row.sample}'"
            def platform = (row.platform ?: '').toUpperCase()
            if (!mm2Presets().containsKey(platform)) error "Unknown platform '${row.platform}' for ${row.sample}; use one of ${mm2Presets().keySet()}"
            def meta = [
                id: row.sample,
                platform: platform,
                preset: mm2Presets()[platform],
                parallax_config: row.parallax_config ?: null,
            ]
            [ meta, file(row.fastq, checkIfExists: true) ]
        }
}

// ─── Workflow ─────────────────────────────────────────────────────────────
workflow {
    if (!params.samplesheet) error "Please specify --samplesheet"
    def refs = []
    if (params.hg38)  refs << [ 'hg38',  params.hg38,  params.hg38_index,  params.strat_hg38 ]
    if (params.chm13) refs << [ 'chm13', params.chm13, params.chm13_index, params.strat_chm13 ]
    if (!refs) error "Please specify --hg38 and/or --chm13"

    SETUP_PYTHON(projectDir)
    python3 = SETUP_PYTHON.out.python3
    // Staged as inputs so that editing a script invalidates cached results
    eval_scripts = channel.fromPath("${projectDir}/scripts/eval/*.py").collect()

    // References: [ref_id, fasta, fai]
    ref_in = channel.fromList(refs).map { id, fa, _idx, _strat -> [ id, file(fa, checkIfExists: true) ] }
        .branch { id, fa ->
            indexed: file("${fa}.fai").exists()
            unindexed: true
        }
    FAIDX(ref_in.unindexed)
    ref_ch = ref_in.indexed.map { id, fa -> [ id, fa, file("${fa}.fai") ] }.mix(FAIDX.out.ref)

    // Strata TSV path (or '') per reference: [ref_id, strata]
    strata_ch = channel.fromList(refs.collect { id, _fa, _idx, strat ->
        [ id, strat ? file(strat, checkIfExists: true).toString() : '' ] })

    // Parallax indexes: [ref_id, index_dir]
    idx_in = channel.fromList(refs).map { id, fa, idx, _strat -> [ id, file(fa), idx ] }
        .branch { _id, _fa, idx ->
            prebuilt: indexExists(idx)
            build: true
        }
    BUILD_INDEX(idx_in.build.map { id, fa, _idx -> [ id, fa ] })
    plx_index_ch = idx_in.prebuilt
        .map { id, _fa, idx -> [ id, file(resolveIndexDir(idx), type: 'dir') ] }
        .mix(BUILD_INDEX.out.index)

    // Reads, chunked: [meta + chunk + variant, fastq]
    samples = parseSamplesheet(params.samplesheet)
    PREP_READS(samples)
    chunks = PREP_READS.out.chunks
        .transpose()
        .map { meta, fq -> [ meta + [ chunk: fq.name.tokenize('.')[1], variant: 'orig' ], fq ] }

    if (params.rc || fragmentLengths()) {
        MAKE_VARIANTS(chunks, python3, eval_scripts)
        variant_reads = MAKE_VARIANTS.out.variants
            .transpose()
            .map { meta, fq -> [ meta + [ variant: fq.name.tokenize('.')[0] ], fq ] }
    } else {
        variant_reads = channel.empty()
    }
    all_reads = chunks.mix(variant_reads)

    // minimap2 indexes, one per reference x preset in use
    presets = samples.map { meta, _fq -> meta.preset }.unique()
    MM2_INDEX(ref_ch.map { id, fa, _fai -> [ id, fa ] }.combine(presets).map { id, fa, preset -> [ id, preset, fa ] })

    // ── Alignment ──
    plx_in = all_reads
        .combine(ref_ch)
        .map { meta, fq, ref_id, fa, fai -> [ ref_id, meta + [ ref: ref_id, aligner: 'parallax' ], fq, fa, fai ] }
        .combine(plx_index_ch, by: 0)
        .map { _ref_id, meta, fq, fa, fai, idx -> [ meta, fq, fa, fai, idx ] }
    ALIGN_PARALLAX(plx_in)

    mm2_in = all_reads
        .combine(ref_ch.map { id, _fa, _fai -> id })
        .map { meta, fq, ref_id -> [ [ ref_id, meta.preset ], meta + [ ref: ref_id, aligner: 'minimap2' ], fq ] }
        .combine(MM2_INDEX.out.mmi.map { id, preset, mmi -> [ [ id, preset ], mmi ] }, by: 0)
        .map { _key, meta, fq, mmi -> [ meta, fq, mmi ] }
    ALIGN_MINIMAP2(mm2_in)

    aligned = ALIGN_PARALLAX.out.bam.mix(ALIGN_MINIMAP2.out.bam)
        .branch { meta, _bam ->
            orig: meta.variant == 'orig'
            variant: true
        }
    ref_files = ref_ch.join(strata_ch)   // [ref_id, fasta, fai, strata]

    // ── Per-alignment statistics (original reads) ──
    ALIGNMENT_STATS(
        aligned.orig
            .map { meta, bam -> [ meta.ref, meta, bam ] }
            .combine(ref_files, by: 0)
            .map { _ref_id, meta, bam, fa, fai, strata -> [ meta, bam, fa, fai, strata ] },
        python3, eval_scripts)

    // ── Invariance (variants vs originals, same aligner/ref/chunk) ──
    INVARIANCE(
        aligned.orig.map { meta, bam -> [ chunkKey(meta), bam ] }
            .combine(aligned.variant.map { meta, bam -> [ chunkKey(meta), meta, bam ] }, by: 0)
            .map { _key, orig_bam, meta, var_bam -> [ meta.ref, meta, orig_bam, var_bam ] }
            .combine(strata_ch, by: 0)
            .map { _ref_id, meta, orig_bam, var_bam, strata -> [ meta, orig_bam, var_bam, strata ] },
        python3, eval_scripts)

    // ── Pileup over merged, sorted original-read alignments ──
    MERGE_SORT(
        aligned.orig
            .map { meta, bam -> [ [ meta.id, meta.ref, meta.aligner ], meta.subMap('id', 'platform', 'preset', 'ref', 'aligner'), bam ] }
            .groupTuple(by: 0)
            .map { _key, metas, bams -> [ metas[0], bams.sort { it.name } ] })

    chroms = ref_ch.flatMap { id, _fa, fai ->
        fai.readLines()
            .collect { it.tokenize('\t')[0] }
            .findAll { it ==~ params.pileup_chroms }
            .collect { [ id, it ] }
    }
    PILEUP(
        MERGE_SORT.out.bam
            .map { meta, bam, bai -> [ meta.ref, meta, bam, bai ] }
            .combine(chroms, by: 0)
            .combine(ref_files, by: 0)
            .map { _ref_id, meta, bam, bai, chrom, fa, fai, strata -> [ meta, bam, bai, chrom, fa, fai, strata ] })

    // ── Report ──
    REPORT(
        ALIGNMENT_STATS.out.files.flatten()
            .mix(INVARIANCE.out.files.flatten(), PILEUP.out.files.flatten())
            .collect(),
        python3, eval_scripts)
}
