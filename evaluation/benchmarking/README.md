# Chaining DP benchmark and gap-alignment study

Tools for measuring parallax's bits-scored chaining DP (`best_chain` in
`src/reads/segments.rs`) offline, on the inputs it sees in real alignment runs.

- `instrumentation.patch` adds timing counters to parallax and two optional
  dumps: every `best_chain` input (`PARALLAX_DUMP_DP`), and every gap alignment
  the DP makes, with its sequences, score and time (`PARALLAX_GAPLOG`). It
  applies to commit 3bceef4; later changes to `segments.rs`, `explanatory2.rs`
  or `reads.rs` may need it refreshed.
- `dpbench` (this crate) replays a dump through several implementations of
  `best_chain`, checks that they agree, and times them. It also runs the
  gap-alignment study.

## Producing a dump

Use a scratch clone: the patch is not meant for production builds.

```sh
git clone <parallax repo> parallax-instr && cd parallax-instr
git apply evaluation/benchmarking/instrumentation.patch
cargo build --release
PARALLAX_DUMP_DP=dp.bin PARALLAX_GAPLOG=gaplog.bin \
  target/release/parallax align2 -x <index> -t 1 --no-secondary -c bits.toml \
    <reference.fa> <reads.fq> out.sam
```

Run with `-t 1`: each record's read index comes from a global counter. At the
end the instrumented build prints a breakdown of time in the bits DP (`chain`,
`mapq`, `best_chain`, gap alignments), with trip-count histograms, to stderr.
1,000 HiFi reads give a dump of about 120 MB and a gap log of a few hundred MB.
The formats are defined by the patch and read by `load` (`main.rs`) and
`load_gaplog` (`gapstudy.rs`).

## Running

```sh
cargo run --release                                  # no dump: SIMD self-test
cargo run --release -- dp.bin [reps] [limit]         # DP variants
RUSTFLAGS="-C target-cpu=x86-64-v3" \
  cargo run --release --target-dir target-v3 -- dp.bin   # with AVX2 enabled
DPBENCH_FX=1 cargo run --release -- dp.bin           # Fx-hashed gap lookups
DPBENCH_GAPSPAN=1 cargo run --release -- dp.bin      # gaps spanning other seeds
DPBENCH_GAPSTUDY=gaplog.bin cargo run --release -- dp.bin   # gap-alignment study
```

The variants:

| Variant | Implementation |
| --- | --- |
| A | copy of `best_chain` (f64) |
| D | A, but predecessors the prefix-max bound rules out are only checked as gap candidates, not costed; identical results |
| B0, B1 | struct-of-arrays, i32 coordinates, f32 scores, blocks of 8 predecessors in branch-free plain Rust (for autovectorisation); B1 adds D's skip per block |
| C0, C1 | B's algorithm with explicit `fearless_simd` `f32x8`/`i32x8` vectors |

All skip gaps that could not win even if they aligned perfectly, as
`best_chain` does. B and C assume chromosomes shorter than 2³¹ bp (as
parallax's seeds already do) and use a polynomial `log2` (error at most
1.8×10⁻⁵ bits), so they can choose differently from A on near-ties.

The gap study reports gap sizes and alignment cost by size, scores cheap
proxies for the alignment (4-mer and 8-mer spectra, order-aware banded 4-mers,
a single-indel ungapped alignment) against the real scores and against
unrelated pairs, then replays the DP with each proxy used as a screen or as
a replacement, counting alignment time saved and emitted records changed.

## Findings so far (November 2026)

Intel i7-9750H (AVX2), hg38, one thread, before commit 3bceef4:

- The DP loop itself is about 3–4% of per-read time on HiFi reads, and about
  9% with 1% substitution errors. Seeding is about 60%; gap alignments called
  from the DP are about 25–40%.
- Against A, D runs 1.5–1.7× faster, B1 2.5–3.2× and C1 2.9–4.0×, depending
  on build flags and gap-lookup hashing; that is about 2–6% of per-read time
  end to end. Building the struct-of-arrays copy costs about 5% of C1's time.
  With f32, 0.67% of calls choose a different chain, all near-ties (at most
  0.011 bits apart).
- On 1,000 PacBio Revio HiFi reads the DP made 435 gap alignments per read,
  27 µs each. Gaps of 512 bp or more were 24% of alignments but 60% of their
  time. 31% of alignments repeated a gap already aligned for the same read,
  because the cache is keyed by seed pair rather than gap coordinates.
- k-mer proxies separate homologous from unrelated gaps well (AUC 0.97–0.99),
  but cannot replace the alignment: the DP's choices are near-ties, and
  replacing alignments changed the emitted records in about 90% of calls.
  The 4-mer spectrum also saturates on long gaps.
- The upper bound now in `best_chain` (`gap_bits_bound`, commit 3bceef4)
  avoids 15.6% of gap alignments and about 20% of their time, with identical
  output.

Wall-clock comparisons of whole runs on a laptop varied by ±20% between
repeats; compare work counts, or repeated interleaved runs on an idle machine.
