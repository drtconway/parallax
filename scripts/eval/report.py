"""Aggregate evaluate.nf outputs into a summary table, plots and an HTML report.

Usage:
    report.py <dir> [--out OUTDIR]

Input files are found in <dir> (searched recursively) by suffix.  Their names
start with <sample>.<ref>.<aligner>., which is how results are grouped:

  *.alignments.tsv.gz *.reads.tsv.gz *.indels.tsv *.errpos.tsv   (alignment_stats.py)
  *.<variant>.invariance.tsv.gz                                 (invariance.py)
  *.hist.tsv *.totals.tsv *.windows.tsv *.sites.bed               (parallax eval-pileup)

Writes OUTDIR/summary.tsv (long: sample, ref, aligner, metric, stratum,
value), OUTDIR/summary_wide.tsv (one column per sample/ref/aligner),
OUTDIR/plots/*.png and OUTDIR/report.html.
"""

import argparse
import html
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ALL = "all"
COLORS = ["#2a6fdb", "#e0672b", "#2f9e5a", "#9b4dca", "#c23b3b", "#7a7a7a"]
STYLES = {0: "-", 1: "--", 2: ":", 3: "-."}
# Invariance classes, best first, with fixed colours so panels are comparable.
INVARIANCE_COLORS = {
    "identical": "#2f9e5a", "agree": "#2f9e5a",
    "near": "#8fd18f",
    "partial": "#f2c14e",
    "disagree_lowmapq": "#e0672b",
    "disjoint": "#c23b3b", "disagree_confident": "#c23b3b",
    "orig_only": "#9b4dca", "unaligned": "#9b4dca",
    "variant_only": "#5b8def", "parent_unaligned": "#5b8def",
    "both_unmapped": "#bbbbbb", "both_unaligned": "#bbbbbb",
}
INVARIANCE_CLASSES = list(INVARIANCE_COLORS)


def key_of(path):
    parts = path.name.split(".")
    return tuple(parts[:3])


def find(root, suffix):
    groups = defaultdict(list)
    for p in sorted(Path(root).rglob(f"*{suffix}")):
        groups[key_of(p)].append(p)
    return groups


def read_concat(paths, **kw):
    frames = [pd.read_csv(p, sep="\t", **kw) for p in paths]
    frames = [f for f in frames if len(f)]
    return pd.concat(frames, ignore_index=True) if frames else None


class Summary:
    def __init__(self):
        self.rows = []

    def add(self, key, metric, value, stratum=ALL):
        if value is None or (isinstance(value, float) and not math.isfinite(value)):
            return
        self.rows.append((*key, metric, stratum, value))

    def frame(self):
        return pd.DataFrame(self.rows, columns=["sample", "ref", "aligner", "metric", "stratum", "value"])


def ratio(a, b):
    return float(a) / float(b) if b else float("nan")


# ---------------------------------------------------------------- metrics


def read_metrics(s, key, reads):
    total_len = reads.read_len.sum()
    mapped = reads[reads.mapped == 1]
    s.add(key, "reads", len(reads))
    s.add(key, "read_bases", int(total_len))
    s.add(key, "frac_reads_mapped", ratio(len(mapped), len(reads)))
    s.add(key, "frac_bases_aligned", ratio(reads.aligned_union.sum(), total_len))
    s.add(key, "frac_bases_multiply_aligned", ratio((reads.aligned_sum - reads.aligned_union).sum(), total_len))
    s.add(key, "frac_bases_interior_gap", ratio(reads.interior_gap.sum(), total_len))
    s.add(key, "frac_bases_end_clipped", ratio(mapped.clip5.sum() + mapped.clip3.sum(), total_len))
    s.add(key, "alignments_per_mapped_read", ratio((mapped.n_primary + mapped.n_supplementary).sum(), len(mapped)))
    s.add(key, "frac_mapped_reads_split", ratio((mapped.n_supplementary > 0).sum(), len(mapped)))
    s.add(key, "secondary_per_mapped_read", ratio(mapped.n_secondary.sum(), len(mapped)))
    s.add(key, "outlier_windows_per_Mb_aligned", 1e6 * ratio(reads.n_outlier_windows.sum(), reads.aligned_union.sum()))
    s.add(key, "frac_reads_with_outlier_window", ratio((mapped.n_outlier_windows > 0).sum(), len(mapped)))
    s.add(key, "split_indels_per_Mb_aligned", 1e6 * ratio(reads.split_events.sum(), reads.aligned_union.sum()))


def alignment_metrics(s, key, alns):
    alns = alns.assign(qlen=alns.qe - alns.qs)
    for cls, g in alns.groupby("class"):
        s.add(key, f"{cls}_count", len(g))
        s.add(key, f"{cls}_median_qlen", float(g.qlen.median()))
        s.add(key, f"{cls}_mean_mapq", float(g.mapq.mean()))
    main = alns[alns["class"] != "secondary"]
    groups = [(ALL, main)] + list(main.groupby("stratum"))
    for stratum, g in groups:
        m, x = g.matches.sum(), g.mismatches.sum()
        ib, db, ie, de = g.ins_bases.sum(), g.del_bases.sum(), g.ins_events.sum(), g.del_events.sum()
        qbases = g.qlen.sum()
        s.add(key, "aligned_query_bases", int(qbases), stratum)
        s.add(key, "identity", ratio(m, m + x + ib + db), stratum)
        s.add(key, "gap_compressed_identity", ratio(m, m + x + ie + de), stratum)
        s.add(key, "mismatches_per_kb", 1e3 * ratio(x, qbases), stratum)
        s.add(key, "ins_events_per_kb", 1e3 * ratio(ie, qbases), stratum)
        s.add(key, "del_events_per_kb", 1e3 * ratio(de, qbases), stratum)
        prim = g[g["class"] == "primary"]
        s.add(key, "primary_frac_mapq0", ratio((prim.mapq == 0).sum(), len(prim)), stratum)
        s.add(key, "primary_frac_mapq_ge20", ratio((prim.mapq >= 20).sum(), len(prim)), stratum)


def indel_metrics(s, key, indels, aligned_bases):
    main = indels[indels["class"] != "secondary"]
    for kind in ("ins", "del"):
        k = main[main.kind == kind]
        big = k[k.band >= 64]
        s.add(key, f"{kind}_ge64_per_Mb", 1e6 * ratio(big.events.sum(), aligned_bases))
        s.add(key, f"{kind}_ge64_frac_as_split", ratio(big[big.source == "split"].events.sum(), big.events.sum()))
        small = k[k.band < 4]
        s.add(key, f"{kind}_lt4_per_kb", 1e3 * ratio(small.events.sum(), aligned_bases))


def errpos_metrics(s, key, err):
    def rate(g):
        return ratio(g.mismatches.sum() + g.ins_events.sum() + g.del_events.sum(), g.positions.sum())

    interior = err[err.dist_bin >= 1024]
    base = rate(interior)
    s.add(key, "error_rate_interior", base)
    for kind, g in err.groupby("end_kind"):
        near = g[g.dist_bin < 32]
        s.add(key, f"error_rate_within32_{kind}", rate(near))
        s.add(key, f"error_excess_within32_{kind}", ratio(rate(near), base))


def invariance_metrics(s, key, variant, inv):
    counts = inv["class"].value_counts()
    total = counts.sum()
    for cls, n in counts.items():
        s.add(key, f"{variant}_frac_{cls}", ratio(n, total))
    scores = pd.to_numeric(inv.score, errors="coerce")
    s.add(key, f"{variant}_mean_score", float(scores.mean()))
    if variant.startswith("frag"):
        for stratum, g in inv.groupby("stratum"):
            s.add(key, f"{variant}_frac_agree", ratio((g["class"] == "agree").sum(), len(g)), stratum)
            s.add(key, f"{variant}_frac_disagree_confident",
                  ratio((g["class"] == "disagree_confident").sum(), len(g)), stratum)


def pileup_metrics(s, key, totals, hist, windows):
    for mq, t in totals.groupby("mapq_class"):
        groups = [(ALL, t)] + list(t.groupby("stratum"))
        for stratum, g in groups:
            g = g.sum(numeric_only=True)
            s.add(key, f"pileup_{mq}_mean_depth", ratio(g.depth_sum, g.positions), stratum)
            s.add(key, f"pileup_{mq}_nonref_rate",
                  ratio(g.mismatch_sum + g.del_sum, g.depth_sum), stratum)
            s.add(key, f"pileup_{mq}_ins_rate", ratio(g.ins_sum, g.depth_sum), stratum)
            s.add(key, f"pileup_{mq}_candidate_sites_per_Mb",
                  1e6 * ratio(g.candidate_sites, g.called_positions), stratum)
            s.add(key, f"pileup_{mq}_major_nonref_per_Mb",
                  1e6 * ratio(g.major_nonref, g.called_positions), stratum)
            s.add(key, f"pileup_{mq}_frac_called", ratio(g.called_positions, g.positions), stratum)
    depth = hist[(hist.kind == "depth") & (hist.mapq_class == "all")]
    if len(depth):
        d = depth.groupby("bin")["count"].sum().sort_index()
        cum = d.cumsum() / d.sum()
        median = float(cum[cum >= 0.5].index[0])
        s.add(key, "pileup_median_depth", median)
        s.add(key, "pileup_frac_depth0", ratio(d.get(0, 0), d.sum()))
        low = d[d.index < 0.2 * median].sum()
        s.add(key, "pileup_frac_depth_lt_0.2x_median", ratio(low, d.sum()))
    if windows is not None and len(windows):
        col = [c for c in windows.columns if c.startswith("mean_depth_")][0]
        w = windows[windows.positions > 0][col]
        med = w.median()
        if med > 0:
            s.add(key, "window_depth_cv", float(w.std() / w.mean()))
            s.add(key, "frac_windows_gt_3x_median", ratio((w > 3 * med).sum(), len(w)))
            s.add(key, "frac_windows_lt_0.2x_median", ratio((w < 0.2 * med).sum(), len(w)))


# ---------------------------------------------------------------- plots


class Plotter:
    def __init__(self, outdir):
        self.dir = Path(outdir) / "plots"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.files = []

    def style(self, key, keys):
        """Colour by aligner, line style by reference."""
        aligners = sorted({k[2] for k in keys})
        refs = sorted({k[1] for k in keys})
        return dict(color=COLORS[aligners.index(key[2]) % len(COLORS)],
                    linestyle=STYLES[refs.index(key[1]) % len(STYLES)],
                    label=f"{key[2]} / {key[1]}")

    def save(self, fig, name, title, caption):
        path = self.dir / f"{name}.png"
        fig.tight_layout()
        fig.savefig(path, dpi=110)
        plt.close(fig)
        self.files.append((f"plots/{name}.png", title, caption))


def plot_lengths(p, sample, data):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
    bins = np.logspace(1, 6, 60)
    for ax, cls in zip(axes, ["primary", "supplementary", "secondary"]):
        for key, alns in data.items():
            q = (alns.qe - alns.qs)[alns["class"] == cls]
            if len(q):
                ax.hist(q, bins=bins, histtype="step", **p.style(key, data))
        ax.set_xscale("log")
        ax.set_title(cls)
        ax.set_xlabel("aligned query length (bp)")
    axes[0].set_ylabel("alignments")
    axes[0].legend(fontsize=8)
    p.save(fig, f"{sample}.lengths", f"{sample}: alignment lengths",
           "Aligned query length by alignment class.")


def plot_identity_curve(p, sample, data, read_bases):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for key, alns in data.items():
        main = alns[alns["class"] != "secondary"]
        ident = pd.to_numeric(main.identity, errors="coerce")
        order = np.argsort(-ident.to_numpy())
        ids = ident.to_numpy()[order]
        cum = np.cumsum((main.qe - main.qs).to_numpy()[order]) / read_bases[key]
        ax.plot(ids, cum, **p.style(key, data))
    ax.set_xlim(1.0, 0.7)
    ax.set_xlabel("alignment identity (at least)")
    ax.set_ylabel("fraction of read bases aligned")
    ax.legend(fontsize=8)
    p.save(fig, f"{sample}.identity_curve", f"{sample}: aligned bases vs identity",
           "Cumulative read bases in primary/supplementary alignments of at least the given identity. "
           "Higher and further left is better; a curve that only rises at low identity is aligning junk.")


def plot_indels(p, sample, data, aligned):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    for ax, kind in zip(axes, ["ins", "del"]):
        for key, ind in data.items():
            k = ind[(ind.kind == kind) & (ind["class"] != "secondary")]
            st = p.style(key, data)
            both = k.groupby("band").events.sum()
            ax.plot(both.index, 1e6 * both / aligned[key], marker="o", ms=3, **st)
            split = k[k.source == "split"].groupby("band").events.sum()
            if len(split):
                ax.plot(split.index, 1e6 * split / aligned[key], marker="x", ms=4, alpha=0.6,
                        color=st["color"], linestyle=":")
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_title(f"{kind} (o = CIGAR + split, x = split only)")
        ax.set_xlabel("event size band (bp, lower bound)")
    axes[0].set_ylabel("events per Mb aligned")
    axes[0].legend(fontsize=8)
    p.save(fig, f"{sample}.indels", f"{sample}: indel size spectrum",
           "Indel events per Mb of aligned query by power-of-two size band. Split = implied by "
           "colinear supplementary alignments, so the two aligners' representations are comparable.")


def plot_errpos(p, sample, data):
    kinds = ["read_end", "junction", "clipped"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for ax, kind in zip(axes, kinds):
        for key, err in data.items():
            g = err[err.end_kind == kind].sort_values("dist_bin")
            if not len(g):
                continue
            rate = (g.mismatches + g.ins_events + g.del_events) / g.positions.replace(0, np.nan)
            rate = rate.replace(0, np.nan)  # log axis: leave gaps rather than cliffs
            ax.plot(g.dist_bin + 1, rate, marker="o", ms=3, **p.style(key, data))
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_title(f"nearest end: {kind}")
        ax.set_xlabel("distance from alignment end + 1 (bp)")
    axes[0].set_ylabel("error events per aligned base")
    axes[0].legend(fontsize=8)
    p.save(fig, f"{sample}.errpos", f"{sample}: error rate vs distance from alignment end",
           "A good aligner's curve is flat: raised error near ends (especially at split junctions or "
           "clips) means alignments are being extended through sequence that does not belong.")


def plot_af(p, sample, data, mq):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, kind in zip(axes, ["major_af", "minor_af", "ins_af"]):
        for key, hist in data.items():
            h = hist[(hist.kind == kind) & (hist.mapq_class == mq)]
            if not len(h):
                continue
            h = h.groupby("bin")["count"].sum()
            ax.plot(h.index, h / h.sum(), **p.style(key, data))
        ax.set_yscale("log")
        ax.set_title(f"{kind} ({mq}, depth >= min)")
        ax.set_xlabel("allele fraction")
    axes[0].set_ylabel("fraction of called positions")
    axes[0].legend(fontsize=8)
    p.save(fig, f"{sample}.af_{mq}", f"{sample}: pileup allele fractions ({mq})",
           "On a haploid sample with a matching reference (CHM13), the major allele fraction should sit "
           "near 1 - error everywhere; mass at 0.5 or a heavy minor-AF tail means misplaced reads.")


def plot_depth(p, sample, data):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for key, hist in data.items():
        h = hist[(hist.kind == "depth") & (hist.mapq_class == "all")].groupby("bin")["count"].sum()
        if len(h):
            ax.plot(h.index, h / h.sum(), **p.style(key, data))
    ax.set_yscale("log")
    ax.set_xlim(0, 200)
    ax.set_xlabel("depth")
    ax.set_ylabel("fraction of positions")
    ax.legend(fontsize=8)
    p.save(fig, f"{sample}.depth", f"{sample}: depth distribution",
           "Spikes well above the mode are collapsed repeats attracting reads; mass near zero is dropout.")


def plot_invariance(p, sample, data):
    variants = sorted({v for (_, v) in data})
    keys = sorted({k for (k, _) in data})
    if not variants:
        return
    fig, axes = plt.subplots(1, len(variants), figsize=(5 * len(variants), 4.5), squeeze=False)
    for ax, variant in zip(axes[0], variants):
        fracs = {}
        for key in keys:
            inv = data.get((key, variant))
            if inv is not None:
                fracs[f"{key[2]}\n{key[1]}"] = inv["class"].value_counts(normalize=True)
        frame = pd.DataFrame(fracs).fillna(0).T
        order = [c for c in INVARIANCE_CLASSES if c in frame.columns]
        order += [c for c in frame.columns if c not in order]
        frame = frame[order]
        colors = [INVARIANCE_COLORS.get(c, "#999999") for c in order]
        frame.plot.bar(stacked=True, ax=ax, legend=False, color=colors, rot=0)
        ax.set_title(variant)
        ax.set_ylim(0, 1)
        ax.legend(fontsize=7, loc="lower right")
    axes[0][0].set_ylabel("fraction of reads / tiles")
    p.save(fig, f"{sample}.invariance", f"{sample}: invariance",
           "rc: reads vs their reverse complements. fragN: N bp tiles vs the parent read's alignment "
           "over the same bases. disagree_confident tiles are placed elsewhere with MAPQ >= 20.")


# ---------------------------------------------------------------- html


def write_html(outdir, wide, plots):
    def table(df):
        return df.to_html(float_format=lambda v: f"{v:.4g}", na_rep="", border=0, classes="t")

    overall = wide[wide.index.get_level_values("stratum") == ALL].droplevel("stratum")
    strat = wide[wide.index.get_level_values("stratum") != ALL]
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>Aligner evaluation</title>",
        "<style>body{font-family:system-ui,sans-serif;margin:24px;max-width:1500px}"
        "table.t{border-collapse:collapse;font-size:12px}"
        "table.t td,table.t th{padding:2px 8px;border-bottom:1px solid #ddd;text-align:right}"
        "table.t th{text-align:left}figure{margin:24px 0}img{max-width:100%}"
        "figcaption{color:#555;font-size:13px;max-width:900px}</style></head><body>",
        "<h1>Aligner evaluation</h1>",
        "<p>Truth-free metrics. Columns are sample / reference / aligner. "
        "Full data in <code>summary.tsv</code>.</p>",
        "<h2>Overall</h2>", table(overall),
    ]
    parts.append("<h2>Plots</h2>")
    for path, title, caption in plots:
        parts.append(f"<figure><h3>{html.escape(title)}</h3><img src='{path}'>"
                     f"<figcaption>{html.escape(caption)}</figcaption></figure>")
    if len(strat):
        parts += ["<h2>By stratum</h2>", table(strat)]
    parts.append("</body></html>")
    (Path(outdir) / "report.html").write_text("\n".join(parts))


# ---------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dir")
    ap.add_argument("--out", default=".")
    args = ap.parse_args()
    Path(args.out).mkdir(parents=True, exist_ok=True)

    s = Summary()
    plots = Plotter(args.out)
    aln_data, indel_data, err_data, hist_data = {}, {}, {}, {}
    inv_data = {}
    read_bases, aligned = {}, {}

    for key, paths in find(args.dir, ".reads.tsv.gz").items():
        reads = read_concat(paths)
        if reads is not None:
            read_metrics(s, key, reads)
            read_bases[key] = reads.read_len.sum()
    for key, paths in find(args.dir, ".alignments.tsv.gz").items():
        alns = read_concat(paths)
        if alns is not None:
            alignment_metrics(s, key, alns)
            aln_data[key] = alns
            main_ = alns[alns["class"] != "secondary"]
            aligned[key] = (main_.qe - main_.qs).sum()
    for key, paths in find(args.dir, ".indels.tsv").items():
        ind = read_concat(paths)
        if ind is not None and aligned.get(key):
            ind = ind.groupby(["class", "stratum", "kind", "source", "band"], as_index=False).sum()
            indel_metrics(s, key, ind, aligned[key])
            indel_data[key] = ind
    for key, paths in find(args.dir, ".errpos.tsv").items():
        err = read_concat(paths)
        if err is not None:
            err = err.groupby(["end_kind", "dist_bin"], as_index=False).sum()
            errpos_metrics(s, key, err)
            err_data[key] = err
    for key, paths in find(args.dir, ".invariance.tsv.gz").items():
        by_variant = defaultdict(list)
        for p in paths:
            by_variant[p.name.split(".")[-4]].append(p)
        for variant, vpaths in by_variant.items():
            inv = read_concat(vpaths)
            if inv is not None:
                invariance_metrics(s, key, variant, inv)
                inv_data[(key, variant)] = inv
    totals = find(args.dir, ".totals.tsv")
    hists = find(args.dir, ".hist.tsv")
    windows = find(args.dir, ".windows.tsv")
    for key in totals:
        t = read_concat(totals[key])
        h = read_concat(hists.get(key, []))
        w = read_concat(windows.get(key, [])) if key in windows else None
        if t is not None and h is not None:
            pileup_metrics(s, key, t, h, w)
            hist_data[key] = h

    summary = s.frame()
    summary.to_csv(Path(args.out) / "summary.tsv", sep="\t", index=False)
    if len(summary):
        wide = summary.pivot_table(index=["metric", "stratum"], columns=["sample", "ref", "aligner"],
                                   values="value", aggfunc="first")
        wide.columns = ["/".join(c) for c in wide.columns]
        metric_order = {m: i for i, m in enumerate(dict.fromkeys(summary.metric))}
        wide = wide.sort_index(key=lambda idx: idx.map(metric_order) if idx.name == "metric" else idx)
    else:
        wide = pd.DataFrame(index=pd.MultiIndex.from_tuples([], names=["metric", "stratum"]))
    wide.to_csv(Path(args.out) / "summary_wide.tsv", sep="\t")

    samples = sorted({k[0] for k in list(aln_data) + list(hist_data)})
    for sample in samples:
        def mine(d):
            return {k: v for k, v in d.items() if k[0] == sample}
        if mine(aln_data):
            plot_lengths(plots, sample, mine(aln_data))
            plot_identity_curve(plots, sample, mine(aln_data), read_bases)
        if mine(indel_data):
            plot_indels(plots, sample, mine(indel_data), aligned)
        if mine(err_data):
            plot_errpos(plots, sample, mine(err_data))
        if mine(hist_data):
            for mq in sorted({m for h in mine(hist_data).values() for m in h.mapq_class.unique()}):
                plot_af(plots, sample, mine(hist_data), mq)
            plot_depth(plots, sample, mine(hist_data))
        inv = {(k, v): d for (k, v), d in inv_data.items() if k[0] == sample}
        if inv:
            plot_invariance(plots, sample, inv)

    write_html(args.out, wide, plots.files)


if __name__ == "__main__":
    main()
