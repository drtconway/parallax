"""Generate the conventional-alignment poster figure as a standalone SVG.

Uses the same example read as global_chaining.py, so the two figures compare
directly.

Regenerate with: python3 poster/conventional_alignment.py poster/conventional-alignment.svg
"""
import sys

W, H = 1600, 1070
FONT = "Helvetica, Arial, sans-serif"

PLUS = "#1f6fb2"     # + strand seeds and records
MINUS = "#d9730d"    # - strand seeds and records
GREY = "#9aa0a6"     # repeat-copy and unchained seeds
INK = "#222222"
MUTED = "#5f6368"
EDGE = "#333333"
CLASH = "#c0392b"    # read bases claimed by two chains
PANEL = "#f4f6f9"


def X(r):            # read position (kb) -> x
    return 130 + 70 * r


def Y7(f):           # chr7 reference position (kb) -> y
    return 700 - 40 * f


def Y3(f):           # chr3 reference position (kb) -> y
    return 242 - 28 * f


out = []


def add(s):
    out.append(s)


def text(x, y, s, size=18, anchor="start", weight="normal", fill=INK, style="normal"):
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{FONT}" font-size="{size}" '
        f'font-weight="{weight}" font-style="{style}" fill="{fill}" text-anchor="{anchor}">{s}</text>')


def lines(x, y, rows, size=18, step=None, **kw):
    step = step or size * 1.35
    for i, row in enumerate(rows):
        text(x, y + i * step, row, size=size, **kw)


def line(x1, y1, x2, y2, stroke=EDGE, width=2.5, dash=None, marker=False, cap="round"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    m = ' marker-end="url(#arrow)"' if marker else ""
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" '
        f'stroke-width="{width}" stroke-linecap="{cap}"{d}{m}/>')


def rect(x, y, w, h, fill, **kw):
    extra = "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items())
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{fill}"{extra}/>')


# ── Seeds: (read_start, read_end, ref_at_start, ref_at_end, colour, band).
# The chr7 seeds are identical to global_chaining.py; the chr3 repeat copy
# here matches a little more of the read.
seeds = {
    "s1": (0.4, 1.0, 0.4, 1.0, PLUS, 7),
    "s2": (1.5, 2.1, 1.5, 2.1, PLUS, 7),
    "s3": (2.5, 3.1, 2.8, 3.4, PLUS, 7),
    "s4": (4.4, 5.0, 4.7, 5.3, PLUS, 7),
    "s5": (5.2, 5.8, 6.7, 6.1, MINUS, 7),
    "s6": (6.0, 6.5, 5.9, 5.4, MINUS, 7),
    "s7": (6.9, 7.5, 7.2, 7.8, PLUS, 7),
    "s8": (8.3, 8.9, 8.6, 9.2, PLUS, 7),
    "s9": (9.3, 9.9, 9.6, 10.2, PLUS, 7),
    "d1": (1.5, 2.1, 0.5, 1.1, GREY, 3),
    "d2": (2.5, 3.1, 1.5, 2.1, GREY, 3),
    "d3": (3.6, 4.2, 2.6, 3.2, GREY, 3),
    "n1": (7.6, 8.0, 2.2, 2.6, GREY, 7),
}

# Local chains: (name, members, colour, read interval after extension,
#                extensions as (read, ref) endpoints beyond the outer seeds)
chains = [
    ("A", ["s1", "s2", "s3", "s4"], PLUS, (0.4, 5.5), [((5.0, 5.3), (5.5, 5.8))]),
    ("B", ["s5", "s6"], MINUS, (4.8, 6.8), [((5.2, 6.7), (4.8, 7.1)), ((6.5, 5.4), (6.8, 5.1))]),
    ("C", ["s7", "s8", "s9"], PLUS, (6.6, 9.9), [((6.9, 7.2), (6.6, 6.9))]),
    ("D", ["d1", "d2", "d3"], GREY, (1.5, 4.2), []),
]


def ends(name):
    r0, r1, f0, f1, _, band = seeds[name]
    Y = Y7 if band == 7 else Y3
    return (X(r0), Y(f0)), (X(r1), Y(f1))


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'role="img" aria-label="Conventional seed-chain-extend alignment: seeds are clustered into '
    f'local chains that are extended and scored one at a time, so chains overlap on the read and '
    f'offer competing explanations of the same read bases, which overlap rules settle only after '
    f'the fact, leaving some read bases in two records.">')
add('<defs><marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" '
    f'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{EDGE}"/>'
    '</marker></defs>')
rect(0, 0, W, H, "#ffffff")

# ── Title
text(60, 66, "Conventional alignment: chains that compete for the read", 34, weight="bold")
text(60, 100, "Seed–chain–extend aligners build and extend each chain on its own, so chains "
     "overlap on the read and offer", 19, fill=MUTED)
text(60, 124, "competing explanations of the same bases, settled afterwards by overlap rules. "
     "The same read as in the global-chaining figure.", 19, fill=MUTED)

# ── Dot plot frame
rect(X(0), Y3(3.4), X(10) - X(0), Y3(0) - Y3(3.4), "#f5f5f5")
rect(X(0), Y7(10.5), X(10) - X(0), Y7(0) - Y7(10.5), "#f6f8fb")
text(X(0) - 12, (Y3(0) + Y3(3.4)) / 2 + 6, "chr3", 18, anchor="end", weight="bold", fill=MUTED)
text(X(0) - 12, (Y7(0) + Y7(10.5)) / 2 + 6, "chr7", 18, anchor="end", weight="bold", fill=MUTED)
add(f'<text transform="translate(40 {(Y3(3.4) + Y7(0)) / 2:.1f}) rotate(-90)" font-family="{FONT}" '
    f'font-size="18" fill="{MUTED}" text-anchor="middle">reference position</text>')

ay = Y7(0) + 12
line(X(0), ay, X(10), ay, stroke=MUTED, width=1.5, cap="butt")
for kb in range(0, 11, 2):
    line(X(kb), ay, X(kb), ay + 7, stroke=MUTED, width=1.5, cap="butt")
    text(X(kb), ay + 27, f"{kb} kb" if kb in (0, 10) else str(kb), 16, anchor="middle", fill=MUTED)
text((X(0) + X(10)) / 2, ay + 52, "read position", 18, anchor="middle", fill=MUTED)

# ── Read intervals two chains claim, as columns from the top of the plot
# down through the records panel, behind everything else.
ry = 800
ROW = 28
NRECORDS = 4
panel_top, panel_bottom = ry - 8, ry + ROW * NRECORDS - 2
clashes = [(1.5, 4.2), (4.8, 5.5), (6.6, 6.8)]
col_top = Y3(3.4)
for a, b in clashes:
    rect(X(a), col_top, X(b) - X(a), panel_bottom - col_top, CLASH, fill_opacity="0.10")
    # Edges only in the records panel, where they don't cross plot labels.
    line(X(a), panel_top, X(a), panel_bottom, stroke=CLASH, width=1, dash="3 3", cap="butt")
    line(X(b), panel_top, X(b), panel_bottom, stroke=CLASH, width=1, dash="3 3", cap="butt")

# ── Chain hulls: a wide translucent stroke through each chain's seeds
for name, members, colour, _, _ in chains:
    pts = []
    for m in members:
        p, q = ends(m)
        pts += [p, q]
    d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    add(f'<polyline points="{d}" fill="none" stroke="{colour}" stroke-opacity="0.16" '
        f'stroke-width="34" stroke-linecap="round" stroke-linejoin="round"/>')

# Each chain's alignment: one solid line over its whole aligned extent,
# including the extensions past its outer seeds.
for name, members, colour, _, exts in chains:
    Y = Y7 if seeds[members[0]][5] == 7 else Y3
    pts = []
    for m in members:
        p, q = ends(m)
        pts += [p, q]
    for (r0, f0), (r1, f1) in exts:
        outward = (X(r1), Y(f1))
        if r0 == seeds[members[0]][0]:
            pts.insert(0, outward)
        else:
            pts.append(outward)
    d = " ".join(f"{x:.1f},{y:.1f}" for x, y in pts)
    add(f'<polyline points="{d}" fill="none" stroke="{colour}" stroke-width="4" '
        f'stroke-linecap="round" stroke-linejoin="round"/>')

# steps between chains: not scored
for a, b in [("s4", "s5"), ("s6", "s7")]:
    (_, ae), (bs, _) = ends(a), ends(b)
    line(ae[0] + 4, ae[1] - 4, bs[0] - 4, bs[1] + 4, stroke=MUTED, width=1.5, dash="2 5")

# seeds
for name, (r0, r1, f0, f1, colour, band) in seeds.items():
    (p, q) = ends(name)
    line(p[0], p[1], q[0], q[1], stroke=colour, width=8)

# chain labels and callouts
text(X(2.3), Y7(1.15), "chain A", 17, weight="bold", fill=PLUS)
text(X(5.45), Y7(4.4), "chain B", 17, weight="bold", fill=MINUS)
text(X(6.95), Y7(9.15), "chain C", 17, weight="bold", fill=PLUS)
text(X(4.45), Y3(3.4) + 24, "chain D: a repeat copy of chain A’s bases", 15, weight="bold",
     fill=MUTED)
lines(X(2.0), Y7(7.6), ["alignments extend past the outer seeds,",
                        "into the neighbouring chain’s bases"], 15, fill=MUTED, style="italic")
line(X(4.4), Y7(7.05), X(5.2), Y7(6.0), stroke=MUTED, width=1, cap="butt")
(n1s, n1e) = ends("n1")
lines(n1e[0] + 12, n1e[1] + 4, ["too few seeds:", "no chain"], 15, fill=MUTED, style="italic")

# ── Records, one row each, with the read bases two records claim
records = [
    ("A", PLUS, "primary", ["s1", "s2", "s3", "s4"]),
    ("B", MINUS, "supplementary", ["s5", "s6"]),
    ("C", PLUS, "supplementary", ["s7", "s8", "s9"]),
    ("D", GREY, "secondary", ["d1", "d2", "d3"]),
]
spans = {name: span for name, _, _, span, _ in chains}
assert len(records) == NRECORDS

for i, (name, colour, kind, members) in enumerate(records):
    y = ry + i * ROW
    r0, r1 = spans[name]
    rect(X(r0), y, X(r1) - X(r0), 18, colour, fill_opacity="0.30", stroke=colour,
         stroke_width="1.5", rx="3")
    for m in members:
        a, b = seeds[m][0], seeds[m][1]
        rect(X(a), y, X(b) - X(a), 18, colour)
    text(X(0) - 8, y + 14, f"{name}: {kind}", 14, anchor="end", weight="bold",
         fill=colour if colour != GREY else MUTED)

ly = panel_bottom + 26
text((X(1.5) + X(4.2)) / 2, ly, "two loci compete for these bases:", 15, anchor="middle",
     fill=CLASH, weight="bold")
text((X(1.5) + X(4.2)) / 2, ly + 20, "D is demoted to secondary", 15, anchor="middle", fill=CLASH)
text((X(4.8) + X(6.8)) / 2 + 30, ly, "these bases appear", 15, anchor="middle", fill=CLASH,
     weight="bold")
text((X(4.8) + X(6.8)) / 2 + 30, ly + 20, "in two records", 15, anchor="middle", fill=CLASH)

# legend
gy = 1012
lx = X(0)
legend = [("seed", PLUS, "+ strand seed"), ("seed", MINUS, "− strand seed"),
          ("seed", GREY, "repeat / unchained seed"), ("aln", PLUS, "alignment"),
          ("clash", CLASH, "claimed twice")]
for kind, colour, label in legend:
    if kind == "seed":
        line(lx, gy, lx + 28, gy - 14, stroke=colour, width=8)
    elif kind == "aln":
        line(lx, gy - 7, lx + 28, gy - 7, stroke=colour, width=4)
    else:
        rect(lx, gy - 15, 28, 16, colour, fill_opacity="0.18", stroke=colour,
             stroke_width="1", stroke_dasharray="3 3")
    text(lx + 38, gy - 1, label, 15, fill=MUTED)
    lx += 52 + 7.6 * len(label)
text(X(0), gy + 34, "Schematic, not to scale; describes seed–chain–extend aligners in general.",
     14, fill=MUTED, style="italic")

# ── Right panel
px = 930
pw = 620


def heading(y, num, s):
    add(f'<circle cx="{px + 15}" cy="{y - 7}" r="15" fill="{INK}"/>')
    text(px + 15, y - 1, num, 17, anchor="middle", weight="bold", fill="#ffffff")
    text(px + 42, y, s, 22, weight="bold")


heading(166, "1", "Cluster seeds into local chains")
lines(px, 202, ["Seeds on one chromosome and strand, close to one diagonal, are",
                "chained within a maximum gap. A strand change, another",
                "chromosome or a large jump ends the chain."], 16, fill=MUTED)

heading(310, "2", "Extend and score each chain alone")
lines(px, 346, ["Each chain is extended past its outer seeds and scored on its",
                "own. Nothing stops one chain’s extension running into the bases",
                "its neighbour explains, and nothing scores the step from one",
                "chain to the next."], 16, fill=MUTED)

heading(476, "3", "Chains overlap on the read")
lines(px, 512, ["So two chains can explain the same read bases: neighbours",
                "overlap at their ends, and a repeat copy (chain D) competes",
                "with the true locus (chain A) for a whole stretch of the read.",
                "Each chain alone looks like a good explanation."], 16, fill=MUTED)

heading(642, "4", "Overlap rules settle it afterwards")
lines(px, 678, ["A chain mostly covered by a better one becomes secondary;",
                "smaller overlaps survive, and those read bases appear in two",
                "records. Which explanation wins depends on thresholds, not on",
                "an account of the whole read."], 16, fill=MUTED)

add(f'<rect x="{px}" y="796" width="{pw}" height="140" rx="6" fill="{PANEL}"/>')
text(px + 20, 830, "Parallax instead", 19, weight="bold")
lines(px + 20, 860, ["builds one explanation of the whole read. Each read base is",
                     "credited at most once, so pieces cannot overlap, and a repeat",
                     "copy competes with the true locus inside the same chain, with",
                     "every step between pieces costed in bits."],
      16, fill=INK)

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out) + "\n")
