"""Generate the global-chaining poster figure as a standalone SVG.

Regenerate with: python3 poster/global_chaining.py poster/global-chaining.svg
"""
import sys

W, H = 1600, 1000
FONT = "Helvetica, Arial, sans-serif"

PLUS = "#1f6fb2"     # + strand seeds and records
MINUS = "#d9730d"    # - strand seeds and records
GREY = "#9aa0a6"     # unchosen seeds
INK = "#222222"
MUTED = "#5f6368"
EDGE = "#333333"
REJECT = "#c0392b"
PANEL = "#f4f6f9"


def X(r):            # read position (kb) -> x
    return 130 + 70 * r


def Y7(f):           # chr7 reference position (kb) -> y
    return 700 - 40 * f


def Y3(f):           # chr3 reference position (kb) -> y
    return 232 - 36 * f


out = []


def add(s):
    out.append(s)


def text(x, y, s, size=18, anchor="start", weight="normal", fill=INK, style="normal"):
    if "₂" in s:
        head, tail = s.split("₂", 1)
        d = size * 0.3
        s = (f'{head}<tspan dy="{d:.1f}" font-size="{size * 0.7:.1f}">2</tspan>'
             f'<tspan dy="{-d:.1f}">{tail}</tspan>')
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


# ── Seeds: (name, read_start, read_end, ref_at_start, ref_at_end, colour, band)
seeds = {
    "s1": (0.4, 1.0, 0.4, 1.0, PLUS, 7),
    "s2": (1.5, 2.1, 1.5, 2.1, PLUS, 7),
    "s3": (2.5, 3.1, 2.8, 3.4, PLUS, 7),     # after an indel shift
    "s4": (4.4, 5.0, 4.7, 5.3, PLUS, 7),     # after a seedless gap
    "s5": (5.2, 5.8, 6.7, 6.1, MINUS, 7),    # inverted segment
    "s6": (6.0, 6.5, 5.9, 5.4, MINUS, 7),
    "s7": (6.9, 7.5, 7.2, 7.8, PLUS, 7),
    "s8": (8.3, 8.9, 8.6, 9.2, PLUS, 7),
    "s9": (9.3, 9.9, 9.6, 10.2, PLUS, 7),
    "d1": (1.5, 2.1, 0.5, 1.1, GREY, 3),     # repeat copies on chr3
    "d2": (2.5, 3.1, 1.5, 2.1, GREY, 3),
    "n1": (7.6, 8.0, 2.2, 2.6, GREY, 7),     # isolated off-diagonal hit
}


def ends(name):
    r0, r1, f0, f1, _, band = seeds[name]
    Y = Y7 if band == 7 else Y3
    return (X(r0), Y(f0)), (X(r1), Y(f1))


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'role="img" aria-label="Parallax global chaining: one dynamic-programming pass over every '
    f'seed on every chromosome and strand, costing each move between seeds in bits, so indels, '
    f'structural breaks and repeats compete on one scale.">')
add('<defs><marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" '
    f'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{EDGE}"/>'
    '</marker>'
    '<marker id="arrow-red" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" '
    f'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{REJECT}"/>'
    '</marker></defs>')
add(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')

# ── Title
text(60, 66, "Global chaining: the read’s best explanation, scored in bits", 34, weight="bold")
text(60, 100, "One dynamic-programming pass over every seed on every chromosome and strand.", 19,
     fill=MUTED)
text(60, 124, "Each move between seeds costs −log₂ of its probability, so indels, structural "
     "breaks and repeats compete on one scale.", 19, fill=MUTED)

# ── Dot plot frame
add(f'<rect x="{X(0)}" y="{Y3(2.4):.1f}" width="{X(10) - X(0)}" height="{Y3(0) - Y3(2.4):.1f}" fill="#f5f5f5"/>')
add(f'<rect x="{X(0)}" y="{Y7(10.5):.1f}" width="{X(10) - X(0)}" height="{Y7(0) - Y7(10.5):.1f}" fill="#f6f8fb"/>')
text(X(0) - 12, (Y3(0) + Y3(2.4)) / 2 + 6, "chr3", 18, anchor="end", weight="bold", fill=MUTED)
text(X(0) - 12, (Y7(0) + Y7(10.5)) / 2 + 6, "chr7", 18, anchor="end", weight="bold", fill=MUTED)
add(f'<text transform="translate(40 {(Y3(2.4) + Y7(0)) / 2:.1f}) rotate(-90)" font-family="{FONT}" '
    f'font-size="18" fill="{MUTED}" text-anchor="middle">reference position</text>')

# read axis
ay = Y7(0) + 12
line(X(0), ay, X(10), ay, stroke=MUTED, width=1.5, cap="butt")
for kb in range(0, 11, 2):
    line(X(kb), ay, X(kb), ay + 7, stroke=MUTED, width=1.5, cap="butt")
    text(X(kb), ay + 27, f"{kb} kb" if kb in (0, 10) else str(kb), 16, anchor="middle", fill=MUTED)
text((X(0) + X(10)) / 2, ay + 52, "read position", 18, anchor="middle", fill=MUTED)

# ── Rejected move: jumping to the chr3 repeat copy
(_, s1e) = ends("s1")
(d1s, _) = ends("d1")
line(s1e[0] + 4, s1e[1] - 8, d1s[0] - 6, d1s[1] + 10, stroke=REJECT, width=2.5, dash="7 6")
t = 0.62
mx = s1e[0] + 4 + t * (d1s[0] - 6 - s1e[0] - 4)
my = s1e[1] - 8 + t * (d1s[1] + 10 - s1e[1] + 8)
line(mx - 9, my - 9, mx + 9, my + 9, stroke=REJECT, width=3)
line(mx - 9, my + 9, mx + 9, my - 9, stroke=REJECT, width=3)
lines(mx + 20, my - 12, ["other chromosome: −66.5 bits", "+ second locus: −40 bits"], 16, fill=REJECT)
text(X(3.35), Y3(2.4) + 26, "repeat copies: matched, not chosen", 15, fill=MUTED, style="italic")

# ── Chosen moves
def move(a, b, dash=None):
    (_, ae) = ends(a)
    (bs, _) = ends(b)
    line(ae[0] + 5, ae[1], bs[0] - 5, bs[1], dash=dash, marker=dash is None)
    return ae, bs


ae, bs = move("s1", "s2")
text((ae[0] + bs[0]) / 2 + 10, (ae[1] + bs[1]) / 2 + 24, "same diagonal: −8", 15)
ae, bs = move("s2", "s3")
text((ae[0] + bs[0]) / 2 + 12, (ae[1] + bs[1]) / 2 + 26, "indel shift Δ: −(6 + Δ)", 15)
ae, bs = move("s3", "s4", dash="3 6")
lines((ae[0] + bs[0]) / 2 + 14, (ae[1] + bs[1]) / 2 + 22,
      ["seedless gap, aligned:", "+2 per match, −8 per mismatch"], 15)
ae, bs = move("s4", "s5")
lines(ae[0] - 112, ae[1] - 40, ["strand flip:", "−(36 + log₂ d)"], 15)
ae, bs = move("s5", "s6")
ae, bs = move("s6", "s7")
lines(bs[0] + 10, bs[1] + 34, ["strand flip:", "−(36 + log₂ d)"], 15)
move("s7", "s8")
move("s8", "s9")

# seeds on top of the moves
for name, (r0, r1, f0, f1, colour, band) in seeds.items():
    (p, q) = ends(name)
    line(p[0], p[1], q[0], q[1], stroke=colour, width=8)

(s9s, _) = ends("s9")
text(s9s[0] - 14, s9s[1] - 14, "each seed earns +2 bits per new read base", 15, anchor="end",
     fill=PLUS)
(n1s, n1e) = ends("n1")
lines(n1e[0] + 12, n1e[1] + 4, ["isolated hit:", "left unexplained"], 15, fill=MUTED, style="italic")

# ── Output records along the read axis
ry = 800
text(X(0) - 12, ry + 17, "records", 18, anchor="end", weight="bold", fill=MUTED)


def record(r0, r1, colour, label, members):
    add(f'<rect x="{X(r0):.1f}" y="{ry}" width="{X(r1) - X(r0):.1f}" height="22" rx="3" '
        f'fill="{colour}" fill-opacity="0.28" stroke="{colour}" stroke-width="1.5"/>')
    for m in members:
        a, b = seeds[m][0], seeds[m][1]
        add(f'<rect x="{X(a):.1f}" y="{ry}" width="{X(b) - X(a):.1f}" height="22" fill="{colour}"/>')
    lines((X(r0) + X(r1)) / 2, ry + 46, label, 15, anchor="middle")


record(0.4, 5.0, PLUS, ["chr7 + primary", "gaps aligned in one record"], ["s1", "s2", "s3", "s4"])
record(5.2, 6.5, MINUS, ["chr7 −", "supplementary"], ["s5", "s6"])
record(6.9, 9.9, PLUS, ["chr7 + supplementary"], ["s7", "s8", "s9"])

# legend
ly = 925
lx = X(0)
for colour, label, kind in [(PLUS, "+ strand seed", "seed"), (MINUS, "− strand seed", "seed"),
                            (GREY, "unchosen seed", "seed"), (EDGE, "chosen move", "edge"),
                            (REJECT, "rejected move", "reject")]:
    if kind == "seed":
        line(lx, ly, lx + 28, ly - 14, stroke=colour, width=8)
    elif kind == "edge":
        line(lx, ly - 7, lx + 28, ly - 7, marker=True)
    else:
        line(lx, ly - 7, lx + 28, ly - 7, stroke=REJECT, dash="7 6")
    text(lx + 38, ly - 1, label, 15, fill=MUTED)
    lx += 150
text(X(0), ly + 34, "Schematic, not to scale. Default parameters; bits = −log₂ P.", 14, fill=MUTED,
     style="italic")

# ── Right panel: how a chain is scored
px = 930
pw = 620


def heading(y, num, s):
    add(f'<circle cx="{px + 15}" cy="{y - 7}" r="15" fill="{INK}"/>')
    text(px + 15, y - 1, num, 17, anchor="middle", weight="bold", fill="#ffffff")
    text(px + 42, y, s, 22, weight="bold")


heading(160, "1", "Extend the chain one seed at a time")
add(f'<rect x="{px}" y="180" width="{pw}" height="100" rx="6" fill="{PANEL}"/>')
# Cases layout: the max's two alternatives stacked after a drawn brace.
case1, case2 = 216, 256
mid = (case1 + case2) / 2 + 1
text(px + 20, mid, "score(i) = max", 21)
bx, bw = px + 170, 8
yt, yb, ym = case1 - 24, case2 + 10, mid - 7
add(f'<path d="M{bx + 2 * bw},{yt} q{-bw},0 {-bw},{bw} V{ym - bw} q0,{bw} {-bw},{bw} '
    f'q{bw},0 {bw},{bw} V{yb - bw} q0,{bw} {bw},{bw}" fill="none" stroke="{INK}" '
    f'stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/>')
text(bx + 3 * bw + 6, case1, "2·len(i)", 21)
text(bx + 3 * bw + 6, case2, "score(j) + 2·new(j, i) − cost(j → i)", 21)
lines(px, 314, ["The first case starts a chain at seed i; the second extends the best",
                "chain ending at an earlier seed j. new(j, i) counts the read bases",
                "i covers past the end of j, so no base is credited twice. Seeds are",
                "sorted by read end and scanned once, across all chromosomes and strands."],
      16, fill=MUTED)

y2 = 430
heading(y2, "2", "Cost every move in bits")
# Each group is one shaded band; its first row names the move, and indented
# rows that follow give per-operation values.
groups = [
    [("same diagonal, read gap (≥ 1 mismatch)", "8")],
    [("indel: diagonal shift Δ", "6 + Δ")],
    [("seedless gap 20–2000 bp, aligned:", ""),
     ("per matched base (a gain)", "−2"),
     ("per mismatched base", "8"),
     ("per indel of n bp", "6 + n")],
    [("break, same strand: jump Δ", "35 + log₂(1+Δ)")],
    [("break, strand flip: distance d", "36 + log₂(1+d)")],
    [("break, other chromosome", "35 + 31.5")],
]
ROW = 34
SUB = 26
ty = y2 + 38
for g, group in enumerate(groups):
    height = ROW + SUB * (len(group) - 1)
    if g % 2 == 0:
        add(f'<rect x="{px}" y="{ty - 24}" width="{pw}" height="{height}" fill="{PANEL}"/>')
    for r, (label, bits) in enumerate(group):
        indent = 34 if r > 0 else 14
        text(px + indent, ty, label, 17 if r == 0 else 16)
        if bits:
            text(px + pw - 190, ty, bits, 17, weight="bold")
        ty += SUB if r < len(group) - 1 else ROW
text(px + 14, ty + 6, "An indel is costed as a break whenever that is cheaper.", 15, fill=MUTED,
     style="italic")

y3 = ty + 48
heading(y3, "3", "Prefer few loci")
lines(px, y3 + 36, ["The k-th distinct locus costs 40·(k − 1) bits; a locus is the seeds",
                "within 20 kb of each other on one chromosome. Loci are dropped one",
                "at a time, re-running the pass, while the score improves."], 16, fill=MUTED)

y4 = y3 + 122
heading(y4, "4", "Trace back and emit records")
lines(px, y4 + 36, ["Moves within one locus (≤ 5 kb) stay in one record and their gaps",
                "are aligned; any other move starts a supplementary record."], 16, fill=MUTED)

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out) + "\n")
