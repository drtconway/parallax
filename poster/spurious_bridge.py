"""Generate the spurious-bridge poster figure as a standalone SVG.

Zooms in on the inversion in the example read from global_chaining.py.

Regenerate with: python3 poster/spurious_bridge.py poster/spurious-bridge.svg
"""
import math
import sys

W, H = 1600, 1100
FONT = "Helvetica, Arial, sans-serif"
MONO = "Menlo, Consolas, 'DejaVu Sans Mono', monospace"

PLUS = "#1f6fb2"
MINUS = "#d9730d"
BRIDGE = "#c0392b"
INK = "#222222"
MUTED = "#5f6368"
EDGE = "#333333"
PANEL = "#f4f6f9"
INVERTED = "#fdebd8"


def X(r):            # read position (kb) -> x
    return 130 + 118 * (r - 3.6)


def Y(f):            # chr7 reference position (kb) -> y
    return 700 - 100 * (f - 3.6)


out = []


def add(s):
    out.append(s)


def text(x, y, s, size=18, anchor="start", weight="normal", fill=INK, style="normal", font=FONT):
    if "₂" in s:
        head, tail = s.split("₂", 1)
        d = size * 0.3
        s = (f'{head}<tspan dy="{d:.1f}" font-size="{size * 0.7:.1f}">2</tspan>'
             f'<tspan dy="{-d:.1f}">{tail}</tspan>')
    add(f'<text x="{x:.1f}" y="{y:.1f}" font-family="{font}" font-size="{size}" '
        f'font-weight="{weight}" font-style="{style}" fill="{fill}" text-anchor="{anchor}" '
        f'xml:space="preserve">{s}</text>')


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


# ── A real alignment of two motif-rich sequences, for the base-level view.
# Unrelated stretches built from the same short motifs: an aligner finds a
# path through them because the motifs repeat, not because they share
# ancestry.
REF_SEQ = "GGAATGGAATGGAGTGGAATCGAATGGAAAGGAATGGAATGGAGTGGAGTGGAATGGA"
READ_SEQ = "GGATTGGAATGAATGGATTGGAGGAATGAATTGGAAGGATTGAATGGAGTGGATTGG"


def align(a, b, match=2, mismatch=-4, gap_open=-5, gap_extend=-2):
    """Global alignment with affine gaps (Gotoh); returns aligned strings."""
    n, m = len(a), len(b)
    NEG = -10 ** 9
    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    Ia = [[NEG] * (m + 1) for _ in range(n + 1)]   # gap in b
    Ib = [[NEG] * (m + 1) for _ in range(n + 1)]   # gap in a
    M[0][0] = 0
    for i in range(1, n + 1):
        Ia[i][0] = gap_open + gap_extend * i
    for j in range(1, m + 1):
        Ib[0][j] = gap_open + gap_extend * j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if a[i - 1] == b[j - 1] else mismatch
            M[i][j] = max(M[i - 1][j - 1], Ia[i - 1][j - 1], Ib[i - 1][j - 1]) + s
            Ia[i][j] = max(M[i - 1][j] + gap_open + gap_extend, Ia[i - 1][j] + gap_extend)
            Ib[i][j] = max(M[i][j - 1] + gap_open + gap_extend, Ib[i][j - 1] + gap_extend)
    i, j = n, m
    state = max((M[i][j], 0), (Ia[i][j], 1), (Ib[i][j], 2))[1]
    ra, rb = [], []
    while i > 0 or j > 0:
        if state == 0:
            s = match if a[i - 1] == b[j - 1] else mismatch
            prev = M[i][j] - s
            ra.append(a[i - 1]); rb.append(b[j - 1])
            i, j = i - 1, j - 1
            state = 0 if M[i][j] == prev else (1 if Ia[i][j] == prev else 2)
        elif state == 1:
            ra.append(a[i - 1]); rb.append("-")
            state = 0 if M[i - 1][j] + gap_open + gap_extend == Ia[i][j] else 1
            i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1])
            state = 0 if M[i][j - 1] + gap_open + gap_extend == Ib[i][j] else 2
            j -= 1
    return "".join(reversed(ra)), "".join(reversed(rb))


aln_read, aln_ref = align(READ_SEQ, REF_SEQ)

# ── Bits for the whole 1.9 kb gap (illustrative counts, real costs)
GAP_BP = 1900
B_BP = 1800
MATCHES, MISMATCHES, INDELS, INDEL_BP = 1300, 360, 90, 240
bridge_cost = -2 * MATCHES + 8 * MISMATCHES + (6 * INDELS + INDEL_BP)
flip_a = 36 + math.log2(1 + 1700)
flip_b = 36 + math.log2(1 + 1900)
flips = flip_a + flip_b
inversion_cost = -2 * B_BP + flips

add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'role="img" aria-label="A spurious bridge: chains either side of an inversion share a '
    f'diagonal, so a conventional aligner can join them and force an alignment through the '
    f'inverted sequence, which repeated motifs make possible; scored in bits, the bridge costs far '
    f'more than explaining the inversion.">')
add('<defs><marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" '
    f'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{EDGE}"/>'
    '</marker></defs>')
rect(0, 0, W, H, "#ffffff")

# ── Title
text(60, 66, "Aligning through an inversion: a bridge built from repeated motifs", 34,
     weight="bold")
text(60, 100, "Chains either side of an inversion share a diagonal, so a conventional aligner "
     "can join them and align the gap by force.", 19, fill=MUTED)
text(60, 124, "Repeated motifs let the alignment struggle through, with many indels, although "
     "the sequences are not related.", 19, fill=MUTED)

# ── Dot plot, zoomed on the inversion
rect(X(3.6), Y(8.4), X(8.0) - X(3.6), Y(3.6) - Y(8.4), "#f6f8fb")
rect(X(3.6), Y(7.0), X(8.0) - X(3.6), Y(5.3) - Y(7.0), INVERTED)
text(X(3.6) + 10, Y(7.0) + 22, "inverted in the read", 15, fill="#a0521a", weight="bold")
add(f'<text transform="translate(52 {(Y(3.6) + Y(8.4)) / 2:.1f}) rotate(-90)" font-family="{FONT}" '
    f'font-size="18" fill="{MUTED}" text-anchor="middle">chr7 reference position</text>')

ay = Y(3.6) + 12
line(X(3.6), ay, X(8.0), ay, stroke=MUTED, width=1.5, cap="butt")
for kb in range(4, 9):
    line(X(kb), ay, X(kb), ay + 7, stroke=MUTED, width=1.5, cap="butt")
    text(X(kb), ay + 27, f"{kb} kb", 16, anchor="middle", fill=MUTED)
text((X(3.6) + X(8.0)) / 2, ay + 52, "read position", 18, anchor="middle", fill=MUTED)

# the shared diagonal of chains A and C
line(X(3.6), Y(3.9), X(8.0), Y(8.3), stroke=MUTED, width=1, dash="4 6")

# the bridge: a wobbly path along the diagonal, with indel steps
pts = [(5.0, 5.3)]
steps = 38
for k in range(1, steps):
    r = 5.0 + 1.9 * k / steps
    wobble = 0.07 * math.sin(k * 1.9) + 0.04 * math.sin(k * 4.3)
    pts.append((r, r + 0.3 + wobble))
pts.append((6.9, 7.2))
d = " ".join(f"{X(r):.1f},{Y(f):.1f}" for r, f in pts)
add(f'<polyline points="{d}" fill="none" stroke="{BRIDGE}" stroke-width="3" '
    f'stroke-linejoin="round"/>')

# chains A and C, and B's inverted seeds
seeds = [
    ((3.7, 4.0), (4.1, 4.4), PLUS), ((4.4, 4.7), (5.0, 5.3), PLUS),      # chain A
    ((6.9, 7.2), (7.5, 7.8), PLUS), ((7.7, 8.0), (8.0, 8.3), PLUS),      # chain C
    ((5.2, 6.7), (5.8, 6.1), MINUS), ((6.0, 5.9), (6.5, 5.4), MINUS),    # B, inverted
]
line(X(4.0), Y(4.3), X(4.4), Y(4.7), stroke=PLUS, width=4)
line(X(5.05), Y(6.85), X(6.75), Y(5.15), stroke=MINUS, width=4)
line(X(7.5), Y(7.8), X(7.7), Y(8.0), stroke=PLUS, width=4)
for (r0, f0), (r1, f1), colour in seeds:
    line(X(r0), Y(f0), X(r1), Y(f1), stroke=colour, width=8)

text(X(3.75), Y(4.55), "chain A", 17, weight="bold", fill=PLUS)
text(X(7.0), Y(8.05), "chain C", 17, weight="bold", fill=PLUS)
lines(X(3.7), Y(6.55), ["chain B, aligned on", "the anti-diagonal"], 15, fill=MINUS,
      weight="bold")
lines(X(6.2), Y(4.95), ["bridge: A joined to C", "through the inverted", "sequence"], 15,
      fill=BRIDGE, weight="bold")

# ── Records along the read
ry = 806


def bar(y, r0, r1, colour, label, ticks=None):
    rect(X(r0), y, X(r1) - X(r0), 22, colour, rx="3")
    if ticks:
        for k in range(ticks[2]):
            x = X(ticks[0]) + (X(ticks[1]) - X(ticks[0])) * (k + 0.5) / ticks[2]
            line(x, y + 2, x, y + 20, stroke="#ffffff", width=1.5, cap="butt")
    text(X(3.6) - 12, y + 16, label, 15, anchor="end", weight="bold", fill=MUTED)


# The inverted bases, explained twice in the conventional output.
rect(X(5.05), ry - 8, X(6.75) - X(5.05), 22 + 30 + 16, BRIDGE, fill_opacity="0.12")
line(X(5.05), ry - 8, X(5.05), ry + 60, stroke=BRIDGE, width=1, dash="3 3", cap="butt")
line(X(6.75), ry - 8, X(6.75), ry + 60, stroke=BRIDGE, width=1, dash="3 3", cap="butt")
rect(X(5.0), ry, X(6.9) - X(5.0), 22, BRIDGE, rx="0")
bar(ry, 3.6, 5.0, PLUS, "conventional")
bar(ry, 6.9, 8.0, PLUS, "")
for k in range(24):
    x = X(5.0) + (X(6.9) - X(5.0)) * (k + 0.5) / 24
    line(x, ry + 3, x, ry + 19, stroke="#ffffff", width=1.5, cap="butt")
rect(X(5.05), ry + 30, X(6.75) - X(5.05), 22, MINUS, rx="3")
text((X(5.05) + X(6.75)) / 2, ry + 46, "chr7 −", 14, anchor="middle", fill="#ffffff",
     weight="bold")
lines(X(6.75) + 12, ry + 34, ["secondary +", "supplementary"], 14, fill=MINUS, weight="bold",
      step=16)
lines((X(3.6) + X(8.0)) / 2, ry + 82,
      ["the bridged chr7 record, with 1.9 kb of spurious indels and mismatches,",
       "and B as well: the inverted bases are explained twice, in contradiction"],
      15, anchor="middle", fill=BRIDGE)

ry2 = ry + 136
bar(ry2, 3.6, 5.0, PLUS, "Parallax")
rect(X(5.0), ry2, X(6.9) - X(5.0), 22, MINUS, rx="3")
rect(X(6.9), ry2, X(8.0) - X(6.9), 22, PLUS, rx="3")
text((X(3.6) + X(5.0)) / 2, ry2 + 16, "chr7 +", 14, anchor="middle", fill="#ffffff", weight="bold")
text((X(5.0) + X(6.9)) / 2, ry2 + 16, "chr7 −", 14, anchor="middle", fill="#ffffff", weight="bold")
text((X(6.9) + X(8.0)) / 2, ry2 + 16, "chr7 +", 14, anchor="middle", fill="#ffffff", weight="bold")
text((X(3.6) + X(8.0)) / 2, ry2 + 44, "three records: the inversion is reported as such", 15,
     anchor="middle", fill=MUTED)

text(X(3.6), H - 26, "Schematic. The sequences and the bridge’s counts are illustrative; costs "
     "are Parallax’s defaults.", 14, fill=MUTED, style="italic")

# ── Right panel
px, pw = 780, 770


def heading(y, num, s):
    add(f'<circle cx="{px + 15}" cy="{y - 7}" r="15" fill="{INK}"/>')
    text(px + 15, y - 1, num, 17, anchor="middle", weight="bold", fill="#ffffff")
    text(px + 42, y, s, 22, weight="bold")


heading(170, "1", "The aligner fills the gap by force")
lines(px, 206, ["Chains A and C lie on one diagonal, so they chain together and the",
                "1.9 kb between them is aligned. The read there is the inverted",
                "segment, unrelated to the forward reference, but both are built",
                "from the same short motifs, so the alignment finds a path:"],
      16, fill=MUTED)

# base-level snippet
sx, sy = px, 296
cw = 10.2
cols = len(aln_read)
rect(sx, sy, pw, 112, PANEL, rx="6")
text(sx + 16, sy + 34, "read", 14, fill=MUTED, weight="bold")
text(sx + 16, sy + 88, "ref", 14, fill=MUTED, weight="bold")
x0 = sx + 70
n_match = n_mismatch = n_gap = 0
for c in range(cols):
    a, b = aln_read[c], aln_ref[c]
    x = x0 + c * cw
    if a == "-" or b == "-":
        colour = BRIDGE
        n_gap += 1
    elif a != b:
        colour = MINUS
        n_mismatch += 1
    else:
        colour = INK
        n_match += 1
    text(x + cw / 2, sy + 34, a, 16, anchor="middle", fill=colour, font=MONO)
    text(x + cw / 2, sy + 88, b, 16, anchor="middle", fill=colour, font=MONO)
    if a == b:
        line(x + cw / 2, sy + 44, x + cw / 2, sy + 70, stroke=MUTED, width=1, cap="butt")
    elif a != "-" and b != "-":
        text(x + cw / 2, sy + 62, "·", 16, anchor="middle", fill=MINUS, font=MONO)
text(sx + pw - 14, sy + 134, f"{n_match} matches, {n_mismatch} mismatches, {n_gap} gap columns "
     "in this stretch", 14, anchor="end", fill=MUTED, style="italic")

lines(px, 470, ["Chain B is still reported as well, as a secondary + supplementary",
                "record, so the output explains the inverted bases two ways."], 16,
      fill=INK, weight="bold")

heading(546, "2", "The same gap, scored in bits")
text(px + pw - 14, 576, "cost in bits", 14, anchor="end", fill=MUTED, style="italic")
rows = [
    ("bridge: align the 1.9 kb gap", None, True),
    (f"{MATCHES:,} matched bases", f"−{2 * MATCHES:,}", False),
    (f"{MISMATCHES} mismatches", f"{8 * MISMATCHES:,}", False),
    (f"{INDELS} indels, {INDEL_BP} bp in all", f"{6 * INDELS + INDEL_BP:,}", False),
    ("total", f"{bridge_cost:,}", True),
    ("inversion: chain through B’s seeds", None, True),
    (f"{B_BP:,} bases explained by B", f"−{2 * B_BP:,}", False),
    ("two strand flips", f"≈ {flips:.0f}", False),
    ("total", f"≈ −{-inversion_cost:,.0f}", True),
]
ty = 604
for label, bits, strong in rows:
    if label == "total":
        rect(px, ty - 22, pw, 30, PANEL)
    indent = 14 if strong else 34
    text(px + indent, ty, label, 17 if strong else 16, weight="bold" if strong else "normal")
    if bits:
        text(px + pw - 170, ty, bits, 17, weight="bold" if strong else "normal",
             fill=BRIDGE if label == "total" and bits.lstrip("≈ ")[0] != "−" else INK)
    ty += 30 if not (label == "total") else 40
margin = bridge_cost - inversion_cost
lines(px, ty + 2, [f"The inversion costs about {round(margin, -2):,.0f} bits less, so Parallax reports it."],
      16, fill=INK, weight="bold")

heading(ty + 64, "3", "Parallax can bridge too, but must earn it")
lines(px, ty + 100, ["Parallax aligns seedless gaps of 20–2,000 bp and would try this one.",
                     "It keeps the bridge only if that costs less than the alternatives,",
                     "on the same scale. With no seeds in the inverted piece at all,",
                     "a bridge could still win."], 16, fill=MUTED)

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out) + "\n")
