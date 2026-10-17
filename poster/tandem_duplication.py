"""Generate the tandem-duplication poster figure as a standalone SVG.

Regenerate with: python3 poster/tandem_duplication.py poster/tandem-duplication.svg
"""
import math
import sys

W, H = 1600, 960
FONT = "Helvetica, Arial, sans-serif"

REC1 = "#1f6fb2"     # record 1: A + first copy of B
REC2 = "#7b4fa0"     # record 2: second copy of B + C
UNIT = "#fff1cc"     # the duplicated unit B
INK = "#222222"
MUTED = "#5f6368"
EDGE = "#333333"
PANEL = "#f4f6f9"

D = 6  # unit length, kb


def X(r):            # read position (kb) -> x
    return 130 + 41 * r


def Y(f):            # chr7 reference position (kb) -> y
    return 700 - 42 * f


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


def rect(x, y, w, h, fill, **kw):
    extra = "".join(f' {k.replace("_", "-")}="{v}"' for k, v in kw.items())
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{fill}"{extra}/>')


# Read layout (kb): A = 0-2, B = 2-8, B again = 8-14, C = 14-16.
# Reference (kb):   A = 1-3, B = 3-9, C = 9-11.
copy1 = [(0.3, 0.9), (1.3, 1.9), (2.5, 3.1), (3.7, 4.3), (4.9, 5.5), (6.1, 6.7), (7.2, 7.8)]
copy2 = [(8.2, 8.8), (9.3, 9.9), (10.5, 11.1), (11.7, 12.3), (12.9, 13.5), (14.1, 14.7),
         (15.1, 15.7)]


def ref1(r):         # first pass: A then B
    return r + 1


def ref2(r):         # second pass: B again then C
    return r - 5


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'role="img" aria-label="A tandem duplication in Parallax output: the read passes through the '
    f'duplicated unit twice, the chain joins the passes with one backward reference jump, and the '
    f'read is reported as two records adjacent on the read that overlap on the reference by the '
    f'unit length.">')
add('<defs><marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" '
    f'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{EDGE}"/>'
    '</marker></defs>')
rect(0, 0, W, H, "#ffffff")

# ── Title
text(60, 66, "A tandem duplication: two records that overlap on the reference", 34, weight="bold")
text(60, 100, "The read passes through the duplicated unit twice. The chain joins the passes with "
     "one backward jump,", 19, fill=MUTED)
text(60, 124, "and the jump starts a new record, so the two records overlap on the reference "
     "by the unit length.", 19, fill=MUTED)

# ── Dot plot
top, bottom = Y(12), Y(0)
rect(X(0), top, X(16) - X(0), bottom - top, "#f6f8fb")
# the duplicated unit B on the reference
rect(X(0), Y(9), X(16) - X(0), Y(3) - Y(9), UNIT)
text(X(0) + 12, Y(9) + 24, f"duplicated unit B ({D} kb)", 16, fill="#8a6d1f", weight="bold")

# reference labels
add(f'<text transform="translate(44 {(top + bottom) / 2:.1f}) rotate(-90)" font-family="{FONT}" '
    f'font-size="18" fill="{MUTED}" text-anchor="middle">chr7 reference position</text>')
for name, lo, hi in [("A", 1, 3), ("B", 3, 9), ("C", 9, 11)]:
    line(X(0) - 18, Y(lo) - 3, X(0) - 18, Y(hi) + 3, stroke=MUTED, width=2, cap="butt")
    text(X(0) - 28, (Y(lo) + Y(hi)) / 2 + 6, name, 17, anchor="end", weight="bold", fill=MUTED)

# read axis
ay = bottom + 12
line(X(0), ay, X(16), ay, stroke=MUTED, width=1.5, cap="butt")
for kb in range(0, 17, 4):
    line(X(kb), ay, X(kb), ay + 7, stroke=MUTED, width=1.5, cap="butt")
    text(X(kb), ay + 27, f"{kb} kb" if kb in (0, 16) else str(kb), 16, anchor="middle", fill=MUTED)
text((X(0) + X(16)) / 2, ay + 52, "read position", 18, anchor="middle", fill=MUTED)

# read composition strip
sy = ay + 70
for name, lo, hi, fill in [("A", 0, 2, "#e8edf3"), ("B", 2, 8, UNIT), ("B", 8, 14, UNIT),
                           ("C", 14, 16, "#e8edf3")]:
    rect(X(lo) + 1, sy, X(hi) - X(lo) - 2, 24, fill, stroke="#c9ced6", stroke_width="1")
    text((X(lo) + X(hi)) / 2, sy + 18, name, 16, anchor="middle", weight="bold", fill=MUTED)
text(X(0) - 12, sy + 18, "read", 16, anchor="end", weight="bold", fill=MUTED)

# moves along each pass
for seeds, ref in [(copy1, ref1), (copy2, ref2)]:
    for (a0, a1), (b0, b1) in zip(seeds, seeds[1:]):
        line(X(a1) + 4, Y(ref(a1)) - 3, X(b0) - 4, Y(ref(b0)) + 3, width=2, marker=True)

# the backward jump between the passes
j0 = (X(copy1[-1][1]), Y(ref1(copy1[-1][1])))
j1 = (X(copy2[0][0]), Y(ref2(copy2[0][0])))
line(j0[0] + 3, j0[1] + 6, j1[0] - 2, j1[1] - 8, width=2.5, marker=True)
jump_bits = 35 + math.log2(1 + D * 1000)
# Left of the jump, below the first pass (left-aligned: subscripts misplace
# in right-aligned text in some renderers).
lines(j0[0] - 150, j1[1] - 30,
      ["backward jump of D:", f"−(35 + log₂(1 + D))", f"= −{jump_bits:.1f} bits"], 15)
# the second pass repays it
lines(X(11.6), Y(2.2), ["second copy earns", f"+2 per base: +{2 * D * 1000:,} bits"], 15,
      fill=REC2)

# seeds
for seeds, ref, colour in [(copy1, ref1, REC1), (copy2, ref2, REC2)]:
    for a, b in seeds:
        line(X(a), Y(ref(a)), X(b), Y(ref(b)), stroke=colour, width=8)

# ── Records, in read coordinates
ry = sy + 56
text(X(0) - 12, ry + 17, "records", 16, anchor="end", weight="bold", fill=MUTED)
for lo, hi, colour, label in [(0, 8, REC1, "record 1: primary"),
                              (8, 16, REC2, "record 2: supplementary")]:
    rect(X(lo) + 1, ry, X(hi) - X(lo) - 2, 24, colour, rx="3")
    text((X(lo) + X(hi)) / 2, ry + 17, label, 15, anchor="middle", fill="#ffffff", weight="bold")
text((X(0) + X(16)) / 2, ry + 50, "adjacent on the read: no read base is reported twice", 15,
     anchor="middle", fill=MUTED, style="italic")

# ── Records, in reference coordinates (beside the plot, sharing its y axis)
bx = X(16) + 22
rect(bx - 6, Y(9), 64, Y(3) - Y(9), UNIT)
for i, (lo, hi, colour) in enumerate([(1, 9, REC1), (3, 11, REC2)]):
    x = bx + i * 30
    rect(x, Y(hi), 22, Y(lo) - Y(hi), colour, rx="3")
    text(x + 11, Y(hi) - 8, str(i + 1), 15, anchor="middle", weight="bold", fill=colour)
lines(bx + 26, Y(1) + 24, ["records on", "the reference"], 14, anchor="middle", fill=MUTED,
      style="italic")
# overlap bracket
ox = bx + 62
line(ox, Y(9), ox + 8, Y(9), stroke=INK, width=2, cap="butt")
line(ox + 8, Y(9), ox + 8, Y(3), stroke=INK, width=2, cap="butt")
line(ox, Y(3), ox + 8, Y(3), stroke=INK, width=2, cap="butt")
lines(ox + 16, (Y(3) + Y(9)) / 2 - 4, ["overlap", f"= D = {D} kb"], 15, weight="bold")

# ── Right panel
px, pw = 1010, 540


def heading(y, num, s):
    add(f'<circle cx="{px + 15}" cy="{y - 7}" r="15" fill="{INK}"/>')
    text(px + 15, y - 1, num, 17, anchor="middle", weight="bold", fill="#ffffff")
    text(px + 42, y, s, 22, weight="bold")


heading(170, "1", "The read passes through B twice")
lines(px, 206, ["Its second copy of B maps to the same reference as the",
                "first, so those seeds sit on a parallel diagonal, shifted",
                "down by the unit length D."], 16, fill=MUTED)

heading(316, "2", "One backward jump joins the passes")
lines(px, 352, ["Same chromosome and strand, diagonal shift D: the jump",
                "costs min(6 + D, 35 + log₂(1 + D)) bits. The second",
                "copy’s bases earn 2 bits each, so the chain takes the jump",
                "rather than leave them unexplained."], 16, fill=MUTED)

heading(484, "3", "The jump starts a new record")
lines(px, 520, ["It steps back further than the next seed’s length, so it is",
                "a repeat edge, not a continuation. Record 1 (A + B) and",
                "record 2 (B + C) meet on the read and overlap on the",
                "reference by exactly the duplicated unit."], 16, fill=MUTED)

heading(652, "4", "Shorter copies stay in one record")
lines(px, 688, ["If the copy steps back at most 5 kb and its flanks can be",
                "joined, its seeds are dropped and the gap aligner reports",
                "the second copy as a D-bp insertion:"], 16, fill=MUTED)
iy = 768
rect(px, iy, 250, 24, REC1, rx="3")
rect(px + 250, iy - 6, 4, 36, "#c0392b")
text(px + 252, iy + 50, "insertion of D bp", 15, anchor="middle", fill="#c0392b", weight="bold")
lines(px, iy + 84, ["one record across A, B and C"], 15, fill=MUTED, style="italic")
rect(px + 254, iy, 160, 24, REC1, rx="3")

text(X(0), H - 24, "Schematic, not to scale. Default parameters; bits = −log₂ P.", 14,
     fill=MUTED, style="italic")

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out) + "\n")
