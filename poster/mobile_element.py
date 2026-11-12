"""Generate the mobile-element-insertion poster figure as a standalone SVG.

Regenerate with: python3 poster/mobile_element.py poster/mobile-element.svg
"""
import math
import sys

W, H = 1600, 980
FONT = "Helvetica, Arial, sans-serif"

FLANK = "#1f6fb2"    # chr7 records (the flanks)
INSERT = "#138d75"   # the inserted SVA_E, placed on its source copy
GREY = "#9aa0a6"     # other SVA_E copies
INK = "#222222"
MUTED = "#5f6368"
EDGE = "#333333"
REJECT = "#c0392b"
PANEL = "#f4f6f9"
ELEMENT = "#dff3ee"  # the source SVA_E copy on chr10

L = 2.0              # insert length, kb


def X(r):            # read position (kb) -> x
    return 130 + 62 * r


def Y10(f):          # chr10 reference position (kb) -> y
    return 330 - 30 * f


def Y7(f):           # chr7 reference position (kb) -> y
    return 710 - 30 * f


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


# Read layout (kb): left flank 0-4, SVA_E insert 4-6, right flank 6-10.
# chr7: left flank at 1-5, right flank at 5-9 (insertion site at 5).
# chr10: the source SVA_E copy at 1-3.
left = [(0.3, 0.8), (1.2, 1.7), (2.1, 2.6), (3.1, 3.7)]
insert = [(4.15, 4.6), (4.85, 5.3), (5.45, 5.85)]
right = [(6.2, 6.7), (7.1, 7.6), (8.0, 8.5), (9.0, 9.6)]


def p_left(r):
    return X(r), Y7(r + 1)


def p_insert(r):
    return X(r), Y10(r - 3)


def p_right(r):
    return X(r), Y7(r - 1)


add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" '
    f'role="img" aria-label="A mobile element insertion in Parallax output: the chain explains '
    f'the inserted SVA_E sequence by jumping to a matching copy on another chromosome and back, '
    f'which outscores skipping the insert, so the read is reported as three records: the two '
    f'chr7 flanks, which meet at the insertion site, and the insert placed on its source copy.">')
add('<defs><marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="6" '
    f'markerHeight="6" orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="{EDGE}"/>'
    '</marker></defs>')
rect(0, 0, W, H, "#ffffff")

# ── Title
text(60, 66, "A mobile element insertion: the insert placed on its source copy", 34,
     weight="bold")
text(60, 100, "Seeds in an inserted SVA_E match a copy on another chromosome. Jumping there and "
     "back costs far less than", 19, fill=MUTED)
text(60, 124, "the inserted bases earn, so the read is reported as three records: two chr7 "
     "flanks and the insert.", 19, fill=MUTED)

# ── Dot plot bands
right_edge = X(10)
rect(X(0), Y10(4), right_edge - X(0), Y10(0) - Y10(4), "#f6f8fb")
rect(X(0), Y10(3), right_edge - X(0), Y10(1) - Y10(3), ELEMENT)
text(X(0) + 10, Y10(3) + 20, "SVA_E copy", 15, fill=INSERT, weight="bold")
grey_top, grey_bottom = 352, 398
rect(X(0), grey_top, right_edge - X(0), grey_bottom - grey_top, "#f5f5f5")
rect(X(0), Y7(10), right_edge - X(0), Y7(0) - Y7(10), "#f6f8fb")

text(X(0) - 12, (Y10(0) + Y10(4)) / 2 + 6, "chr10", 18, anchor="end", weight="bold", fill=MUTED)
lines(X(0) - 12, (grey_top + grey_bottom) / 2 - 3, ["other", "copies"], 14, anchor="end",
      fill=MUTED, step=16)
text(X(0) - 12, (Y7(0) + Y7(10)) / 2 + 6, "chr7", 18, anchor="end", weight="bold", fill=MUTED)
add(f'<text transform="translate(30 {(Y10(4) + Y7(0)) / 2:.1f}) rotate(-90)" font-family="{FONT}" '
    f'font-size="18" fill="{MUTED}" text-anchor="middle">reference position</text>')

# insertion site on chr7
line(X(0), Y7(5), right_edge, Y7(5), stroke=MUTED, width=1.5, dash="2 5")
text(right_edge - 8, Y7(5) + 20, "insertion site", 15, anchor="end", fill=MUTED, style="italic")

# read axis
ay = Y7(0) + 12
line(X(0), ay, X(10), ay, stroke=MUTED, width=1.5, cap="butt")
for kb in range(0, 11, 2):
    line(X(kb), ay, X(kb), ay + 7, stroke=MUTED, width=1.5, cap="butt")
    text(X(kb), ay + 27, f"{kb} kb" if kb in (0, 10) else str(kb), 16, anchor="middle", fill=MUTED)
text((X(0) + X(10)) / 2, ay + 52, "read position", 18, anchor="middle", fill=MUTED)

# read composition strip
sy = ay + 70
for name, lo, hi, fill in [("chr7 flank", 0, 4, "#e8edf3"), ("SVA_E", 4, 6, ELEMENT),
                           ("chr7 flank", 6, 10, "#e8edf3")]:
    rect(X(lo) + 1, sy, X(hi) - X(lo) - 2, 24, fill, stroke="#c9ced6", stroke_width="1")
    text((X(lo) + X(hi)) / 2, sy + 18, name, 15, anchor="middle", weight="bold", fill=MUTED)
text(X(0) - 12, sy + 18, "read", 16, anchor="end", weight="bold", fill=MUTED)

# ── The rejected alternative: skip the insert, joining the flanks directly
a = p_left(left[-1][1])
b = p_right(right[0][0])
line(a[0] + 6, a[1] + 4, b[0] - 6, b[1] + 4, stroke=REJECT, width=2.5, dash="7 6")
mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + 4
line(mx - 8, my - 8, mx + 8, my + 8, stroke=REJECT, width=3)
line(mx - 8, my + 8, mx + 8, my - 8, stroke=REJECT, width=3)
skip_bits = 35 + math.log2(1 + L * 1000)
lines(mx - 70, my + 34, [f"skip the insert: −{skip_bits:.0f} bits,",
                         f"and its {L * 1000:,.0f} bases earn nothing"], 15, fill=REJECT)

# ── Other SVA_E copies: shorter matches, not chosen
for (r0, r1), dy in zip(insert, (0, 18, 6)):
    for k, frac in enumerate((0.55, 0.35)):
        x0 = X(r0) + k * 6
        x1 = x0 + (X(r1) - X(r0)) * frac
        y0 = grey_bottom - 10 - dy - k * 12
        line(x0, y0, x1, y0 - (x1 - x0) * 0.4, stroke=GREY, width=6)

# ── Chosen moves
def chain(seeds, at):
    for (a0, a1), (b0, b1) in zip(seeds, seeds[1:]):
        p, q = at(a1), at(b0)
        line(p[0] + 4, p[1] - 2, q[0] - 4, q[1] + 2, width=2, marker=True)


chain(left, p_left)
chain(insert, p_insert)
chain(right, p_right)

out_of = p_left(left[-1][1])
into = p_insert(insert[0][0])
line(out_of[0] + 3, out_of[1] - 6, into[0] - 2, into[1] + 8, width=2.5, marker=True)
back = p_insert(insert[-1][1])
home = p_right(right[0][0])
line(back[0] + 3, back[1] + 6, home[0] - 2, home[1] - 8, width=2.5, marker=True)
lines(out_of[0] - 186, (out_of[1] + into[1]) / 2 + 30,
      ["other chromosome: −66.5", "second locus: −40"], 15)
lines(home[0] + 14, (back[1] + home[1]) / 2 + 30, ["other chromosome:", "−66.5"], 15)
lines(X(6.1), Y10(3.3), [f"inserted bases earn +2 each:", f"≈ +{2 * L * 1000:,.0f} bits"], 15,
      fill=INSERT)

# seeds
for seeds, at, colour in [(left, p_left, FLANK), (insert, p_insert, INSERT),
                          (right, p_right, FLANK)]:
    for r0, r1 in seeds:
        p, q = at(r0), at(r1)
        line(p[0], p[1], q[0], q[1], stroke=colour, width=8)

# ── Records, in read coordinates
ry = sy + 56
text(X(0) - 12, ry + 17, "records", 16, anchor="end", weight="bold", fill=MUTED)
for lo, hi, colour, label in [(0, 4, FLANK, "1: chr7 primary"), (4, 6, INSERT, "2: chr10"),
                              (6, 10, FLANK, "3: chr7 supplementary")]:
    rect(X(lo) + 1, ry, X(hi) - X(lo) - 2, 24, colour, rx="3")
    text((X(lo) + X(hi)) / 2, ry + 17, label, 14, anchor="middle", fill="#ffffff", weight="bold")
text((X(0) + X(10)) / 2, ry + 50, "adjacent on the read, linked by SA tags", 15,
     anchor="middle", fill=MUTED, style="italic")

# ── Records, in reference coordinates (beside the plot, sharing its y axes)
bx = right_edge + 24
for num, lo, hi in [(1, 1, 5), (3, 5, 9)]:
    rect(bx, Y7(hi) + 1, 22, Y7(lo) - Y7(hi) - 2, FLANK, rx="3")
    text(bx + 34, (Y7(lo) + Y7(hi)) / 2 + 5, str(num), 15, weight="bold", fill=FLANK)
rect(bx, Y10(3), 22, Y10(1) - Y10(3), INSERT, rx="3")
text(bx + 34, (Y10(1) + Y10(3)) / 2 + 5, "2", 15, weight="bold", fill=INSERT)
line(bx - 6, Y7(5), bx + 28, Y7(5), stroke=INK, width=2.5, cap="butt")
lines(bx + 50, Y7(5) - 4, ["1 and 3 meet", "at the site"], 14, weight="bold")
lines(bx + 11, Y10(4) - 26, ["records on", "the reference"], 14, anchor="middle", fill=MUTED,
      style="italic")

# ── Right panel
px, pw = 1010, 540


def heading(y, num, s):
    add(f'<circle cx="{px + 15}" cy="{y - 7}" r="15" fill="{INK}"/>')
    text(px + 15, y - 1, num, 17, anchor="middle", weight="bold", fill="#ffffff")
    text(px + 42, y, s, 22, weight="bold")


heading(170, "1", "The insert matches SVA_E copies")
lines(px, 206, ["Seeds in the inserted 2 kb hit SVA_E copies across the",
                "genome. The copy on chr10 explains the most read bases,",
                "so the chain uses it; shorter matches elsewhere lose."], 16, fill=MUTED)

heading(312, "2", "Explaining the insert beats skipping it")
rows = [
    ("jump to chr10 and back", "2 × 66.5"),
    ("second locus", "40"),
    ("2 kb of insert explained", "−2 × 2,000"),
    ("explain via chr10: total", f"≈ −{2 * L * 1000 - 133 - 40:,.0f}"),
    ("skip the insert: diagonal shift", f"≈ {skip_bits:.0f}"),
]
text(px + pw - 14, 344, "cost in bits", 14, anchor="end", fill=MUTED, style="italic")
ty = 372
for i, (label, bits) in enumerate(rows):
    if i in (3, 4):
        add(f'<rect x="{px}" y="{ty - 24}" width="{pw}" height="34" fill="{PANEL}"/>')
    weight = "bold" if i in (3, 4) else "normal"
    text(px + 14 + (0 if i >= 3 else 20), ty, label, 17 if i >= 3 else 16, weight=weight)
    text(px + pw - 150, ty, bits, 17, weight="bold")
    ty += 34 if i >= 2 else 28
lines(px, ty + 4, ["Explaining wins once about 64 inserted bases match:",
                   "2m − 173 > −46."], 16, fill=MUTED)

heading(ty + 92, "3", "Each chromosome change starts a record")
lines(px, ty + 128, ["Records 1 and 3 hold the chr7 flanks, which meet at the",
                     "insertion site. Record 2 places the inserted sequence on",
                     "its source copy, naming the element that was inserted."],
      16, fill=MUTED)

heading(ty + 234, "4", "Weak matches stay in one record")
lines(px, ty + 270, ["If too little of the insert matches, skipping wins: the",
                     "flanks are joined in one chr7 record, with the insert",
                     "reported as an insertion."], 16, fill=MUTED)

text(X(0), H - 24, "Schematic, not to scale; chromosomes and coordinates are illustrative. "
     "Default parameters; bits = −log₂ P.", 14, fill=MUTED, style="italic")

add("</svg>")
open(sys.argv[1], "w").write("\n".join(out) + "\n")
