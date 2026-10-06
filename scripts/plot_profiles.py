#!/usr/bin/env python3
"""Plot density, velocity and pressure profiles from a solver CSV file as an SVG image.

Usage:
    python3 scripts/plot_profiles.py sod_first_order.csv -o sod_first_order.svg \
        --title "Sod shock tube, first order (Rusanov), 100 cells, t = 0.2"

The CSV must have the header "x,rho,u,p". Only the Python standard library is used,
so no plotting package needs to be installed.
"""

import argparse
import csv
import math

# Colours (light surface): one data hue, text in ink tokens, recessive grid and axes.
SURFACE = "#fcfcfb"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID = "#e4e3df"
AXIS = "#b9b8b2"
SERIES = "#2a78d6"

PANEL_W, PANEL_H = 640, 180
MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 72, 24, 56, 40
GAP = 68  # vertical space between panels (x tick labels + next panel title)


def read_profiles(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {key: [float(r[key]) for r in rows] for key in ("x", "rho", "u", "p")}


def nice_ticks(lo, hi, count=5):
    """Round tick positions covering [lo, hi]."""
    if hi - lo < 1e-12:
        lo, hi = lo - 0.5, hi + 0.5
    raw = (hi - lo) / count
    magnitude = 10 ** math.floor(math.log10(raw))
    step = min((s * magnitude for s in (1, 2, 2.5, 5, 10) if s * magnitude >= raw))
    start = math.floor(lo / step) * step
    ticks = []
    t = start
    while t <= hi + step * 1e-9:
        ticks.append(round(t, 10))
        t += step
    if ticks[-1] < hi:
        ticks.append(round(ticks[-1] + step, 10))
    return ticks


def fmt(v):
    return f"{v:g}"


def panel(x, y, label, top):
    """SVG elements for one profile panel whose plot area starts at y = top."""
    left, right = MARGIN_L, MARGIN_L + PANEL_W
    bottom = top + PANEL_H
    xt = nice_ticks(min(x), max(x))
    yt = nice_ticks(min(y), max(y))
    x0, x1, y0, y1 = xt[0], xt[-1], yt[0], yt[-1]

    def sx(v):
        return left + (v - x0) / (x1 - x0) * PANEL_W

    def sy(v):
        return bottom - (v - y0) / (y1 - y0) * PANEL_H

    out = [f'<text x="{left}" y="{top - 12}" font-size="14" font-weight="600" '
           f'fill="{TEXT_PRIMARY}">{label}</text>']
    for t in yt:
        out.append(f'<line x1="{left}" x2="{right}" y1="{sy(t):.1f}" y2="{sy(t):.1f}" '
                   f'stroke="{GRID}" stroke-width="1"/>')
        out.append(f'<text x="{left - 8}" y="{sy(t) + 4:.1f}" font-size="12" text-anchor="end" '
                   f'fill="{TEXT_SECONDARY}">{fmt(t)}</text>')
    for t in xt:
        out.append(f'<text x="{sx(t):.1f}" y="{bottom + 18}" font-size="12" text-anchor="middle" '
                   f'fill="{TEXT_SECONDARY}">{fmt(t)}</text>')
    out.append(f'<line x1="{left}" x2="{right}" y1="{bottom}" y2="{bottom}" '
               f'stroke="{AXIS}" stroke-width="1"/>')
    points = " ".join(f"{sx(a):.2f},{sy(b):.2f}" for a, b in zip(x, y))
    out.append(f'<polyline points="{points}" fill="none" stroke="{SERIES}" stroke-width="2" '
               f'stroke-linejoin="round" stroke-linecap="round"/>')
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", help="solver output with columns x,rho,u,p")
    parser.add_argument("-o", "--output", default="profiles.svg", help="SVG file to write")
    parser.add_argument("--title", default="", help="title shown above the panels")
    args = parser.parse_args()

    data = read_profiles(args.csv)
    panels = [("rho", "Density ρ"), ("u", "Velocity u"), ("p", "Pressure p")]
    width = MARGIN_L + PANEL_W + MARGIN_R
    height = MARGIN_T + len(panels) * PANEL_H + (len(panels) - 1) * GAP + MARGIN_B + 20

    body = []
    if args.title:
        body.append(f'<text x="{MARGIN_L}" y="24" font-size="16" font-weight="600" '
                    f'fill="{TEXT_PRIMARY}">{args.title}</text>')
    for k, (key, label) in enumerate(panels):
        top = MARGIN_T + 12 + k * (PANEL_H + GAP)
        body += panel(data["x"], data[key], label, top)
    body.append(f'<text x="{MARGIN_L + PANEL_W / 2}" y="{height - 8}" font-size="12" '
                f'text-anchor="middle" fill="{TEXT_SECONDARY}">x</text>')

    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
           f'viewBox="0 0 {width} {height}" font-family="system-ui, sans-serif">\n'
           f'<rect width="100%" height="100%" fill="{SURFACE}"/>\n' + "\n".join(body) + "\n</svg>\n")
    with open(args.output, "w") as f:
        f.write(svg)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
