#!/usr/bin/env python3
"""Gera diagrama esquemático circular do genoma mitocondrial de vertebrados."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np

# ---------- dados ----------
genes = [
    ("D-loop",    0,    1000, "H", "dloop"),
    ("12S rRNA",  1000, 1970, "H", "rrna"),
    ("16S rRNA",  2000, 3600, "H", "rrna"),
    ("ND1",       3700, 4700, "H", "ci"),
    ("ND2",       4900, 5950, "H", "ci"),
    ("COX1",      6500, 8050, "H", "civ"),
    ("COX2",      8050, 8750, "H", "civ"),
    ("ATP8",      8750, 9000, "H", "cv"),
    ("ATP6",      9000, 9700, "H", "cv"),
    ("COX3",      9700, 10500, "H", "civ"),
    ("ND3",       10500, 10850, "H", "ci"),
    ("ND4L",      10850, 11150, "H", "ci"),
    ("ND4",       11150, 12500, "H", "ci"),
    ("ND5",       12500, 14350, "H", "ci"),
    ("ND6",       14350, 14900, "L", "ci"),
    ("CYTB",      14900, 16050, "H", "ciii"),
]

trna_positions = [
    (1970, 2000), (3600, 3700), (4700, 4900),
    (5950, 6050), (6050, 6150), (6150, 6250),
    (6250, 6350), (6350, 6450), (6450, 6500),
    (8700, 8750), (10450, 10500), (11100, 11150),
    (12450, 12500), (14300, 14350), (14850, 14900),
    (16050, 16150), (16150, 16250), (16250, 16350),
    (16350, 16450), (16450, 16500),
]

genome_size = 16500

colors = {
    "ci":    "#E74C3C",
    "ciii":  "#8E44AD",
    "civ":   "#3498DB",
    "cv":    "#F39C12",
    "rrna":  "#27AE60",
    "trna":  "#BDC3C7",
    "dloop": "#95A5A6",
}


def bp_to_deg(bp):
    """0 bp = topo (90°), sentido horário."""
    return 90 - (bp / genome_size) * 360


def draw_arc_cart(ax, start_bp, end_bp, r_inner, r_outer, color):
    a_start = np.deg2rad(bp_to_deg(start_bp))
    a_end   = np.deg2rad(bp_to_deg(end_bp))
    angles = np.linspace(a_start, a_end, 80)
    x_out = r_outer * np.cos(angles)
    y_out = r_outer * np.sin(angles)
    x_in  = r_inner * np.cos(angles[::-1])
    y_in  = r_inner * np.sin(angles[::-1])
    xs = np.concatenate([x_out, x_in, [x_out[0]]])
    ys = np.concatenate([y_out, y_in, [y_out[0]]])
    verts = list(zip(xs, ys))
    codes = [1] + [2] * (len(verts) - 2) + [79]
    patch = mpatches.PathPatch(Path(verts, codes),
                               facecolor=color, edgecolor='white',
                               linewidth=1.2, zorder=2)
    ax.add_patch(patch)


# ---------- figura ----------
fig, ax = plt.subplots(figsize=(12, 14))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.set_aspect('equal')
ax.axis('off')

ring_inner = 3.5
ring_outer = 5.0
label_r = 5.6

# tRNAs
for s, e in trna_positions:
    draw_arc_cart(ax, s, e, ring_inner, ring_outer, colors["trna"])

# Genes
for name, s, e, strand, cplx in genes:
    draw_arc_cart(ax, s, e, ring_inner, ring_outer, colors[cplx])

# ---------- labels radiais ----------
for name, s, e, strand, cplx in genes:
    mid_bp = (s + e) / 2
    angle_deg = bp_to_deg(mid_bp)
    angle_rad = np.deg2rad(angle_deg)

    lx = label_r * np.cos(angle_rad)
    ly = label_r * np.sin(angle_rad)

    # Normalizar ângulo para [0, 360)
    a_norm = angle_deg % 360

    # Orientação radial: texto aponta para fora do centro
    # No lado esquerdo (90° < a < 270°), inverter 180° para não ficar de cabeça para baixo
    if 90 < a_norm < 270:
        rot = angle_deg - 180
        ha = 'right'
    else:
        rot = angle_deg
        ha = 'left'

    fontsize = 12
    if name in ("ATP8", "ND3", "ND4L"):
        fontsize = 10

    ax.text(lx, ly, name,
            ha=ha, va='center',
            fontsize=fontsize, fontweight='bold',
            rotation=rot, rotation_mode='anchor',
            color='#2C3E50', zorder=5)

# ---------- texto central ----------
ax.text(0, 0.6, "Genoma\nMitocondrial",
        ha='center', va='center',
        fontsize=22, fontweight='bold', color='#2C3E50')
ax.text(0, -0.5, "Vertebrado típico\n~16.500 bp",
        ha='center', va='center', fontsize=14, color='#555555')
ax.text(0, -1.4, "37 genes: 13 CDS, 22 tRNA, 2 rRNA\n+ região controle (D-loop)",
        ha='center', va='center', fontsize=11, color='#777777')
ax.text(0, -2.2, "Fita H (pesada, exterior) — maioria dos genes\nFita L (leve, interior) — ND6 + 8 tRNAs",
        ha='center', va='center', fontsize=10, color='#999999')

# ---------- título ----------
fig.text(0.5, 0.96,
         "Organização típica do genoma mitocondrial de vertebrados",
         ha='center', va='center', fontsize=18, fontweight='bold', color='#2C3E50')

# ---------- legenda na faixa inferior ----------
legend_items = [
    ("Complexo I (ND1–ND6, ND4L)", colors["ci"]),
    ("Complexo III (CYTB)",         colors["ciii"]),
    ("Complexo IV (COX1–COX3)",     colors["civ"]),
    ("Complexo V (ATP6, ATP8)",     colors["cv"]),
    ("rRNA (12S, 16S)",             colors["rrna"]),
    ("tRNA (22 genes)",             colors["trna"]),
    ("D-loop (região controle)",    colors["dloop"]),
]

handles = [mpatches.Patch(facecolor=c, edgecolor='#666666', linewidth=0.5, label=l)
           for l, c in legend_items]

fig.legend(handles=handles, loc='lower center',
           bbox_to_anchor=(0.5, 0.01),
           ncol=3, fontsize=11, frameon=False,
           handlelength=1.5, handleheight=1.2,
           labelspacing=0.6, columnspacing=1.5)

# limites
m = 1.8
ax.set_xlim(-ring_outer - m, ring_outer + m)
ax.set_ylim(-ring_outer - m - 0.5, ring_outer + m + 0.5)

plt.savefig('/home/matheus/mitogenome-pipeline/latex/imagens/mitogenoma_vertebrado_esquematico.png',
            dpi=200, bbox_inches='tight', facecolor='white', pad_inches=0.3)
print("OK")
