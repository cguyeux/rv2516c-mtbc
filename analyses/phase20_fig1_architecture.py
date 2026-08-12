#!/usr/bin/env python3
"""phase20_fig1_architecture.py -- Figure 1 du manuscrit : architecture à trois domaines et gain
apporté par l'interrogation par domaine (P7.6.1 / P2.1).

Panneau A : diagramme linéaire des 267 résidus, trois unités structurales colorées (contact-map,
P7.6.1) avec la fenêtre 111-121 (charge/conservation, P2.3/P3.1) repérée dans l'unité B.
Panneau B : gain d'E-value HHpred en isolant l'unité B (P2.1) -- barres log, protéine entière vs
unité B seule, avec le facteur de gain annoté.

Chiffres repris de main.tex (vérifiés par recalcul direct depuis les fichiers de résultats bruts) :
E=5.9e-7 (protéine entière), E=1.3e-10 (unité B seule).

Sorties : article/figures/{figure1_architecture.pdf,.png}
Run: python analyses/phase20_fig1_architecture.py

Dépend d'une bibliothèque de style de figures externe à ce dépôt (module `figstyle`, presets
journal) ; définir SCI_FIGURE_SCRIPTS pour pointer vers son répertoire si elle n'est pas au
même endroit que sur la machine d'origine.
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, os.environ.get(
    "SCI_FIGURE_SCRIPTS",
    str(Path.home() / "docs/codes/claude_plugins/bio_population_genetics/skills/sci-figure/scripts"),
))
import figstyle as fs
import matplotlib.patches as mpatches

OUT = ROOT / "article" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

PROT_LEN = 267
UNITS = [
    ("Unit A\nferredoxin-like\n(DUF8830)", 1, 87, fs.PALETTE_CATEGORICAL[0]),
    ("Unit B\nwHTH AlpA/excisionase", 88, 147, fs.PALETTE_CATEGORICAL[1]),
    ("Linker", 148, 177, "#cccccc"),
    ("Unit C\nIg-like $\\beta$-sandwich", 178, 267, fs.PALETTE_CATEGORICAL[2]),
]
RECOGNITION_WINDOW = (111, 121)  # tour + helice de reconnaissance, P2.3/P3.1

EVALUES = [
    ("Full-length\nprotein", 5.9e-7, "grey"),
    ("Unit B\nalone", 1.3e-10, fs.PALETTE_CATEGORICAL[1]),
]


def panel_a(ax):
    ax.set_xlim(0, PROT_LEN)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    bar_h = 0.5
    y0 = 0.25
    for label, start, end, color in UNITS:
        ax.add_patch(mpatches.Rectangle((start, y0), end - start, bar_h,
                                         facecolor=color, edgecolor="black", linewidth=0.8))
        ax.text((start + end) / 2, y0 + bar_h / 2, label, ha="center", va="center",
                fontsize=6.5, color="black" if color == "#cccccc" else "white")
    # recognition window bracket, above the domain bar
    rs, re_ = RECOGNITION_WINDOW
    ax.annotate("", xy=(re_, y0 + bar_h + 0.22), xytext=(rs, y0 + bar_h + 0.22),
                arrowprops=dict(arrowstyle="-", color="black", lw=1.2))
    ax.plot([rs, rs], [y0 + bar_h, y0 + bar_h + 0.22], color="black", lw=1.0)
    ax.plot([re_, re_], [y0 + bar_h, y0 + bar_h + 0.22], color="black", lw=1.0)
    ax.text((rs + re_) / 2, y0 + bar_h + 0.30, f"candidate DNA-binding\nwindow ({rs}--{re_})",
            ha="center", va="bottom", fontsize=6)
    # Bornes adjacentes (87/88, 147/148, 177/178) fusionnées en UN seul repère label combiné,
    # sinon les deux nombres, tournés à 90 degres, se recouvrent visuellement (un seul pixel d'ecart).
    boundary_labels = [(1, "1"), (87.5, "87/88"), (147.5, "147/148"),
                        (177.5, "177/178"), (267, "267")]
    for x, label in boundary_labels:
        ax.plot([x, x], [y0 - 0.06, y0], color="black", lw=0.6)
        ax.text(x, y0 - 0.10, label, ha="center", va="top", fontsize=6, rotation=90)
    ax.set_xlabel("Residue position")
    ax.spines[["top", "right", "left"]].set_visible(False)


def panel_b(ax):
    labels = [e[0] for e in EVALUES]
    values = [e[1] for e in EVALUES]
    colors = [e[2] for e in EVALUES]
    xpos = [0, 1]
    ax.set_yscale("log")
    ax.plot(xpos, values, color="black", lw=1.2, zorder=1, linestyle="--")
    ax.scatter(xpos, values, s=140, color=colors, edgecolor="black", linewidth=1.0, zorder=2)
    for x, val in zip(xpos, values):
        ax.annotate(f"$E={val:.1e}$".replace("e-0", "e-"), (x, val),
                    textcoords="offset points", xytext=(0, 10), ha="center", fontsize=7)
    ax.set_xticks(xpos, labels)
    ax.set_xlim(-0.4, 1.4)
    ax.set_ylim(max(values) * 20, min(values) / 20)  # inverted: "stronger" (lower E) goes up
    ax.set_ylabel("HHpred $E$-value (log scale; higher on plot = stronger)")
    mid_y = (values[0] * values[1]) ** 0.5
    ax.text(0.5, mid_y, ">4,500-fold\ngain", ha="center", va="center", fontsize=7.5,
            style="italic", rotation=-28)
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    fig, axes = fs.panel_grid("nature_double", nrows=1, ncols=2, width_ratios=[1.6, 1])
    panel_a(axes[0])
    panel_b(axes[1])
    fs.panel_labels(axes)
    fig.suptitle("")
    out = fs.save(fig, OUT / "figure1_architecture", "nature_double", raster=True)
    print("Ecrit:", *out, sep="\n  ")


if __name__ == "__main__":
    main()
