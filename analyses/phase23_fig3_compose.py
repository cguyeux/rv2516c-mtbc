#!/usr/bin/env python3
"""phase23_fig3_compose.py -- assemble les deux rendus PyMOL (phase22) en une seule figure 3
a deux panneaux (a/b), au format et a la resolution du preset nature_double, comme les figures 1
et 2 du manuscrit.

Sorties : article/figures/{figure3_site.pdf,.png}
Run: python analyses/phase23_fig3_compose.py (apres phase22_fig3_site_figure.pml)

Depend d'une bibliotheque de style de figures externe a ce depot (module `figstyle`, presets
journal) ; definir SCI_FIGURE_SCRIPTS pour pointer vers son repertoire si elle n'est pas au
meme endroit que sur la machine d'origine.
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
import matplotlib.image as mpimg

OUT = ROOT / "article" / "figures"


def main() -> None:
    fig, axes = fs.panel_grid("nature_double", nrows=1, ncols=2, width_ratios=[1, 1], height_mm=95)
    for ax, fname in zip(axes, ["figure3a_complex_overview.png", "figure3b_site_closeup.png"]):
        img = mpimg.imread(OUT / fname)
        ax.imshow(img)
        ax.axis("off")
    fs.panel_labels(list(axes))
    out = fs.save(fig, OUT / "figure3_site", "nature_double", raster=True)
    print("Ecrit:", *out, sep="\n  ")


if __name__ == "__main__":
    main()
