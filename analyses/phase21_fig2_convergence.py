#!/usr/bin/env python3
"""phase21_fig2_convergence.py -- Figure 2 du manuscrit : convergence de quatre lignes de preuve
indépendantes sur le même site candidat de liaison à l'ADN (module HTH, unité B).

Panneau A : fraction de résidus basiques exposés (P2.3, SASA Shrake-Rupley) -- fenêtre de
reconnaissance (111-121) vs protéine entière vs unité A (contrôle de spécificité) vs hélice
adjacente (101-110).
Panneau B : conservation par position sur la fenêtre HTH (98-121), 4 lignées NTM phylogénétiquement
indépendantes après correction de la pseudo-réplication du clade M. avium (P3.1), bande 111-118
invariante mise en évidence.
Panneau C : distances protéine-ADN au contact le plus proche (J1, sauvage) pour R113/H115/R118 vs
la médiane sur les 60 résidus de l'unité B (P2.6, corrigé cette session : 7,19 Å, pas 7,24).
Panneau D : réponse du modèle à la mutation R113A/H115A/R118A -- delta de distance par résidu,
restreint aux résidus réellement en contact avec l'ADN au départ (P2.6).

Toutes les valeurs recalculées directement depuis les fichiers de résultats bruts, pas recopiées
de la prose.

Sorties : article/figures/{figure2_convergence.pdf,.png}
Run: python analyses/phase21_fig2_convergence.py

Dépend d'une bibliothèque de style de figures externe à ce dépôt (module `figstyle`, presets
journal) ; définir SCI_FIGURE_SCRIPTS pour pointer vers son répertoire si elle n'est pas au
même endroit que sur la machine d'origine.
"""
from __future__ import annotations
import csv
import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, os.environ.get(
    "SCI_FIGURE_SCRIPTS",
    str(Path.home() / "docs/codes/claude_plugins/bio_population_genetics/skills/sci-figure/scripts"),
))
import figstyle as fs

OUT = ROOT / "article" / "figures"
OUT.mkdir(parents=True, exist_ok=True)
RES = ROOT / "résultats"

C_TEST, C_CTRL, C_NEG = fs.PALETTE_CATEGORICAL[1], fs.PALETTE_CATEGORICAL[0], "#999999"


def load_charge():
    d = json.load(open(RES / "p2_3_charge_hth" / "resume.json"))["sasa"]
    return {
        "Recognition\nwindow\n(111-121)": d["tour_plus_helice"]["frac_basic_of_exposed"],
        "Whole\nprotein": d["proteine_entiere"]["frac_basic_of_exposed"],
        "Unit A\n(specificity\nctrl)": d["unite_A"]["frac_basic_of_exposed"],
        "Adjacent helix\n(101-110)": d["helice2"]["frac_basic_of_exposed"],
    }


def load_conservation():
    d = json.load(open(RES / "p3_1_ntm_orthologs" / "identite_par_position.json"))
    lineages = ["M_avium", "M_celatum", "M_shinjukuense", "M_simiae"]
    out = {}
    for pos in range(98, 122):
        vals = [d[sp].get(str(pos)) for sp in lineages if d[sp].get(str(pos)) is not None]
        out[pos] = sum(vals) / len(vals) if vals else None
    return out


def load_contacts():
    def load(path):
        with open(path) as f:
            return {int(r["seq_id"]): float(r["min_dist_any_atom"])
                    for r in csv.DictReader(f, delimiter="\t")}
    j1 = load(RES / "p2_6_interface_contacts" / "j1_contacts.tsv")
    j2 = load(RES / "p2_6_interface_contacts" / "j2_contacts.tsv")
    residues = {113: 26, 115: 28, 118: 31}  # H37Rv residue -> seq_id (offset 87)
    dist_wt = {r: j1[s] for r, s in residues.items()}
    median_all = sorted(j1.values())[len(j1) // 2 - 1:len(j1) // 2 + 1]
    median_all = sum(median_all) / 2
    deltas = {r: j2[s] - j1[s] for r, s in residues.items()}
    contact_deltas = sorted(((j2[s] - j1[s]) for s, d in j1.items() if d < 5.0), reverse=True)
    next_best = contact_deltas[3] if len(contact_deltas) > 3 else None  # 4th largest among contacts
    return dist_wt, median_all, deltas, next_best


def panel_a(ax):
    data = load_charge()
    labels = list(data.keys())
    values = list(data.values())
    colors = [C_TEST, C_NEG, C_NEG, C_NEG]
    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.8)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01, f"{val:.2f}",
                ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("Fraction of exposed residues\nthat are basic")
    ax.set_ylim(0, 0.55)
    ax.tick_params(axis="x", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)


def panel_b(ax):
    cons = load_conservation()
    positions = sorted(cons)
    values = [cons[p] * 100 for p in positions]
    colors = [C_TEST if 111 <= p <= 118 else C_NEG for p in positions]
    ax.bar(positions, values, color=colors, edgecolor="black", linewidth=0.4, width=0.8)
    ax.axvspan(110.5, 118.5, color=C_TEST, alpha=0.12, zorder=0)
    ax.text(114.5, 104, "invariant\n(111-118)", ha="center", va="bottom", fontsize=6.5)
    ax.set_xlabel("Residue position (HTH window)")
    ax.set_ylabel("Identity across 4 independent\nNTM lineages (%)")
    ax.set_ylim(0, 118)
    ax.set_xticks([98, 105, 111, 118, 121])
    ax.spines[["top", "right"]].set_visible(False)


def panel_c(ax):
    dist_wt, median_all, _, _ = load_contacts()
    labels = [f"Arg{r}" if r != 115 else f"His{r}" for r in [113, 115, 118]]
    values = [dist_wt[r] for r in [113, 115, 118]]
    bars = ax.bar(labels, values, color=C_TEST, edgecolor="black", linewidth=0.8)
    ax.axhline(median_all, color="black", linestyle="--", linewidth=1.0)
    ax.text(2.05, median_all, f"median over 60\nresidues: {median_all:.2f} Å",
            ha="left", va="center", fontsize=6.5)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.1, f"{val:.2f}",
                ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("Closest distance to DNA\nin wild-type model (Å)")
    ax.set_xlim(-0.6, 3.0)
    ax.set_ylim(0, 8.5)
    ax.spines[["top", "right"]].set_visible(False)


def panel_d(ax):
    _, _, deltas, next_best = load_contacts()
    labels = ["Arg113A", "His115A", "Arg118A", "next-best\n(contact res.)"]
    values = [deltas[113], deltas[115], deltas[118], next_best]
    colors = [C_TEST, C_TEST, C_TEST, C_NEG]
    bars = ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.8)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.05, f"+{val:.2f}",
                ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("Distance increase to DNA,\nwild-type $\\to$ mutant (Å)")
    ax.set_ylim(0, 4.0)
    ax.tick_params(axis="x", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    fig, axes = fs.panel_grid("nature_double", nrows=2, ncols=2, height_mm=150)
    panel_a(axes[0, 0])
    panel_b(axes[0, 1])
    panel_c(axes[1, 0])
    panel_d(axes[1, 1])
    fs.panel_labels(list(axes.flat))
    out = fs.save(fig, OUT / "figure2_convergence", "nature_double", raster=True)
    print("Ecrit:", *out, sep="\n  ")


if __name__ == "__main__":
    main()
