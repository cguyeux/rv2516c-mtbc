#!/usr/bin/env python3
"""P10.5 -- PAE inter-domaines du modele AlphaFold DB de Rv2516c.

Question : les unites A (1-87, ferredoxin-like DUF8830), B (88-147, wHTH AlpA) et
C (178-267, beta-sandwich Ig-like) sont-elles rigidement empaquetees les unes contre
les autres, ou orientees librement ? Enjeu : si l'unite A ou C peut occulter la face
de liaison a l'ADN de l'unite B, cela change l'interpretation du site candidat
R111/R113/H115/R118 et le design d'un futur construct experimental.

GARDE-FOU CENTRAL, ecrit dans l'enonce de la piste AVANT le calcul : un PAE
inter-domaine eleve ne distingue PAS une vraie flexibilite biologique d'une simple
absence de signal de co-evolution dans le MSA. Sans temoin, le chiffre est
ininterpretable. Ce script fournit donc la calibration manquante : quatre proteines
de reference passees par le MEME pipeline AFDB, deux a arrangement inter-domaine
experimentalement RIGIDE et deux a flexibilite inter-domaine ETABLIE.

  RIGIDES
    P9WIE5  KatG de M. tuberculosis (740 aa). Domaines N-terminal et C-terminal
            resolus ensemble comme un seul corps rigide (PDB 1SJ2), interface
            etendue. Meme organisme, meme pipeline : le temoin le mieux apparie.
    P0A786  Chaine catalytique de l'aspartate transcarbamylase d'E. coli (311 aa).
            Domaines polaire et equatorial a large interface, cristallises en un
            corps unique.
  FLEXIBLES
    P0DP23  Calmoduline humaine (149 aa). Deux lobes EF-hand relies par un linker
            central dont la flexibilite est etablie par RMN (les lobes culbutent
            independamment en solution). Cas d'ecole.
    P03023  Represseur LacI d'E. coli (360 aa). Tete de lecture HTH reliee au coeur
            par une charniere mobile en l'absence d'operateur -- analogie
            structurale directe avec la question posee ici (un HTH relie a d'autres
            domaines).

Sortie : resultats/p10_5_pae_interdomaine/
"""
import json
import statistics
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "résultats" / "p10_5_pae_interdomaine"

# (accession, etiquette, {nom_domaine: (debut, fin)} en numerotation UniProt 1-based)
PANEL = [
    ("I6YDM0", "Rv2516c (SUJET)", None, {
        "A_ferredoxin": (1, 87),
        "B_wHTH": (88, 147),
        "C_Ig_like": (178, 267),
    }),
    ("P9WIE5", "KatG M. tuberculosis", "RIGIDE", {
        "N_term": (1, 430),
        "C_term": (445, 740),
    }),
    ("P0A786", "ATCase chaine catalytique", "RIGIDE", {
        "polaire": (1, 133),
        "equatorial": (134, 284),
    }),
    ("P0DP23", "Calmoduline humaine", "FLEXIBLE", {
        "lobe_N": (5, 75),
        "lobe_C": (82, 146),
    }),
    ("P03023", "LacI E. coli", "FLEXIBLE", {
        "tete_HTH": (1, 49),
        "coeur": (62, 330),
    }),
]

FILES = {
    "I6YDM0": "AF-I6YDM0-F1-pae_v6.json",
    "P9WIE5": "pae_P9WIE5.json",
    "P0A786": "pae_P0A786.json",
    "P0DP23": "pae_P0DP23.json",
    "P03023": "pae_P03023.json",
}


def load_pae(path):
    d = json.loads(Path(path).read_text())
    e = d[0] if isinstance(d, list) else d
    return e["predicted_aligned_error"], e.get("max_predicted_aligned_error")


def block(pae, r1, r2):
    """Valeurs PAE du bloc lignes r1 x colonnes r2 (bornes UniProt 1-based inclusives)."""
    a1, b1 = r1
    a2, b2 = r2
    return [pae[i][j] for i in range(a1 - 1, b1) for j in range(a2 - 1, b2)]


def med(xs):
    return statistics.median(xs) if xs else float("nan")


def analyse(pae, domains):
    """Rend (intra par domaine, inter par paire non ordonnee symetrisee)."""
    intra = {n: med(block(pae, r, r)) for n, r in domains.items()}
    names = list(domains)
    inter = {}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            n1, n2 = names[i], names[j]
            vals = block(pae, domains[n1], domains[n2]) + block(pae, domains[n2], domains[n1])
            inter[(n1, n2)] = med(vals)
    return intra, inter


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    lines = []
    summary_rows = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("P10.5 -- PAE inter-domaines, sujet vs panel de calibration")
    emit("=" * 78)
    emit("PAE = erreur d'alignement predite (Angstroms) : confiance du modele sur la")
    emit("POSITION RELATIVE de deux residus. Faible = arrangement contraint ; eleve =")
    emit("le modele ne sait pas comment les deux parties sont orientees l'une par")
    emit("rapport a l'autre. Mediane par bloc, symetrisee sur les deux triangles.")
    emit("")

    for acc, label, expected, domains in PANEL:
        pae, maxpae = load_pae(OUT / FILES[acc])
        intra, inter = analyse(pae, domains)
        intra_ref = med(list(intra.values()))
        emit(f"--- {label}  [{acc}]" + (f"   attendu : {expected}" if expected else ""))
        emit(f"    PAE max du modele : {maxpae} A")
        for n, v in intra.items():
            emit(f"    intra {n:<16s} mediane {v:6.2f} A")
        for (n1, n2), v in inter.items():
            ratio = v / intra_ref if intra_ref else float("nan")
            emit(f"    INTER {n1} <-> {n2}   mediane {v:6.2f} A   "
                 f"(x{ratio:.1f} le PAE intra-domaine)")
            summary_rows.append({
                "accession": acc, "proteine": label, "attendu": expected or "SUJET",
                "paire": f"{n1}|{n2}", "pae_inter_median": round(v, 2),
                "pae_intra_median_ref": round(intra_ref, 2),
                "ratio_inter_sur_intra": round(ratio, 2),
                "pct_du_max": round(100 * v / maxpae, 1) if maxpae else None,
            })
        emit("")

    # --- lecture calibree -------------------------------------------------
    emit("=" * 78)
    emit("LECTURE CALIBREE")
    emit("")
    rig = [r for r in summary_rows if r["attendu"] == "RIGIDE"]
    flx = [r for r in summary_rows if r["attendu"] == "FLEXIBLE"]
    suj = [r for r in summary_rows if r["attendu"] == "SUJET"]
    rig_max = max(r["pae_inter_median"] for r in rig)
    flx_min = min(r["pae_inter_median"] for r in flx)
    emit(f"  Temoins RIGIDES   : PAE inter-domaine mediane de "
         f"{min(r['pae_inter_median'] for r in rig):.2f} a {rig_max:.2f} A")
    emit(f"  Temoins FLEXIBLES : PAE inter-domaine mediane de "
         f"{flx_min:.2f} a {max(r['pae_inter_median'] for r in flx):.2f} A")
    emit(f"  Les deux classes {'SE SEPARENT' if rig_max < flx_min else 'SE CHEVAUCHENT'} "
         f"(rigide max {rig_max:.2f} vs flexible min {flx_min:.2f} A)")
    emit("")
    for r in suj:
        v = r["pae_inter_median"]
        if rig_max < flx_min:
            if v <= rig_max:
                verdict = "dans la plage RIGIDE"
            elif v >= flx_min:
                verdict = "dans la plage FLEXIBLE"
            else:
                verdict = "ENTRE les deux plages, non classable"
        else:
            verdict = "panel non discriminant, non classable"
        emit(f"  {r['paire']:<28s} {v:6.2f} A  -> {verdict}")

    (OUT / "pae_interdomaine.tsv").write_text(
        "\t".join(summary_rows[0].keys()) + "\n"
        + "\n".join("\t".join(str(v) for v in r.values()) for r in summary_rows) + "\n")
    (OUT / "rapport.txt").write_text("\n".join(lines) + "\n")
    emit("")
    emit(f"Ecrit : {OUT/'pae_interdomaine.tsv'}")
    emit(f"Ecrit : {OUT/'rapport.txt'}")


if __name__ == "__main__":
    main()
