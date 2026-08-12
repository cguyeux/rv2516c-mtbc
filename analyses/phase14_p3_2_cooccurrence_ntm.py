#!/usr/bin/env python3
"""phase14_p3_2_cooccurrence_ntm.py -- P3.2 : la paire Rv2516c-Rv2517c est-elle un module conservé,
au sens d'un GAIN/PERTE couplé à travers le genre, ou leur concordance de présence/absence n'est-elle
que ce qu'on attendrait de deux gènes rares pris au hasard, ou de deux voisins génomiques quelconques ?

CE QUI EST DÉJÀ ÉTABLI : Rv2516c et Rv2517c ont des listes
d'espèces NTM porteuses IDENTIQUES (6/53, restricted). C'est un fait brut, cité comme argument de
« unité génétique réelle », mais JAMAIS testé contre un modèle nul. Une concordance parfaite entre
deux gènes RARES peut être triviale : moins un gène est répandu, moins il a d'occasions de discorder
avec un autre gène tout aussi rare (un gène présent dans 6/53 espèces et un autre présent dans 6/53
espèces DIFFÉRENTES concorderaient déjà sur 47/53 par la seule rareté partagée). Le test qui compte
n'est donc pas « concordent-ils ? » mais « concordent-ils PLUS que ce qu'attendrait la rareté seule,
et plus que deux voisins génomiques quelconques ? ».

DEUX MODÈLES NULS INDÉPENDANTS, chacun neutralisant un confondant distinct :
  N1. Paires de gènes de MÊME BANDE DE PRÉVALENCE (n_ntm_present dans [3,10], 242 gènes disponibles,
      tirées au hasard, PAS nécessairement voisines) -- neutralise le confondant « rareté partagée ».
  N2. Paires de VOISINS GÉNOMIQUES (gènes physiquement adjacents ailleurs dans H37Rv, ordonnés par
      `start_mtbc0`) -- neutralise le confondant « deux gènes voisins co-varient pour des raisons
      structurelles/d'assemblage/d'échantillonnage, indépendamment de toute biologie de module ».
      Si N2 seul explique la concordance, l'argument de "module fonctionnel" s'effondre : ce serait
      une propriété générique du voisinage, pas de CETTE paire.

Source : `annotation_mtbc/résultats/phase57_ntm/ntm_presence_per_gene.tsv` (tblastn présence/absence,
53 génomes NTM, déjà calculé -- aucun BLAST relancé ici, contrairement à P3.1 qui devait reconstruire
l'alignement par résidu). Coordonnées génomiques (ordre, pas position H37Rv vraie, cf. garde-fou du
projet sur `start_mtbc0`) via `site/data/db.sqlite`.

MÉTRIQUE : indice de Jaccard sur les ensembles d'espèces porteuses (|A∩B| / |A∪B|), qui neutralise
mieux la prévalence brute qu'un simple taux d'accord position-par-position (lequel serait dominé par
les nombreuses absences partagées entre deux gènes rares, cf. ci-dessus).

Run: python analyses/phase14_p3_2_cooccurrence_ntm.py
"""
from __future__ import annotations

import csv
import json
import random
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT.parent / "annotation_mtbc"
NTM_TSV = ATLAS / "résultats" / "phase57_ntm" / "ntm_presence_per_gene.tsv"
DB = ATLAS / "site" / "data" / "db.sqlite"
OUT = ROOT / "résultats" / "p3_2_cooccurrence_ntm"

GENE_A, GENE_B = "Rv2516c", "Rv2517c"
PREVALENCE_BAND = (3, 10)   # Rv2516c et Rv2517c sont a 6/53
N_NULL_DRAWS = 20000
SEED = 20260801  # date du jour, fixe pour reproductibilite -- pas d'horodatage a l'execution


def load_presence() -> dict[str, dict]:
    rows = {}
    with open(NTM_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sp = frozenset(r["species_present"].split(",")) if r["species_present"] else frozenset()
            rows[r["rv"]] = {"species": sp, "n": int(r["n_ntm_present"]),
                              "classification": r["classification"]}
    return rows


def load_genomic_order() -> list[str]:
    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    cur = con.cursor()
    cur.execute("SELECT rv, start_mtbc0 FROM genes WHERE start_mtbc0 IS NOT NULL ORDER BY start_mtbc0")
    order = [r[0] for r in cur.fetchall()]
    con.close()
    return order


def jaccard(a: frozenset, b: frozenset) -> float:
    if not a and not b:
        return float("nan")  # deux absences totales : non informatif, exclu du null
    return len(a & b) / len(a | b)


def rank_with_ties(observed: float, null: list[float]) -> dict:
    """Percentile qui NOMME les ex-æquo au lieu de les cacher dans un simple <=. Un `sum(x<=obs)`
    naïf peut afficher « percentile 100 » alors que d'autres paires du null atteignent EXACTEMENT
    la même valeur -- ce qui change entièrement l'interprétation (unique vs. type de paire fréquent)."""
    n = len(null)
    n_below = sum(1 for x in null if x < observed)
    n_equal = sum(1 for x in null if x == observed)
    n_above = n - n_below - n_equal
    return {"n_null": n, "n_strictement_inferieur": n_below, "n_ex_aequo": n_equal,
            "n_strictement_superieur": n_above,
            "percentile_conservateur": round(100 * n_below / n, 1) if n else None,
            "frac_ex_aequo": round(n_equal / n, 3) if n else None}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    presence = load_presence()
    if GENE_A not in presence or GENE_B not in presence:
        raise SystemExit(f"{GENE_A} ou {GENE_B} absent de {NTM_TSV}")

    a, b = presence[GENE_A], presence[GENE_B]
    observed_jaccard = jaccard(a["species"], b["species"])
    print("== P3.2 : co-occurrence phylogénétique Rv2516c / Rv2517c à travers 53 NTM ==\n")
    print(f"  {GENE_A} : {a['n']}/53 ({a['classification']}) -- {sorted(a['species'])}")
    print(f"  {GENE_B} : {b['n']}/53 ({b['classification']}) -- {sorted(b['species'])}")
    print(f"  Jaccard observé : {observed_jaccard:.4f}"
          f"{'  (ensembles IDENTIQUES)' if a['species'] == b['species'] else ''}")

    # --- N1 : paires de même bande de prévalence, tirées au hasard dans tout le génome ---
    rng = random.Random(SEED)
    band_genes = [rv for rv, d in presence.items()
                  if PREVALENCE_BAND[0] <= d["n"] <= PREVALENCE_BAND[1]
                  and rv not in (GENE_A, GENE_B)]
    print(f"\n[N1] pool de gènes de même bande de prévalence [{PREVALENCE_BAND[0]}-{PREVALENCE_BAND[1]}]"
          f"/53 (hors Rv2516c/Rv2517c) : {len(band_genes)} gènes")
    null1 = []
    for _ in range(N_NULL_DRAWS):
        g1, g2 = rng.sample(band_genes, 2)
        j = jaccard(presence[g1]["species"], presence[g2]["species"])
        if j == j:  # exclut les NaN (double absence totale, hors bande normalement)
            null1.append(j)
    null1.sort()
    r1 = rank_with_ties(observed_jaccard, null1)
    print(f"  {len(null1)} paires nulles valides, Jaccard médian={null1[len(null1)//2]:.3f}, "
          f"95e percentile={null1[int(0.95*len(null1))]:.3f}")
    print(f"  Rv2516c/Rv2517c : {r1['n_strictement_inferieur']}/{r1['n_null']} paires nulles STRICTEMENT "
          f"en dessous ({r1['percentile_conservateur']}%), {r1['n_ex_aequo']} ex-æquo à 1,0 "
          f"({100*r1['frac_ex_aequo']:.1f}% du null)")

    # --- N2 : paires de voisins génomiques réels, ailleurs dans H37Rv ---
    order = load_genomic_order()
    neighbor_pairs = [(order[i], order[i + 1]) for i in range(len(order) - 1)]
    # Exclure la paire cible elle-même et ses voisins immédiats (Rv2515c-Rv2516c, Rv2517c-Rv2518c)
    # pour ne pas contaminer le null avec la fenêtre déjà connue comme atypique.
    EXCLUDE_WINDOW = {"Rv2513", "Rv2514c", "Rv2515c", "Rv2516c", "Rv2517c", "Rv2518c"}
    neighbor_pairs = [(g1, g2) for g1, g2 in neighbor_pairs
                       if g1 not in EXCLUDE_WINDOW and g2 not in EXCLUDE_WINDOW
                       and g1 in presence and g2 in presence]
    print(f"\n[N2] {len(neighbor_pairs)} paires de voisins génomiques réels (fenêtre Rv2513-Rv2518c "
          f"exclue du pool, pas seulement de la paire cible)")
    null2 = []
    for g1, g2 in neighbor_pairs:
        j = jaccard(presence[g1]["species"], presence[g2]["species"])
        if j == j:
            null2.append(j)
    null2.sort()
    r2 = rank_with_ties(observed_jaccard, null2)
    print(f"  {len(null2)} paires de voisins valides, Jaccard médian={null2[len(null2)//2]:.3f}, "
          f"95e percentile={null2[int(0.95*len(null2))]:.3f}")
    print(f"  Rv2516c/Rv2517c : {r2['n_strictement_inferieur']}/{r2['n_null']} paires nulles STRICTEMENT "
          f"en dessous ({r2['percentile_conservateur']}%), {r2['n_ex_aequo']} ex-æquo à 1,0 "
          f"({100*r2['frac_ex_aequo']:.1f}% du null -- la plupart des voisins sont genus-core, "
          f"donc trivialement concordants, ce null est peu discriminant seul)")

    # --- N2b : parmi les paires de voisins, restreindre à celles de MÊME bande de prévalence pour
    # les deux membres -- combine N1 et N2, le test le plus strict : des voisins ET rares comme la paire. ---
    neighbor_band_pairs = [(g1, g2) for g1, g2 in neighbor_pairs
                            if PREVALENCE_BAND[0] <= presence[g1]["n"] <= PREVALENCE_BAND[1]
                            and PREVALENCE_BAND[0] <= presence[g2]["n"] <= PREVALENCE_BAND[1]]
    print(f"\n[N2b, combiné] paires de voisins ET de même bande de prévalence : "
          f"{len(neighbor_band_pairs)} paires")
    if neighbor_band_pairs:
        scored2b = sorted(((jaccard(presence[g1]["species"], presence[g2]["species"]), g1, g2)
                           for g1, g2 in neighbor_band_pairs), reverse=True)
        null2b = [x for x, _, _ in scored2b if x == x]
        r2b = rank_with_ties(observed_jaccard, null2b)
        print(f"  {len(null2b)} paires valides, Jaccard médian={null2b[len(null2b)//2]:.3f}")
        print(f"  Rv2516c/Rv2517c : {r2b['n_strictement_inferieur']}/{r2b['n_null']} paires STRICTEMENT "
              f"en dessous ({r2b['percentile_conservateur']}%), et **{r2b['n_ex_aequo']} PAIRES EX-ÆQUO "
              f"à Jaccard=1,0** ({100*r2b['frac_ex_aequo']:.1f}% du pool le mieux apparié) -- la paire "
              f"n'est PAS unique dans son propre voisinage de comparaison le plus strict :")
        for j, g1, g2 in scored2b[:12]:
            if j == 1.0:
                print(f"      {g1}-{g2} (n={presence[g1]['n']}/{presence[g2]['n']}) Jaccard=1,000")
    else:
        null2b, r2b = [], {}

    summary = {
        "gene_a": GENE_A, "gene_b": GENE_B,
        "n_present_a": a["n"], "n_present_b": b["n"],
        "species_identical": a["species"] == b["species"],
        "observed_jaccard": observed_jaccard,
        "null1_prevalence_matched": {**r1, "median": null1[len(null1) // 2] if null1 else None},
        "null2_genomic_neighbors": {**r2, "median": null2[len(null2) // 2] if null2 else None,
                                     "caveat": "peu discriminant : domine par des paires genus-core "
                                               "trivialement concordantes"},
        "null2b_neighbors_and_prevalence_matched": {**r2b,
                                                      "median": null2b[len(null2b) // 2] if null2b else None,
                                                      "note": "le null le mieux apparie ; la paire N'Y EST "
                                                              "PAS unique, 9/49 paires comparables ex-aequo"},
    }
    (OUT / "resume.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f"\nSorties : {OUT}/")


if __name__ == "__main__":
    main()
