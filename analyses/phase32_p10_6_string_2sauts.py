#!/usr/bin/env python3
"""phase32_p10_6_string_2sauts.py -- P10.6 : le module Rv2516c-Rv2517c est-il topologiquement
ISOLÉ dans le réseau STRING de H37Rv, ou branché sur quelque chose que le top-15 tronqué cache ?

CE QUE LA FICHE NE PEUT PAS DIRE. Le champ `string.json` des fiches d'atlas est tronqué au TOP-15
partenaires par gène. Toute affirmation d'« isolement » tirée
de cette fiche est donc une déduction sur données tronquées, pas un résultat de graphe. Le fichier
BULK complet est sur disque (`annotation_mtbc/résultats/phase2h_string/`, 13 Mo, ~1,4 M arêtes) :
il n'y a aucune raison de continuer à raisonner sur la troncature.

LE CONTRE-ARGUMENT DE LA PISTE EST LE CŒUR DU DISPOSITIF, PAS UNE NOTE DE BAS DE PAGE.
La piste prévoit que « les 3 partenaires immédiats sont uniquement `context_driven` par le canal de
voisinage génomique, donc une expansion à 2 sauts risque de n'ajouter que du bruit de proximité en
cascade ». Si l'on se contente d'expanser le graphe combiné, on ne pourra pas distinguer un vrai
module fonctionnel d'une chaîne d'adjacence chromosomique : le canal `neighborhood` de STRING relie
mécaniquement les gènes voisins, donc un 2-sauts sur ce canal REDÉCOUVRE le chromosome. Le script
construit donc DEUX graphes et les compare :

  - graphe CONTEXTE-INCLUS  : noisy-OR des 6 canaux sans le textmining (`combined_no_tm`), c'est-à-dire
    voisinage + fusion + cooccurrence + coexpression + expérimental + base de données ;
  - graphe SANS CONTEXTE   : noisy-OR de coexpression + expérimental + base de données SEULEMENT.
    Le voisinage génomique, la fusion et la cooccurrence sont retirés parce qu'ils sont des
    prédicteurs DE POSITION, pas des mesures d'interaction. Ce que Rv2516c garde ici est ce qui reste
    quand on interdit à la proximité chromosomique de parler.

CALIBRATION, SANS QUOI « ISOLÉ » N'EST QU'UNE IMPRESSION. Un voisinage à 2 sauts de 40 protéines est
petit ou grand selon le graphe. On rapporte donc systématiquement le PERCENTILE de Rv2516c parmi les
3 906 protéines, et surtout un nul APPARIÉ SUR LE DEGRÉ : parmi les protéines de même degré ±1, quelle
est la taille médiane du 2-sauts ? Un gène peu connecté a mécaniquement un petit 2-sauts ; la seule
question intéressante est s'il est plus petit que ce que son degré impose.

QUESTION SUBSIDIAIRE POSÉE PAR LA PISTE : existe-t-il ailleurs dans le protéome un paralogue de repli
AlpA/excisionase mieux branché, dont le voisinage nommerait une fonction ? On interroge db.sqlite
(couches `domains` et `pfam_tentative`, seuil ET sous-seuil) sur les quatre familles Pfam de la
super-famille identifiée en P8.4, puis on lit leur connectivité.

Entrées : annotation_mtbc/résultats/phase2h_string/83332.protein.links.detailed.v12.0.txt.gz
          annotation_mtbc/site/data/db.sqlite
Sorties : résultats/p10_6_string_2sauts/{voisinages.tsv, calibration.tsv, paralogues.tsv, resume.md}
Run: python analyses/phase32_p10_6_string_2sauts.py
"""
from __future__ import annotations

import gzip
import json
import random
import sqlite3
import statistics as st
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT.parent / "annotation_mtbc"
LINKS = ATLAS / "résultats" / "phase2h_string" / "83332.protein.links.detailed.v12.0.txt.gz"
DB = ATLAS / "site" / "data" / "db.sqlite"
OUT = ROOT / "résultats" / "p10_6_string_2sauts"

PRIOR = 0.041
MEDIUM, HIGH = 400, 700
FOCUS = ["Rv2516c", "Rv2517c"]
CASSETTE = ["Rv2512c", "Rv2513", "Rv2514c", "Rv2515c", "Rv2516c", "Rv2517c", "Rv2518c"]
# familles Pfam de la super-famille AlpA/excisionase recensées en P8.4
ALPA_PFAM = {"PF05930": "Phage_AlpA", "PF06806": "DUF1233", "PF11112": "PyocinActivator",
             "PF09035": "Tn916-Xis", "PF31169": "DUF8830"}
SEED = 20260810


def combine(scores) -> int:
    """noisy-OR de STRING sur une liste de canaux (0..1000) -> 0..1000."""
    prod = 1.0
    for s in scores:
        sp = (s / 1000.0 - PRIOR) / (1.0 - PRIOR)
        if sp > 0:
            prod *= (1.0 - sp)
    tot = 1.0 - prod
    return round((tot + PRIOR * (1.0 - tot)) * 1000)


def load_edges():
    """(u, v) -> (score_contexte_inclus, score_sans_contexte)."""
    edges: dict[tuple[str, str], tuple[int, int]] = {}
    with gzip.open(LINKS, "rt") as fh:
        header = fh.readline().split()
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            p = line.split()
            u, v = p[idx["protein1"]].split(".")[-1], p[idx["protein2"]].split(".")[-1]
            if u >= v:
                continue                      # le fichier est symétrique : une arête suffit
            ch = {c: float(p[idx[c]]) for c in
                  ("neighborhood", "fusion", "cooccurence", "coexpression",
                   "experimental", "database")}
            with_ctx = combine(list(ch.values()))
            no_ctx = combine([ch["coexpression"], ch["experimental"], ch["database"]])
            edges[(u, v)] = (with_ctx, no_ctx)
    return edges


def graph(edges, which: int, thr: int) -> dict[str, set[str]]:
    g: dict[str, set[str]] = defaultdict(set)
    for (u, v), sc in edges.items():
        if sc[which] >= thr:
            g[u].add(v)
            g[v].add(u)
    return g


def hops(g, seeds: list[str]) -> tuple[set[str], set[str]]:
    one = set()
    for s in seeds:
        one |= g.get(s, set())
    one -= set(seeds)
    two = set()
    for n in one:
        two |= g.get(n, set())
    two -= one | set(seeds)
    return one, two


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    edges = load_edges()
    print(f"{len(edges)} arêtes STRING (non redondantes) chargées")

    variants = {
        "contexte_inclus_400": (0, MEDIUM), "contexte_inclus_700": (0, HIGH),
        "sans_contexte_400": (1, MEDIUM), "sans_contexte_700": (1, HIGH),
    }
    rng = random.Random(SEED)
    resume = ["# P10.6 — expansion STRING à 2 sauts depuis Rv2516c/Rv2517c", ""]
    vlines = ["variante\tgraine\tdegre\tn_1saut\tn_2sauts\tpartenaires_1saut"]
    clines = ["variante\tgraine\tdegre\tpct_degre\tn_2sauts\tpct_2sauts\t"
              "median_2sauts_degre_apparie\tn_apparies"]

    for vname, (which, thr) in variants.items():
        g = graph(edges, which, thr)
        deg = {n: len(v) for n, v in g.items()}
        alldeg = sorted(deg.values())
        two_sizes = {}
        for n in g:
            _, t2 = hops(g, [n])
            two_sizes[n] = len(t2)
        all2 = sorted(two_sizes.values())
        print(f"\n== {vname} : {len(g)} noeuds, {sum(deg.values())//2} arêtes, "
              f"degré médian {st.median(alldeg):.0f}")

        for seed in FOCUS:
            d = deg.get(seed, 0)
            one, two = hops(g, [seed])
            vlines.append(f"{vname}\t{seed}\t{d}\t{len(one)}\t{len(two)}\t"
                          f"{','.join(sorted(one))}")
            pct_d = 100.0 * sum(1 for x in alldeg if x < d) / len(alldeg) if alldeg else float("nan")
            pct_2 = 100.0 * sum(1 for x in all2 if x < len(two)) / len(all2) if all2 else float("nan")
            matched = [two_sizes[n] for n in g if abs(deg[n] - d) <= 1 and n != seed]
            clines.append(f"{vname}\t{seed}\t{d}\t{pct_d:.1f}\t{len(two)}\t{pct_2:.1f}\t"
                          f"{st.median(matched) if matched else float('nan')}\t{len(matched)}")
            print(f"   {seed}: degré {d} (percentile {pct_d:.0f}), 1-saut {len(one)}, "
                  f"2-sauts {len(two)} (percentile {pct_2:.0f}, nul apparié médian "
                  f"{st.median(matched) if matched else float('nan')})")
            print(f"      1-saut = {sorted(one)}")

        # la cassette entière comme graine : le module a-t-il une porte de sortie collective ?
        one, two = hops(g, [c for c in CASSETTE if c in g])
        vlines.append(f"{vname}\tCASSETTE\t-\t{len(one)}\t{len(two)}\t{','.join(sorted(one))}")
        print(f"   cassette Rv2512c-Rv2518c : 1-saut {len(one)}, 2-sauts {len(two)}")
        if one:
            print(f"      1-saut = {sorted(one)}")

    (OUT / "voisinages.tsv").write_text("\n".join(vlines) + "\n")
    (OUT / "calibration.tsv").write_text("\n".join(clines) + "\n")

    # ---------- paralogues de repli AlpA/excisionase dans le protéome -----------------------
    con = sqlite3.connect(DB)
    g400 = graph(edges, 0, MEDIUM)
    plines = ["rv\tsource\tpfam\tproduit\tdegre_string_400\tpartenaires_top"]
    found = 0
    for rv, domains, tentative, product in con.execute(
            "select rv, domains, pfam_tentative, product_h37rv from genes"):
        hitfams = []
        for src, blob in (("domains", domains), ("pfam_tentative", tentative)):
            if not blob:
                continue
            for acc, name in ALPA_PFAM.items():
                if acc in blob or name in blob:
                    hitfams.append(f"{src}:{name}")
        if not hitfams:
            continue
        found += 1
        parts = sorted(g400.get(rv, set()))
        plines.append(f"{rv}\t{'|'.join(hitfams)}\t-\t{(product or '')[:60]}\t"
                      f"{len(parts)}\t{','.join(parts[:12])}")
    (OUT / "paralogues.tsv").write_text("\n".join(plines) + "\n")
    print(f"\nprotéines portant un Pfam de la super-famille AlpA/excisionase "
          f"(seuil ou sous-seuil) : {found}")
    for l in plines[1:]:
        print("   " + l.replace("\t", "  "))

    resume.append(f"{len(edges)} arêtes STRING non redondantes ; voir voisinages.tsv, "
                  f"calibration.tsv, paralogues.tsv.")
    (OUT / "resume.md").write_text("\n".join(resume) + "\n")
    print(f"\nsorties dans {OUT}")


if __name__ == "__main__":
    main()
