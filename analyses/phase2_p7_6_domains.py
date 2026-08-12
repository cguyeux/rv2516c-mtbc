#!/usr/bin/env python3
"""phase2_p7_6_domains.py -- P7.6 : que contient le segment C-terminal orphelin ?

P7.1 a montré que le repli de type S6 ne couvre que les résidus ~2-83 et que la
fenêtre HTH prédite s'arrête à 121 : plus de la moitié des 267 résidus n'a aucun
candidat, tout en étant repliée (pLDDT 78,4). Cette phase demande ce que contient
ce segment, en découpant la protéine et en interrogeant chaque morceau séparément
-- le signal de la protéine entière étant manifestement dominé par son domaine
N-terminal, c'est probablement pour cela que le Foldseek aveugle de l'atlas n'a
jamais rendu que des hits S6.

Jeu de référence ÉLARGI par rapport à P7.1, à partir des 45 hits Foldseek complets
de l'atlas (la fiche n'en montrait que 8, ce qui masquait deux thèmes) :
  - S6            : permutants circulaires de la protéine ribosomale S6
  - ACT_ferredoxin: domaines ACT / RAM / PII, qui partagent le repli ferredoxin-like
                    (βαββαβ) avec S6 -- « type S6 » et « type ACT » sont le MÊME repli
  - Lrp_AsnC      : régulateurs Lrp/AsnC entiers = HTH + domaine RAM/ACT. Inclut
                    Rv3291c (2w29) et le régulateur de réponse à la leucine de Mtb
                    (2qz8), c'est-à-dire le paralogue MTBC de cette famille.
  - sigma_r4, MerR: les deux classes de liaison à l'ADN déjà testées en P7.1.

GARDE-FOU PARALOGUE-DE-FOLD (KB) : la fonction « régulateur Lrp/AsnC » est DÉJÀ
attribuée dans ce génome, à Rv3291c. Un score élevé contre cette famille ne peut
donc PAS conduire à écrire que Rv2516c est un Lrp ; au mieux qu'il en partage
l'architecture. C'est écrit ici pour que la conclusion ne dérive pas.

Frontières de domaines : déterminées OBJECTIVEMENT par balayage de la carte de
contacts (le point de coupure qui minimise les contacts inter-segments), pas à
l'œil ni sur les bornes Pfam -- sans quoi le découpage préjugerait du résultat.

Entrées : résultats/p7_1/structures/, sigma_meta.json
Sorties : résultats/p7_6/{domain_split.tsv, tmalign_segments.tsv, baseline.tsv}
Run: python analyses/phase2_p7_6_domains.py
"""
from __future__ import annotations
import json
import math
import sys
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase1_p7_1_tmalign import (  # noqa: E402
    parse_pdb_chains, write_chain, slice_residues, tmalign, aligned_span_on_second, AA3,
)

STRUCT = ROOT / "résultats" / "p7_1" / "structures"
OUT = ROOT / "résultats" / "p7_6"
SEG = OUT / "segments"

# Références ajoutées en P7.6, tirées des 45 hits Foldseek complets de l'atlas.
NEW_REFS = {
    "2f1f": ("A", "ACT_ferredoxin", "acetohydroxyacid synthase regulatory subunit (ACT)"),
    "3tvi": ("H", "ACT_ferredoxin", "aspartate kinase, C. acetobutylicum (ACT)"),
    "4c3l": ("A", "ACT_ferredoxin", "PII signal transduction protein, S. elongatus"),
    "2qz8": ("B", "Lrp_AsnC", "M. tuberculosis leucine response regulator (Lrp/AsnC)"),
    "2w29": ("B", "Lrp_AsnC", "Rv3291c G102T (M. tuberculosis Lrp/AsnC)"),
    "3i4p": ("A", "Lrp_AsnC", "AsnC family transcriptional regulator, Agrobacterium"),
}
OLD_REFS = {
    "3hh0": ("A", "MerR", "MerR family regulator, B. cereus"),
    "5d8c": ("A", "MerR", "HiNmlR MerR-family, DNA-bound"),
    "5d90": ("C", "MerR", "HiNmlR MerR-family, DNA-bound"),
    "7b90": ("E", "S6", "circular permutant of ribosomal protein S6"),
    "7bff": ("E", "S6", "circular permutant of ribosomal protein S6"),
    "7bfd": ("K", "S6", "circular permutant of ribosomal protein S6"),
}


def ca_coords(path: Path) -> dict[int, tuple[float, float, float]]:
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[17:20].strip() in AA3:
            out[int(line[22:26])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out


def domain_split(ca: dict, min_seg: int = 40, cutoff: float = 8.0) -> list[tuple[int, float]]:
    """Balayage du point de coupure minimisant les contacts inter-segments.

    Pour chaque frontière b, compte les paires de Cα à moins de `cutoff` situées de
    part et d'autre, en excluant les voisins de chaîne (|i-j| <= 4) qui donneraient
    un minimum trivial. Normalise par min(taille des deux segments) pour ne pas
    favoriser mécaniquement les coupures près des extrémités.
    """
    res = sorted(ca)
    n = len(res)
    scores = []
    for k in range(min_seg, n - min_seg):
        b = res[k]
        left, right = res[:k], res[k:]
        c = 0
        for i in left:
            xi, yi, zi = ca[i]
            for j in right:
                if j - i <= 4:
                    continue
                xj, yj, zj = ca[j]
                if (xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2 < cutoff * cutoff:
                    c += 1
        scores.append((b, c / min(len(left), len(right))))
    return scores


def main() -> None:
    print("== P7.6 : le segment C-terminal orphelin ==")
    SEG.mkdir(parents=True, exist_ok=True)
    query = STRUCT / "query_Rv2516c.pdb"
    ca = ca_coords(query)
    print(f"  modèle : {len(ca)} résidus")

    # -- 1. frontières objectives --
    scores = domain_split(ca)
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "domain_split.tsv", "w") as fh:
        fh.write("frontiere\tcontacts_inter_normalises\n")
        for b, s in scores:
            fh.write(f"{b}\t{s:.4f}\n")
    ordered = sorted(scores, key=lambda t: t[1])
    print("\n-- Frontières candidates (contacts inter-segments les plus faibles) --")
    for b, s in ordered[:8]:
        print(f"     résidu {b:3d} : {s:.3f}")
    best = ordered[0][0]
    print(f"  -> coupure principale retenue : {best}")

    # Segments testés. On garde la coupure objective ET les bornes issues de P7.1,
    # pour que le résultat ne dépende pas d'un seul découpage.
    segments = {
        "full_1-267": (1, 267),
        f"Nterm_1-{best-1}": (1, best - 1),
        f"Cterm_{best}-267": (best, 267),
        "Nterm_S6span_1-95": (1, 95),
        "Cterm_afterHTH_122-267": (122, 267),
        "HTHwindow_85-135": (85, 135),
    }
    seg_paths = {}
    for name, (lo, hi) in segments.items():
        dest = SEG / f"{name}.pdb"
        n = slice_residues(query, lo, hi, dest)
        if n >= 20:
            seg_paths[name] = dest
        print(f"  segment {name:24} {n:3d} résidus")

    # -- 2. jeu de références --
    meta = json.load(open(ROOT / "résultats" / "p7_1" / "sigma_meta.json"))
    entries = []
    for rv, m in sorted(meta.items()):
        src = STRUCT / f"sigma_{rv}.pdb"
        if not src.exists():
            continue
        dest = SEG / f"r4_{rv}.pdb"
        if slice_residues(src, m["r4_from"], m["r4_to"], dest) >= 20:
            entries.append((f"{m['gene']}_r4", "sigma_r4", dest))
    for pid, (chain, klass, _desc) in {**OLD_REFS, **NEW_REFS}.items():
        src = STRUCT / f"pdb_{pid}.pdb"
        if not src.exists():
            print(f"  [absent] {pid}")
            continue
        chains = parse_pdb_chains(src)
        use = chain if chain in chains else (sorted(chains)[0] if chains else None)
        if use is None:
            continue
        dest = SEG / f"{pid}_{use}.pdb"
        n = write_chain(chains[use], dest)
        entries.append((f"{pid}_{use}", klass, dest))
    from collections import Counter
    print(f"\n  {len(entries)} références : {dict(Counter(k for _, k, _ in entries))}")

    # -- 3. chaque segment contre chaque référence --
    rows = []
    for sname, spath in seg_paths.items():
        for rname, klass, rpath in entries:
            r = tmalign(rpath, spath)
            if not r:
                continue
            span = aligned_span_on_second(r["raw"])
            rows.append({"segment": sname, "reference": rname, "classe": klass,
                         "tm": r["tm_ref"], "aligned": r.get("aligned"), "rmsd": r.get("rmsd"),
                         "span": f"{span[0]}-{span[1]}" if span else "?"})
    with open(OUT / "tmalign_segments.tsv", "w") as fh:
        fh.write("segment\treference\tclasse\ttm_norm_ref\taligned\trmsd\tspan_sur_segment\n")
        for d in rows:
            fh.write(f"{d['segment']}\t{d['reference']}\t{d['classe']}\t{d['tm']:.4f}\t"
                     f"{d['aligned']}\t{d['rmsd']:.2f}\t{d['span']}\n")

    # -- 4. ligne de base entre références (indispensable pour lire l'échelle) --
    base = []
    for (n1, k1, p1), (n2, k2, p2) in combinations(entries, 2):
        r = tmalign(p1, p2)
        if r:
            base.append((k1, k2, r["tm_ref"]))
    with open(OUT / "baseline.tsv", "w") as fh:
        fh.write("classe_a\tclasse_b\ttm\n")
        for k1, k2, t in base:
            fh.write(f"{k1}\t{k2}\t{t:.4f}\n")

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else float("nan")

    classes = sorted({k for _, k, _ in entries})
    print("\n-- Ligne de base : références entre elles (médiane TM) --")
    hdr = "            " + "".join(f"{c[:11]:>13}" for c in classes)
    print(hdr)
    for ka in classes:
        line = f"{ka[:11]:12}"
        for kb in classes:
            v = [t for k1, k2, t in base if {k1, k2} == {ka, kb} or (ka == kb and k1 == k2 == ka)]
            line += f"{med(v):13.3f}" if v else f"{'-':>13}"
        print(line)

    print("\n-- Chaque segment de Rv2516c, par classe (médiane TM) --")
    hdr = f"{'segment':26}" + "".join(f"{c[:11]:>13}" for c in classes)
    print(hdr)
    for sname in seg_paths:
        line = f"{sname:26}"
        for c in classes:
            v = [d["tm"] for d in rows if d["segment"] == sname and d["classe"] == c]
            line += f"{med(v):13.3f}" if v else f"{'-':>13}"
        print(line)

    print("\n-- Meilleurs appariements du segment C-terminal --")
    ct = sorted([d for d in rows if d["segment"].startswith("Cterm")], key=lambda d: -d["tm"])[:8]
    for d in ct:
        print(f"  {d['segment']:24} {d['reference']:12} {d['classe']:15} TM={d['tm']:.3f} "
              f"RMSD={d['rmsd']:.2f} alig={d['aligned']}")

    print(f"\nÉcrit {OUT}/")


if __name__ == "__main__":
    main()
