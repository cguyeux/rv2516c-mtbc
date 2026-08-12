#!/usr/bin/env python3
"""phase4_p5_2_1_conservation_positionnelle.py -- P5.2.1 : contrainte sélective par unité.

Question. Les unités de Rv2516c sont-elles sous des contraintes sélectives différentes ?
En particulier le module HTH AlpA/excisionase (unité B, 88-147), dont HHpred donne la
famille avec confiance, est-il plus contraint que le reste du gène ?

Intérêt : c'est le seul argument de fonctionnalité INDÉPENDANT de la structure et de
l'homologie. P7.1, P7.6.1 et P2.1 reposent tous sur de la comparaison à des bases ; ici
on ne regarde que la variation observée dans la population, donc rien n'est circulaire
vis-à-vis des assignations de repli.

Méthode, et ses garde-fous (KB, pathogène clonal massivement échantillonné).
  - Témoin interne = la classe SYNONYME du MÊME gène. Elle subit la même mutabilité, le
    même échantillonnage et la même structure clonale que la classe non-synonyme ; c'est
    la seule référence honnête. Comparer un taux brut au génome ne veut rien dire ici.
  - Stratification par FRÉQUENCE : à 145 000 souches, les singletons sont dominés par le
    bruit de séquençage et par la structure clonale. On rapporte séparément toutes les
    positions et celles au-dessus d'un plancher.
  - On compte des SITES, pas des fréquences d'allèles : un marqueur de lignée à forte
    fréquence compterait sinon pour des milliers d'observations (piège documenté).
  - Le gène est sur le brin MOINS : l'allèle SPDI est donné sur le brin plus et doit être
    complémenté avant de reconstituer le codon muté.

Entrées : annotation_mtbc/résultats/phase2f_conservation/snp_counts.tsv (comptes par
          position, génome entier), NC_000962.3.fasta, coordonnées H37Rv du GFF3.
Sorties : résultats/p5_2_1/{sites.tsv, resume.md}
Run: python analyses/phase4_p5_2_1_conservation_positionnelle.py
"""
from __future__ import annotations
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
GENOME = MTBC / "investigate_phylo" / "resources" / "NC_000962.3.fasta"
GFF = MTBC / "investigate_phylo" / "resources" / "NC_000962.3.gff3"
SNP = MTBC / "annotation_mtbc" / "résultats" / "phase2f_conservation" / "snp_counts.tsv"
OUT = ROOT / "résultats" / "p5_2_1"

RV = "Rv2516c"
UNITS = [("A", 1, 87), ("B", 88, 147), ("L", 148, 177), ("C", 178, 267)]
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
BASES = "TCAG"
AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AAS[i] for i, (a, b, c) in
        enumerate((x, y, z) for x in BASES for y in BASES for z in BASES)}


def unit_of(res: int) -> str:
    return next((n for n, a, b in UNITS if a <= res <= b), "?")


def main() -> None:
    print("== P5.2.1 : contrainte sélective par unité ==")
    seq = "".join(l.strip() for l in GENOME.read_text().splitlines() if not l.startswith(">"))
    start = end = None
    for line in GFF.read_text().splitlines():
        if "\tgene\t" in line and f"locus_tag={RV}" in line:
            f = line.split("\t")
            start, end = int(f[3]), int(f[4])
            break
    if start is None:
        raise SystemExit("coordonnées introuvables dans le GFF3")
    print(f"  {RV} H37Rv {start:,}-{end:,} (brin −, {end-start+1} pb)")

    # CDS sur le brin codant, et contrôle de traduction (le garde-fou qui a rattrapé
    # l'erreur de coordonnées MTBC0 : sans lui on analyse la mauvaise région en silence)
    cds = "".join(COMP[c] for c in reversed(seq[start - 1:end]))
    prot = "".join(CODE[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))
    assert prot[-1] == "*", "le dernier codon n'est pas un stop : coordonnées suspectes"
    print(f"  contrôle : {len(prot)-1} résidus, dernier codon = stop  OK")

    rows = []
    with open(SNP) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pos = int(f[0])
            if not (start <= pos <= end):
                continue
            alt_plus, count = f[1], int(f[2])
            if alt_plus not in COMP:
                continue
            off = end - pos                      # 0-based depuis le début du CDS
            ci, cp = off // 3, off % 3           # index de codon, position dans le codon
            res = ci + 1
            if res > 267:
                continue                          # codon stop
            codon = cds[ci * 3:ci * 3 + 3]
            alt = COMP[alt_plus]                 # brin −
            mut = codon[:cp] + alt + codon[cp + 1:]
            if codon[cp] == alt:
                continue                          # allèle identique à la référence
            aa_ref, aa_alt = CODE[codon], CODE[mut]
            eff = "syn" if aa_ref == aa_alt else ("nonsense" if aa_alt == "*" else "missense")
            rows.append({"pos": pos, "res": res, "unit": unit_of(res), "count": count,
                         "codon": codon, "mut": mut, "aa": f"{aa_ref}{res}{aa_alt}", "eff": eff})

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "sites.tsv", "w") as fh:
        fh.write("pos\tresidu\tunite\tcount\tcodon\tcodon_mut\taa\teffet\n")
        for r in sorted(rows, key=lambda d: d["res"]):
            fh.write(f"{r['pos']}\t{r['res']}\t{r['unit']}\t{r['count']}\t{r['codon']}\t"
                     f"{r['mut']}\t{r['aa']}\t{r['eff']}\n")
    print(f"  {len(rows)} sites variables dans le CDS")

    def table(sel, label):
        print(f"\n-- {label} ({len(sel)} sites) --")
        print(f"{'unité':7} {'aa':>4} {'syn':>5} {'mis':>5} {'non-sens':>9} {'mis/syn':>8} {'sites/100aa':>12}")
        out = []
        for name, a, b in UNITS:
            sub = [r for r in sel if r["unit"] == name]
            c = Counter(r["eff"] for r in sub)
            n_aa = b - a + 1
            ratio = c["missense"] / c["syn"] if c["syn"] else float("nan")
            print(f"{name:7} {n_aa:>4} {c['syn']:>5} {c['missense']:>5} {c['nonsense']:>9} "
                  f"{ratio:>8.2f} {100*len(sub)/n_aa:>12.1f}")
            out.append((name, n_aa, c["syn"], c["missense"], c["nonsense"], ratio, len(sub)))
        tot = Counter(r["eff"] for r in sel)
        r_tot = tot["missense"] / tot["syn"] if tot["syn"] else float("nan")
        print(f"{'TOTAL':7} {267:>4} {tot['syn']:>5} {tot['missense']:>5} {tot['nonsense']:>9} "
              f"{r_tot:>8.2f} {100*len(sel)/267:>12.1f}")
        return out, tot

    all_rows = rows
    floor = [r for r in rows if r["count"] >= 145]     # ~0,1 % de 145 000 souches (règle KB)
    t_all, tot_all = table(all_rows, "TOUTES les positions variables")
    t_flo, tot_flo = table(floor, "positions au-dessus du plancher de fréquence (count >= 145)")

    print("\n-- Sites au-dessus du plancher, en détail --")
    for r in sorted(floor, key=lambda d: -d["count"]):
        print(f"   {r['pos']}  res {r['res']:3d} ({r['unit']})  {r['aa']:9} {r['eff']:9} n={r['count']:,}")

    with open(OUT / "resume.md", "w") as fh:
        fh.write(f"# P5.2.1 -- contrainte selective par unite ({RV})\n\n")
        fh.write("Temoin interne = classe synonyme du meme gene. Sites comptes, pas frequences.\n\n")
        for label, t, tot in (("toutes positions", t_all, tot_all), ("count>=145", t_flo, tot_flo)):
            fh.write(f"\n## {label}\n\n| unite | aa | syn | mis | non-sens | mis/syn | sites/100aa |\n")
            fh.write("|---|---|---|---|---|---|---|\n")
            for n, naa, s, m, ns, ratio, nsit in t:
                fh.write(f"| {n} | {naa} | {s} | {m} | {ns} | {ratio:.2f} | {100*nsit/naa:.1f} |\n")
    print(f"\nEcrit {OUT}/")


if __name__ == "__main__":
    main()
