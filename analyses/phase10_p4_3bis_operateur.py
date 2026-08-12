#!/usr/bin/env python3
"""phase10_p4_3bis_operateur.py -- P4.3bis, volet séquence : y a-t-il un OPÉRATEUR candidat
dans les régions intergéniques que Rv2516c pourrait commander ?

POURQUOI CE TEST PLUTÔT QU'UN CO-REPLIEMENT.
P4.3 voulait des partenaires PROTÉIQUES (β-flap de RpoB, sigma, anti-sigma) : liste morte avec la
piste région-4. Sa reformulation P4.3bis proposait la machinerie d'élément mobile. Mais le mécanisme
de la famille assignée dit autre chose : **AlpA, chez *E. coli*, est un ACTIVATEUR TRANSCRIPTIONNEL**
du gène d'excision du prophage CP4-57 — il se fixe à l'ADN en amont d'un gène, il n'est pas un
partenaire protéique de la recombinase. Si Rv2516c suit ce mécanisme, le test direct est la recherche
d'un site opérateur, pas un co-repliement. C'est aussi le test le moins cher du dossier.

LA CARTE GÉNOMIQUE FIXE LES CANDIDATS, et elle corrige la priorité annoncée dans P4.3bis.
Tout le voisinage est sur le brin MOINS, donc les promoteurs sont du côté des COORDONNÉES HAUTES.

    Rv2518c  2834109-2835335 (-)
      [ R1 : 2833762-2834108, 347 pb ]  <- promoteur de l'opéron Rv2517c-Rv2516c (AUTORÉGULATION ?)
    Rv2517c  2833510-2833761 (-)
    Rv2516c  2832710-2833513 (-)   chevauchement de 4 pb avec Rv2517c
      [ R3 : 2832592-2832709, 118 pb ]  <- promoteur de l'opéron TA Rv2515c-Rv2514c
    Rv2515c  2831344-2832591 (-)   antitoxine
    Rv2514c  2830877-2831338 (-)   toxine, 5 pb seulement -> co-transcrite avec l'antitoxine
      [ R4 : 2830584-2830876, 293 pb ]
    Rv2513   2830161-2830583 (+)
    Rv2512a  2829954-2830139 (+)
      [ R2 : 2829804-2829953, 150 pb ]  <- promoteur de Rv2512c, transposase IS1081
    Rv2512c  2828556-2829803 (-)

**R3 devient le candidat le plus serré**, et ce n'est pas ce que P4.3bis avait annoncé : le promoteur
du couple toxine-antitoxine est à 118 pb de Rv2516c, alors que la transposase est à ~3 kb et séparée
par deux gènes sur l'autre brin. Un régulateur commandant la TA adjacente est une hypothèse plus
économique que le même régulateur commandant un élément mobile distant.

CE QUE CHERCHE CE SCRIPT. Les winged HTH de cette famille se fixent en dimère, donc sur des
**répétitions inversées** (palindromes imparfaits). On énumère, dans chaque région, les répétitions
inversées de bras 6 à 12 pb, avec 0 à 2 mésappariements et un espaceur de 0 à 8 pb.

LE MODÈLE NUL EST L'ESSENTIEL, ET IL EST DOUBLE.
H37Rv est à 65 % GC : les palindromes GC-riches y sont mécaniquement fréquents, et « j'ai trouvé une
répétition inversée » ne veut rien dire sans dénominateur. Donc :
  (a) **nul par permutation locale** : mêmes régions, séquences mélangées en préservant la composition
      en DINUCLÉOTIDES (préserver seulement le %GC ne suffit pas, la structure dinucléotidique crée
      déjà des palindromes) ;
  (b) **nul empirique génomique** : N régions intergéniques tirées au hasard dans H37Rv, de même
      longueur, qui donnent la distribution réelle du meilleur palindrome par région.
Un candidat n'est retenu que s'il dépasse les DEUX nuls.

LIMITE À ÉNONCER, et elle est sévère. Aucun site de liaison n'est connu pour Rv2516c ni pour aucun
homologue proche ; il n'y a donc pas de matrice à chercher, seulement une propriété générique de
symétrie. Un résultat positif dirait « il existe ici une répétition inversée inhabituelle », ce qui
est un INDICE de site, pas un site. Un résultat négatif, lui, est informatif : il retirerait
l'argument le plus simple en faveur d'une autorégulation ou d'un contrôle de la TA voisine.

Entrées : ../investigate_phylo/resources/{NC_000962.3.fasta, NC_000962.3.gff3}
Sorties : résultats/p4_3bis_operateur/{candidats.tsv, resume.md}
Run: python analyses/phase10_p4_3bis_operateur.py
"""
from __future__ import annotations
import random, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT.parent / "investigate_phylo" / "resources"
OUT = ROOT / "résultats" / "p4_3bis_operateur"

REGIONS = {
    "R1_prom_operon_2516_2517": (2833762, 2834108, "promoteur de l'opéron Rv2517c-Rv2516c (autorégulation ?)"),
    "R3_prom_TA_2515_2514":     (2832592, 2832709, "promoteur de l'opéron toxine-antitoxine, 118 pb de Rv2516c"),
    "R4_inter_2513_2514c":      (2830584, 2830876, "intergénique Rv2513 / Rv2514c"),
    "R2_prom_transposase":      (2829804, 2829953, "promoteur de Rv2512c, transposase IS1081"),
}
COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
ARM = range(6, 13)          # bras de la répétition inversée
SPACER = range(0, 9)        # espaceur central
MAXMM = 2                   # mésappariements tolérés dans le bras
N_NULL = 400                # régions intergéniques tirées, et permutations par région


def rc(s: str) -> str:
    return "".join(COMP.get(c, "N") for c in reversed(s))


def best_ir(seq: str) -> tuple[int, int, int, int, str]:
    """Meilleure répétition inversée : (score, bras, espaceur, position, motif).
    Score = longueur du bras moins les mésappariements ; on privilégie le bras long."""
    best = (0, 0, 0, -1, "")
    n = len(seq)
    for arm in ARM:
        for sp in SPACER:
            span = 2 * arm + sp
            if span > n:
                continue
            for i in range(n - span + 1):
                left = seq[i:i + arm]
                right = seq[i + arm + sp:i + span]
                if "N" in left or "N" in right:
                    continue
                mm = sum(1 for a, b in zip(left, rc(right)) if a != b)
                if mm > MAXMM:
                    continue
                score = arm - mm
                if score > best[0]:
                    best = (score, arm, sp, i, seq[i:i + span])
    return best


# dinuc_shuffle vit désormais dans lib_sequence_nulls.py (reforge du 2026-08-10).
# La version locale était une marche eulérienne NAÏVE : sur impasse elle sautait vers un
# autre nucléotide, inventant un dinucléotide absent de l'original — le contrat annoncé
# « préserve la composition en dinucléotides » était donc faux par intermittence. La
# version partagée REJETTE le tirage au lieu de le rapiécer. Voir l'en-tête de la lib
# pour la portée sur les résultats déjà acquis : conclusions NÉGATIVES de P4.3bis-a/-c
# inchangées (le biais était dilué sur 400 tirages et conservatif), mais les p-values
# ne seront pas reproduites au chiffre près si on relance.
from lib_sequence_nulls import dinuc_shuffle  # noqa: E402


def main() -> None:
    print("== P4.3bis, volet séquence : opérateur candidat dans le voisinage de Rv2516c ==")
    genome = "".join(l.strip() for l in (RES / "NC_000962.3.fasta").read_text().splitlines()
                     if not l.startswith(">")).upper()
    print(f"  génome {len(genome):,} pb, GC {100*sum(c in 'GC' for c in genome)/len(genome):.1f} %\n")

    # ── nul (b) : distribution du meilleur palindrome dans des intergéniques réelles ──
    genes = []
    for line in (RES / "NC_000962.3.gff3").read_text().splitlines():
        if "\tgene\t" in line:
            f = line.split("\t")
            genes.append((int(f[3]), int(f[4])))
    genes.sort()
    inter = [(genes[i][1] + 1, genes[i + 1][0] - 1) for i in range(len(genes) - 1)
             if genes[i + 1][0] - genes[i][1] - 1 >= 60]
    print(f"  {len(inter)} régions intergéniques >= 60 pb dans H37Rv (nul génomique)\n")
    rng = random.Random(0)

    rows = []
    for name, (a, b, desc) in REGIONS.items():
        seq = genome[a - 1:b]
        L = len(seq)
        gc = 100 * sum(c in "GC" for c in seq) / L
        obs = best_ir(seq)
        # nul (a) : permutations dinucléotidiques de CETTE région
        perm = [best_ir(dinuc_shuffle(seq, rng))[0] for _ in range(N_NULL)]
        p_perm = (sum(1 for s in perm if s >= obs[0]) + 1) / (N_NULL + 1)
        # nul (b) : intergéniques réelles de longueur comparable (+/- 30 %)
        pool = [r for r in inter if 0.7 * L <= r[1] - r[0] + 1 <= 1.3 * L]
        samp = rng.sample(pool, min(N_NULL, len(pool)))
        gnull = [best_ir(genome[s - 1:e][:400])[0] for s, e in samp]
        p_gen = (sum(1 for s in gnull if s >= obs[0]) + 1) / (len(gnull) + 1)
        rows.append({"name": name, "a": a, "b": b, "L": L, "gc": gc, "obs": obs,
                     "p_perm": p_perm, "p_gen": p_gen,
                     "mu_perm": st.mean(perm), "mu_gen": st.mean(gnull) if gnull else float("nan"),
                     "n_gen": len(gnull), "desc": desc})
        print(f"-- {name} ({a}-{b}, {L} pb, GC {gc:.0f} %) --")
        print(f"   {desc}")
        print(f"   meilleure répétition inversée : bras {obs[1]} pb, espaceur {obs[2]}, "
              f"score {obs[0]}")
        print(f"     motif : {obs[4]}")
        print(f"   nul par permutation dinucléotidique : moyenne {st.mean(perm):.2f}, p = {p_perm:.3f}")
        print(f"   nul intergénique génomique (n={len(gnull)}) : moyenne "
              f"{st.mean(gnull) if gnull else float('nan'):.2f}, p = {p_gen:.3f}")
        verdict = ("CANDIDAT : dépasse les DEUX nuls" if p_perm < 0.05 and p_gen < 0.05
                   else "non retenu : ne dépasse pas les deux nuls")
        print(f"   -> {verdict}\n")

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "candidats.tsv", "w") as fh:
        fh.write("region\tdebut\tfin\tlongueur\tgc\tbras\tespaceur\tscore\tmotif\t"
                 "moy_nul_permut\tp_permut\tmoy_nul_genomique\tp_genomique\tretenu\tdescription\n")
        for r in rows:
            o = r["obs"]
            keep = int(r["p_perm"] < 0.05 and r["p_gen"] < 0.05)
            fh.write(f"{r['name']}\t{r['a']}\t{r['b']}\t{r['L']}\t{r['gc']:.1f}\t{o[1]}\t{o[2]}\t"
                     f"{o[0]}\t{o[4]}\t{r['mu_perm']:.2f}\t{r['p_perm']:.3f}\t{r['mu_gen']:.2f}\t"
                     f"{r['p_gen']:.3f}\t{keep}\t{r['desc']}\n")
    kept = [r for r in rows if r["p_perm"] < 0.05 and r["p_gen"] < 0.05]
    print("-- BILAN --")
    if kept:
        print(f"   {len(kept)} région(s) portent une répétition inversée dépassant les deux nuls :")
        for r in kept:
            print(f"     {r['name']} -> {r['obs'][4]}")
        print("   À lire comme un INDICE de site, jamais comme un site : aucune matrice de liaison")
        print("   n'est connue pour cette famille chez les actinobactéries.")
    else:
        print("   AUCUNE région ne dépasse les deux nuls. C'est un négatif informatif : il retire")
        print("   l'argument le plus simple en faveur d'une autorégulation de l'opéron ou d'un")
        print("   contrôle direct de la TA voisine par une symétrie d'opérateur détectable.")
    print(f"\nÉcrit {OUT}/candidats.tsv")


if __name__ == "__main__":
    main()
