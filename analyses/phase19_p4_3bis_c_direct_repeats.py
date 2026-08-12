#!/usr/bin/env python3
"""phase19_p4_3bis_c_direct_repeats.py -- P4.3bis-c, volet répétitions DIRECTES.

POURQUOI CE SCRIPT, ET CE QU'IL AJOUTE À PHASE10.
P4.3bis-a (phase10) n'a testé qu'une symétrie GÉNÉRIQUE, les répétitions inversées, parce qu'un wHTH
dimérique lit un palindrome. Négatif sur les deux nuls dans les 4 régions. P4.3bis-b (dimérisation
apo du module HTH) est ressorti INCONCLUANT le 2026-08-03 pour deux raisons indépendantes (biais
d'empilement générique de Boltz sur petits domaines + dimérisation de BldC probablement couplée à
l'ADN, invisible en test apo) -- il ne tranche donc PAS entre "l'unité B dimérise" et "elle ne
dimérise pas". Beaucoup de régulateurs bactériens à un seul domaine de liaison lisent des
**répétitions DIRECTES** (deux copies du même demi-site, même orientation, souvent en tandem) plutôt
qu'un palindrome -- ce mode ne suppose PAS de dimère structural obligatoire, une protéine peut lier
deux copies successives du site en tandem sans jamais se dimériser en solution. C'est donc un test
complémentaire de phase10, pas conditionné à son résultat.

CE QUI EST PRÉ-ENREGISTRÉ AVANT LE CALCUL (contre le dragage après un négatif) :
  - MÊMES 4 régions que phase10 (R1, R3, R4, R2), aucune région ajoutée ni retirée.
  - MÊME fenêtre de bras (6-12 pb), MÊME espaceur (0-8 pb), MÊME tolérance (2 mésappariements) --
    changer uniquement la fonction de comparaison (identité directe des deux bras, pas complémentarité
    inversée), rien d'autre.
  - MÊME double nul (permutation dinucléotidique + intergéniques génomiques réelles).
  - Seuil de significativité CORRIGÉ pour tenir compte des tests déjà faits en phase10 : 4 régions x
    2 formes de motif (inversé + direct) = 8 tests indépendants sur ce voisinage. Bonferroni :
    alpha = 0,05 / 8 = 0,00625, appliqué aux DEUX p-values (permutation ET génomique) de CE script.

CE QUE CE SCRIPT NE FAIT PAS, ET POURQUOI -- "sites asymétriques".
La piste P4.3bis-c demandait aussi de chercher des "sites asymétriques" (un demi-site UNIQUE, sans
répétition ni symétrie). Ce test n'est PAS implémenté ici, délibérément : sans consensus connu ni
donnée de footprinting/ChIP pour cette famille (AlpA/excisionase) chez les actinobactéries, un demi-
site asymétrique unique n'a AUCUN dénominateur statistique -- n'importe quel k-mer de 6-8 pb apparaît
une fois par hasard dans une région de quelques centaines de pb, et "trouver" un tel k-mer ne serait
rien d'autre que du bruit habillé en résultat. Un test honnête de ce mode nécessiterait soit une
matrice de la famille (alignement curé de sites connus, absent ici), soit une donnée expérimentale
(footprinting). Consigné comme limite plutôt que comme test fabriqué.

Entrées : ../investigate_phylo/resources/{NC_000962.3.fasta, NC_000962.3.gff3}
Sorties : résultats/p4_3bis_c_direct_repeats/{candidats.tsv, resume.md}
Run: python analyses/phase19_p4_3bis_c_direct_repeats.py
"""
from __future__ import annotations
import random, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT.parent / "investigate_phylo" / "resources"
OUT = ROOT / "résultats" / "p4_3bis_c_direct_repeats"

# Régions identiques à phase10_p4_3bis_operateur.py, ne pas en ajouter/retirer (pré-enregistrement).
REGIONS = {
    "R1_prom_operon_2516_2517": (2833762, 2834108, "promoteur de l'opéron Rv2517c-Rv2516c (autorégulation ?)"),
    "R3_prom_TA_2515_2514":     (2832592, 2832709, "promoteur de l'opéron toxine-antitoxine, 118 pb de Rv2516c"),
    "R4_inter_2513_2514c":      (2830584, 2830876, "intergénique Rv2513 / Rv2514c"),
    "R2_prom_transposase":      (2829804, 2829953, "promoteur de Rv2512c, transposase IS1081"),
}
ARM = range(6, 13)          # bras de la répétition, identique à phase10
SPACER = range(0, 9)        # espaceur central, identique à phase10
MAXMM = 2                   # mésappariements tolérés dans le bras, identique à phase10
N_NULL = 400                # régions/permutations tirées, identique à phase10
N_TESTS_TOTAL = 8           # 4 régions x (inversé [phase10] + direct [ici]) -- Bonferroni
ALPHA_CORRIGE = 0.05 / N_TESTS_TOTAL


def best_dr(seq: str) -> tuple[int, int, int, int, str]:
    """Meilleure répétition DIRECTE (mêmes bras dans le MÊME sens, pas de complémentarité inversée).
    Score = longueur du bras moins les mésappariements, comme phase10.best_ir."""
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
                mm = sum(1 for a, b in zip(left, right) if a != b)   # identité directe, pas rc()
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
    print("== P4.3bis-c : répétitions DIRECTES dans le voisinage de Rv2516c ==")
    print(f"   (répliqué de phase10_p4_3bis_operateur.py, seule la comparaison de bras change ;")
    print(f"   seuil corrigé Bonferroni sur {N_TESTS_TOTAL} tests : alpha = {ALPHA_CORRIGE:.5f})\n")
    genome = "".join(l.strip() for l in (RES / "NC_000962.3.fasta").read_text().splitlines()
                     if not l.startswith(">")).upper()
    print(f"  génome {len(genome):,} pb, GC {100*sum(c in 'GC' for c in genome)/len(genome):.1f} %\n")

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
        obs = best_dr(seq)
        perm = [best_dr(dinuc_shuffle(seq, rng))[0] for _ in range(N_NULL)]
        p_perm = (sum(1 for s in perm if s >= obs[0]) + 1) / (N_NULL + 1)
        pool = [r for r in inter if 0.7 * L <= r[1] - r[0] + 1 <= 1.3 * L]
        samp = rng.sample(pool, min(N_NULL, len(pool)))
        gnull = [best_dr(genome[s - 1:e][:400])[0] for s, e in samp]
        p_gen = (sum(1 for s in gnull if s >= obs[0]) + 1) / (len(gnull) + 1)
        rows.append({"name": name, "a": a, "b": b, "L": L, "gc": gc, "obs": obs,
                     "p_perm": p_perm, "p_gen": p_gen,
                     "mu_perm": st.mean(perm), "mu_gen": st.mean(gnull) if gnull else float("nan"),
                     "n_gen": len(gnull), "desc": desc})
        print(f"-- {name} ({a}-{b}, {L} pb, GC {gc:.0f} %) --")
        print(f"   {desc}")
        print(f"   meilleure répétition directe : bras {obs[1]} pb, espaceur {obs[2]}, score {obs[0]}")
        print(f"     motif : {obs[4]}")
        print(f"   nul par permutation dinucléotidique : moyenne {st.mean(perm):.2f}, p = {p_perm:.3f}")
        print(f"   nul intergénique génomique (n={len(gnull)}) : moyenne "
              f"{st.mean(gnull) if gnull else float('nan'):.2f}, p = {p_gen:.3f}")
        verdict = ("CANDIDAT : dépasse les DEUX nuls, au seuil CORRIGÉ"
                   if p_perm < ALPHA_CORRIGE and p_gen < ALPHA_CORRIGE
                   else "non retenu : ne dépasse pas les deux nuls au seuil corrigé")
        print(f"   -> {verdict}\n")

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "candidats.tsv", "w") as fh:
        fh.write("region\tdebut\tfin\tlongueur\tgc\tbras\tespaceur\tscore\tmotif\t"
                 "moy_nul_permut\tp_permut\tmoy_nul_genomique\tp_genomique\tretenu_seuil_corrige\tdescription\n")
        for r in rows:
            o = r["obs"]
            keep = int(r["p_perm"] < ALPHA_CORRIGE and r["p_gen"] < ALPHA_CORRIGE)
            fh.write(f"{r['name']}\t{r['a']}\t{r['b']}\t{r['L']}\t{r['gc']:.1f}\t{o[1]}\t{o[2]}\t"
                     f"{o[0]}\t{o[4]}\t{r['mu_perm']:.2f}\t{r['p_perm']:.3f}\t{r['mu_gen']:.2f}\t"
                     f"{r['p_gen']:.3f}\t{keep}\t{r['desc']}\n")
    kept = [r for r in rows if r["p_perm"] < ALPHA_CORRIGE and r["p_gen"] < ALPHA_CORRIGE]
    print("-- BILAN --")
    if kept:
        print(f"   {len(kept)} région(s) portent une répétition directe dépassant les deux nuls "
              f"au seuil corrigé ({ALPHA_CORRIGE:.5f}) :")
        for r in kept:
            print(f"     {r['name']} -> {r['obs'][4]}")
        print("   À lire comme un INDICE de site, jamais comme un site : aucune matrice de liaison")
        print("   n'est connue pour cette famille chez les actinobactéries.")
    else:
        print("   AUCUNE région ne dépasse les deux nuls au seuil corrigé. Combiné au négatif de")
        print("   phase10 (répétitions inversées), aucune des deux symétries génériques testées ne")
        print("   ressort dans ce voisinage. Les sites asymétriques restent non testables faute de")
        print("   consensus ou de donnée de footprinting (cf. docstring) -- limite, pas un négatif.")
    print(f"\nÉcrit {OUT}/candidats.tsv")


if __name__ == "__main__":
    main()
