#!/usr/bin/env python3
"""phase33_p10_3_coexpression_genome.py -- P10.3 : scan de co-expression GÉNOME ENTIER, au-delà
des candidats choisis d'avance.

CE QUE P7.4 A FAIT, ET SA LIMITE. P7.4 a mesuré la co-expression de Rv2516c avec HUIT gènes choisis
a priori (Rv2517c, recA, dnaE2, hspX, Rv2623, Rv2514c, Rv2518c, sigA). Un test sur huit candidats ne
peut jamais découvrir un partenaire inattendu : il ne peut que confirmer ou infirmer une liste déjà
écrite. Les archives brutes contiennent l'expression de TOUT le génome et dormaient sur disque.

DEUX JEUX INDÉPENDANTS, ET C'EST LE POINT MÉTHODOLOGIQUE CENTRAL.
  - GSE166501 : RNA-seq, 27 échantillons, moxifloxacine 0 / 0,12x / 0,3x CMI x 4 / 24 / 72 h.
    Un seul antibiotique : la structure de covariance y est dominée par le TEMPS et la DOSE.
  - GSE71200 : puces bicolores TIGR, 106 lames, ~30 composés à 2x / 4x / 8x CMI, 16 h (Ma et al.
    2015, PMID 26618656). Chaque valeur est déjà un log-ratio traité/contrôle : la structure de
    covariance y est dominée par la NATURE DU COMPOSÉ.
Un partenaire retrouvé dans les DEUX est un partenaire qui survit à deux structures expérimentales
sans rapport. Un partenaire retrouvé dans un seul est un candidat, et rien de plus.

LE NUL, ET POURQUOI LE NUL NAÏF EST FAUX ICI. Tester r contre zéro serait absurde : dans un jeu où
tous les gènes répondent à un gradient commun (dose, temps, toxicité), la corrélation de deux gènes
tirés au hasard n'est PAS centrée sur zéro et sa dispersion est large. Le nul empirique est donc la
distribution des corrélations de 200 000 PAIRES DE GÈNES TIRÉES AU HASARD dans le même jeu. Une
corrélation ne compte que par son rang dans cette distribution.

LE TÉMOIN POSITIF, SANS QUOI UN NÉGATIF NE VEUT RIEN DIRE. Quelle corrélation atteint, dans CE jeu,
une paire de gènes dont on SAIT qu'ils sont co-transcrits ? On construit le témoin depuis
l'annotation : toutes les paires de gènes ADJACENTS, MÊME BRIN, séparées de moins de 20 pb
(proxy d'opéron largement utilisé chez les procaryotes). Si la médiane de ces paires est basse, le
jeu n'a pas le pouvoir de détecter un opéron et un négatif sur Rv2516c ne prouve rien. C'est le même
raisonnement que le contrôle recA de P7.4, transposé au génome entier.

Corrections : Spearman (robuste aux non-linéarités et aux valeurs extrêmes des puces),
Benjamini-Hochberg sur les ~3 900 tests par jeu.

Entrées : data/gse166501/GSE166501_RAW.tar (déjà extrait), data/gse71200/series_matrix.txt.gz
          + GPL5774_family.soft.gz, investigate_phylo/resources/NC_000962.3.gff3,
          annotation_mtbc/site/data/db.sqlite (produits, pour lire les têtes de liste).
Sorties : résultats/p10_3_coexpression/{gse166501_correlations.tsv, gse71200_correlations.tsv,
          convergence.tsv, resume.md}
Run: python analyses/phase33_p10_3_coexpression_genome.py
"""
from __future__ import annotations

import gzip
import random
import re
import sqlite3
import statistics as st
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
OUT = ROOT / "résultats" / "p10_3_coexpression"
GSE166 = ROOT / "data" / "gse166501"
GSE712 = ROOT / "data" / "gse71200"
GFF = MTBC / "investigate_phylo" / "resources" / "NC_000962.3.gff3"
DB = MTBC / "annotation_mtbc" / "site" / "data" / "db.sqlite"

FOCUS = ["Rv2516c", "Rv2517c"]
N_NULL_PAIRS = 200_000
MAX_INTERGENIC = 20      # pb : seuil de co-transcription présumée pour le témoin positif
SEED = 20260810


# ------------------------------------------------------------------ jeux de données
def load_gse166501() -> tuple[list[str], np.ndarray, list[str]]:
    """RNA-seq RPKM : un fichier par échantillon, colonnes (gène, ..., RPKM)."""
    files = sorted(GSE166.glob("GSM*.MTb.Transcript.txt.gz"))
    if not files:
        with tarfile.open(GSE166 / "GSE166501_RAW.tar") as tf:
            tf.extractall(GSE166)
        files = sorted(GSE166.glob("GSM*.MTb.Transcript.txt.gz"))
    per_sample: list[dict[str, float]] = []
    names = []
    for f in files:
        d: dict[str, float] = {}
        with gzip.open(f, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            # colonne d'expression : on prend la dernière colonne numérique nommée RPKM/expression
            icol = next((i for i, c in enumerate(head) if "rpkm" in c.lower()), len(head) - 1)
            igene = next((i for i, c in enumerate(head)
                          if c.lower() in ("gene", "gene_id", "locus", "name")), 0)
            for line in fh:
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(icol, igene):
                    continue
                m = re.match(r"(Rv\d+[A-Za-z]?)", p[igene])
                if not m:
                    continue
                try:
                    d[m.group(1)] = float(p[icol])
                except ValueError:
                    continue
        if d:
            per_sample.append(d)
            names.append(f.name.split("_")[0] + "_" + f.name.split("_")[1])
    genes = sorted(set.intersection(*[set(d) for d in per_sample]))
    mat = np.array([[per_sample[j][g] for j in range(len(per_sample))] for g in genes])
    mat = np.log2(mat + 1.0)                      # RPKM -> log, la corrélation de rang s'en moque
    return genes, mat, names


def load_gse71200() -> tuple[list[str], np.ndarray, list[str]]:
    """Puces bicolores : la matrice de série porte déjà des log-ratios traité/contrôle."""
    probe2rv: dict[str, str] = {}
    with gzip.open(GSE712 / "GPL5774_family.soft.gz", "rt", errors="replace") as fh:
        in_tab = False
        cols: list[str] = []
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_tab = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_tab and not cols:
                cols = line.rstrip("\n").split("\t")
                continue
            if in_tab:
                p = line.rstrip("\n").split("\t")
                row = dict(zip(cols, p))
                cand = row.get("Strain H37Rv") or row.get("ORF") or ""
                m = re.match(r"(Rv\d+[A-Za-z]?)", cand)
                if m:
                    probe2rv[row[cols[0]]] = m.group(1)

    with gzip.open(GSE712 / "series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
        header = [c.strip('"') for c in fh.readline().rstrip("\n").split("\t")]
        samples = header[1:]
        vals: dict[str, list[list[float]]] = defaultdict(list)
        for line in fh:
            if line.startswith("!series_matrix_table_end"):
                break
            p = line.rstrip("\n").split("\t")
            rv = probe2rv.get(p[0].strip('"'))
            if not rv:
                continue
            row = []
            for x in p[1:]:
                try:
                    row.append(float(x))
                except ValueError:
                    row.append(np.nan)
            vals[rv].append(row)
    genes = sorted(vals)
    mat = np.array([np.nanmean(np.array(vals[g], dtype=float), axis=0) for g in genes])
    keep = ~np.all(np.isnan(mat), axis=1)
    genes = [g for g, k in zip(genes, keep) if k]
    mat = mat[keep]
    # les valeurs manquantes résiduelles sont remplacées par la médiane du GÈNE (pas de la lame :
    # ce serait injecter la structure inter-lames que l'on cherche justement à mesurer)
    for i in range(mat.shape[0]):
        m = np.isnan(mat[i])
        if m.any():
            mat[i, m] = np.nanmedian(mat[i]) if not np.all(m) else 0.0
    return genes, mat, samples


# ------------------------------------------------------------------ témoin positif
def operon_pairs(genes: set[str]) -> list[tuple[str, str]]:
    """Paires de gènes adjacents, même brin, < MAX_INTERGENIC pb : proxy d'opéron."""
    feats = []
    for line in GFF.read_text(errors="replace").splitlines():
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        m = re.search(r"locus_tag=([^;]+)", p[8])
        if not m:
            continue
        feats.append((int(p[3]), int(p[4]), p[6], m.group(1)))
    feats.sort()
    out = []
    for (s1, e1, st1, g1), (s2, e2, st2, g2) in zip(feats, feats[1:]):
        if st1 == st2 and (s2 - e1) <= MAX_INTERGENIC and g1 in genes and g2 in genes:
            out.append((g1, g2))
    return out


def bh(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((pvals[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def analyse(tag: str, genes: list[str], mat: np.ndarray, samples: list[str],
            products: dict[str, str]) -> dict[str, dict[str, float]]:
    idx = {g: i for i, g in enumerate(genes)}
    print(f"\n=== {tag} : {len(genes)} gènes x {len(samples)} échantillons ===")
    ranks = np.apply_along_axis(stats.rankdata, 1, mat)       # Spearman = Pearson sur les rangs
    rz = (ranks - ranks.mean(axis=1, keepdims=True))
    rz /= np.linalg.norm(rz, axis=1, keepdims=True) + 1e-12

    rng = random.Random(SEED)
    null = np.array([float(rz[rng.randrange(len(genes))] @ rz[rng.randrange(len(genes))])
                     for _ in range(N_NULL_PAIRS)])
    print(f"nul empirique (200k paires au hasard) : médiane {np.median(null):+.3f}, "
          f"IQR [{np.percentile(null,25):+.3f}, {np.percentile(null,75):+.3f}], "
          f"p95 {np.percentile(null,95):+.3f}, p99.9 {np.percentile(null,99.9):+.3f}")

    op = operon_pairs(set(genes))
    opr = np.array([float(rz[idx[a]] @ rz[idx[b]]) for a, b in op])
    print(f"témoin positif : {len(op)} paires adjacentes même brin <{MAX_INTERGENIC} pb — "
          f"r médian {np.median(opr):+.3f}, quartiles [{np.percentile(opr,25):+.3f}, "
          f"{np.percentile(opr,75):+.3f}] ; percentile de ce médian dans le nul : "
          f"{100*np.mean(null < np.median(opr)):.1f}")

    res: dict[str, dict[str, float]] = {}
    lines = ["focus\tpartenaire\tspearman\tp\tq_BH\tpercentile_nul\tproduit"]
    for f in FOCUS:
        if f not in idx:
            print(f"{f} absent de ce jeu")
            continue
        r = rz @ rz[idx[f]]
        n = mat.shape[1]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = r * np.sqrt((n - 2) / np.clip(1 - r ** 2, 1e-12, None))
        p = 2 * stats.t.sf(np.abs(t), n - 2)
        p[idx[f]] = 1.0
        q = bh(p)
        order = np.argsort(-r)
        res[f] = {genes[i]: float(r[i]) for i in range(len(genes))}
        print(f"\n{f} : r(Rv2517c)={r[idx['Rv2517c']] if 'Rv2517c' in idx else float('nan'):+.3f}"
              f"  |  top 12 partenaires :")
        for i in order[:12]:
            if genes[i] == f:
                continue
            pct = 100 * float(np.mean(null < r[i]))
            print(f"    {genes[i]:9s} r={r[i]:+.3f} q={q[i]:.2e} percentile_nul={pct:6.3f}  "
                  f"{products.get(genes[i],'')[:52]}")
        for i in order:
            if genes[i] == f:
                continue
            lines.append(f"{f}\t{genes[i]}\t{r[i]:.4f}\t{p[i]:.3e}\t{q[i]:.3e}\t"
                         f"{100*float(np.mean(null < r[i])):.3f}\t{products.get(genes[i],'')[:60]}")
        n_sig = int(np.sum(q < 0.05))
        print(f"    -> {n_sig} gènes à q<0,05 ; percentile de Rv2517c dans le nul : "
              f"{100*float(np.mean(null < r[idx['Rv2517c']])):.2f}" if "Rv2517c" in idx else "")
    (OUT / f"{tag}_correlations.tsv").write_text("\n".join(lines) + "\n")
    return res


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(DB)
    products = {rv: (p or "") for rv, p in
                con.execute("select rv, product_h37rv from genes")}

    g1, m1, s1 = load_gse166501()
    r1 = analyse("gse166501", g1, m1, s1, products)
    g2, m2, s2 = load_gse71200()
    r2 = analyse("gse71200", g2, m2, s2, products)

    # ---------- convergence entre deux jeux sans rapport ------------------------------------
    lines = ["focus\tgene\tr_gse166501\tr_gse71200\trang_166501\trang_71200\tproduit"]
    print("\n=== CONVERGENCE ENTRE LES DEUX JEUX ===")
    for f in FOCUS:
        if f not in r1 or f not in r2:
            continue
        common = sorted(set(r1[f]) & set(r2[f]))
        o1 = {g: i for i, g in enumerate(sorted(common, key=lambda x: -r1[f][x]))}
        o2 = {g: i for i, g in enumerate(sorted(common, key=lambda x: -r2[f][x]))}
        both = sorted(common, key=lambda g: o1[g] + o2[g])
        rho = stats.spearmanr([r1[f][g] for g in common], [r2[f][g] for g in common])
        print(f"\n{f} : {len(common)} gènes communs ; corrélation des deux profils de "
              f"co-expression rho={rho.statistic:+.3f} (p={rho.pvalue:.1e})")
        print("  meilleurs par rang CUMULÉ dans les deux jeux :")
        for g in both[:15]:
            if g == f:
                continue
            print(f"    {g:9s} r1={r1[f][g]:+.3f} (rang {o1[g]:4d})  "
                  f"r2={r2[f][g]:+.3f} (rang {o2[g]:4d})  {products.get(g,'')[:46]}")
            lines.append(f"{f}\t{g}\t{r1[f][g]:.4f}\t{r2[f][g]:.4f}\t{o1[g]}\t{o2[g]}\t"
                         f"{products.get(g,'')[:60]}")
        if "Rv2517c" in o1:
            print(f"  Rv2517c : rang {o1['Rv2517c']} / {len(common)} dans GSE166501, "
                  f"rang {o2['Rv2517c']} dans GSE71200")
    (OUT / "convergence.tsv").write_text("\n".join(lines) + "\n")
    print(f"\nsorties dans {OUT}")


if __name__ == "__main__":
    main()
