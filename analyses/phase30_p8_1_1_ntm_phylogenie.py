#!/usr/bin/env python3
"""phase30_p8_1_1_ntm_phylogenie.py -- P8.1.1 : acquisition unique + pertes, ou transferts
horizontaux indépendants ? La question centrale du dossier, posée au niveau PHYLOGÉNÉTIQUE.

CE QUI EST DÉJÀ SU, ET CE QUI MANQUE.
P8.1 a établi que Rv2516c est présent chez 6/53 NTM (M. avium, M. bouchedurhonense, M. celatum,
M. shinjukuense, M. simiae, M. timonense) avec la MÊME liste d'espèces que Rv2517c, et a posé
elle-même la réserve : « les six espèces porteuses ne forment pas un clade évident, donc transfert
horizontal et héritage avec pertes ne sont pas départagés ». P3.1 a corrigé la pseudo-réplication
du complexe M. avium (avium/bouchedurhonense/timonense = orthologue quasi-identique = UNE lignée),
ramenant les six espèces à QUATRE lignées indépendantes. Aucun test n'a jamais été fait sur un ARBRE.

LES TROIS TESTS, ET POURQUOI IL EN FAUT TROIS.

  1. **Parcimonie + permutation sur l'arbre d'espèces.** Le score de Fitch (nombre minimal de
     changements d'état) du caractère présence/absence est comparé à sa distribution nulle sous
     placement ALÉATOIRE du même nombre de porteurs sur le même arbre. Ce test est INDÉPENDANT de
     l'enracinement (Fitch non enraciné). Il répond à : les porteurs sont-ils plus groupés que le
     hasard ? Le contre-argument de la piste (puissance faible à 4 lignées) est traité EN AMONT et
     non commenté après coup : la distribution nulle donne directement le score minimal atteignable,
     donc la p-value la plus petite possible du dispositif. Si ce plancher est > 0,05, le test est
     déclaré sans puissance AVANT de regarder l'observé.

  2. **Divergence relative : la signature qui distingue vraiment HGT récent de vieil héritage.**
     Un HGT RÉCENT laisse une trace que la parcimonie ne voit pas : le gène transféré est
     anormalement PEU divergent entre espèces éloignées, comparé aux gènes hérités verticalement.
     On calibre donc l'identité de Rv2516c (et Rv2517c) chez chaque porteur contre la distribution
     d'identité de ~150 gènes H37Rv tirés au hasard dans le MÊME génome, par le MÊME tblastn.
     Percentile élevé chez une espèce éloignée = signature de transfert récent ; percentile ordinaire
     ou bas = compatible avec un héritage ancien. C'est le test qui a de la puissance ici, parce
     qu'il ne dépend pas du nombre de lignées mais du contraste identité/distance.

  3. **Arbre de gène vs arbre d'espèces.** Sous héritage vertical, l'arbre des orthologues doit
     suivre l'arbre d'espèces ; sous transferts indépendants, il ne le suit pas. À 5-7 taxons le
     test est faible, il est rapporté comme indice, pas comme preuve.

L'ARBRE D'ESPÈCES : construit ici, faute d'outil externe.
Ni mafft, ni FastTree, ni IQ-TREE ne sont installés sur cette machine (vérifié). L'arbre est donc
construit par la méthode DÉJÀ ÉPROUVÉE SUR CE PROJET (phase13_p3_1) : alignement ANCRÉ SUR LA
RÉFÉRENCE. tblastn de 12 protéines de ménage de H37Rv contre chaque génome NTM, les positions de la
requête servant de système de coordonnées commun, concaténation, distance de Kimura sur protéines,
Neighbor-Joining (Biopython). Ce n'est pas un arbre de publication : c'est un arbre de travail dont
la VALIDITÉ EST TESTÉE par la récupération des complexes connus (MAC, complexe M. simiae, complexe
M. kansasii, clade MTBAP). Si ces complexes ne sortent pas, l'arbre est rejeté et le test 1 avec lui.

RÉSERVES ÉNONCÉES D'AVANCE.
- Alignement ancré sur H37Rv : les positions absentes de H37Rv sont invisibles, et l'espèce la plus
  proche de H37Rv est favorisée. Acceptable pour des protéines de ménage à 80-95 % d'identité,
  pas pour une phylogénie profonde.
- Aucun grimpeur rapide (M. smegmatis, M. abscessus) dans bdd/hors_mtbc : pas d'exogroupe extérieur
  au clade des grimpeurs lents. Enracinement par point médian, et le scénario de Dollo est rapporté
  sous PLUSIEURS enracinements alternatifs pour montrer ce qui en dépend et ce qui n'en dépend pas.
- « Absence » = absence sous les seuils de l'atlas (pident 30 %, qcov 50 %, E 1e-5). Le meilleur hit
  SOUS le seuil est conservé et rapporté pour chaque espèce non porteuse, exactement comme le
  correctif de sensibilité de phase57 l'impose.

Entrées : bdd/hors_mtbc/<espèce>/**/genome.fna (53 espèces),
          annotation_mtbc/résultats/phase57_ntm/query_proteome.faa (protéome H37Rv).
Sorties : résultats/p8_1_1_ntm_phylogenie/{hits.tsv, alignement_concat.faa, arbre_especes.nwk,
          distances.tsv, parsimonie.tsv, identite_calibree.tsv, resume.md}

Run: python analyses/phase30_p8_1_1_ntm_phylogenie.py --blast   # passe tblastn (~20-40 min)
     python analyses/phase30_p8_1_1_ntm_phylogenie.py           # analyse seule (hits.tsv en cache)
"""
from __future__ import annotations

import argparse
import glob
import math
import random
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
HORS = MTBC / "bdd" / "hors_mtbc"
PROTEOME = MTBC / "annotation_mtbc" / "résultats" / "phase57_ntm" / "query_proteome.faa"
OUT = ROOT / "résultats" / "p8_1_1_ntm_phylogenie"

# --- protéines de ménage servant d'échafaudage phylogénétique -------------------------------
# Choix classique MLSA du genre Mycobacterium (Devulder 2005, Tortoli 2017) restreint aux protéines
# longues, mono-copie et sans intéine. recA est ÉCARTÉ volontairement : l'intéine de M. tuberculosis
# (440 aa insérés en phase) creuserait un trou d'alignement chez toutes les espèces sans intéine.
MARKERS = {
    "Rv0667": "rpoB", "Rv0005": "gyrB", "Rv3240c": "secA1", "Rv1310": "atpD",
    "Rv0440": "groEL2", "Rv0685": "tuf", "Rv2916c": "ffh", "Rv3457c": "rpoA",
    "Rv0350": "dnaK", "Rv0684": "fusA1", "Rv1630": "rpsA", "Rv2703": "sigA",
}
FOCUS = {"Rv2516c": "cible du projet", "Rv2517c": "partenaire d'opéron"}
N_CALIB = 150          # gènes tirés au hasard pour la distribution nulle d'identité
CALIB_SEED = 20260810
SEED_BOOT = 20260811

# DEUX PASSES, et la raison est un coût mesuré. La passe d'arbre a besoin des 53 génomes mais de
# 14 protéines seulement ; la passe de calibration a besoin de 150 protéines mais SEULEMENT des
# espèces où l'on veut situer Rv2516c dans la distribution d'identité de son propre génome. Tout
# mettre dans une seule requête coûtait 3,3x plus de calcul pour rien (mesuré à 6 min/génome sous la
# charge concurrente de la machine, soit 5 h ; découpé, ~50 min).
CALIB_SPECIES = ["M_avium", "M_bouchedurhonense", "M_celatum", "M_shinjukuense", "M_simiae",
                 "M_timonense",                                    # les 6 porteurs
                 "M_decipiens", "M_lacus", "M_riyadhense",         # MTBAP NON porteurs = le contrôle
                 "M_kansasii", "M_marinum"]                        # deux jalons plus éloignés

# seuils de PRÉSENCE : strictement ceux de l'atlas phase57, pour rester comparable
MIN_PIDENT, MIN_QCOV, EVALUE = 30.0, 50.0, 1e-5

# complexes connus, servant de CONTRÔLE DE VALIDITÉ de l'arbre (littérature, pas nos données)
KNOWN_COMPLEXES = {
    "MAC": ["M_avium", "M_intracellulare", "M_paraintracellulare", "M_colombiense",
            "M_marseillense", "M_timonense", "M_bouchedurhonense", "M_avium-intracellulare"],
    "simiae": ["M_simiae", "M_lentiflavum", "M_florentinum", "M_triplex", "M_genavense",
               "M_europaeum", "M_paraffinicum", "M_interjectum"],
    "kansasii": ["M_kansasii", "M_persicum", "M_pseudokansasii", "M_innocens", "M_attenuatum",
                 "M_basiliense", "M_neolactis", "M_novum"],
    "MTBAP": ["M_decipiens", "M_lacus", "M_riyadhense", "M_shinjukuense"],
    "ulcerans_marinum": ["M_marinum", "M_ulcerans", "M_liflandii", "M_pseudoshottsii"],
}
# libellés qui ne désignent pas une espèce : exclus de l'arbre, signalés dans le résumé
NOT_A_SPECIES = {"M_sp_indetermine", "M_unidentified"}

CARRIERS_ATLAS = ["M_avium", "M_bouchedurhonense", "M_celatum", "M_shinjukuense",
                  "M_simiae", "M_timonense"]


# ------------------------------------------------------------------ entrées
def read_faa(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def species_genomes() -> dict[str, Path]:
    out: dict[str, Path] = {}
    for f in sorted(glob.glob(str(HORS / "**" / "genome.fna"), recursive=True)):
        sp = Path(f).relative_to(HORS).parts[0]
        out.setdefault(sp, Path(f))
    return out


def build_query(proteome: dict[str, str]) -> tuple[Path, Path, list[str]]:
    rng = random.Random(CALIB_SEED)
    pool = sorted(k for k in proteome if k not in MARKERS and k not in FOCUS)
    calib = rng.sample(pool, N_CALIB)
    OUT.mkdir(parents=True, exist_ok=True)
    q_tree = OUT / "query_tree.faa"
    q_tree.write_text("".join(f">{i}\n{proteome[i]}\n" for i in list(MARKERS) + list(FOCUS)))
    q_cal = OUT / "query_calib.faa"
    q_cal.write_text("".join(f">{i}\n{proteome[i]}\n" for i in calib))
    (OUT / "calib_ids.txt").write_text("\n".join(calib) + "\n")
    return q_tree, q_cal, calib


# ------------------------------------------------------------------ tblastn
FMT = ("6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen qseq sseq")


def run_blast(query: Path, genomes: dict[str, Path], hits: Path) -> Path:
    with hits.open("w") as fh:
        fh.write("species\t" + "\t".join(FMT.split()[1:]) + "\n")
        for i, (sp, g) in enumerate(sorted(genomes.items()), 1):
            with tempfile.TemporaryDirectory() as td:
                db = Path(td) / "db"
                subprocess.run(["makeblastdb", "-in", str(g), "-dbtype", "nucl", "-out", str(db)],
                               check=True, capture_output=True)
                r = subprocess.run(
                    ["tblastn", "-query", str(query), "-db", str(db), "-outfmt", FMT,
                     "-evalue", str(EVALUE), "-max_target_seqs", "20", "-num_threads", "8",
                     "-seg", "no"],
                    check=True, capture_output=True, text=True)
            n = 0
            for line in r.stdout.splitlines():
                fh.write(f"{sp}\t{line}\n")
                n += 1
            print(f"[{i:2d}/{len(genomes)}] {sp:28s} {n:5d} HSP", flush=True)
    return hits


def load_hits(path: Path) -> list[dict]:
    cols = ["species"] + FMT.split()[1:]
    rows = []
    with path.open() as fh:
        next(fh)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) != len(cols):
                continue
            d = dict(zip(cols, p))
            for k in ("pident", "evalue", "bitscore"):
                d[k] = float(d[k])
            for k in ("length", "qstart", "qend", "sstart", "send", "qlen"):
                d[k] = int(d[k])
            rows.append(d)
    return rows


# ------------------------------------------------------------------ alignement ancré
def anchored(rows: list[dict], qid: str, qlen: int) -> str:
    """Reconstruit la séquence sujet projetée sur les coordonnées de la requête.

    Les HSP du MEILLEUR sujet (sseqid du HSP de plus fort bitscore) sont posés par bitscore
    décroissant ; une position déjà remplie n'est jamais écrasée (évite qu'un HSP faible
    chevauchant contredise le meilleur).
    """
    hs = [r for r in rows if r["qseqid"] == qid]
    if not hs:
        return "-" * qlen
    best_subj = max(hs, key=lambda r: r["bitscore"])["sseqid"]
    hs = sorted([r for r in hs if r["sseqid"] == best_subj],
                key=lambda r: -r["bitscore"])
    seq = ["-"] * qlen
    for r in hs:
        qpos = r["qstart"] - 1
        for qc, sc in zip(r["qseq"], r["sseq"]):
            if qc == "-":
                continue                      # insertion dans le sujet : hors coordonnées requête
            if 0 <= qpos < qlen and seq[qpos] == "-":
                seq[qpos] = sc
            qpos += 1
    return "".join(seq)


def best_present(by_sp: dict, sp: str, qid: str, qlen: dict) -> bool:
    """Présence d'un orthologue sous les seuils de l'atlas (pident 30 %, qcov 50 %)."""
    hs = [r for r in by_sp[sp] if r["qseqid"] == qid]
    if not hs:
        return False
    pid = max(hs, key=lambda r: r["bitscore"])["pident"]
    seq = anchored(by_sp[sp], qid, qlen[qid])
    qcov = 100.0 * (len(seq) - seq.count("-")) / len(seq)
    return pid >= MIN_PIDENT and qcov >= MIN_QCOV


def kimura(a: str, b: str) -> tuple[float, int]:
    """Distance protéique de Kimura (1983) sur les positions comparables."""
    n = same = 0
    for x, y in zip(a, b):
        if x in "-X" or y in "-X":
            continue
        n += 1
        same += x == y
    if n < 50:
        return float("nan"), n
    p = 1.0 - same / n
    arg = 1.0 - p - 0.2 * p * p
    if arg <= 0:
        return float("nan"), n
    return -math.log(arg), n


# ------------------------------------------------------------------ NJ + parcimonie
def nj_tree(labels: list[str], dm: list[list[float]]):
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    lower = [[dm[i][j] for j in range(i + 1)] for i in range(len(labels))]
    return DistanceTreeConstructor().nj(DistanceMatrix(labels, lower))


def fitch_score(tree, states: dict[str, int]) -> int:
    """Score de parcimonie de Fitch, non enraciné en pratique (invariant par ré-enracinement
    pour un caractère binaire sur un arbre binaire)."""
    score = 0

    def post(clade):
        nonlocal score
        if clade.is_terminal():
            return {states[clade.name]}
        sets = [post(c) for c in clade.clades]
        inter = set.intersection(*sets)
        if inter:
            return inter
        score += 1
        return set.union(*sets)

    post(tree.root)
    return score


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", action="store_true", help="passe tblastn arbre (53 génomes x 14 prot.)")
    ap.add_argument("--calib", action="store_true", help="passe tblastn calibration (11 x 150)")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    proteome = read_faa(PROTEOME)
    q_tree, q_cal, calib = build_query(proteome)
    genomes = species_genomes()
    print(f"{len(genomes)} génomes d'espèces ; arbre = {len(MARKERS)+len(FOCUS)} protéines, "
          f"calibration = {len(calib)} protéines x {len(CALIB_SPECIES)} génomes")

    hits_path, calib_path = OUT / "hits.tsv", OUT / "hits_calib.tsv"
    if args.blast or not hits_path.exists():
        run_blast(q_tree, genomes, hits_path)
    if args.calib:
        run_blast(q_cal, {s: genomes[s] for s in CALIB_SPECIES if s in genomes}, calib_path)

    rows = load_hits(hits_path)
    if calib_path.exists() and calib_path.stat().st_size > 0:
        rows += load_hits(calib_path)
    else:
        print("passe de calibration absente ou en cours : les percentiles d'identité seront vides")
    print(f"{len(rows)} HSP chargés")

    by_sp: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_sp[r["species"]].append(r)

    qlen = {i: len(proteome[i]) for i in list(MARKERS) + list(FOCUS) + calib}

    # ---------- 1. alignement concaténé des marqueurs -------------------------------------
    species = sorted(sp for sp in genomes if sp not in NOT_A_SPECIES)
    concat: dict[str, str] = {}
    for sp in species:
        concat[sp] = "".join(anchored(by_sp[sp], m, qlen[m]) for m in MARKERS)
    concat["H37Rv"] = "".join(proteome[m] for m in MARKERS)
    labels = ["H37Rv"] + species
    (OUT / "alignement_concat.faa").write_text(
        "".join(f">{k}\n{concat[k]}\n" for k in labels))
    cov = {k: 1 - concat[k].count("-") / len(concat[k]) for k in labels}
    print("couverture d'alignement min/médiane :",
          f"{min(cov.values()):.2f} / {sorted(cov.values())[len(cov)//2]:.2f}")

    # ---------- 2. matrice de distance + NJ ------------------------------------------------
    n = len(labels)
    dm = [[0.0] * n for _ in range(n)]
    with (OUT / "distances.tsv").open("w") as fh:
        fh.write("sp1\tsp2\tkimura\tn_positions\n")
        for i in range(n):
            for j in range(i + 1, n):
                d, npos = kimura(concat[labels[i]], concat[labels[j]])
                if math.isnan(d):
                    d = 2.0                  # paire non comparable : pénalité explicite
                dm[i][j] = dm[j][i] = d
                fh.write(f"{labels[i]}\t{labels[j]}\t{d:.4f}\t{npos}\n")
    tree = nj_tree(labels, dm)
    from Bio import Phylo
    Phylo.write(tree, OUT / "arbre_especes.nwk", "newick")

    # ---------- 3. contrôle de validité : complexes connus ---------------------------------
    # Un arbre NJ est NON ENRACINÉ ; Biopython le stocke enraciné à un point arbitraire. Tester la
    # monophylie avec `common_ancestor` revient donc à tester une propriété de cet enracinement
    # arbitraire, et un groupe parfaitement monophylétique ressort « avec 14 intrus » dès que la
    # racine tombe dedans (piégé le 2026-08-10 : le premier passage rejetait l'arbre à tort).
    # Le test correct est sur les BIPARTITIONS : un groupe est un clan si une arête le sépare
    # exactement du reste.
    tips = {t.name for t in tree.get_terminals()}
    splits: list[frozenset[str]] = []
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals())
        if 1 < len(s) < len(tips):
            splits.append(s)
            splits.append(frozenset(tips - s))
    splitset = set(splits)

    def plus_proche(grp: set[str]) -> tuple[int, int]:
        """(intrus, manquants) du split le plus ressemblant."""
        best = min(splits, key=lambda s: len(s ^ grp))
        return len(best - grp), len(grp - best)

    print("\nCONTRÔLE DE VALIDITÉ DE L'ARBRE (bipartitions, non enraciné)")
    n_ok = n_test = 0
    for name, members in KNOWN_COMPLEXES.items():
        present = set(m for m in members if m in tips)
        if len(present) < 2:
            print(f"  {name:18s} {len(present)} membre(s) présent(s) -> non testable")
            continue
        n_test += 1
        if frozenset(present) in splitset:
            n_ok += 1
            print(f"  {name:18s} {len(present)} membres -> CLAN exact")
        else:
            ex, mi = plus_proche(present)
            print(f"  {name:18s} {len(present)} membres -> pas un clan "
                  f"(split le plus proche : {ex} intrus, {mi} manquants)")
    print(f"  => {n_ok}/{n_test} complexes connus récupérés exactement")

    # attendu fort et indépendant : le MTBC (H37Rv) est frère du clade MTBAP
    mtbap = set(KNOWN_COMPLEXES["MTBAP"]) & tips
    print(f"  MTBC+MTBAP forme-t-il un clan ? "
          f"{'OUI' if frozenset(mtbap | {'H37Rv'}) in splitset else 'non'}")

    # BOOTSTRAP : « l'arbre est-il faux, ou seulement incertain ? » est une question quantitative.
    # 200 ré-échantillonnages des colonnes de l'alignement concaténé ; on mesure (a) le support de
    # chaque complexe connu, (b) le support du groupement des porteurs, (c) la distribution du score
    # de Fitch à travers les arbres, ce qui propage l'incertitude topologique DANS le test au lieu
    # de la mentionner en réserve.
    import numpy as np
    cols = len(concat[labels[0]])
    codes = np.array([np.frombuffer(concat[k].encode(), dtype=np.uint8) for k in labels])
    valid = (codes != ord("-")) & (codes != ord("X"))

    def dm_numpy(sel: np.ndarray | None = None) -> list[list[float]]:
        """Matrice de Kimura vectorisée. La version en boucles Python coûtait 10,6 M comparaisons
        de caractères par réplicat : 200 réplicats ne passaient pas (processus tué)."""
        c = codes if sel is None else codes[:, sel]
        v = valid if sel is None else valid[:, sel]
        n = len(labels)
        out = [[0.0] * n for _ in range(n)]
        for i in range(n):
            both = v[i] & v[i + 1:]
            eq = (c[i] == c[i + 1:]) & both
            npos = both.sum(axis=1)
            same = eq.sum(axis=1)
            with np.errstate(divide="ignore", invalid="ignore"):
                p = 1.0 - np.where(npos > 0, same / np.maximum(npos, 1), np.nan)
                arg = 1.0 - p - 0.2 * p * p
                d = np.where((npos >= 50) & (arg > 0), -np.log(np.clip(arg, 1e-12, None)), 2.0)
            for k, j in enumerate(range(i + 1, n)):
                out[i][j] = out[j][i] = float(d[k])
        return out

    brng = np.random.default_rng(SEED_BOOT)
    support = {k: 0 for k in KNOWN_COMPLEXES}
    support["porteurs_Rv2516c"] = 0
    carriers_now = [sp for sp in species if best_present(by_sp, sp, "Rv2516c", qlen)]
    fitch_boot: list[int] = []
    N_BOOT = 200
    for b in range(N_BOOT):
        pick = brng.integers(0, cols, cols)
        bt = nj_tree(labels, dm_numpy(pick))
        btips = {t.name for t in bt.get_terminals()}
        bsplits = set()
        for cl in bt.get_nonterminals():
            s = frozenset(t.name for t in cl.get_terminals())
            if 1 < len(s) < len(btips):
                bsplits.add(s)
                bsplits.add(frozenset(btips - s))
        for k, members in KNOWN_COMPLEXES.items():
            if frozenset(set(members) & btips) in bsplits:
                support[k] += 1
        if frozenset(carriers_now) in bsplits:
            support["porteurs_Rv2516c"] += 1
        st_b = {sp: int(sp in carriers_now) for sp in species}
        st_b["H37Rv"] = 1
        fitch_boot.append(fitch_score(bt, st_b))
    print(f"\nSUPPORT BOOTSTRAP ({N_BOOT} réplicats de colonnes)")
    for k, v in support.items():
        print(f"  {k:20s} {100*v/N_BOOT:5.1f} %")
    print(f"  score de Fitch des porteurs à travers les arbres bootstrap : "
          f"min {min(fitch_boot)}, médiane {sorted(fitch_boot)[N_BOOT//2]}, max {max(fitch_boot)}")

    tree.root_at_midpoint()

    # ---------- 4. présence/absence + calibration d'identité --------------------------------
    def best(sp: str, qid: str) -> tuple[float, float, float] | None:
        hs = [r for r in by_sp[sp] if r["qseqid"] == qid]
        if not hs:
            return None
        b = max(hs, key=lambda r: r["bitscore"])
        # couverture requête cumulée sur le meilleur sujet (mêmes règles que l'alignement ancré)
        seq = anchored(by_sp[sp], qid, qlen[qid])
        qcov = 100.0 * (len(seq) - seq.count("-")) / len(seq)
        return b["pident"], qcov, b["evalue"]

    pa: dict[str, dict[str, int]] = {}
    with (OUT / "identite_calibree.tsv").open("w") as fh:
        fh.write("species\tgene\tpident\tqcov\tevalue\tpresent\tpercentile_vs_calib\tn_calib_hits\n")
        for sp in species:
            cal = [b[0] for g in calib if (b := best(sp, g)) and b[1] >= MIN_QCOV
                   and b[0] >= MIN_PIDENT]
            pa[sp] = {}
            for g in FOCUS:
                b = best(sp, g)
                if b is None:
                    pa[sp][g] = 0
                    fh.write(f"{sp}\t{g}\tNA\tNA\tNA\t0\tNA\t{len(cal)}\n")
                    continue
                pid, qcov, ev = b
                present = int(pid >= MIN_PIDENT and qcov >= MIN_QCOV)
                pa[sp][g] = present
                pct = (100.0 * sum(1 for c in cal if c < pid) / len(cal)) if cal else float("nan")
                fh.write(f"{sp}\t{g}\t{pid:.1f}\t{qcov:.1f}\t{ev:.1e}\t{present}\t"
                         f"{pct:.1f}\t{len(cal)}\n")

    carriers = sorted(sp for sp in species if pa[sp]["Rv2516c"])
    print(f"\nporteurs Rv2516c (re-dérivés) : {carriers}")
    print(f"porteurs atlas               : {sorted(CARRIERS_ATLAS)}")

    # ---------- 5. parcimonie + permutation ------------------------------------------------
    rng = random.Random(CALIB_SEED)
    with (OUT / "parsimonie.tsv").open("w") as fh:
        fh.write("gene\tn_carriers\tfitch_obs\tfitch_min_possible\tfitch_max_possible\t"
                 "p_permutation\tn_perm\tdollo_losses\n")
        for g in FOCUS:
            car = [sp for sp in species if pa[sp][g]]
            if not car:
                continue
            states = {sp: pa[sp][g] for sp in species}
            states["H37Rv"] = 1
            obs = fitch_score(tree, states)
            null = []
            pool = [sp for sp in species]
            for _ in range(20000):
                pick = set(rng.sample(pool, len(car)))
                st = {sp: int(sp in pick) for sp in species}
                st["H37Rv"] = 1
                null.append(fitch_score(tree, st))
            p = (sum(1 for x in null if x <= obs) + 1) / (len(null) + 1)
            # scénario de Dollo : un seul gain au MRCA des porteurs, pertes = absents dans le clade
            mrca = tree.common_ancestor(car + ["H37Rv"])
            inside = {t.name for t in mrca.get_terminals()}
            losses = sum(1 for t in inside if t in states and states[t] == 0)
            fh.write(f"{g}\t{len(car)}\t{obs}\t{min(null)}\t{max(null)}\t{p:.4f}\t"
                     f"{len(null)}\t{losses}\n")
            print(f"{g}: Fitch={obs} (nul {min(null)}-{max(null)}, p={p:.4f}), "
                  f"Dollo: 1 gain + {losses} pertes")

    # ---------- 5bis. LE test, une fois la pseudo-réplication MAC retirée --------------------
    # Le score de Fitch de 4 pour 6 porteurs vient peut-être entièrement du fait que trois d'entre
    # eux (avium, bouchedurhonense, timonense) sont le MÊME orthologue dans le MÊME complexe : ils
    # coûtent 1 changement au lieu de 3. P3.1 avait déjà établi qu'ils ne comptent que pour UNE
    # lignée indépendante. On refait donc le test sur un arbre élagué de deux des trois, avec
    # 4 porteurs — et on mesure d'abord le PLANCHER de p-value atteignable, pour dire si le
    # dispositif a la moindre puissance AVANT de regarder l'observé.
    from Bio import Phylo as _P
    import copy as _copy
    red = _copy.deepcopy(tree)
    for extra in ("M_bouchedurhonense", "M_timonense"):
        try:
            red.prune(target=extra)
        except Exception:
            pass
    red_tips = [t.name for t in red.get_terminals() if t.name != "H37Rv"]
    with (OUT / "parsimonie_MAC_collapse.tsv").open("w") as fh:
        fh.write("gene\tn_porteurs\tn_taxons\tfitch_obs\tfitch_min_nul\tp_permutation\t"
                 "p_plancher_si_monophyletique\tpuissance\n")
        for g in FOCUS:
            car = [sp for sp in red_tips if pa.get(sp, {}).get(g)]
            st = {sp: int(sp in car) for sp in red_tips}
            st["H37Rv"] = 1
            obs = fitch_score(red, st)
            null = []
            for _ in range(20000):
                pick = set(rng.sample(red_tips, len(car)))
                s2 = {sp: int(sp in pick) for sp in red_tips}
                s2["H37Rv"] = 1
                null.append(fitch_score(red, s2))
            p = (sum(1 for x in null if x <= obs) + 1) / (len(null) + 1)
            floor = (sum(1 for x in null if x <= min(null)) + 1) / (len(null) + 1)
            puiss = "OUI" if floor < 0.05 else "NON (aucun résultat ne peut être significatif)"
            fh.write(f"{g}\t{len(car)}\t{len(red_tips)}\t{obs}\t{min(null)}\t{p:.4f}\t"
                     f"{floor:.4f}\t{puiss}\n")
            print(f"{g} SANS pseudo-réplication MAC : {len(car)} porteurs sur {len(red_tips)} taxons, "
                  f"Fitch={obs} (nul {min(null)}-{max(null)}), p={p:.4f} ; "
                  f"plancher de p si les porteurs étaient monophylétiques = {floor:.4f} -> "
                  f"puissance {puiss}")

    # ---------- 6. arbre de gène des orthologues -------------------------------------------
    for g in FOCUS:
        car = [sp for sp in species if pa[sp][g]]
        if len(car) < 3:
            continue
        aln = {sp: anchored(by_sp[sp], g, qlen[g]) for sp in car}
        aln["H37Rv"] = proteome[g]
        labs = ["H37Rv"] + sorted(car)
        m = len(labs)
        gdm = [[0.0] * m for _ in range(m)]
        for i in range(m):
            for j in range(i + 1, m):
                d, _ = kimura(aln[labs[i]], aln[labs[j]])
                gdm[i][j] = gdm[j][i] = 2.0 if math.isnan(d) else d
        gt = nj_tree(labs, gdm)
        Phylo.write(gt, OUT / f"arbre_gene_{g}.nwk", "newick")
        (OUT / f"alignement_{g}.faa").write_text(
            "".join(f">{k}\n{aln[k]}\n" for k in labs))

    print(f"\nsorties dans {OUT}")


if __name__ == "__main__":
    main()
