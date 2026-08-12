#!/usr/bin/env python3
"""P10.3.1 -- test d'enrichissement FORMEL du voisinage de co-expression de Rv2516c
contre des regulons de stress/persistance curés, plutot qu'une lecture qualitative
des rangs (P10.3).

Regulons :
  - DosR   : Rustad et al. 2008 (PLOS ONE), Table S1, 49 genes -- fichier deja sur
             disque dans le projet frere Rv2699c (data/dosr_regulon/rustad2008_tableS1_dosr.tsv).
  - SigH   : reseau TFOE Minch/Rustad (ISB MTB Network Portal), genes regules par
             Rv3223c (sigH) dans le champ `regulation` de l'atlas annotation_mtbc.
  - SigF   : meme reseau, genes regules par Rv3286c (sigF).
  - PERSIST: meme reseau, union des genes regules par l'un des cinq regulateurs
             connus de Rv2516c lui-meme (Rv0081, Rv0324, ramB/Rv0465c, trcR/Rv1033c,
             lsr2/Rv3597c) -- teste si le voisinage de co-expression partage un
             regulateur avec Rv2516c, plutot qu'un regulon publie generique.

Univers et listes-requete : les tables de correlation deja calculees par P10.3
(resultats/p10_3_coexpression/{gse}_correlations.tsv), focus=Rv2516c, top-50 et
top-100 par |spearman| decroissant (deja trie par construction de P10.3, on retrie
par surete). Univers = tous les partenaires presents dans la table de CE jeu (le
gene set est restreint a cet univers avant le test).

Test : hypergeometrique (scipy.stats.hypergeom), correction Benjamini-Hochberg sur
l'ensemble des tests (2 jeux x 2 tailles de liste x 4 regulons = 16 tests).
Lecture disciplinee : un regulon n'est retenu comme resultat que s'il est
significatif (q < 0.05) DANS LES DEUX JEUX independamment, meme discipline que le
recouvrement top-100 de P10.3.
"""
import csv
import json
import os
import sqlite3

from scipy.stats import hypergeom

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
COEXPR_DIR = os.path.join(ROOT, "résultats", "p10_3_coexpression")
OUT_DIR = os.path.join(ROOT, "résultats", "p10_3_1_enrichment")
os.makedirs(OUT_DIR, exist_ok=True)

ATLAS_DB = os.path.expanduser(
    "~/docs/codes/mtbc/annotation_mtbc/site/data/db.sqlite"
)
RUSTAD_DOSR = os.path.expanduser(
    "~/docs/codes/mtbc/Rv2699c/data/dosr_regulon/rustad2008_tableS1_dosr.tsv"
)

FOCUS = "Rv2516c"
DATASETS = {
    "GSE166501": os.path.join(COEXPR_DIR, "gse166501_correlations.tsv"),
    "GSE71200": os.path.join(COEXPR_DIR, "gse71200_correlations.tsv"),
}
TOP_NS = [50, 100]

RV2516C_REGULATORS = {
    "Rv0081": None,
    "Rv0324": None,
    "Rv0465c": "ramB",
    "Rv1033c": "trcR",
    "Rv3597c": "lsr2",
}


def load_dosr_regulon():
    genes = set()
    with open(RUSTAD_DOSR) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            rv = row.get("Rv#") or row.get("Rv")
            if rv:
                genes.add(rv.strip())
    return genes


def load_atlas_regulon(tf_rv_list):
    """Genes dont `regulation.regulated_by` contient un TF de tf_rv_list."""
    db = sqlite3.connect(ATLAS_DB)
    cur = db.cursor()
    cur.execute("SELECT rv, regulation FROM genes WHERE regulation IS NOT NULL")
    members = set()
    for rv, reg_json in cur.fetchall():
        if not reg_json:
            continue
        try:
            reg = json.loads(reg_json)
        except (json.JSONDecodeError, TypeError):
            continue
        for edge in reg.get("regulated_by", []) or []:
            if edge.get("tf") in tf_rv_list:
                members.add(rv)
    db.close()
    return members


def load_dataset(path):
    rows = []
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row["focus"] != FOCUS:
                continue
            rows.append(
                {
                    "partenaire": row["partenaire"],
                    "spearman": float(row["spearman"]),
                }
            )
    rows.sort(key=lambda x: -abs(x["spearman"]))
    return rows


def hypergeom_test(universe_size, set_size, sample_size, overlap):
    """P(X >= overlap) sous hypergeometrique -- test unilateral d'ENRICHISSEMENT."""
    if overlap == 0:
        return 1.0
    # sf(k-1) = P(X >= k)
    return float(hypergeom.sf(overlap - 1, universe_size, set_size, sample_size))


def bh_correct(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    running_min = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        val = pvals[idx] * n / rank
        running_min = min(running_min, val)
        q[idx] = running_min
    return q


def main():
    print("Chargement des regulons...")
    regulons_raw = {
        "DosR (Rustad 2008, 49g)": load_dosr_regulon(),
        "SigH (TFOE Minch/Rustad)": load_atlas_regulon({"Rv3223c"}),
        "SigF (TFOE Minch/Rustad)": load_atlas_regulon({"Rv3286c"}),
        "PERSIST (5 regulateurs connus de Rv2516c)": load_atlas_regulon(
            set(RV2516C_REGULATORS)
        ),
    }
    for name, s in regulons_raw.items():
        print(f"  {name}: {len(s)} genes (avant restriction a l'univers)")

    all_tests = []
    for ds_name, ds_path in DATASETS.items():
        rows = load_dataset(ds_path)
        universe = {r["partenaire"] for r in rows}
        print(f"\n{ds_name}: univers = {len(universe)} partenaires de {FOCUS}")
        for regulon_name, regulon_full in regulons_raw.items():
            regulon_in_universe = regulon_full & universe
            for top_n in TOP_NS:
                top_genes = {r["partenaire"] for r in rows[:top_n]}
                overlap = len(top_genes & regulon_in_universe)
                expected = top_n * len(regulon_in_universe) / len(universe)
                p = hypergeom_test(
                    len(universe), len(regulon_in_universe), top_n, overlap
                )
                all_tests.append(
                    {
                        "dataset": ds_name,
                        "regulon": regulon_name,
                        "regulon_size_universe": len(regulon_in_universe),
                        "regulon_size_total": len(regulon_full),
                        "top_n": top_n,
                        "universe_size": len(universe),
                        "overlap": overlap,
                        "expected": round(expected, 2),
                        "fold": round(overlap / expected, 2) if expected > 0 else None,
                        "p": p,
                        "genes_overlap": ";".join(
                            sorted(top_genes & regulon_in_universe)
                        ),
                    }
                )

    pvals = [t["p"] for t in all_tests]
    qvals = bh_correct(pvals)
    for t, q in zip(all_tests, qvals):
        t["q_BH"] = q

    out_tsv = os.path.join(OUT_DIR, "enrichment.tsv")
    fields = [
        "dataset",
        "regulon",
        "top_n",
        "universe_size",
        "regulon_size_universe",
        "regulon_size_total",
        "overlap",
        "expected",
        "fold",
        "p",
        "q_BH",
        "genes_overlap",
    ]
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for t in sorted(all_tests, key=lambda x: x["p"]):
            w.writerow({k: t[k] for k in fields})
    print(f"\nEcrit : {out_tsv}")

    # Croisement obligatoire entre les deux jeux avant de citer un resultat
    print("\n=== Lecture croisee (discipline : significatif q<0.05 dans les DEUX jeux) ===")
    by_regulon_topn = {}
    for t in all_tests:
        key = (t["regulon"], t["top_n"])
        by_regulon_topn.setdefault(key, {})[t["dataset"]] = t
    cross_validated = []
    for (regulon, top_n), by_ds in by_regulon_topn.items():
        if len(by_ds) < 2:
            continue
        sig_both = all(by_ds[ds]["q_BH"] < 0.05 for ds in DATASETS)
        if sig_both:
            cross_validated.append((regulon, top_n, by_ds))
        print(
            f"{regulon:45s} top{top_n:<4d} "
            + " | ".join(
                f"{ds}: overlap={by_ds[ds]['overlap']}/{by_ds[ds]['regulon_size_universe']} "
                f"fold={by_ds[ds]['fold']} q={by_ds[ds]['q_BH']:.3g}"
                for ds in DATASETS
            )
            + ("  <-- CROISE" if sig_both else "")
        )

    with open(os.path.join(OUT_DIR, "resume.md"), "w") as fh:
        fh.write("# P10.3.1 -- enrichissement formel du voisinage de co-expression\n\n")
        fh.write(f"Script : `analyses/phase35_p10_3_1_enrichment.py`.\n\n")
        fh.write("## Regulons testes\n\n")
        for name, s in regulons_raw.items():
            fh.write(f"- **{name}** : {len(s)} genes au total.\n")
        fh.write(
            "\nSource DosR : Rustad et al. 2008, PLOS ONE, Table S1 (fichier deja sur "
            "disque, projet frere `Rv2699c`). Source SigH/SigF/PERSIST : reseau "
            "TFOE Minch et al. 2015 / Rustad et al. 2014 (champ `regulation` de "
            "l'atlas `annotation_mtbc`), condition-aveugle par construction (voir "
            "caveat de l'atlas lui-meme).\n\n"
        )
        fh.write("## Tous les tests (tries par p croissante)\n\n")
        fh.write(
            "| jeu | regulon | top-N | overlap/regulon(univers) | attendu | fold | p | q_BH |\n"
        )
        fh.write("|---|---|---|---|---|---|---|---|\n")
        for t in sorted(all_tests, key=lambda x: x["p"]):
            fh.write(
                f"| {t['dataset']} | {t['regulon']} | {t['top_n']} | "
                f"{t['overlap']}/{t['regulon_size_universe']} | {t['expected']} | "
                f"{t['fold']} | {t['p']:.3g} | {t['q_BH']:.3g} |\n"
            )
        fh.write("\n## Resultats croises (q<0.05 dans les DEUX jeux)\n\n")
        if cross_validated:
            for regulon, top_n, by_ds in cross_validated:
                fh.write(f"- **{regulon}**, top-{top_n} : ")
                fh.write(
                    " ; ".join(
                        f"{ds} overlap={by_ds[ds]['overlap']}/{by_ds[ds]['regulon_size_universe']} "
                        f"(fold {by_ds[ds]['fold']}, q={by_ds[ds]['q_BH']:.3g})"
                        for ds in DATASETS
                    )
                )
                fh.write("\n")
        else:
            fh.write("**Aucun regulon ne croise les deux jeux a q<0.05.**\n")

    print(f"Ecrit : {os.path.join(OUT_DIR, 'resume.md')}")
    return cross_validated


if __name__ == "__main__":
    main()
