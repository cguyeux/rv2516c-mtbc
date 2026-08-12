#!/usr/bin/env python3
"""phase7_p7_4_specificite_antibiotiques.py -- P7.4, second volet : l'induction est-elle
SPÉCIFIQUE des dommages à l'ADN, ou une réponse générale au stress ?

CE QUE CE SCRIPT TESTE, ET POURQUOI IL EST NÉCESSAIRE.
phase6 a établi sur GSE166501 que Rv2516c est induit par la moxifloxacine (RPKM 57 -> 138 à
24 h / 0,3x CMI, percentile 93 du génome), en même temps que son partenaire d'opéron Rv2517c,
et avec une cinétique RETARDÉE par rapport au coeur SOS (rien à 4 h, alors que recA est déjà
à +2,3). Le contre-argument le plus fort contre une lecture « réponse aux dommages de l'ADN »
est immédiat : **à dose sub-inhibitrice pendant 24-72 h, presque n'importe quel antibiotique
finit par induire un programme de stress et de persistance**. Une induction sous une seule
molécule ne distingue pas « répond aux cassures d'ADN » de « répond au stress ».

Le test qui tranche est un panel MULTI-CLASSES. GSE71200 (même institut, ISB) expose H37Rv à
neuf antibiotiques nommés couvrant six mécanismes d'action distincts, dont un seul endommage
l'ADN :

    moxifloxacine        ADN gyrase, coupure double brin       -> DOMMAGES À L'ADN
    rifampicine          ARN polymérase                        -> transcription
    isoniazide           InhA, acides mycoliques               -> paroi
    éthionamide          InhA (pro-drogue)                     -> paroi
    cyclosérine          D-Ala-D-Ala ligase, peptidoglycane    -> paroi
    streptomycine        ribosome 30S                          -> traduction
    kanamycine           ribosome 30S                          -> traduction
    capréomycine         ribosome, interface 30S-50S           -> traduction
    PA-824 (prétomanide) nitroimidazole, respiration / NO      -> respiration

PRÉDICTIONS OPPOSÉES, donc test réellement falsifiant :
  - si Rv2516c suit les dommages à l'ADN, il monte sous moxifloxacine et PAS sous les autres,
    et il suit recA classe par classe ;
  - s'il s'agit d'une réponse générale au stress, il monte sous la plupart des classes, et sa
    corrélation avec recA à travers le panel est faible.

GARDE-FOUS.
  1. **recA est le contrôle de spécificité de l'essai lui-même.** Il DOIT être haut sous
     moxifloxacine et bas ailleurs. S'il ne l'est pas, le panel ne mesure pas ce qu'on croit
     et aucune conclusion sur Rv2516c n'est recevable.
  2. **Percentile intra-lame, pas log-ratio brut.** Ce sont des puces deux couleurs, une lame
     par condition, avec des amplitudes très différentes. Comparer des log-ratios bruts entre
     lames compare des normalisations. Le rang du gène DANS SA PROPRE lame est la seule
     quantité comparable (même discipline qu'en P8.3 et phase6).
  3. **Échantillons signalés défectueux par les déposants** ((BAD), (dead), (poor), « ? »)
     écartés, et la liste des exclusions est imprimée : on ne filtre pas en silence.
  4. Les composés propriétaires IMTB* (mécanisme non public) ne servent PAS à conclure ; ils
     ne servent qu'à élargir la distribution nulle du comportement de Rv2516c sous molécule
     quelconque, ce qui est précisément la question « stress général ».

Entrées : data/gse71200/{series_matrix.txt.gz, GPL5774_family.soft.gz}
Sorties : résultats/p7_4_moxi/{specificite.tsv, specificite.md}
Run: python analyses/phase7_p7_4_specificite_antibiotiques.py
"""
from __future__ import annotations
import bisect, gzip, random, re, statistics as st
from math import sqrt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "gse71200"
OUT = ROOT / "résultats" / "p7_4_moxi"

CLASSES = {
    "Moxifloxacin":  ("ADN gyrase", "DOMMAGES ADN"),
    "Rifampicin":    ("ARN polymerase", "transcription"),
    "INH":           ("InhA / acides mycoliques", "paroi"),
    "Ethionamide":   ("InhA (pro-drogue)", "paroi"),
    "Cycloserine":   ("D-Ala-D-Ala ligase", "paroi"),
    "Streptomycin":  ("ribosome 30S", "traduction"),
    "Kanamycin":     ("ribosome 30S", "traduction"),
    "Capreomycin":   ("ribosome 30S-50S", "traduction"),
    "PA-824":        ("nitroimidazole / respiration", "respiration"),
}
WATCH = {
    "Rv2516c": "cible du projet",
    "Rv2517c": "partenaire d'operon",
    "Rv2514c": "toxine TA",
    "Rv2515c": "antitoxine TA",
    "Rv2737c": "recA -- CONTROLE DE SPECIFICITE",
    "Rv2720":  "lexA",
    "Rv3370c": "dnaE2",
    "Rv2031c": "hspX (DosR)",
    "Rv2703":  "sigA (stable)",
}
BAD = re.compile(r"\(BAD\)|\(dead\)|\(poor\)|\?")
RVRE = re.compile(r"^(Rv\d{4}[A-Za-z]?)")


def load_platform() -> dict[str, str]:
    """sonde -> locus Rv (colonne ORF, sinon la correspondance H37Rv des sondes CDC1551)."""
    probe2rv, inside, cols = {}, False, None
    with gzip.open(DATA / "GPL5774_family.soft.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                inside = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not inside:
                continue
            f = line.rstrip("\n").split("\t")
            if cols is None:
                cols = {n: i for i, n in enumerate(f)}
                continue
            pid = f[0]
            rv = ""
            o = f[cols["ORF"]] if cols.get("ORF", 99) < len(f) else ""
            m = RVRE.match(o.strip())
            if m:
                rv = m.group(1)
            else:                                    # sonde CDC1551 : lire l'homologue H37Rv
                h = f[cols["Strain H37Rv"]] if cols.get("Strain H37Rv", 99) < len(f) else ""
                m = RVRE.match(h.strip())
                if m:
                    rv = m.group(1)
            if rv:
                probe2rv[pid] = rv
    return probe2rv


def load_matrix() -> tuple[list[str], list[str], dict[str, list[float]]]:
    titles, gsms, table, inside = [], [], {}, False
    with gzip.open(DATA / "series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif inside:
                f = line.rstrip("\n").split("\t")
                if not gsms:
                    gsms = [x.strip('"') for x in f[1:]]
                    continue
                vals = []
                for x in f[1:]:
                    try:
                        vals.append(float(x))
                    except ValueError:
                        vals.append(float("nan"))
                table[f[0].strip('"')] = vals
    return titles, gsms, table


def mean_pct(sel: list[dict], rv: str) -> float:
    """Percentile intra-lame moyen d'un locus sur un sous-ensemble de lames."""
    v = [r[rv] for r in sel if r.get(rv) == r.get(rv)]
    return st.mean(v) if v else float("nan")


def main() -> None:
    print("== P7.4 volet 2 : spécificité de classe (GSE71200, 106 lames, ISB) ==")
    probe2rv = load_platform()
    titles, gsms, table = load_matrix()
    print(f"  plateforme GPL5774 : {len(probe2rv)} sondes mappées sur un locus Rv")
    print(f"  matrice : {len(table)} sondes x {len(titles)} lames\n")

    excluded = [t for t in titles if BAD.search(t)]
    print(f"  écartés (signalés défectueux par les déposants) : {len(excluded)} lames")
    print(f"    {', '.join(excluded)}\n")

    # valeur par gène et par lame = médiane des sondes du gène (log-ratio traité/témoin)
    per_gene: dict[str, list[float]] = {}
    for probe, vals in table.items():
        rv = probe2rv.get(probe)
        if rv:
            per_gene.setdefault(rv, []).append(vals)
    gene_vals = {rv: [st.median([v[i] for v in mat if v[i] == v[i]] or [float("nan")])
                      for i in range(len(titles))] for rv, mat in per_gene.items()}
    print(f"  {len(gene_vals)} loci Rv quantifiés\n")

    # percentile intra-lame : le rang dans SA lame, seule quantite comparable entre lames.
    # Colonnes triees une fois pour toutes puis recherche dichotomique : sans cela le calcul
    # est quadratique en nombre de genes (3916^2 x 106 comparaisons).
    sorted_cols = [sorted(g[i] for g in gene_vals.values() if g[i] == g[i])
                   for i in range(len(titles))]

    def pct(rv: str, i: int) -> float:
        x = gene_vals[rv][i]
        col = sorted_cols[i]
        if x != x or not col:
            return float("nan")
        return 100.0 * bisect.bisect_left(col, x) / len(col)

    rows = []
    for i, t in enumerate(titles):
        if BAD.search(t):
            continue
        drug = t.split("-")[0] if not t.startswith("PA-824") else "PA-824"
        rows.append({"lame": gsms[i], "titre": t, "drug": drug, "i": i,
                     **{rv: pct(rv, i) for rv in WATCH if rv in gene_vals}})

    named = [r for r in rows if r["drug"] in CLASSES]
    imtb = [r for r in rows if r["drug"].startswith("IMTB")]
    print(f"  {len(named)} lames de molécules nommées, {len(imtb)} lames IMTB (mécanisme non public)\n")

    print("-- CONTRÔLE DE SPÉCIFICITÉ : recA distingue-t-il bien les classes ? --")
    print(f"   {'molecule':14} {'classe':16} {'n':>2}  {'recA':>7} {'dnaE2':>7} {'lexA':>7}")
    order = sorted(CLASSES, key=lambda d: -st.mean([r["Rv2737c"] for r in named if r["drug"] == d] or [0]))
    for d in order:
        sel = [r for r in named if r["drug"] == d]
        if not sel:
            continue
        print(f"   {d:14} {CLASSES[d][1]:16} {len(sel):>2}  "
              f"{mean_pct(sel,'Rv2737c'):>7.0f} {mean_pct(sel,'Rv3370c'):>7.0f} "
              f"{mean_pct(sel,'Rv2720'):>7.0f}")

    print("\n-- RÉSULTAT : percentile intra-lame de la cassette Rv2514c-Rv2517c --")
    print(f"   {'molecule':14} {'classe':16} {'Rv2516c':>8} {'Rv2517c':>8} {'Rv2514c':>8} {'Rv2515c':>8} {'sigA':>6}")
    for d in order:
        sel = [r for r in named if r["drug"] == d]
        if not sel:
            continue
        print(f"   {d:14} {CLASSES[d][1]:16} {mean_pct(sel,'Rv2516c'):>8.0f} "
              f"{mean_pct(sel,'Rv2517c'):>8.0f} {mean_pct(sel,'Rv2514c'):>8.0f} "
              f"{mean_pct(sel,'Rv2515c'):>8.0f} {mean_pct(sel,'Rv2703'):>6.0f}")

    if imtb:
        print(f"   {'IMTB* (nul)':14} {'non public':16} {mean_pct(imtb,'Rv2516c'):>8.0f} "
              f"{mean_pct(imtb,'Rv2517c'):>8.0f} {mean_pct(imtb,'Rv2514c'):>8.0f} "
              f"{mean_pct(imtb,'Rv2515c'):>8.0f} {mean_pct(imtb,'Rv2703'):>6.0f}")

    # corrélation Rv2516c vs recA à travers TOUTES les lames retenues : le discriminant
    def corr(a: list[float], b: list[float]) -> float:
        p = [(x, y) for x, y in zip(a, b) if x == x and y == y]
        if len(p) < 3:
            return float("nan")
        xs, ys = [x for x, _ in p], [y for _, y in p]
        mx, my = st.mean(xs), st.mean(ys)
        num = sum((x - mx) * (y - my) for x, y in p)
        den = sqrt(sum((x - mx) ** 2 for x in xs)) * sqrt(sum((y - my) ** 2 for y in ys))
        return num / den if den else float("nan")

    keep = [r["i"] for r in rows]
    print(f"\n-- Corrélation à travers les {len(keep)} lames retenues (percentiles intra-lame) --")
    rec = [r["Rv2737c"] for r in rows]
    for rv, lab in (("Rv2516c", "cible"), ("Rv2517c", "operon"), ("Rv2514c", "toxine"),
                    ("Rv2515c", "antitoxine"), ("Rv3370c", "dnaE2"), ("Rv2703", "sigA")):
        if rv in gene_vals:
            print(f"   r(recA, {rv:8}) = {corr(rec, [r[rv] for r in rows]):+.3f}   {lab}")

    # ── calibration du couplage d'opéron contre une distribution nulle ──
    # +0,36 n'est ni fort ni faible dans l'absolu : il l'est par rapport à ce qu'obtiennent
    # deux gènes tirés au hasard sur les mêmes lames. Et le comparateur qui compte est la
    # paire toxine-antitoxine Rv2514c-Rv2515c, seul opéron du voisinage validé
    # EXPÉRIMENTALEMENT (Tandon 2019) : elle donne le plafond réaliste de la méthode.
    print("\n-- Le couplage d'opéron est-il réel ? calibration contre 4000 paires au hasard --")
    allpct = {rv: [pct(rv, i) for i in keep] for rv in gene_vals}
    usable = [rv for rv, v in allpct.items() if sum(1 for x in v if x == x) >= 0.8 * len(keep)]
    rng = random.Random(0)
    null = sorted(c for c in (corr(allpct[a], allpct[b])
                              for a, b in (rng.sample(usable, 2) for _ in range(4000))) if c == c)
    mu, sd = st.mean(null), st.pstdev(null)
    print(f"   nulle : {len(null)} paires sur {len(usable)} loci | moyenne {mu:+.3f} | "
          f"écart-type {sd:.3f} | 95e pct {null[int(.95*len(null))]:+.3f}")
    print(f"   {'paire':28} {'r':>7} {'pct nulle':>10} {'z':>7}")
    for a, b, lab in (("Rv2516c", "Rv2517c", "l'opéron testé"),
                      ("Rv2514c", "Rv2515c", "TA validé (Tandon 2019) -- PLAFOND"),
                      ("Rv2516c", "Rv2514c", "cassette, opérons différents"),
                      ("Rv2516c", "Rv2518c", "voisin HORS cassette -- plancher"),
                      ("Rv2737c", "Rv3370c", "recA vs dnaE2 (régulon SOS)")):
        if a in allpct and b in allpct:
            c = corr(allpct[a], allpct[b])
            p = 100.0 * bisect.bisect_left(null, c) / len(null)
            print(f"   {a+' / '+b:28} {c:>+7.3f} {p:>10.1f} {(c-mu)/sd:>+7.2f}   {lab}")

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "specificite.md", "w") as fh:
        fh.write("# P7.4 volet 2 -- specificite de classe (GSE71200, GPL5774, ISB)\n\n")
        fh.write(f"{len(rows)} lames retenues sur {len(titles)} ({len(excluded)} ecartees comme "
                 "defectueuses par les deposants). Valeurs = percentile intra-lame du log-ratio\n"
                 "traite/temoin, seule quantite comparable entre lames deux couleurs.\n\n")
        fh.write("## recA valide l'essai\n\n| molecule | classe | n | recA | dnaE2 | lexA |\n|---|---|---|---|---|---|\n")
        for d in order:
            sel = [r for r in named if r["drug"] == d]
            if not sel:
                continue
            fh.write(f"| {d} | {CLASSES[d][1]} | {len(sel)} | {mean_pct(sel,'Rv2737c'):.0f} | "
                     f"{mean_pct(sel,'Rv3370c'):.0f} | {mean_pct(sel,'Rv2720'):.0f} |\n")
        fh.write("\n## La cassette Rv2514c-Rv2517c par classe\n\n"
                 "| molecule | classe | Rv2516c | Rv2517c | Rv2514c | Rv2515c | sigA |\n|---|---|---|---|---|---|---|\n")
        for d in order:
            sel = [r for r in named if r["drug"] == d]
            if not sel:
                continue
            fh.write(f"| {d} | {CLASSES[d][1]} | {mean_pct(sel,'Rv2516c'):.0f} | "
                     f"{mean_pct(sel,'Rv2517c'):.0f} | {mean_pct(sel,'Rv2514c'):.0f} | "
                     f"{mean_pct(sel,'Rv2515c'):.0f} | {mean_pct(sel,'Rv2703'):.0f} |\n")
        fh.write(f"\n## Couplage d'operon, calibre sur {len(null)} paires au hasard\n\n"
                 f"Nulle : moyenne {mu:+.3f}, ecart-type {sd:.3f}, 95e percentile {null[int(.95*len(null))]:+.3f}.\n\n"
                 "| paire | r | pct nulle | z | lecture |\n|---|---|---|---|---|\n")
        for a, b, lab in (("Rv2516c", "Rv2517c", "l'operon teste"),
                          ("Rv2514c", "Rv2515c", "TA valide (Tandon 2019) -- plafond"),
                          ("Rv2516c", "Rv2514c", "cassette, operons differents"),
                          ("Rv2516c", "Rv2518c", "voisin hors cassette -- plancher"),
                          ("Rv2737c", "Rv3370c", "recA vs dnaE2 (regulon SOS)")):
            if a in allpct and b in allpct:
                c = corr(allpct[a], allpct[b])
                fh.write(f"| {a} / {b} | {c:+.3f} | {100.0*bisect.bisect_left(null, c)/len(null):.1f} | "
                         f"{(c-mu)/sd:+.2f} | {lab} |\n")
    with open(OUT / "specificite.tsv", "w") as fh:
        cols = [rv for rv in WATCH if rv in gene_vals]
        fh.write("gsm\ttitre\tmolecule\t" + "\t".join(cols) + "\n")
        for r in rows:
            fh.write(f"{r['lame']}\t{r['titre']}\t{r['drug']}\t" +
                     "\t".join(f"{r[c]:.1f}" if r.get(c) == r.get(c) else "NA" for c in cols) + "\n")
    print(f"\nÉcrit {OUT}/specificite.tsv")


if __name__ == "__main__":
    main()
