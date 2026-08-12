#!/usr/bin/env python3
"""phase6_p7_4_geo_moxifloxacine.py -- P7.4 : Rv2516c répond-il à une fluoroquinolone ?

POURQUOI CE SCRIPT EXISTE, ET POURQUOI IL REMPLACE UN HAND-OFF.
P1.4 avait identifié un verrou : le papier AAC 2024 (doi 10.1128/aac.01185-23, 99 profils
sur 17 antibiotiques, Rv2516c dedans) est sous paywall, et la piste demandait à l'auteur un
téléchargement institutionnel. Mais ce papier est une MÉTA-ANALYSE de jeux publics : son
résumé dit « by pooling and analyzing Mtb microarray and RNA-seq data sets ». Les données
primaires sont donc dans GEO, en accès libre. Leçon déjà payée deux fois aujourd'hui :
**un hôte bloqué n'est pas une ressource bloquée.**

Jeu utilisé : GSE166501, « M. tuberculosis transcriptional response to Moxifloxacin »
(Institute for Systems Biology, labo Baliga — le même groupe que le réseau de TF Minch 2015
déjà cité dans la fiche atlas). Publications associées : Immanuel et al. 2021, npj Syst Biol
Appl, doi 10.1038/s41540-021-00205-6 ; Srinivas et al. 2021, Cell Rep Methods,
doi 10.1016/j.crmeth.2021.100123.

Design : H37Rv, moxifloxacine à 0x / 0,12x / 0,3x CMI, temps 4 / 24 / 72 h, triplicats
biologiques = 27 échantillons. Un fichier RPKM par gène et par échantillon.

CE QUE LA QUESTION EXIGE, ET SES GARDE-FOUS.

  1. **Contrôle interne de déclenchement de la réponse.** Les doses sont SUB-inhibitrices.
     Si la réponse SOS n'est pas déclenchée du tout, un résultat nul sur Rv2516c ne veut
     rien dire. On mesure donc d'abord recA, lexA et dnaE2 : ce sont EUX qui disent si
     l'expérience a le pouvoir de répondre à la question. Sans ce contrôle, le script
     produirait une absence de preuve maquillée en preuve d'absence.

  2. **Calibration génome entier, pas seuil arbitraire.** Un log2FC de 0,8 n'est ni grand
     ni petit dans l'absolu : il l'est par rapport à la distribution des ~4000 gènes du
     même contraste. On rapporte le percentile et le z-score empiriques (méthode validée
     en P8.3, où un écart cru « substantiel » s'est révélé au 4e percentile).

  3. **Le test le plus informatif est le CONTRASTE Rv2516c vs Rv2517c.** Les deux gènes
     partagent un opéron avec chevauchement de 4 pb. Si Rv2517c bouge et pas Rv2516c, le
     couplage transcriptionnel supposé tombe, et avec lui l'argument « Rv2516c est dans la
     réponse aux dommages de l'ADN par contiguïté ». Si les deux bougent ensemble, l'opéron
     est confirmé fonctionnellement — ce que ni la synténie ni STRING ne peuvent établir.

  4. n = 3 vs 3 : le t de Welch est faible. Il est rapporté, mais c'est la **cohérence des
     trois réplicats et la reproductibilité à travers doses et temps** qui décide, pas p.

Entrée  : archive GEO GSE166501_RAW.tar (2,6 Mo), téléchargée par --fetch ou déjà extraite.
Sorties : résultats/p7_4_moxi/{par_gene.tsv, genome_log2fc.tsv, resume.md}
Run: python analyses/phase6_p7_4_geo_moxifloxacine.py --fetch
"""
from __future__ import annotations
import argparse, gzip, math, re, statistics as st, subprocess, tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p7_4_moxi"
RAW = ROOT / "data" / "gse166501"
URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166501/suppl/"
       "GSE166501_RAW.tar")
TAR_SIZE = 2754560

# Gènes suivis, groupés par ce qu'ils servent à démontrer.
PANEL = {
    "cible": [
        ("Rv2516c", "la protéine du projet (wHTH AlpA + ferredoxin-like + Ig-like)"),
        ("Rv2517c", "partenaire d'opéron, chevauchement 4 pb ; dit induit en SOS (Iacobino 2021)"),
    ],
    "voisinage": [
        ("Rv2514c", "toxine du couple TA validé (Tandon 2019)"),
        ("Rv2515c", "antitoxine du couple TA (HTH_3 Xre + Peptidase_M78)"),
        ("Rv2518c", "ldtB, bordure de la cassette"),
    ],
    "controle_SOS": [
        ("Rv2737c", "recA — LE contrôle de déclenchement"),
        ("Rv2720", "lexA — répresseur SOS, auto-induit"),
        ("Rv3370c", "dnaE2 — polymérase translésionnelle, SOS tardif"),
    ],
    "controle_dormance": [
        ("Rv2031c", "hspX/acr — régulon DosR, marqueur de dormance"),
        ("Rv2623", "USP DosR-dépendant"),
    ],
    "controle_stable": [
        ("Rv2703", "sigA — facteur sigma domestique, ne doit pas bouger"),
    ],
}
WATCH = {rv: (grp, desc) for grp, lst in PANEL.items() for rv, desc in lst}
SAMPLE = re.compile(r"GSM\d+_H37Rv_T(\d+)_MXF_([\d.x]+)_([ABC])")
TIMES, DOSES, REPS = ["4", "24", "72"], ["0.12", "0.3"], ["A", "B", "C"]


def fetch() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    tar = RAW / "GSE166501_RAW.tar"
    if not (tar.exists() and tar.stat().st_size == TAR_SIZE):
        print(f"  téléchargement de {URL} ({TAR_SIZE/1e6:.1f} Mo)…")
        subprocess.run(["curl", "-sS", "--retry", "3", "-o", str(tar), URL], check=True)
        got = tar.stat().st_size
        # comparer à content-length, jamais faire confiance à un download muet (KB boltz)
        assert got == TAR_SIZE, f"archive tronquée : {got}/{TAR_SIZE} o"
    with tarfile.open(tar) as t:
        t.extractall(RAW, filter="data")
    print(f"  {len(list(RAW.glob('*.txt.gz')))} échantillons dans {RAW}")


def load() -> dict:
    """{(temps, dose, rep): {rv: rpkm}} sur l'ensemble du génome."""
    data = {}
    for f in sorted(RAW.glob("*.MTb.Transcript.txt.gz")):
        m = SAMPLE.match(f.name)
        if not m:
            continue
        expr = {}
        with gzip.open(f, "rt") as fh:
            next(fh)
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 4 and c[1]:
                    try:
                        expr[c[1]] = float(c[3])
                    except ValueError:
                        pass
        data[(m.group(1), m.group(2), m.group(3))] = expr
    return data


def l2(x: float) -> float:
    return math.log2(x + 1.0)


def welch_p(a: list[float], b: list[float]) -> float:
    """t de Welch bilatéral, approximation normale (n=3 : indicatif seulement)."""
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    va, vb = st.variance(a), st.variance(b)
    se = math.sqrt(va / len(a) + vb / len(b))
    if se == 0:
        return float("nan")
    t = abs(st.mean(a) - st.mean(b)) / se
    return math.erfc(t / math.sqrt(2))          # borne haute, conservatrice


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="télécharger l'archive GEO")
    a = ap.parse_args()
    print("== P7.4 : réponse de Rv2516c à la moxifloxacine (GSE166501, ISB) ==")
    if a.fetch or not list(RAW.glob("*.txt.gz")):
        fetch()
    data = load()
    print(f"  {len(data)} échantillons chargés, {len(next(iter(data.values())))} gènes\n")

    genes = sorted(set(next(iter(data.values()))))
    rows, genome = [], {}

    for t in TIMES:
        ctrl = {r: data[(t, "0x", r)] for r in REPS if (t, "0x", r) in data}
        for d in DOSES:
            trt = {r: data[(t, d, r)] for r in REPS if (t, d, r) in data}
            if not ctrl or not trt:
                continue
            # distribution génome entier du même contraste -> calibration
            fc_all = {}
            for g in genes:
                c = [l2(s.get(g, 0.0)) for s in ctrl.values()]
                x = [l2(s.get(g, 0.0)) for s in trt.values()]
                if max(max(c), max(x)) < l2(5.0):     # gène quasi éteint : bruit pur
                    continue
                fc_all[g] = st.mean(x) - st.mean(c)
            vals = sorted(fc_all.values())
            mu, sd = st.mean(vals), st.pstdev(vals)
            genome[(t, d)] = (len(vals), mu, sd)

            for rv in WATCH:
                if rv not in fc_all:
                    continue
                c = [l2(s.get(rv, 0.0)) for s in ctrl.values()]
                x = [l2(s.get(rv, 0.0)) for s in trt.values()]
                fc = fc_all[rv]
                below = sum(1 for v in vals if v < fc)
                rows.append({
                    "rv": rv, "groupe": WATCH[rv][0], "temps": t, "dose": d,
                    "rpkm_ctrl": st.mean([s.get(rv, 0.0) for s in ctrl.values()]),
                    "rpkm_trt": st.mean([s.get(rv, 0.0) for s in trt.values()]),
                    "log2fc": fc, "pct": 100.0 * below / len(vals),
                    "z": (fc - mu) / sd if sd else float("nan"),
                    "p": welch_p(x, c),
                    "reps_ctrl": c, "reps_trt": x,
                })

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "par_gene.tsv", "w") as fh:
        fh.write("rv\tgroupe\ttemps_h\tdose_xCMI\trpkm_ctrl\trpkm_traite\tlog2FC\t"
                 "percentile_genome\tz_genome\tp_welch\n")
        for r in rows:
            fh.write(f"{r['rv']}\t{r['groupe']}\t{r['temps']}\t{r['dose']}\t"
                     f"{r['rpkm_ctrl']:.1f}\t{r['rpkm_trt']:.1f}\t{r['log2fc']:+.3f}\t"
                     f"{r['pct']:.1f}\t{r['z']:+.2f}\t{r['p']:.3g}\n")

    print("-- Calibration : distribution génome entier du log2FC par contraste --")
    print(f"{'temps':>6} {'dose':>6} {'n genes':>8} {'moyenne':>9} {'ecart-type':>11}")
    for (t, d), (n, mu, sd) in sorted(genome.items(), key=lambda kv: (int(kv[0][0]), kv[0][1])):
        print(f"{t+'h':>6} {d+'x':>6} {n:>8} {mu:>+9.3f} {sd:>11.3f}")

    print("\n-- CONTRÔLE DE DÉCLENCHEMENT : la moxifloxacine a-t-elle induit SOS ? --")
    print("   (si recA ne bouge pas, aucun résultat nul n'est interprétable)")
    for rv in ("Rv2737c", "Rv2720", "Rv3370c"):
        line = [f"{r['log2fc']:+.2f}" for r in rows if r["rv"] == rv]
        pct = [f"{r['pct']:.0f}" for r in rows if r["rv"] == rv]
        nm = WATCH[rv][1].split(" — ")[0]
        print(f"   {rv:9} {nm:8} log2FC {' '.join(f'{v:>6}' for v in line)}")
        print(f"   {'':9} {'':8} pct    {' '.join(f'{v:>6}' for v in pct)}")

    print("\n-- RÉSULTAT : la cible et son opéron --")
    hdr = [f"{t}h/{d}x" for t in TIMES for d in DOSES]
    print(f"   {'gene':9} " + " ".join(f"{h:>9}" for h in hdr))
    for grp in ("cible", "voisinage", "controle_dormance", "controle_stable"):
        for rv, _ in PANEL[grp]:
            sel = {(r["temps"], r["dose"]): r for r in rows if r["rv"] == rv}
            cells = []
            for t in TIMES:
                for d in DOSES:
                    r = sel.get((t, d))
                    cells.append(f"{r['log2fc']:+.2f}" if r else "  ·  ")
            print(f"   {rv:9} " + " ".join(f"{c:>9}" for c in cells))
        print()

    with open(OUT / "resume.md", "w") as fh:
        fh.write("# P7.4 -- reponse de Rv2516c a la moxifloxacine (GSE166501)\n\n")
        fh.write("H37Rv, moxifloxacine 0.12x et 0.3x CMI vs non traite, 4/24/72 h, "
                 "triplicats. RPKM GEO, log2FC sur log2(RPKM+1).\n\n")
        fh.write("| gene | groupe | t (h) | dose | RPKM ctrl | RPKM trt | log2FC | pct genome | z | p Welch |\n")
        fh.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write(f"| {r['rv']} | {r['groupe']} | {r['temps']} | {r['dose']}x | "
                     f"{r['rpkm_ctrl']:.1f} | {r['rpkm_trt']:.1f} | {r['log2fc']:+.3f} | "
                     f"{r['pct']:.1f} | {r['z']:+.2f} | {r['p']:.3g} |\n")
    print(f"Écrit {OUT}/par_gene.tsv et resume.md")


if __name__ == "__main__":
    main()
