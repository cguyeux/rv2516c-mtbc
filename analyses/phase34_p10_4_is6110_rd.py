#!/usr/bin/env python3
"""phase34_p10_4_is6110_rd.py -- P10.4 : le locus Rv2516c-Rv2517c est-il une cible d'insertion
IS6110 ou une borne de délétion structurale à travers les lignées du MTBC ?

L'ANGLE, ET EN QUOI IL DIFFÈRE DE P8.3. P8.3 a cherché les CICATRICES d'un événement mobile ANCIEN
dans la fenêtre de 9 kb (bornes attL/attR, biais de GC, répétitions directes) : négatif sur les trois
signatures. La question posée ici est l'autre moitié : le locus est-il, ou a-t-il RÉCEMMENT été, une
cible ACTIVE de transposition ou de délétion dans une ou plusieurs lignées vivantes ? Ce sont deux
questions différentes, et la seconde se lit dans des données que P8.3 n'a jamais ouvertes : les
positions d'insertion IS6110 par souche et le catalogue de délétions structurales.

TROIS SOURCES, TOUTES DÉJÀ SUR DISQUE.
  1. `insertion_sequences` des report.json de bdd/actuelle : chaque souche porte la liste de ses IS
     avec position et un drapeau `is_reference` (présente chez H37Rv) / non-référence (insertion
     ACQUISE par rapport à H37Rv). C'est exactement la donnée de transposition récente.
  2. Le catalogue local complet de gaps structuraux CUS_GS du projet voisin
     `structural_variation_is_rd` (67 694 gaps avec bornes et nombre de porteurs).
  3. Le catalogue canonique de 187 RD (Bespiatykh, via le classeur Sola) du même projet.

CALIBRATION, PARCE QU'UN COMPTE BRUT NE VEUT RIEN DIRE. « 3 insertions dans la fenêtre » est
ininterprétable sans savoir combien en porte une fenêtre quelconque. On construit donc l'histogramme
GÉNOME ENTIER des insertions non-référence par fenêtre de 1 kb, et le locus est situé par son
PERCENTILE dans cette distribution, avec les vrais points chauds affichés comme échelle. Même
discipline que P8.3 sur le GC, dont la leçon était précisément qu'un écart « substantiel » comparé à
une moyenne se révélait au 4e percentile une fois la dispersion regardée.

ÉCHANTILLONNAGE. Lire les 146 000 report.json coûterait ~100 Go d'E/S pour une question à laquelle
un échantillon stratifié répond aussi bien. On tire au plus N_PAR_CLADE souches par répertoire-clade
(631 clades), ce qui couvre toute la diversité de lignées sans sur-représenter L2/L4. Chaque fichier
est lu PARTIELLEMENT : `insertion_sequences` se situe vers 45 % du fichier, avant le gros champ
`snp` ; on lit par blocs jusqu'à fermeture du tableau.

CE QUE LA PISTE PRÉVOIT, ET POURQUOI LE TEST GARDE DE LA VALEUR SI ELLE A RAISON. Le contre-argument
posé était : « le locus est quasi-invariant à travers 145k génomes, un chevauchement avec un point
chaud actif est a priori peu probable ». Un négatif ici n'est donc pas une déception : c'est la
première mesure DIRECTE de la stabilité structurale du locus, là où le dossier ne disposait que d'un
argument par absence de RD connue. Un locus qui n'est ni inséré ni délété dans un échantillon large,
alors qu'il porte un module d'origine mobile, est un argument POSITIF de domestication achevée.

Sorties : résultats/p10_4_is6110_rd/{is_positions.tsv, fenetres_1kb.tsv, cusgs_locus.tsv,
          rd_locus.tsv, resume.md}
Run: python analyses/phase34_p10_4_is6110_rd.py
"""
from __future__ import annotations

import csv
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
except Exception:
    pass

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
BDD = MTBC / "bdd" / "actuelle"
VOISIN = MTBC / "structural_variation_is_rd"
OUT = ROOT / "résultats" / "p10_4_is6110_rd"

RV2516C = (2_832_710, 2_833_513)
RV2517C = (2_833_510, 2_833_761)
LOCUS = (RV2516C[0], RV2517C[1])          # la paire, 1 052 pb
FENETRE_LARGE = (LOCUS[0] - 10_000, LOCUS[1] + 10_000)
CASSETTE = (2_827_000, 2_836_000)         # Rv2512c-Rv2518c, bornes larges (cf. P8.1)
GENOME = 4_411_532
BIN = 1000
N_PAR_CLADE = 12
SEED = 20260810
MARKER = b'"insertion_sequences"'


DECODER = json.JSONDecoder()


def read_is(path: Path) -> list[dict] | None:
    """Lecture PARTIELLE : on s'arrête dès que le tableau insertion_sequences est fermé.

    PIÈGE RENCONTRÉ ET CORRIGÉ (2026-08-10). La première version cherchait le premier `]` après
    le `[` ouvrant. Or les entrées NON-RÉFÉRENCE — les seules qui nous intéressent — portent un
    champ imbriqué `"is_range":[3,1356]` que les entrées de référence n'ont pas. Le tableau était
    donc tronqué EXACTEMENT à la première insertion acquise : le lecteur était aveugle à la donnée
    mesurée, et le biais allait dans le sens du résultat attendu (aucune insertion nulle part).
    `raw_decode` fait l'analyse syntaxique réelle et rend la position de fin ; s'il manque des
    octets, on lit un bloc de plus au lieu de conclure.
    """
    try:
        with path.open("rb") as fh:
            buf = b""
            while True:
                chunk = fh.read(262_144)
                buf += chunk
                i = buf.find(MARKER)
                if i >= 0:
                    j = buf.find(b"[", i)
                    if j >= 0:
                        try:
                            obj, _ = DECODER.raw_decode(buf[j:].decode("utf-8", "replace"))
                            return obj
                        except ValueError:
                            pass          # tableau incomplet : il faut lire davantage
                if not chunk:
                    return None           # EOF : marqueur absent (report.json de schéma réduit)
                if len(buf) > 8_000_000:
                    return None
    except Exception:
        return None


def sample_strains() -> list[tuple[str, Path]]:
    rng = random.Random(SEED)
    out = []
    for clade in sorted(p for p in BDD.iterdir() if p.is_dir()):
        strains = [d for d in clade.iterdir() if d.is_dir()]
        rng.shuffle(strains)
        for s in strains[:N_PAR_CLADE]:
            r = s / "NC_000962.3" / "report.json"
            if r.exists():
                out.append((clade.name, r))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # ---------- 1. paysage d'insertion IS -------------------------------------------------
    strains = sample_strains()
    print(f"{len(strains)} souches échantillonnées sur {len(set(c for c, _ in strains))} clades")
    novel: list[tuple[str, str, str, int]] = []     # clade, souche, IS, position
    n_ok = n_fail = 0
    for i, (clade, path) in enumerate(strains, 1):
        rec = read_is(path)
        if rec is None:
            n_fail += 1
            continue
        n_ok += 1
        for x in rec:
            if not x.get("is_reference", True):
                novel.append((clade, path.parts[-3], x.get("name", "?"), int(x["position"])))
        if i % 1000 == 0:
            print(f"  {i}/{len(strains)} lus, {len(novel)} insertions non-référence")
    print(f"{n_ok} report.json lus ({n_fail} illisibles), "
          f"{len(novel)} insertions NON-RÉFÉRENCE au total")

    with (OUT / "is_positions.tsv").open("w") as fh:
        fh.write("clade\tsouche\tIS\tposition\n")
        for c, s, n, p in novel:
            fh.write(f"{c}\t{s}\t{n}\t{p}\n")

    par_is = Counter(n for _, _, n, _ in novel)
    print("familles d'IS mobilisées :", par_is.most_common(8))

    is6110 = [p for _, _, n, p in novel if n == "IS6110"]
    bins = Counter(p // BIN for p in is6110)
    all_bins = [bins.get(b, 0) for b in range(GENOME // BIN + 1)]
    with (OUT / "fenetres_1kb.tsv").open("w") as fh:
        fh.write("bin_1kb\tdebut\tfin\tn_insertions_IS6110\n")
        for b, n in sorted(bins.items(), key=lambda x: -x[1]):
            fh.write(f"{b}\t{b*BIN}\t{b*BIN+BIN-1}\t{n}\n")

    def compte(a: int, b: int) -> int:
        return sum(1 for p in is6110 if a <= p <= b)

    n_locus = compte(*LOCUS)
    n_cassette = compte(*CASSETTE)
    n_large = compte(*FENETRE_LARGE)
    occupes = [n for n in all_bins if n > 0]
    print(f"\nIS6110 non-référence : {len(is6110)} insertions, "
          f"{len(occupes)} fenêtres de 1 kb occupées sur {len(all_bins)}")
    print(f"  médiane des fenêtres occupées : {sorted(occupes)[len(occupes)//2] if occupes else 0} ; "
          f"max : {max(all_bins)}")
    print(f"  PAIRE Rv2516c-Rv2517c ({LOCUS[0]}-{LOCUS[1]}, 1 052 pb) : {n_locus} insertion(s)")
    print(f"  cassette Rv2512c-Rv2518c (9 kb)                        : {n_cassette}")
    print(f"  fenêtre ±10 kb                                         : {n_large}")
    pct = 100.0 * sum(1 for n in all_bins if n <= n_locus) / len(all_bins)
    print(f"  percentile du locus parmi les fenêtres de 1 kb : {pct:.1f}")
    print("\n  points chauds réels (échelle de comparaison) :")
    for b, n in bins.most_common(12):
        print(f"    {b*BIN:>9}-{b*BIN+BIN-1:<9} {n:5d} insertions")
    if is6110:
        proche = min(is6110, key=lambda p: min(abs(p - LOCUS[0]), abs(p - LOCUS[1])))
        d = min(abs(proche - LOCUS[0]), abs(proche - LOCUS[1]))
        print(f"\n  insertion non-référence la plus PROCHE du locus : {proche} (à {d:,} pb)")

    # ---------- 2. gaps structuraux CUS_GS ------------------------------------------------
    cus = VOISIN / "résultats" / "cusgs_local_lineage_stats.csv"
    lines = ["cus_gs\tstart\tend\tlength_bp\tn_local_carriers\tdominant_lineage\tf_in\trelation"]
    n_over = 0
    if cus.exists():
        with cus.open() as fh:
            for r in csv.DictReader(fh):
                s, e = int(r["start"]), int(r["end"])
                if e < FENETRE_LARGE[0] or s > FENETRE_LARGE[1]:
                    continue
                rel = ("COUVRE la paire" if s <= LOCUS[0] and e >= LOCUS[1]
                       else "chevauche la paire" if s <= LOCUS[1] and e >= LOCUS[0]
                       else "dans ±10 kb")
                n_over += 1
                lines.append(f"{r['cus_gs']}\t{s}\t{e}\t{r['length_bp']}\t"
                             f"{r['n_local_carriers']}\t{r['dominant_lineage']}\t{r['f_in']}\t{rel}")
        print(f"\ngaps structuraux CUS_GS dans la fenêtre ±10 kb : {n_over}")
        for l in lines[1:]:
            print("   " + l.replace("\t", "  "))
    else:
        print("\ncatalogue CUS_GS introuvable")
    (OUT / "cusgs_locus.tsv").write_text("\n".join(lines) + "\n")

    # ---------- 3. RD canoniques ------------------------------------------------------------
    rdf = VOISIN / "data" / "bespiatykh_canonical_rd.csv"
    rlines = ["rd_name\tstart\tlength\tend\tgenes\trelation"]
    n_rd = 0
    if rdf.exists():
        with rdf.open() as fh:
            for r in csv.DictReader(fh):
                if not r["start"]:
                    continue
                s = int(r["start"])
                e = int(r["end"]) if r.get("end") else (
                    s + int(r["length"]) if r.get("length") else s)
                if e < FENETRE_LARGE[0] or s > FENETRE_LARGE[1]:
                    continue
                n_rd += 1
                rlines.append(f"{r['rd_name']}\t{s}\t{r.get('length','')}\t{e}\t"
                              f"{r.get('genes','')}\tdans ±10 kb")
        print(f"\nRD canoniques (187) dans la fenêtre ±10 kb : {n_rd}")
        for l in rlines[1:]:
            print("   " + l.replace("\t", "  "))
    (OUT / "rd_locus.tsv").write_text("\n".join(rlines) + "\n")

    (OUT / "resume.md").write_text(
        f"# P10.4 — IS6110 et délétions structurales au locus Rv2516c-Rv2517c\n\n"
        f"- {n_ok} souches lues sur {len(set(c for c, _ in strains))} clades ; "
        f"{len(is6110)} insertions IS6110 non-référence.\n"
        f"- paire Rv2516c-Rv2517c : {n_locus} insertion(s) ; cassette 9 kb : {n_cassette} ; "
        f"±10 kb : {n_large} ; percentile {pct:.1f}.\n"
        f"- gaps CUS_GS dans ±10 kb : {n_over} ; RD canoniques : {n_rd}.\n")
    print(f"\nsorties dans {OUT}")


if __name__ == "__main__":
    main()
