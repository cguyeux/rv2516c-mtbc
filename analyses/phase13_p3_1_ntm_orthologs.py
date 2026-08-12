#!/usr/bin/env python3
"""phase13_p3_1_ntm_orthologs.py -- P3.1 : conservation par residu des orthologues NTM, mappee sur
l'architecture a trois unites.

QUESTION. La couche NTM de l'atlas (data/Rv2516c_ntm_conservation.json) donne une identite MOYENNE de
59,1% sur 6/53 NTM porteurs -- un chiffre GLOBAL, pas positionnel. La divergence est une chance :
si le module HTH (98-121, et plus specifiquement le tour+helice de
reconnaissance 111-121 identifie en P2.3) ressort NETTEMENT plus conserve que la moyenne du gene a
travers 6 especes independantes, c'est un argument de fonctionnalite a une profondeur evolutive tres
superieure a celle du MTBC (9 SNP/145k souches, cf. atlas P2f), et INDEPENDANT de la structure (P2.1/
P2.2) et de la charge de surface (P2.3) deja mobilisees.

CE QUE L'ATLAS N'A PAS. `phase57_ntm_orthology.py` lance tblastn en `-outfmt "6 qseqid pident qcovs
bitscore"` -- identite et couverture GLOBALES du hit, aucune sequence, aucune position. Impossible d'en
tirer un profil par residu. Ce script relance tblastn UNE FOIS sur les 6 genomes NTM porteurs deja
identifies (memes seuils, meme base de genomes bdd/hors_mtbc/), avec qseq/sseq/qstart/qend en plus, pour
reconstruire l'alignement complet et mapper l'identite sur les coordonnees H37Rv de la requete.

GENOMES (bdd/hors_mtbc/<espece>/**/genome.fna, verifie present pour les 6) :
    M_avium, M_bouchedurhonense, M_celatum, M_shinjukuense (seul membre MTBAP porteur, cf. atlas P10.2),
    M_simiae, M_timonense.

GARDE-FOU (deja pose par P7.3) : a 59% d'identite moyenne, l'alignement lui-meme est proche de la zone
d'incertitude -- un HTH peut etre conserve par contrainte de REPLI sans que la specificite de LIAISON le
soit, la conservation ne dit pas QUOI est lie. Et un residu peut etre invariant parce qu'il est enfoui
(contrainte structurale generique) plutot que parce qu'il est fonctionnel en surface -- a croiser avec
P2.3 (charge/exposition) plutot qu'a interpreter seul.

NULL MODEL INTERNE (meme discipline que P2.3) : comparer la fenetre HTH/tour+helice a l'identite MOYENNE
du reste de la meme proteine, dans le MEME alignement -- pas a un seuil generique importe d'ailleurs.

Run: python analyses/phase13_p3_1_ntm_orthologs.py
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BDD_HORS = ROOT.parent / "bdd" / "hors_mtbc"
OUT = ROOT / "résultats" / "p3_1_ntm_orthologs"
QUERY = ROOT / "data" / "Rv2516c.faa"

SPECIES = ["M_avium", "M_bouchedurhonense", "M_celatum", "M_shinjukuense", "M_simiae", "M_timonense"]
MTBAP_MEMBER = "M_shinjukuense"  # cf. atlas phase57 : MTBAP = {decipiens, lacus, riyadhense, shinjukuense}

# DECOUVERT A L'EXECUTION (pas suppose a priori) : M_avium / M_bouchedurhonense / M_timonense rendent
# des statistiques BLAST identiques a la decimale pres, et un controle position-par-position confirme
# un profil d'identite RIGOUREUSEMENT identique aux 168 positions communes -- pas un bug (genomes.fna
# distincts, sha256 differents, tailles differentes), mais la biologie du complexe M. avium (MAC) :
# ces trois especes partagent un orthologue quasi-identique. Compter les 6 especes comme 6 lignees
# INDEPENDANTES gonflerait artificiellement le signal de conservation (pseudo-replication phylogenetique,
# meme piege que l'homoplasie/densite d'echantillonnage documente pour les lignees MTBC). Un representant
# du clade MAC suffit ; PHYLO_GROUPS ci-dessous porte cette correction dans le calcul, pas seulement le
# commentaire -- toute fenetre est rapportee sous les DEUX comptages (brut 6 especes, deduplique 4 lignees).
PHYLO_GROUPS = {
    "MAC_clade_avium_bouchedurhonense_timonense": ["M_avium"],  # representant ; verifie identique aux 2 autres
    "M_celatum": ["M_celatum"],
    "M_shinjukuense_MTBAP": ["M_shinjukuense"],
    "M_simiae": ["M_simiae"],
}

MIN_PIDENT = 30.0   # seuils identiques a l'atlas phase57 (presence deja etablie, on ne re-decide rien)
MIN_QCOV = 50.0
QLEN = 267

# Fenetres de domaine, memes bornes que P2.1/P2.3 :
WINDOWS = {
    "proteine_entiere": (1, 267),
    "unite_A_ferredoxin": (1, 87),
    "unite_B_wHTH": (88, 147),
    "HTH_98_121": (98, 121),
    "helice2_101_110": (101, 110),
    "tour_helice_recon_111_121": (111, 121),   # cf. P2.3, la fenetre qui portait le patch basique
    "linker_148_177": (148, 177),
    "unite_C_Ig_178_267": (178, 267),
}


def find_genome(species: str) -> Path | None:
    hits = sorted(BDD_HORS.glob(f"{species}/**/genome.fna"))
    return hits[0] if hits else None


def run_tblastn(genome: Path, tmp: Path) -> list[dict]:
    db = tmp / "db"
    subprocess.run(["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db)],
                    check=True, capture_output=True)
    out = tmp / "hits.tsv"
    fmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq"
    r = subprocess.run(
        ["tblastn", "-query", str(QUERY), "-db", str(db), "-out", str(out),
         "-evalue", "1e-5", "-max_target_seqs", "5", "-max_hsps", "3", "-outfmt", fmt],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"tblastn a echoue : {r.stderr[:300]}")
    rows = []
    for line in out.read_text().splitlines():
        f = line.split("\t")
        qcov = 100.0 * (int(f[5]) - int(f[4]) + 1) / QLEN
        rows.append({"sseqid": f[1], "pident": float(f[2]), "length": int(f[3]),
                     "qstart": int(f[4]), "qend": int(f[5]), "sstart": int(f[6]), "send": int(f[7]),
                     "evalue": float(f[8]), "bitscore": float(f[9]), "qseq": f[10], "sseq": f[11],
                     "qcov": qcov})
    return rows


def best_hit(rows: list[dict]) -> dict | None:
    """Meilleur HSP par bitscore -- un seul retenu par espece, coherent avec le max_hsps=1 de l'atlas
    pour la classification de presence (on ne re-ouvre pas ce choix, on l'instrumente)."""
    passed = [r for r in rows if r["pident"] >= MIN_PIDENT and r["qcov"] >= MIN_QCOV]
    pool = passed or rows
    return max(pool, key=lambda r: r["bitscore"]) if pool else None


def per_residue_identity(hit: dict) -> dict[int, bool]:
    """qseq/sseq alignes (avec gaps '-') -> {position H37Rv 1-based : identique(bool)}.
    Les positions hors du HSP (non alignees) sont absentes du dict (donnee manquante, pas False)."""
    pos = hit["qstart"]
    ident: dict[int, bool] = {}
    for qc, sc in zip(hit["qseq"], hit["sseq"]):
        if qc != "-":
            if sc != "-":
                ident[pos] = (qc.upper() == sc.upper())
            # qc present, sc gap (deletion chez le NTM) : position query existe mais rien a comparer -> absente.
            pos += 1
        # qc gap (insertion chez le NTM) : ne consomme pas de position query.
    return ident


def window_identity(ident_by_species: dict[str, dict[int, bool]], lo: int, hi: int) -> dict:
    """Identite MOYENNE dans la fenetre [lo,hi], calculee POSITION PAR POSITION puis moyennee sur les
    especes qui couvrent cette position -- pas une simple moyenne de pourcentages globaux, pour que le
    denominateur soit homogene entre fenetres de tailles differentes."""
    per_pos = []
    for p in range(lo, hi + 1):
        vals = [sp[p] for sp in ident_by_species.values() if p in sp]
        if vals:
            per_pos.append(sum(vals) / len(vals))
    if not per_pos:
        return {"n_positions_couvertes": 0, "identite_moyenne": None}
    return {"n_positions_couvertes": len(per_pos), "n_positions_totales": hi - lo + 1,
            "identite_moyenne": round(100 * sum(per_pos) / len(per_pos), 1)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P3.1 : conservation par residu, orthologues NTM (6 especes porteuses) ==\n")

    ident_by_species: dict[str, dict[int, bool]] = {}
    hit_summaries = []
    for sp in SPECIES:
        genome = find_genome(sp)
        if genome is None:
            print(f"  [ABSENT] {sp} : genome.fna introuvable sous bdd/hors_mtbc/{sp}/")
            continue
        with tempfile.TemporaryDirectory() as td:
            rows = run_tblastn(genome, Path(td))
        hit = best_hit(rows)
        if hit is None:
            print(f"  [AUCUN HIT] {sp} -- inattendu, l'atlas le classait present")
            continue
        ident = per_residue_identity(hit)
        ident_by_species[sp] = ident
        n_ident = sum(ident.values())
        marque = " <- membre MTBAP" if sp == MTBAP_MEMBER else ""
        print(f"  {sp:22s}{marque:18s} pident={hit['pident']:5.1f}% qcov={hit['qcov']:5.1f}% "
              f"span={hit['qstart']}-{hit['qend']} ({n_ident}/{len(ident)} positions identiques)")
        hit_summaries.append({"species": sp, "is_mtbap": sp == MTBAP_MEMBER,
                              "pident_blast": hit["pident"], "qcov": round(hit["qcov"], 1),
                              "qstart": hit["qstart"], "qend": hit["qend"],
                              "evalue": hit["evalue"], "n_aligned_positions": len(ident),
                              "n_identical_positions": n_ident})

    if not ident_by_species:
        print("\nAucune espece exploitable -- arret.")
        return

    # Verification EMPIRIQUE (pas supposee) du regroupement phylogenetique : les especes du complexe
    # MAC censees partager un orthologue quasi-identique doivent avoir un profil d'identite identique
    # a TOUTES les positions communes ; sinon le regroupement est faux et il ne faut pas dedupliquer.
    print("[CONTROLE] profils d'identite intra-clade MAC (verifie la pseudo-replication avant de corriger)")
    MAC_TRIO = ["M_avium", "M_bouchedurhonense", "M_timonense"]
    for a_sp, b_sp in [(MAC_TRIO[0], MAC_TRIO[1]), (MAC_TRIO[0], MAC_TRIO[2])]:
        if a_sp in ident_by_species and b_sp in ident_by_species:
            a, b = ident_by_species[a_sp], ident_by_species[b_sp]
            common = set(a) & set(b)
            same = sum(1 for p in common if a[p] == b[p])
            verdict = "IDENTIQUES, regroupement MAC confirme" if same == len(common) else "DIVERGENT, ne pas regrouper"
            print(f"  {a_sp} vs {b_sp} : {same}/{len(common)} positions communes au meme etat -- {verdict}")

    def report_windows(source: dict[str, dict[int, bool]], label_suffix: str) -> dict:
        print(f"\n[FENETRES] identite moyenne position-par-position -- {label_suffix} "
              f"({len(source)} groupes)")
        results = {}
        for label, (lo, hi) in WINDOWS.items():
            w = window_identity(source, lo, hi)
            results[label] = w
            cov = w["n_positions_couvertes"]
            tot = w.get("n_positions_totales", hi - lo + 1)
            idn = w["identite_moyenne"]
            idn_s = f"{idn:5.1f}%" if idn is not None else "  n/a"
            print(f"  {label:28s} [{lo:3d}-{hi:3d}] couverture={cov:3d}/{tot:3d}  identite_moyenne={idn_s}")
        return results

    def report_contrasts(source: dict[str, dict[int, bool]], results: dict, label_suffix: str) -> dict:
        print(f"\n[CONTRASTE] fenetre d'interet vs reste de la proteine -- {label_suffix}")
        contrasts = {}
        whole = set(range(1, 268))
        for label in ("HTH_98_121", "tour_helice_recon_111_121", "unite_B_wHTH"):
            lo, hi = WINDOWS[label]
            rest_positions = sorted(whole - set(range(lo, hi + 1)))
            per_pos_rest = []
            for p in rest_positions:
                vals = [sp[p] for sp in source.values() if p in sp]
                if vals:
                    per_pos_rest.append(sum(vals) / len(vals))
            rest_pct = round(100 * sum(per_pos_rest) / len(per_pos_rest), 1) if per_pos_rest else None
            fen_pct = results[label]["identite_moyenne"]
            ecart = round(fen_pct - rest_pct, 1) if (fen_pct is not None and rest_pct is not None) else None
            contrasts[label] = {"fenetre_pct": fen_pct, "reste_proteine_pct": rest_pct, "ecart_points": ecart}
            print(f"  {label:28s} fenetre={fen_pct}% vs reste={rest_pct}%  ecart={ecart:+.1f} pts"
                  if ecart is not None else f"  {label:28s} donnees insuffisantes")
        return contrasts

    window_results = report_windows(ident_by_species, "BRUT, 6 especes (gonfle par le clade MAC)")
    contrasts = report_contrasts(ident_by_species, window_results, "brut, 6 especes")

    dedup = {g: ident_by_species[reps[0]] for g, reps in PHYLO_GROUPS.items()
              if reps[0] in ident_by_species}
    window_results_dedup = report_windows(dedup, "DEDUPLIQUE, 4 lignees phylogenetiques -- A CITER")
    contrasts_dedup = report_contrasts(dedup, window_results_dedup, "deduplique, 4 lignees -- A CITER")

    # Detail par position de la fenetre d'interet (111-121), espece par espece, pour tracabilite.
    print("\n[DETAIL] tour + helice de reconnaissance (111-121), identite par position et par espece")
    detail_rows: list[dict[str, int | bool | None]] = []
    for p in range(111, 122):
        row: dict[str, int | bool | None] = {"pos": p}
        marks = []
        for sp in SPECIES:
            if sp in ident_by_species and p in ident_by_species[sp]:
                v = ident_by_species[sp][p]
                row[sp] = v
                marks.append("=" if v else "x")
            else:
                row[sp] = None
                marks.append(".")
        detail_rows.append(row)
        print(f"  {p:3d}  " + " ".join(f"{sp[2:6]:6s}" for sp in SPECIES))
        print(f"       " + " ".join(f"{m:^6s}" for m in marks))

    resume = {
        "especes_brut": list(ident_by_species.keys()),
        "hits": hit_summaries,
        "avertissement_pseudo_replication": (
            "M_avium/M_bouchedurhonense/M_timonense (clade MAC) partagent un profil d'identite "
            "IDENTIQUE a toutes les positions communes -- les chiffres 'brut' comptent 6 especes mais "
            "n'ont que 4 degres de liberte phylogenetiques reels. Citer les resultats '_deduplique'."
        ),
        "fenetres_brut_6_especes": window_results,
        "contrastes_brut_6_especes": contrasts,
        "groupes_phylogenetiques_4_lignees": list(PHYLO_GROUPS.keys()),
        "fenetres_deduplique_4_lignees": window_results_dedup,
        "contrastes_deduplique_4_lignees": contrasts_dedup,
    }
    (OUT / "resume.json").write_text(json.dumps(resume, indent=1, ensure_ascii=False))
    (OUT / "identite_par_position.json").write_text(
        json.dumps({sp: {str(k): v for k, v in d.items()} for sp, d in ident_by_species.items()},
                   indent=1, ensure_ascii=False))
    print(f"\nSorties : {OUT}/")


if __name__ == "__main__":
    main()
