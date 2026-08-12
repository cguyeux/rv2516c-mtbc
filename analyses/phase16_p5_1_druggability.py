#!/usr/bin/env python3
"""phase16_p5_1_druggability.py -- P5.1 : poche de liaison / druggabilité sur le modele AF, CALIBREE.

GARDE-FOU CENTRAL (skill `pocket-detection`) : un score de poche brut n'est jamais interpretable sans
temoin positif ET negatif apparies a la population etudiee. Ne PAS lire le score de Rv2516c seul.

REUTILISATION -- rien a recalculer pour la calibration. L'atlas a deja construit et scanne par P2Rank
deux populations de reference sur MEME pipeline / memes modeles AlphaFold :
    annotation_mtbc/résultats/phase72_pockets_af/pockets_af_enzymes.tsv   (n=150, TEMOIN POSITIF,
        enzymes averees EC, plafond atteignable par P2Rank sur un modele predit)
    annotation_mtbc/résultats/phase72_pockets_af/pockets_af_noncat.tsv    (n=150, PLANCHER, proteines
        regulatrices/liaison ADN-ARN SANS EC -- le groupe le plus pertinent pour Rv2516c)
    annotation_mtbc/résultats/phase72_pockets_af/pockets_af.tsv           (contient DEJA Rv2516c,
        calcule par l'atlas sur le MEME modele AF que celui utilise dans tout ce projet)
Ce script ne relance PAS P2Rank : il lit ces trois TSV et situe Rv2516c dans les deux distributions.

SECOND OUTIL, CROISE (skill : "toujours croiser au moins deux outils"). fpocket, geometrique pur,
lance ICI sur le meme modele AF (résultats/p5_1_druggability/fpocket_out/), independant de P2Rank.
Sous-note les sites polaires/metal (garde-fou skill) -- pertinent ici car un score fpocket bas
n'argumente PAS contre une surface de liaison a l'ADN (plate, etendue, chargee), seulement contre une
poche a petite molecule classique (profonde, apolaire).

Run: python analyses/phase16_p5_1_druggability.py
"""
from __future__ import annotations

import csv
import json
import subprocess
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT.parent / "annotation_mtbc"
CALIB = ATLAS / "résultats" / "phase72_pockets_af"
OUT = ROOT / "résultats" / "p5_1_druggability"
MODEL = ROOT / "résultats" / "p7_1" / "structures" / "query_Rv2516c_AF.pdb"
FPOCKET = ATLAS / "tools" / "fpocket" / "bin" / "fpocket"

# Bornes de domaine, memes que P2.1/P2.3/P3.1.
DOMAINS = {"unite_A_ferredoxin": (1, 87), "unite_B_wHTH": (88, 147),
           "linker": (148, 177), "unite_C_Ig": (178, 267)}


def domain_of(pos: int) -> str:
    for name, (lo, hi) in DOMAINS.items():
        if lo <= pos <= hi:
            return name
    return "hors_bornes"


def prob(row: dict) -> float:
    v = row.get("best_probability", "")
    return float(v) if v not in ("", None) else 0.0


def load_tsv(path: Path) -> list[dict]:
    return list(csv.DictReader(open(path), delimiter="\t"))


def p2rank_calibration() -> dict:
    enz = load_tsv(CALIB / "pockets_af_enzymes.tsv")
    noncat = load_tsv(CALIB / "pockets_af_noncat.tsv")
    targets = load_tsv(CALIB / "pockets_af.tsv")
    rv = next(r for r in targets if r["rv"] == "Rv2516c")

    def stats(rows, label):
        probs = sorted(prob(r) for r in rows)
        n_conf = sum(1 for p in probs if p >= 0.5)
        n_any = sum(1 for r in rows if int(r.get("n_pockets") or 0) > 0)
        return {"label": label, "n": len(rows), "median_prob": probs[len(probs) // 2],
                "frac_prob_ge_0_5": n_conf / len(rows), "frac_has_pocket": n_any / len(rows),
                "sorted_probs": probs}

    st_enz = stats(enz, "enzymes (positif)")
    st_noncat = stats(noncat, "non_catalytiques (plancher)")
    obs = prob(rv)

    def percentile(sorted_probs, x):
        n_below = sum(1 for p in sorted_probs if p < x)
        n_equal = sum(1 for p in sorted_probs if p == x)
        return {"pct_conservateur": 100 * n_below / len(sorted_probs),
                "n_ex_aequo": n_equal}

    residues = [int(x) for x in rv["best_residues"].split()] if rv.get("best_residues") else []
    dom_count: dict[str, int] = {}
    for r in residues:
        d = domain_of(r)
        dom_count[d] = dom_count.get(d, 0) + 1

    result = {
        "rv2516c": {"best_probability": obs, "n_pockets": int(rv["n_pockets"]),
                    "best_score": float(rv["best_score"]), "plddt_af": float(rv["plddt_af"]),
                    "best_residues": residues, "best_pocket_domaine": dom_count},
        "enzymes": {k: v for k, v in st_enz.items() if k != "sorted_probs"},
        "non_catalytiques": {k: v for k, v in st_noncat.items() if k != "sorted_probs"},
        "percentile_vs_enzymes": percentile(st_enz["sorted_probs"], obs),
        "percentile_vs_non_catalytiques": percentile(st_noncat["sorted_probs"], obs),
    }
    return result


def run_fpocket() -> list[dict]:
    OUT.mkdir(parents=True, exist_ok=True)
    work = OUT / "fpocket_run"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    local_pdb = work / "query_Rv2516c_AF.pdb"
    shutil.copy(MODEL, local_pdb)
    subprocess.run([str(FPOCKET), "-f", local_pdb.name], cwd=work, check=True,
                   capture_output=True)
    info = work / "query_Rv2516c_AF_out" / "query_Rv2516c_AF_info.txt"
    pockets_dir = work / "query_Rv2516c_AF_out" / "pockets"

    pockets = []
    current: dict = {}
    for line in info.read_text().splitlines():
        line = line.strip()
        if line.startswith("Pocket "):
            if current:
                pockets.append(current)
            current = {"id": line.split()[1].rstrip(":")}
        elif ":" in line:
            k, v = line.split(":", 1)
            current[k.strip()] = v.strip()
    if current:
        pockets.append(current)

    for p in pockets:
        pdb_file = pockets_dir / f"pocket{p['id']}_atm.pdb"
        residues = set()
        if pdb_file.exists():
            for line in pdb_file.read_text().splitlines():
                if line.startswith("ATOM"):
                    residues.add(int(line[22:26]))
        p["residues"] = sorted(residues)
        dom_count: dict[str, int] = {}
        for r in p["residues"]:
            d = domain_of(r)
            dom_count[d] = dom_count.get(d, 0) + 1
        p["domaines"] = dom_count
    return pockets


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P5.1 : druggabilité calibrée sur le modèle AF de Rv2516c ==\n")

    print("[P2RANK, réutilisé depuis l'atlas -- aucun recalcul]")
    cal = p2rank_calibration()
    rv = cal["rv2516c"]
    print(f"  Rv2516c : best_probability={rv['best_probability']}  n_pockets={rv['n_pockets']}  "
          f"best_score={rv['best_score']}  pLDDT_AF={rv['plddt_af']}")
    print(f"  meilleure poche, résidus par domaine : {rv['best_pocket_domaine']}")
    e, n = cal["enzymes"], cal["non_catalytiques"]
    print(f"  ENZYMES (positif)         n={e['n']}  médiane={e['median_prob']:.3f}  "
          f"frac(prob>=0.5)={100*e['frac_prob_ge_0_5']:.1f}%  frac(a une poche)={100*e['frac_has_pocket']:.1f}%")
    print(f"  NON-CATALYTIQUES (plancher) n={n['n']}  médiane={n['median_prob']:.3f}  "
          f"frac(prob>=0.5)={100*n['frac_prob_ge_0_5']:.1f}%  frac(a une poche)={100*n['frac_has_pocket']:.1f}%")
    pe, pn = cal["percentile_vs_enzymes"], cal["percentile_vs_non_catalytiques"]
    print(f"  Rv2516c percentile vs enzymes : {pe['pct_conservateur']:.1f}%")
    print(f"  Rv2516c percentile vs non-catalytiques : {pn['pct_conservateur']:.1f}%")

    print("\n[FPOCKET, lancé ici, second outil indépendant]")
    pockets = run_fpocket()
    for p in pockets[:5]:
        print(f"  Poche {p['id']:>2} : score={p.get('Score','?')} "
              f"druggability={p.get('Druggability Score','?')} "
              f"volume={p.get('Volume','?')} résidus_par_domaine={p['domaines']}")

    resume = {"p2rank": cal, "fpocket_pockets": pockets}
    (OUT / "resume.json").write_text(json.dumps(resume, indent=1, ensure_ascii=False))
    print(f"\nSorties : {OUT}/")


if __name__ == "__main__":
    main()
