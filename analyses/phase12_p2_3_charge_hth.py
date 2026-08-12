#!/usr/bin/env python3
"""phase12_p2_3_charge_hth.py -- P2.3 : charge de surface du modele AF sur la face de l'helice
de reconnaissance (98-121), et sur le reste de l'unite B (88-147) en controle.

QUESTION. Un wHTH liant l'ADN presente typiquement, sur l'helice de reconnaissance, une face
enrichie en residus basiques EXPOSES (contacts avec le squelette phosphate et les bases), tandis
que la face opposee, tournee vers le coeur globulaire, est plutot hydrophobe. C'est un signal
GEOMETRIQUE independant de l'homologie de sequence (HHpred) et de la structure (Foldseek) deja
mobilisees : un troisieme axe de preuve, orthogonal aux deux premiers.

DEUX MESURES INDEPENDANTES, calculees sur le MEME modele AlphaFold (query_Rv2516c_AF.pdb,
identique au fichier canonique annotation_mtbc/data/af_models/Rv2516c.pdb -- verifie par diff) :

  1. SASA (Shrake-Rupley, Bio.PDB, pas de binaire externe requis) : charge nette EXPOSEE de
     l'helice de reconnaissance (114-121) contre le reste de l'unite B (88-147)
     et contre la proteine entiere, en pourcentage de surface exposee par residu chimiquement
     basique/acide.
  2. AZIMUT DE CHAINE LATERALE autour de l'axe local de l'helice : si l'helice de reconnaissance
     est une vraie face de liaison, les CB des residus basiques doivent pointer dans une plage
     angulaire etroite, distincte de celle des residus hydrophobes (qui pointent vers le coeur).
     Axe local estime par un fit lineaire (PCA) des CA de l'helice ; angle = position du CB dans
     le plan perpendiculaire a cet axe, origine arbitraire (seul l'ECART inter-groupe compte).

GARDE-FOU. Ce test ne peut PAS prouver la liaison a l'ADN -- une face basique est necessaire, pas
suffisante, et de nombreuses proteines non-ADN ont des patches basiques (liaison ARN, membrane,
partenaire acide). Il ne fait que dire si la geometrie est COMPATIBLE avec le role predit, et sert
de controle croise avec R113A/H115A/R118A deja choisis pour P2.4 sur une base de sequence seule --
si ce test montre que ces trois residus ne sont PAS co-exposes sur la meme face, le design du
mutant de P2.4 devrait etre reexamine.

Run: python analyses/phase12_p2_3_charge_hth.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parent.parent
PDB = ROOT / "résultats" / "p7_1" / "structures" / "query_Rv2516c_AF.pdb"
OUT = ROOT / "résultats" / "p2_3_charge_hth"

BASIC = {"LYS", "ARG", "HIS"}
ACIDIC = {"ASP", "GLU"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO"}

# Bornes de domaine (verifiees par sequence ci-dessus) :
UNIT_A = (1, 87)
UNIT_B = (88, 147)
HTH_FULL = (98, 121)          # helice 2 + tour + helice de reconnaissance
HELIX2 = (101, 110)           # AEIADELGVS
TURN = (111, 113)             # RQR
REC_HELIX = (114, 121)        # VHQLRSTA, l'helice de reconnaissance
TURN_PLUS_REC = (111, 121)    # RQRVHQLRSTA -- le tour ET l'helice, cf. R111/R113 ci-dessous
WING_CANDIDATE = (122, 147)   # reste de l'unite B, wing putatif du wHTH
MUTANT_SITE = {113, 115, 118}  # R113A/H115A/R118A, deja choisis pour P2.4


def load_structure():
    p = PDBParser(QUIET=True)
    s = p.get_structure("Rv2516c", str(PDB))
    chain = list(s[0].child_list)[0]
    residues = {r.id[1]: r for r in chain if r.id[0] == " "}
    return s, residues


def charge_class(resname: str) -> str:
    if resname in BASIC:
        return "basique"
    if resname in ACIDIC:
        return "acide"
    if resname in HYDROPHOBIC:
        return "hydrophobe"
    return "polaire"


def sasa_analysis(structure, residues, out_lines):
    sr = ShrakeRupley()
    sr.compute(structure, level="R")
    # Reference Tien et al. 2013 (theoretical max ASA, Gly-X-Gly), pour la SASA RELATIVE.
    MAXASA = {
        "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225, "GLU": 223,
        "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236, "MET": 224, "PHE": 240,
        "PRO": 159, "SER": 155, "THR": 172, "TRP": 285, "TYR": 263, "VAL": 174,
    }
    rows = []
    for num, res in sorted(residues.items()):
        rsa = float(res.sasa) / MAXASA[res.resname]
        rows.append({"pos": num, "aa": seq1(res.resname), "resname": res.resname,
                     "sasa": round(float(res.sasa), 1), "rsa": round(rsa, 3),
                     "classe": charge_class(res.resname), "exposed": bool(rsa >= 0.25)})

    def window_stats(lo, hi, label):
        w = [r for r in rows if lo <= r["pos"] <= hi]
        exposed = [r for r in w if r["exposed"]]
        n_basic_exp = sum(1 for r in exposed if r["classe"] == "basique")
        n_acid_exp = sum(1 for r in exposed if r["classe"] == "acide")
        n_hydro_exp = sum(1 for r in exposed if r["classe"] == "hydrophobe")
        charge_nette = n_basic_exp - n_acid_exp
        frac_basic = n_basic_exp / len(exposed) if exposed else float("nan")
        out_lines.append(
            f"  {label:24s} [{lo:3d}-{hi:3d}] n={len(w):3d} exposes={len(exposed):3d} "
            f"basiques_exp={n_basic_exp} acides_exp={n_acid_exp} hydrophobes_exp={n_hydro_exp} "
            f"charge_nette_exposee={charge_nette:+d} frac_basique/expose={frac_basic:.2f}"
        )
        return {"label": label, "lo": lo, "hi": hi, "n": len(w), "n_exposed": len(exposed),
                "n_basic_exp": n_basic_exp, "n_acid_exp": n_acid_exp,
                "n_hydro_exp": n_hydro_exp, "net_charge_exposed": charge_nette,
                "frac_basic_of_exposed": frac_basic}

    out_lines.append("\n[SASA] fenetres (seuil expose : RSA >= 0.25, Tien 2013)")
    summary = {
        "proteine_entiere": window_stats(1, 267, "proteine entiere (reference)"),
        "unite_A": window_stats(*UNIT_A, "unite A (ferredoxin-like)"),
        "unite_B": window_stats(*UNIT_B, "unite B (wHTH complet)"),
        "hth_98_121": window_stats(*HTH_FULL, "HTH complet 98-121"),
        "helice2": window_stats(*HELIX2, "helice 2 (101-110)"),
        "helice_reconnaissance": window_stats(*REC_HELIX, "helice de reconnaissance (114-121)"),
        "tour_plus_helice": window_stats(*TURN_PLUS_REC, "tour + helice recon. (111-121)"),
        "wing_candidat": window_stats(*WING_CANDIDATE, "wing candidat (122-147)"),
    }

    # Detail residu par residu du tour + helice de reconnaissance (111-121) et du site mutant
    # P2.4 : la fenetre 114-121 seule masquait R111 et R113, EXPOSES eux aussi (cf. resume).
    out_lines.append("\n[SASA] detail par residu, tour + helice de reconnaissance (111-121) :")
    for r in rows:
        if TURN_PLUS_REC[0] <= r["pos"] <= TURN_PLUS_REC[1]:
            marque = " <- site MUTANT P2.4" if r["pos"] in MUTANT_SITE else ""
            out_lines.append(f"    {r['aa']}{r['pos']} ({r['classe']:10s}) RSA={r['rsa']:.2f} "
                             f"{'EXPOSE' if r['exposed'] else 'enfoui'}{marque}")
    return summary, rows


def helix_axis(residues, lo, hi):
    """Axe local par PCA des CA -- premiere composante principale, orientee N->C."""
    ca = np.array([residues[i]["CA"].coord for i in range(lo, hi + 1) if i in residues])
    centroid = ca.mean(axis=0)
    _, _, vt = np.linalg.svd(ca - centroid)
    axis = vt[0]
    if np.dot(ca[-1] - ca[0], axis) < 0:
        axis = -axis
    return centroid, axis


def azimuth_analysis(residues_raw, out_lines):
    """Angle azimutal du vecteur CA->CB (ou CA->CG pseudo pour GLY) autour de l'axe local
    de l'helice de reconnaissance. Deux residus a la meme face du cylindre helicoidal ont
    un azimut proche ; deux residus a des faces opposees sont separes de ~180 degres."""
    residues = {num: res for num, res in residues_raw.items()}
    lo, hi = REC_HELIX
    _centroid, axis = helix_axis({n: {"CA": r["CA"]} for n, r in residues.items()}, lo, hi)
    axis = axis / np.linalg.norm(axis)
    # Base orthonormee du plan perpendiculaire a l'axe.
    tmp = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = tmp - np.dot(tmp, axis) * axis
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)

    out_lines.append(f"\n[AZIMUT] axe local de l'helice de reconnaissance (fit PCA sur CA {lo}-{hi})")
    entries = []
    for num in range(lo, hi + 1):
        res = residues.get(num)
        if res is None:
            continue
        ca = res["CA"].coord
        cb_atom = res.child_dict.get("CB")
        if cb_atom is None:
            if res.resname != "GLY":
                continue
            # Gly n'a pas de CB : direction approximee par le plan N-CA-C (bissectrice opposee).
            n_atom, c_atom = res.child_dict.get("N"), res.child_dict.get("C")
            if n_atom is None or c_atom is None:
                continue
            v = -( (n_atom.coord - ca) + (c_atom.coord - ca) )
            side_dir = v
        else:
            side_dir = cb_atom.coord - ca
        # Projection dans le plan perpendiculaire a l'axe.
        proj = side_dir - np.dot(side_dir, axis) * axis
        angle = math.degrees(math.atan2(np.dot(proj, e2), np.dot(proj, e1))) % 360
        entries.append({"pos": num, "aa": seq1(res.resname), "classe": charge_class(res.resname),
                        "azimut_deg": round(angle, 1)})
        marque = " <- site MUTANT P2.4" if num in MUTANT_SITE else ""
        out_lines.append(f"    {seq1(res.resname)}{num} ({charge_class(res.resname):10s}) "
                         f"azimut={angle:6.1f} deg{marque}")

    basiques = [e["azimut_deg"] for e in entries if e["classe"] == "basique"]
    hydro = [e["azimut_deg"] for e in entries if e["classe"] == "hydrophobe"]

    def circ_mean_spread(angles):
        if not angles:
            return None, None
        rad = np.radians(angles)
        x, y = np.cos(rad).mean(), np.sin(rad).mean()
        mean = math.degrees(math.atan2(y, x)) % 360
        r = math.hypot(x, y)  # 1 = tres groupe, 0 = disperse
        return round(mean, 1), round(r, 3)

    mean_b, r_b = circ_mean_spread(basiques)
    mean_h, r_h = circ_mean_spread(hydro)
    ecart = None
    if mean_b is not None and mean_h is not None:
        d = abs(mean_b - mean_h) % 360
        ecart = min(d, 360 - d)

    out_lines.append(f"\n  Azimut moyen (circulaire) basiques : {mean_b} deg (compacite r={r_b}, "
                     f"n={len(basiques)})")
    out_lines.append(f"  Azimut moyen (circulaire) hydrophobes : {mean_h} deg (compacite r={r_h}, "
                     f"n={len(hydro)})")
    out_lines.append(f"  Ecart angulaire basiques vs hydrophobes : {ecart} deg "
                     f"(attendu proche de 180 deg si faces opposees, comme pour un vrai HTH)")
    return {"entries": entries, "mean_basic_deg": mean_b, "spread_basic_r": r_b,
            "mean_hydrophobic_deg": mean_h, "spread_hydrophobic_r": r_h,
            "angular_separation_deg": ecart}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    structure, residues = load_structure()
    out_lines = ["== P2.3 : charge de surface et geometrie de face, helice de reconnaissance 114-121 ==",
                 f"Modele : {PDB.relative_to(ROOT)} (AlphaFold Monomer v2.0, verifie identique au "
                 f"fichier canonique annotation_mtbc/data/af_models/Rv2516c.pdb)"]

    sasa_summary, sasa_rows = sasa_analysis(structure, residues, out_lines)
    azimuth_summary = azimuth_analysis(residues, out_lines)

    report = "\n".join(out_lines)
    print(report)
    (OUT / "rapport.txt").write_text(report + "\n")
    (OUT / "sasa_par_residu.json").write_text(json.dumps(sasa_rows, indent=1, ensure_ascii=False))
    (OUT / "resume.json").write_text(json.dumps(
        {"sasa": sasa_summary, "azimut": azimuth_summary}, indent=1, ensure_ascii=False))
    print(f"\nSorties : {OUT}/")


if __name__ == "__main__":
    main()
