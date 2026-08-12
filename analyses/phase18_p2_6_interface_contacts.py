#!/usr/bin/env python3
"""phase18_p2_6_interface_contacts.py -- P2.6 : recroiser les contacts d'interface RÉELS de P2.4,
pas seulement son score agrégé, pour départager les deux lectures laissées ouvertes par le contrôle
négatif mutant.

RAPPEL DU PROBLÈME. J2 (unité B R113A/H115A/R118A + ADN) est QUASI IDENTIQUE à J1 (unité B sauvage +
ADN) sur l'ipTM et le PAE agrégés -- aucun effet de la mutation. Deux lectures, non départagées par
le score seul : (i) ces résidus ne portent pas la spécificité de liaison prédite (négatif biologique) ;
(ii) le modèle garde la même conformation globale malgré la perte des trois chaînes latérales
(insensibilité de méthode à une mutation ponctuelle). ANALYSE ICI, avec les données déjà produites,
aucun nouveau calcul Boltz.

DEUX TESTS INDÉPENDANTS SUR LES COORDONNÉES DU .cif :

  T1. Les résidus R113/H115/R118 sont-ils RÉELLEMENT au contact de l'ADN dans J1 (avant mutation) ?
      Si NON (ils sont déjà loin de l'ADN dans le modèle sauvage), la piste "ce sont les résidus de
      contact" était déjà faible AVANT même de regarder le mutant -- le négatif de J2 devient
      simplement cohérent, pas surprenant.
      Si OUI (contact net en J1), la question se déplace vers T2.

  T2. Le contact protéine-ADN se redistribue-t-il en J2 (le modèle bouge/compense), ou le domaine
      reste-t-il dans une pose QUASI IDENTIQUE à celle de J1 malgré la perte des trois chaînes
      latérales ? Une pose quasi identique malgré la perte de contacts de résidus qui étaient bien
      au contact en J1 est le signal le plus direct d'une insensibilité de méthode : le modèle n'a
      simplement pas recalculé une interface différente.

CALIBRATION : mêmes calculs sur J4 (BldC, contrôle positif, régulateur RÉEL lié à cet ADN) pour savoir
à quoi ressemble un contact protéine-ADN "réel" dans CE pipeline -- distances, nombre de résidus
d'interface -- avant de juger si J1/J2 s'en approchent ou s'en éloignent.

Numérotation : chaque .cif est un fragment autonome (chaîne A = protéine, auth_seq_id 1-based DANS LE
FRAGMENT). Unité B = résidus H37Rv 88-147 -> position_fragment = résidu_H37Rv - 87 (vérifié : positions
26/28/31 = ARG/HIS/ARG en J1, ALA/ALA/ALA en J2, cf. cahier). BldC (J4) n'a pas cette correspondance
(protéine différente), ses résidus sont numérotés dans son propre référentiel 1-68.

Entrées : résultats/p2_4_boltz/out_{j1_unitB_dna,j2_unitB_mut_dna,j4_bldc_dna}/.../. cif (déjà sur disque)
Sorties : résultats/p2_6_interface_contacts/{j1,j2,j4}_contacts.tsv + resume.md
Run: python analyses/phase18_p2_6_interface_contacts.py
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / "résultats" / "p2_4_boltz"
OUT = ROOT / "résultats" / "p2_6_interface_contacts"

JOBS = {
    "j1": (BASE / "out_j1_unitB_dna" / "boltz_results_j1_unitB_dna" / "predictions" /
           "j1_unitB_dna" / "j1_unitB_dna_model_0.cif", "unité B sauvage (TEST)"),
    "j2": (BASE / "out_j2_unitB_mut_dna" / "boltz_results_j2_unitB_mut_dna" / "predictions" /
           "j2_unitB_mut_dna" / "j2_unitB_mut_dna_model_0.cif", "unité B R113A/H115A/R118A (MUTANT)"),
    "j4": (BASE / "out_j4_bldc_dna" / "boltz_results_j4_bldc_dna" / "predictions" /
           "j4_bldc_dna" / "j4_bldc_dna_model_0.cif", "BldC (CONTRÔLE POSITIF)"),
}

FRAGMENT_TO_H37RV_OFFSET = 87   # position_fragment + 87 = residu H37Rv, unite B (88-147) seulement
MUTANT_SITE_FRAGMENT = {26: "R113A", 28: "H115A", 31: "R118A"}
# chaines laterales : atomes lourds AU-DELA du CB (donc exclut le squelette N/CA/C/O ET le CB commun
# a tout residu y compris Ala) -- teste specifiquement l'engagement de LA chaine laterale, pas la
# simple proximite du squelette qui ne bouge quasiment jamais entre WT et mutant par construction.
BACKBONE_ATOMS = {"N", "CA", "C", "O"}


def parse_cif_atoms(path: Path) -> list[dict]:
    """Parse minimal des lignes ATOM d'un .cif Boltz (format verifie manuellement, cf. docstring)."""
    atoms = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            continue
        f = line.split()
        atoms.append({
            "elem": f[2], "atom_name": f[3], "resname": f[5],
            "seq_id": int(f[7]), "chain": f[15],
            "x": float(f[10]), "y": float(f[11]), "z": float(f[12]),
        })
    return atoms


def dist2(a: dict, b: dict) -> float:
    return (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2


def analyze(path: Path, label: str) -> dict:
    atoms = parse_cif_atoms(path)
    protein = [a for a in atoms if a["chain"] == "A"]
    dna = [a for a in atoms if a["chain"] in ("B", "C")]

    by_res: dict[int, list[dict]] = {}
    for a in protein:
        by_res.setdefault(a["seq_id"], []).append(a)

    rows = []
    for seq_id, res_atoms in sorted(by_res.items()):
        resname = res_atoms[0]["resname"]
        all_min = min((dist2(a, d) for a in res_atoms for d in dna), default=None)
        sc_atoms = [a for a in res_atoms if a["atom_name"] not in BACKBONE_ATOMS]
        sc_min = min((dist2(a, d) for a in sc_atoms for d in dna), default=None)
        rows.append({
            "seq_id": seq_id, "resname": resname,
            "min_dist_any_atom": all_min ** 0.5 if all_min is not None else None,
            "min_dist_sidechain": sc_min ** 0.5 if sc_min is not None else None,
        })

    n_interface_4_5 = sum(1 for r in rows if r["min_dist_any_atom"] is not None and r["min_dist_any_atom"] <= 4.5)
    return {"label": label, "n_protein_res": len(by_res), "n_dna_atoms": len(dna),
            "n_interface_4_5A": n_interface_4_5, "rows": rows}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P2.6 : contacts d'interface réels, J1 (sauvage) vs J2 (mutant) vs J4 (positif) ==\n")

    results = {}
    for tag, (path, label) in JOBS.items():
        if not path.exists():
            print(f"  [ABSENT] {tag} : {path}")
            continue
        r = analyze(path, label)
        results[tag] = r
        print(f"-- {tag.upper()} : {label} --")
        print(f"   {r['n_protein_res']} résidus protéine, {r['n_dna_atoms']} atomes ADN, "
              f"{r['n_interface_4_5A']} résidus à l'interface (<= 4,5 Å, tout atome)")
        with open(OUT / f"{tag}_contacts.tsv", "w") as fh:
            fh.write("seq_id\tresname\tmin_dist_any_atom\tmin_dist_sidechain\n")
            for row in r["rows"]:
                fh.write(f"{row['seq_id']}\t{row['resname']}\t"
                         f"{row['min_dist_any_atom']:.2f}\t"
                         f"{row['min_dist_sidechain'] if row['min_dist_sidechain'] is not None else 'NA'}\n")

    print("\n-- T1 : R113/H115/R118 sont-ils au contact de l'ADN dans J1 (avant mutation) ? --")
    if "j1" in results:
        j1_by_id = {r["seq_id"]: r for r in results["j1"]["rows"]}
        for frag_id, name in MUTANT_SITE_FRAGMENT.items():
            r = j1_by_id.get(frag_id)
            h37rv = frag_id + FRAGMENT_TO_H37RV_OFFSET
            if r:
                print(f"   position {frag_id} (H37Rv {h37rv}, {name[:-1]} en J1) : "
                      f"distance min tout-atome = {r['min_dist_any_atom']:.2f} Å, "
                      f"chaîne latérale = {r['min_dist_sidechain']:.2f} Å" if r['min_dist_sidechain'] is not None
                      else f"   position {frag_id} : pas de chaîne latérale au-delà de CB")
        # contexte : distribution generale en J1 pour comparer
        all_min = sorted(r["min_dist_any_atom"] for r in results["j1"]["rows"] if r["min_dist_any_atom"] is not None)
        print(f"\n   Pour comparaison, distance min tout-atome sur les 60 résidus de J1 : "
              f"médiane {all_min[len(all_min)//2]:.2f} Å, minimum {all_min[0]:.2f} Å, "
              f"{sum(1 for d in all_min if d <= 4.5)} résidus <= 4,5 Å au total.")

    print("\n-- T2 : le contact se redistribue-t-il en J2, ou la pose reste-t-elle quasi identique ? --")
    if "j1" in results and "j2" in results:
        j1_by_id = {r["seq_id"]: r["min_dist_any_atom"] for r in results["j1"]["rows"]}
        j2_by_id = {r["seq_id"]: r["min_dist_any_atom"] for r in results["j2"]["rows"]}
        deltas = []
        for seq_id in sorted(set(j1_by_id) & set(j2_by_id)):
            d1, d2 = j1_by_id[seq_id], j2_by_id[seq_id]
            if d1 is not None and d2 is not None:
                deltas.append((seq_id, d1, d2, d2 - d1))
        rmsd_like = (sum(d[3] ** 2 for d in deltas) / len(deltas)) ** 0.5 if deltas else float("nan")
        print(f"   écart quadratique moyen des distances protéine-ADN par résidu, J1 vs J2 : "
              f"{rmsd_like:.2f} Å (0 = pose protéine-ADN rigoureusement identique résidu par résidu)")
        big_moves = [d for d in deltas if abs(d[3]) > 1.0]
        print(f"   résidus dont la distance à l'ADN change de plus de 1 Å entre J1 et J2 : "
              f"{len(big_moves)} / {len(deltas)}")
        for seq_id, d1, d2, delta in sorted(big_moves, key=lambda x: -abs(x[3]))[:10]:
            h37rv = seq_id + FRAGMENT_TO_H37RV_OFFSET
            marque = f" <- {MUTANT_SITE_FRAGMENT[seq_id]}" if seq_id in MUTANT_SITE_FRAGMENT else ""
            print(f"     position {seq_id} (H37Rv {h37rv}) : J1={d1:.2f} Å -> J2={d2:.2f} Å "
                  f"(Δ={delta:+.2f}){marque}")
        for frag_id, name in MUTANT_SITE_FRAGMENT.items():
            if frag_id in j1_by_id and frag_id in j2_by_id:
                print(f"   {name} : J1={j1_by_id[frag_id]:.2f} Å -> J2={j2_by_id[frag_id]:.2f} Å "
                      f"(Δ={j2_by_id[frag_id]-j1_by_id[frag_id]:+.2f})")

    print("\n-- Calibration : à quoi ressemble un contact réel dans ce pipeline (J4, BldC) ? --")
    if "j4" in results:
        all_min4 = sorted(r["min_dist_any_atom"] for r in results["j4"]["rows"] if r["min_dist_any_atom"] is not None)
        print(f"   J4 (BldC, positif) : médiane {all_min4[len(all_min4)//2]:.2f} Å, "
              f"minimum {all_min4[0]:.2f} Å, {results['j4']['n_interface_4_5A']} résidus <= 4,5 Å "
              f"sur {results['j4']['n_protein_res']} résidus totaux "
              f"({100*results['j4']['n_interface_4_5A']/results['j4']['n_protein_res']:.1f} %)")
        if "j1" in results:
            print(f"   J1 (unité B, test) : {results['j1']['n_interface_4_5A']} résidus <= 4,5 Å "
                  f"sur {results['j1']['n_protein_res']} "
                  f"({100*results['j1']['n_interface_4_5A']/results['j1']['n_protein_res']:.1f} %)")

    summary = {tag: {"label": r["label"], "n_protein_res": r["n_protein_res"],
                     "n_interface_4_5A": r["n_interface_4_5A"]} for tag, r in results.items()}
    (OUT / "resume.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    print(f"\nSorties : {OUT}/")


if __name__ == "__main__":
    main()
