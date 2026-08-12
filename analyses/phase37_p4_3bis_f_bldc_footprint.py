#!/usr/bin/env python3
"""P4.3bis-f, volet 2 (partie 1) -- lecture directe, dans le co-cristal 6AMA, des
contacts BASE-SPECIFIQUES (par opposition aux contacts de squelette phosphate/sucre)
entre BldC (chaine A) et son operateur (chaines N/R), pour construire le demi-site
de reconnaissance REEL de ce wHTH -- pas une hypothese, une lecture directe de
structure resolue.

Limite assumee et explicitee (pas un oubli) : AlpA (8C3T) et RdfS (8DGL), les deux
autres homologues cites par P2.1, sont des structures APO (sans ADN dans le fichier
deposé) -- aucun contact base-residu n'est lisible pour elles avec les donnees
disponibles ici. Elles ne peuvent donc alimenter qu'la comparaison de repli (deja
faite en P2.1), pas un second demi-site independant. Le test qui suit porte donc sur
le demi-site de BldC specifiquement, explicitement qualifie comme tel (« motif de
reconnaissance de type BldC », pas « motif de Rv2516c ») en l'absence d'un
alignement structural valide reliant residu a residu l'unite B a BldC.
"""
import json
import os

import numpy as np
from Bio.PDB import PDBParser

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
STRUCT_DIR = os.path.join(ROOT, "résultats", "p4_3bis_f_orthogonal", "structures")
OUT_DIR = os.path.join(ROOT, "résultats", "p4_3bis_f_orthogonal")

BLDC_PDB = os.path.join(STRUCT_DIR, "6AMA.pdb")
CONTACT_CUTOFF = 4.5

BACKBONE_ATOMS = {
    "P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'",
}
BASE_ATOM_PREFIX = None  # tout ce qui n'est pas backbone/sucre est traite "base"

DNA3_TO_1 = {"DA": "A", "DT": "T", "DC": "C", "DG": "G"}


def main():
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("bldc_full", BLDC_PDB)
    model = s[0]
    protein_chain = model["A"]
    dna_chains = {"N": model["N"], "R": model["R"]}

    dna_atoms = []
    for cid, chain in dna_chains.items():
        for res in chain:
            if res.id[0] != " ":
                continue
            for atom in res:
                if atom.element == "H":
                    continue
                dna_atoms.append(
                    {
                        "chain": cid,
                        "resnum": res.id[1],
                        "resname": res.get_resname(),
                        "atomname": atom.get_name(),
                        "coord": atom.coord,
                        "is_base": atom.get_name() not in BACKBONE_ATOMS,
                    }
                )
    dna_coords = np.array([a["coord"] for a in dna_atoms])

    base_contacts = {}  # resnum(BldC) -> list of dna contacts (base atoms only)
    backbone_only = set()
    for res in protein_chain:
        if res.id[0] != " ":
            continue
        prot_atoms = [a for a in res if a.element != "H"]
        best_base_contact = None
        has_backbone_contact = False
        for atom in prot_atoms:
            d = np.linalg.norm(dna_coords - atom.coord, axis=1)
            imin = int(np.argmin(d))
            if d[imin] <= CONTACT_CUTOFF:
                dna_hit = dna_atoms[imin]
                if dna_hit["is_base"]:
                    if best_base_contact is None or d[imin] < best_base_contact["dist"]:
                        best_base_contact = {
                            "dna_chain": dna_hit["chain"],
                            "dna_resnum": dna_hit["resnum"],
                            "dna_resname": dna_hit["resname"],
                            "dna_atom": dna_hit["atomname"],
                            "prot_atom": atom.get_name(),
                            "dist": float(d[imin]),
                        }
                else:
                    has_backbone_contact = True
        if best_base_contact is not None:
            base_contacts[res.id[1]] = {
                "resname": res.get_resname(),
                **best_base_contact,
            }
        elif has_backbone_contact:
            backbone_only.add(res.id[1])

    print(f"Residus BldC (chaine A) a contact de BASE (<= {CONTACT_CUTOFF} A) : {len(base_contacts)}")
    for rn, info in sorted(base_contacts.items()):
        print(
            f"  {info['resname']}{rn} -- {info['prot_atom']} vs "
            f"{info['dna_resname']}{info['dna_resnum']}({info['dna_chain']}).{info['dna_atom']} "
            f"= {info['dist']:.2f} A"
        )
    print(f"\nResidus a contact BACKBONE seulement (squelette phosphate/sucre) : {sorted(backbone_only)}")

    # Sequence des bases contactees, dans l'ordre de la chaine N (brin de reference)
    chain_n_seq = {}
    for res in dna_chains["N"]:
        if res.id[0] == " ":
            chain_n_seq[res.id[1]] = DNA3_TO_1.get(res.get_resname(), "N")
    contacted_n_positions = sorted(
        {info["dna_resnum"] for info in base_contacts.values() if info["dna_chain"] == "N"}
    )
    if contacted_n_positions:
        lo, hi = min(contacted_n_positions), max(contacted_n_positions)
        window_seq = "".join(chain_n_seq.get(p, "N") for p in range(lo, hi + 1))
        print(f"\nFenetre de bases contactees sur le brin N, positions {lo}-{hi} : 5'-{window_seq}-3'")
        marked = "".join(
            (chain_n_seq.get(p, "N") if p in contacted_n_positions else chain_n_seq.get(p, "N").lower())
            for p in range(lo, hi + 1)
        )
        print(f"(MAJUSCULE = position effectivement contactee par une base laterale) : {marked}")

    with open(os.path.join(OUT_DIR, "bldc_base_footprint.json"), "w") as fh:
        json.dump(
            {
                "base_contacts": {str(k): v for k, v in base_contacts.items()},
                "backbone_only_residues": sorted(backbone_only),
                "contacted_dna_window_chainN": (
                    {"start": lo, "end": hi, "seq": window_seq, "marked": marked}
                    if contacted_n_positions
                    else None
                ),
            },
            fh,
            indent=2,
        )
    print(f"\nEcrit : {os.path.join(OUT_DIR, 'bldc_base_footprint.json')}")


if __name__ == "__main__":
    main()
