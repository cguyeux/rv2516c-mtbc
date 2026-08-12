#!/usr/bin/env python3
"""P4.3bis-f, volet 1 -- complementarite electrostatique (APBS/PDB2PQR) de la face
concave de l'unite B, calibree AVANT application sur un panel temoin (BldC = positif
REEL, avec ses contacts ADN lus directement dans le co-cristal 6AMA ; unite A =
negatif de meme provenance que Rv2516c ; RdfS 8DGL et AlpA 8C3T = comparateurs de
repli supplementaires).

Methode, pour eviter la circularite holo/site (definir un site par la conservation puis demander
a un predicteur d'y loger un ligand) : le potentiel electrostatique de chaque structure est calcule SEUL
(PDB2PQR + APBS, PBE linearisee), sans jamais docker de duplex d'ADN -- ce n'est
donc structurellement PAS un test de co-repliement de type Boltz/AF, et le biais
demontre par P4.3bis-e (E1-E4) ne peut pas s'y appliquer. Le patch candidat de
chaque proteine est defini INDEPENDAMMENT de ce calcul (residus de contact reels
lus dans 6AMA pour BldC ; residus deja etablis par P2.3/P2.6 pour l'unite B, avant
ce script), puis calibre contre la distribution des fenetres glissantes de MEME
TAILLE sur toute la surface exposee de la MEME proteine -- pas un seuil absolu.

Etapes :
  1. Extraction des sous-structures (Biopython) depuis les fichiers deja sur disque
     ou telecharges pour cette piste.
  2. PDB2PQR (pdb2pqr30, force field PARSE, protonation pH 7.0) -> .pqr
  3. APBS (mg-auto, PBE linearisee, solvant implicite, 0.150 M NaCl, 298.15 K) -> .dx
  4. Lecture du potentiel electrostatique (OpenDX) par interpolation trilineaire aux
     positions des atomes de chaque residu.
  5. Calibration : fenetres glissantes de meme taille que le patch candidat sur les
     residus exposes (rSASA >= 25 %, Shrake-Rupley) de CHAQUE proteine ; rang du
     patch candidat dans sa propre distribution + comparaison inter-proteines.
"""
import os
import subprocess
import sys

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "résultats", "p4_3bis_f_orthogonal")
STRUCT_DIR = os.path.join(OUT, "structures")
WORK_DIR = os.path.join(OUT, "electro")
os.makedirs(WORK_DIR, exist_ok=True)

PDB2PQR = os.path.join(ROOT, "tools", "venv_electro", "bin", "pdb2pqr30")
APBS = os.path.join(ROOT, "tools", "apbs-3.4.1", "bin", "apbs")

AF_MODEL = os.path.join(ROOT, "résultats", "p7_1", "structures", "query_Rv2516c_AF.pdb")
BLDC_PDB = os.path.join(STRUCT_DIR, "6AMA.pdb")
ALPA_PDB = os.path.join(STRUCT_DIR, "8C3T.pdb")
RDFS_PDB = os.path.join(STRUCT_DIR, "8DGL.pdb")

DNA_CONTACT_CUTOFF = 4.5  # Angstrom, heavy-atom -- heavy-atom
SASA_EXPOSED_THRESHOLD = 0.25  # fraction de la SASA max theorique (Tien et al. 2013)

parser = PDBParser(QUIET=True)


class ChainSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

    def accept_residue(self, residue):
        return residue.id[0] == " "  # ATOM only, pas les heteroatomes/eau


class ResidueRangeSelect(Select):
    def __init__(self, chain_id, start, end):
        self.chain_id = chain_id
        self.start = start
        self.end = end

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " " and self.start <= residue.id[1] <= self.end


def save(structure, select, path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(path, select)


def extract_unit(af_model_path, start, end, out_path, label):
    s = parser.get_structure(label, af_model_path)
    save(s, ResidueRangeSelect("A", start, end), out_path)
    print(f"  {label} ({start}-{end}) -> {out_path}")
    return out_path


def extract_bldc_with_contacts(bldc_pdb_path, protein_chain, dna_chains, out_protein_path):
    """Extrait UNE copie de BldC (proteine seule, pas d'ADN passe a PDB2PQR/APBS),
    et retourne les numeros de residus dont un atome lourd est a moins de
    DNA_CONTACT_CUTOFF d'un atome lourd de l'ADN (lu sur le co-cristal reel)."""
    s = parser.get_structure("bldc_full", bldc_pdb_path)
    model = s[0]
    protein_atoms = [a for a in model[protein_chain].get_atoms() if a.element != "H"]
    dna_atoms = []
    for cid in dna_chains:
        dna_atoms.extend([a for a in model[cid].get_atoms() if a.element != "H"])
    dna_coords = np.array([a.coord for a in dna_atoms])

    contact_residues = set()
    for atom in protein_atoms:
        d = np.linalg.norm(dna_coords - atom.coord, axis=1)
        if d.min() <= DNA_CONTACT_CUTOFF:
            contact_residues.add(atom.get_parent().id[1])

    save(s, ChainSelect([protein_chain]), out_protein_path)
    print(f"  BldC chaine {protein_chain} -> {out_protein_path}; "
          f"{len(contact_residues)} residus de contact reel avec l'ADN (6AMA, <= {DNA_CONTACT_CUTOFF} A)")
    return out_protein_path, sorted(contact_residues)


def run_pdb2pqr(pdb_path, label):
    pqr_path = os.path.join(WORK_DIR, f"{label}.pqr")
    log_path = os.path.join(WORK_DIR, f"{label}_pdb2pqr.log")
    cmd = [
        PDB2PQR,
        "--ff=PARSE",
        "--with-ph=7.0",
        "--titration-state-method=propka",
        pdb_path,
        pqr_path,
    ]
    with open(log_path, "w") as log:
        r = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if r.returncode != 0 or not os.path.exists(pqr_path):
        raise RuntimeError(f"pdb2pqr30 a echoue pour {label}, voir {log_path}")
    print(f"  PDB2PQR {label} -> {pqr_path}")
    return pqr_path


def pqr_extent(pqr_path):
    coords = []
    with open(pqr_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    coords = np.array(coords)
    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)
    return center, span


def write_apbs_input(pqr_path, label):
    center, span = pqr_extent(pqr_path)
    # grille grossiere (coarse) englobant largement la molecule, grille fine centree,
    # espacement cible ~0.5 A, dimensions valides APBS (2^n c1 + 1)
    fine_dim = span + 20.0  # marge de 10 A de chaque cote
    coarse_dim = fine_dim * 1.7

    def valid_grid(n_target, spacing_target):
        # dimensions APBS : 2^k * c + 1 ; on prend c=1 ou 3, k tel que la dimension
        # couvre au moins n_target/spacing_target points
        import math
        for c in (97, 65, 49, 33):
            for k in range(0, 8):
                d = c * 2 ** k + 1
                if d * spacing_target >= n_target:
                    return d
        return 129

    dims = [valid_grid(fine_dim[i], 0.5) for i in range(3)]
    dime = f"{dims[0]} {dims[1]} {dims[2]}"

    in_path = os.path.join(WORK_DIR, f"{label}.in")
    dx_prefix = os.path.join(WORK_DIR, f"{label}_pot")
    content = f"""read
    mol pqr {pqr_path}
end
elec
    mg-auto
    dime {dime}
    cglen {coarse_dim[0]:.1f} {coarse_dim[1]:.1f} {coarse_dim[2]:.1f}
    fglen {fine_dim[0]:.1f} {fine_dim[1]:.1f} {fine_dim[2]:.1f}
    cgcent {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}
    fgcent {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}
    mol 1
    lpbe
    bcfl sdh
    pdie 4.0
    sdie 78.54
    srfm smol
    chgm spl2
    sdens 10.0
    srad 1.4
    swin 0.3
    temp 298.15
    ion charge 1 conc 0.150 radius 2.0
    ion charge -1 conc 0.150 radius 1.8
    calcenergy no
    calcforce no
    write pot dx {dx_prefix}
end
quit
"""
    with open(in_path, "w") as fh:
        fh.write(content)
    return in_path, dx_prefix + ".dx"


def run_apbs(pqr_path, label):
    in_path, dx_path = write_apbs_input(pqr_path, label)
    log_path = os.path.join(WORK_DIR, f"{label}_apbs.log")
    with open(log_path, "w") as log:
        r = subprocess.run([APBS, in_path], stdout=log, stderr=subprocess.STDOUT, cwd=WORK_DIR)
    if r.returncode != 0 or not os.path.exists(dx_path):
        raise RuntimeError(f"APBS a echoue pour {label}, voir {log_path}")
    print(f"  APBS {label} -> {dx_path}")
    return dx_path


def read_opendx(dx_path):
    """Lecteur minimal OpenDX pour grille reguliere APBS (format 'gridpositions
    regular', donnees en 3 colonnes par ligne, unites kT/e)."""
    with open(dx_path) as fh:
        lines = fh.readlines()
    dims = origin = deltas = None
    data_start = None
    n_expected = None
    deltas = []
    for i, line in enumerate(lines):
        if line.startswith("object 1 class gridpositions counts"):
            dims = tuple(int(x) for x in line.split()[-3:])
        elif line.startswith("origin"):
            origin = np.array([float(x) for x in line.split()[1:4]])
        elif line.startswith("delta"):
            deltas.append(np.array([float(x) for x in line.split()[1:4]]))
        elif line.startswith("object 3 class array"):
            n_expected = dims[0] * dims[1] * dims[2]
            data_start = i + 1
            break
    if dims is None or data_start is None:
        raise RuntimeError(f"Format OpenDX inattendu dans {dx_path}")
    delta = np.array([deltas[0][0], deltas[1][1], deltas[2][2]])
    values = []
    for line in lines[data_start:]:
        parts = line.split()
        if not parts or not parts[0][0].isdigit() and parts[0][0] not in "-+.":
            break
        try:
            values.extend(float(x) for x in parts)
        except ValueError:
            break
        if len(values) >= n_expected:
            break
    values = np.array(values[:n_expected]).reshape(dims)
    return values, origin, delta


def trilinear_interp(grid, origin, delta, point):
    idx = (point - origin) / delta
    i0 = np.floor(idx).astype(int)
    frac = idx - i0
    dims = grid.shape
    val = 0.0
    total_w = 0.0
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                ix, iy, iz = i0[0] + dx, i0[1] + dy, i0[2] + dz
                if 0 <= ix < dims[0] and 0 <= iy < dims[1] and 0 <= iz < dims[2]:
                    w = (
                        (frac[0] if dx else 1 - frac[0])
                        * (frac[1] if dy else 1 - frac[1])
                        * (frac[2] if dz else 1 - frac[2])
                    )
                    val += w * grid[ix, iy, iz]
                    total_w += w
    return val / total_w if total_w > 0 else np.nan


PROBE_DISTANCE = 3.4  # Angstrom -- approche la plus proche d'un oxygene de phosphate ADN


def domain_centroid(pdb_path, chain_id="A"):
    s = parser.get_structure("centroid", pdb_path)
    coords = np.array(
        [a.coord for a in s[0][chain_id].get_atoms() if a.element != "H"]
    )
    return coords.mean(axis=0)


def probe_potential(grid, origin, delta, atom_coord, centroid):
    """Potentiel a un point-sonde deporte de PROBE_DISTANCE vers l'exterieur de la
    proteine (loin du champ proche/singulier au contact immediat de la charge
    ponctuelle), approximant la position d'un atome d'ADN qui s'approcherait de ce
    residu -- pas le potentiel AU NIVEAU de l'atome lui-meme."""
    direction = atom_coord - centroid
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        return trilinear_interp(grid, origin, delta, atom_coord)
    direction = direction / norm
    probe = atom_coord + PROBE_DISTANCE * direction
    return trilinear_interp(grid, origin, delta, probe)


def residue_potential(pdb_for_numbering, grid, origin, delta, resnum, chain_id="A", centroid=None):
    """Potentiel moyen (kT/e) a des points-sondes au large des atomes lourds de
    chaine laterale d'un residu (voir probe_potential)."""
    s = parser.get_structure("x", pdb_for_numbering)
    if centroid is None:
        centroid = domain_centroid(pdb_for_numbering, chain_id)
    try:
        res = s[0][chain_id][(" ", resnum, " ")]
    except KeyError:
        return None
    vals = []
    for atom in res:
        if atom.element == "H" or atom.get_name() in ("N", "CA", "C", "O"):
            continue
        vals.append(probe_potential(grid, origin, delta, atom.coord, centroid))
    if not vals:
        for atom in res:
            if atom.element != "H":
                vals.append(probe_potential(grid, origin, delta, atom.coord, centroid))
    return float(np.mean(vals)) if vals else None


# ASA max theoriques (Tien et al. 2013, PLOS ONE, colonne "theoretical"), Angstrom^2
MAX_ASA = {
    "ALA": 129, "ARG": 274, "ASN": 195, "ASP": 193, "CYS": 167, "GLN": 225,
    "GLU": 223, "GLY": 104, "HIS": 224, "ILE": 197, "LEU": 201, "LYS": 236,
    "MET": 224, "PHE": 240, "PRO": 159, "SER": 155, "THR": 172, "TRP": 285,
    "TYR": 263, "VAL": 174,
}


def exposed_residues(pdb_path, chain_id="A"):
    s = parser.get_structure("sasa", pdb_path)
    sr = ShrakeRupley()
    sr.compute(s, level="R")
    exposed = []
    for res in s[0][chain_id]:
        if res.id[0] != " ":
            continue
        max_asa = MAX_ASA.get(res.get_resname())
        if max_asa is None:
            continue
        rel_sasa = res.sasa / max_asa
        if rel_sasa >= SASA_EXPOSED_THRESHOLD:
            exposed.append((res.id[1], rel_sasa))
    return exposed


def sliding_window_scores(pdb_path, grid, origin, delta, window_len, chain_id="A", centroid=None):
    if centroid is None:
        centroid = domain_centroid(pdb_path, chain_id)
    exposed = exposed_residues(pdb_path, chain_id)
    resnums = sorted(r for r, _ in exposed)
    scores = []
    for i in range(len(resnums) - window_len + 1):
        window = resnums[i : i + window_len]
        if window[-1] - window[0] != window_len - 1:
            continue  # pas contigu en sequence (chaine coupee) -- on saute
        vals = [
            residue_potential(pdb_path, grid, origin, delta, r, chain_id, centroid)
            for r in window
        ]
        vals = [v for v in vals if v is not None]
        if vals:
            scores.append((window[0], window[-1], float(np.mean(vals))))
    return scores


def process(label, pdb_path, patch_residues, chain_id="A", window_len=None):
    print(f"\n=== {label} ===")
    pqr = run_pdb2pqr(pdb_path, label)
    dx = run_apbs(pqr, label)
    grid, origin, delta = read_opendx(dx)
    centroid = domain_centroid(pdb_path, chain_id)

    patch_vals = [
        residue_potential(pdb_path, grid, origin, delta, r, chain_id, centroid)
        for r in patch_residues
    ]
    patch_vals = [v for v in patch_vals if v is not None]
    patch_mean = float(np.mean(patch_vals)) if patch_vals else None

    if window_len is None:
        window_len = len(patch_residues)
    windows = sliding_window_scores(pdb_path, grid, origin, delta, window_len, chain_id, centroid)
    window_scores = [w[2] for w in windows]
    percentile = (
        100.0 * sum(1 for w in window_scores if w <= patch_mean) / len(window_scores)
        if window_scores and patch_mean is not None
        else None
    )

    result = {
        "label": label,
        "patch_residues": patch_residues,
        "patch_mean_potential_kT_e": patch_mean,
        "n_windows_same_size": len(window_scores),
        "window_mean": float(np.mean(window_scores)) if window_scores else None,
        "window_max": float(np.max(window_scores)) if window_scores else None,
        "window_min": float(np.min(window_scores)) if window_scores else None,
        "patch_percentile_in_own_windows": percentile,
    }
    patch_str = f"{patch_mean:.2f} kT/e" if patch_mean is not None else "aucun patch candidat (calibration seule)"
    pct_str = f"{percentile:.1f}%" if percentile is not None else "n/a"
    print(
        f"  patch {patch_residues}: potentiel moyen = {patch_str} | "
        f"fenetres ({window_len} res, meme proteine) : moy={result['window_mean']:.2f} "
        f"min={result['window_min']:.2f} max={result['window_max']:.2f} "
        f"(n={result['n_windows_same_size']}) | percentile du patch = {pct_str}"
    )
    return result


def main():
    os.makedirs(STRUCT_DIR, exist_ok=True)
    print("Extraction des sous-structures...")
    unit_b_path = extract_unit(AF_MODEL, 88, 147, os.path.join(STRUCT_DIR, "unitB.pdb"), "unitB")
    unit_a_path = extract_unit(AF_MODEL, 1, 87, os.path.join(STRUCT_DIR, "unitA.pdb"), "unitA")
    bldc_path, bldc_contacts = extract_bldc_with_contacts(
        BLDC_PDB, "A", ["N", "R"], os.path.join(STRUCT_DIR, "bldc_chainA.pdb")
    )

    results = []
    # Patch candidat de l'unite B : residus de contact reel deja etablis par P2.6
    # (R113/H115/R118 parmi les plus proches de l'ADN dans le test Boltz) + R111
    # (charge de surface P2.3), AVANT et INDEPENDAMMENT de ce calcul.
    results.append(process("unitB", unit_b_path, [111, 113, 115, 118], window_len=4))
    # Unite A : aucun patch candidat connu -- calibration seule (fenetres glissantes)
    results.append(process("unitA", unit_a_path, [], window_len=4))
    # BldC : patch = contacts REELS lus dans le co-cristal 6AMA (pas une hypothese)
    if len(bldc_contacts) >= 2:
        results.append(process("bldc", bldc_path, bldc_contacts[:4], window_len=4))
        results.append(
            {"label": "bldc_all_contacts", "patch_residues": bldc_contacts, "note": "voir bldc_contacts.txt"}
        )
    with open(os.path.join(OUT, "bldc_contacts.txt"), "w") as fh:
        fh.write("Residus BldC (chaine A, 6AMA) a <= %.1f A d'un atome ADN (chaines N,R):\n" % DNA_CONTACT_CUTOFF)
        fh.write(",".join(str(r) for r in bldc_contacts) + "\n")

    import json

    with open(os.path.join(OUT, "electrostatics_results.json"), "w") as fh:
        json.dump(results, fh, indent=2, default=str)
    print(f"\nEcrit : {os.path.join(OUT, 'electrostatics_results.json')}")


if __name__ == "__main__":
    sys.exit(main())
