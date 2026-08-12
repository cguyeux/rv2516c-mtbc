#!/usr/bin/env python3
"""phase3_p7_6_1_foldseek_units.py -- P7.6.1 : recherche NON BIAISÉE, unité par unité.

Pourquoi interroger par MORCEAUX. Le Foldseek aveugle de l'atlas, lancé sur la
protéine entière, n'a jamais rendu que des hits vers la protéine ribosomale S6 :
le signal du modèle complet est capté par le domaine N-terminal, et tout le reste
de la protéine n'a jamais été interrogé. P7.6 a montré que plus de la moitié des
267 résidus n'a aucun candidat de repli. On redécoupe donc, et on interroge chaque
unité SÉPARÉMENT contre la base PDB complète.

Découpage OBJECTIF, et corrigé par rapport à P7.6. Le profil des contacts
inter-segments sur le modèle AlphaFold montre non pas deux mais TROIS unités,
séparées par deux zones de contacts strictement nuls (coupes 88-97 et 148-178) :
    A = 1-87     domaine globulaire N-terminal (repli ferredoxin-like, cf. P7.1)
    B = 88-147   unité portant le HTH prédit (98-121)
    L = 148-177  jonction : contacts nuls des deux côtés ET creux de pLDDT
                 (161-180 à 55,7 sur ESMFold) -- probablement un linker flexible
    C = 178-267  unité globulaire C-terminale, celle qui n'a JAMAIS été interrogée
Le « segment orphelin 122-267 » de P7.6 mélangeait donc la queue de B, le linker
et C : c'est ce mauvais découpage qu'on corrige ici.

Modèle : AlphaFold (v6, du cache reconstitué en P18.11) et non plus ESMFold. Les
deux prédicteurs donnent INDÉPENDAMMENT la même frontière principale (résidu 86) et
le même profil de pLDDT, ce qui lève le caveat de méthode noté dans l'état.

Seuils Foldseek (KB) : E < 1e-3 significatif, ~1e-2 suggestif, > 1 bruit.

Entrées : data/af_models/Rv2516c.pdb (via résultats/p7_1/structures/)
Sorties : résultats/p7_6_1/{units/*.pdb, foldseek_<unit>.m8, rapport.md}
Run: python analyses/phase3_p7_6_1_foldseek_units.py
"""
from __future__ import annotations
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))
from phase1_p7_1_tmalign import slice_residues  # noqa: E402

ATLAS = ROOT.parent / "annotation_mtbc"
FOLDSEEK = ATLAS / "tools" / "foldseek" / "bin" / "foldseek"
FSDB = ATLAS / "tools" / "foldseek_db" / "pdb"
QUERY = ROOT / "résultats" / "p7_1" / "structures" / "query_Rv2516c_AF.pdb"
OUT = ROOT / "résultats" / "p7_6_1"
UNITS_DIR = OUT / "units"

UNITS = {
    "full_1-267": (1, 267),
    "A_1-87": (1, 87),
    "B_88-147": (88, 147),
    "C_178-267": (178, 267),
    "BLC_88-267": (88, 267),      # tout sauf A : contrôle, montre ce que A masquait
    "Corphan_122-267": (122, 267),  # le découpage de P7.6, pour comparaison directe
}

FMT = "query,target,fident,alnlen,qstart,qend,tstart,tend,evalue,bits,prob,alntmscore"


def run_foldseek(pdb: Path, tag: str) -> list[dict]:
    m8 = OUT / f"foldseek_{tag}.m8"
    tmp = OUT / "_tmp"
    cmd = [str(FOLDSEEK), "easy-search", str(pdb), str(FSDB), str(m8), str(tmp),
           "--format-output", FMT, "-e", "10", "--max-seqs", "2000", "-v", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [échec foldseek] {tag} : {r.stderr.strip()[:200]}")
        return []
    rows = []
    for line in m8.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        rows.append({"target": f[1], "alnlen": int(f[3]), "qstart": int(f[4]), "qend": int(f[5]),
                     "evalue": float(f[8]), "bits": float(f[9]), "prob": float(f[10]),
                     "tm": float(f[11])})
    rows.sort(key=lambda d: d["evalue"])
    return rows


def verdict(e: float) -> str:
    return "SIGNIFICATIF" if e < 1e-3 else "suggestif" if e < 1e-2 else "faible" if e < 1 else "bruit"


def main() -> None:
    print("== P7.6.1 : Foldseek local, unité par unité ==")
    if not FOLDSEEK.exists() or not Path(str(FSDB) + "_ca.index").exists():
        raise SystemExit(f"Foldseek ou sa base introuvable : {FOLDSEEK} / {FSDB}")
    UNITS_DIR.mkdir(parents=True, exist_ok=True)

    results = {}
    for tag, (lo, hi) in UNITS.items():
        dest = UNITS_DIR / f"{tag}.pdb"
        n = slice_residues(QUERY, lo, hi, dest)
        if n < 20:
            print(f"  {tag} : {n} résidus, ignoré")
            continue
        rows = run_foldseek(dest, tag)
        results[tag] = rows
        best = rows[0] if rows else None
        print(f"\n-- {tag} ({n} résidus) : {len(rows)} hits --")
        if not best:
            print("   AUCUN hit, même au seuil permissif E<10")
            continue
        for d in rows[:6]:
            print(f"   E={d['evalue']:.2e} [{verdict(d['evalue']):12}] TM={d['tm']:.3f} "
                  f"q{d['qstart']}-{d['qend']}  {d['target'][:52]}")

    # -- rapport --
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "rapport.md", "w") as fh:
        fh.write("# P7.6.1 -- Foldseek local, unite par unite (modele AlphaFold v6)\n\n")
        fh.write("Seuils : E < 1e-3 significatif, ~1e-2 suggestif, > 1 bruit.\n\n")
        fh.write("| unite | n hits | meilleur E | verdict | TM | span requete | cible |\n")
        fh.write("|---|---|---|---|---|---|---|\n")
        for tag, rows in results.items():
            if not rows:
                fh.write(f"| {tag} | 0 | - | aucun hit | - | - | - |\n")
                continue
            d = rows[0]
            fh.write(f"| {tag} | {len(rows)} | {d['evalue']:.2e} | {verdict(d['evalue'])} | "
                     f"{d['tm']:.3f} | {d['qstart']}-{d['qend']} | {d['target'][:60]} |\n")
    print(f"\nÉcrit {OUT}/rapport.md")


if __name__ == "__main__":
    main()
