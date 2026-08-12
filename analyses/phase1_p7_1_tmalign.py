#!/usr/bin/env python3
"""phase1_p7_1_tmalign.py -- P7.1 : comparaison structurale CIBLEE de Rv2516c.

Question. La fiche atlas porte deux pistes de repli CONCURRENTES, toutes deux non
significatives : (i) liaison a l'ADN via un HTH de type Sigma70 region 4 (Pfam
tentatif sous le seuil de gathering, res. 98-121 ; nom UniProt ; regulateurs MerR
lies a leur promoteur parmi les hits Foldseek) ; (ii) repli de proteine ribosomale
S6, porte par les seuls meilleurs hits Foldseek en E-value. Un Foldseek AVEUGLE
rend un classement de TM-scores ~0,6-0,7 tous non significatifs, donc peu
interpretable. On teste ici la question directement : le modele de Rv2516c
ressemble-t-il davantage aux regions 4 des facteurs sigma de M. tuberculosis, aux
HTH de la famille MerR, ou aux permutants circulaires de S6 ?

Design. Deux garde-fous, sans lesquels le resultat ne veut rien dire.
  1. MATRICE COMPLETE, pas seulement requete-vs-references. Les trois classes de
     reference sont aussi comparees ENTRE ELLES. Si les regions 4 et les MerR se
     ressemblent deja fortement (ce sont tous des HTH), alors un score eleve de
     Rv2516c contre les deux dit "c'est un HTH", PAS "c'est un sigma". La ligne de
     base inter-classes est la seule facon de lire l'echelle.
  2. ALIGNEMENT NON CIRCULAIRE. On aligne le modele ENTIER (267 aa), jamais le
     fragment 98-121 pre-decoupe : ou tombe l'alignement sur la requete est une
     SORTIE du test, pas une entree. Si les references HTH se posent spontanement
     sur 98-121 et les S6 ailleurs, l'information est reelle ; si on avait
     charcute la requete a la fenetre Pfam, la conclusion aurait ete acquise
     d'avance.

TM-score normalise par la chaine de REFERENCE (la plus courte) : c'est la mesure
pertinente quand on cherche un petit domaine dans une proteine plus grande.

Entrees  : résultats/p7_1/structures/ (cf. phase1_p7_1_fetch_structures.sh), sigma_meta.json
Sorties  : résultats/p7_1/{tmalign_matrix.tsv, tmalign_query.tsv, alignment_spans.tsv}
Run: python analyses/phase1_p7_1_tmalign.py
"""
from __future__ import annotations
import json
import subprocess
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STRUCT = ROOT / "résultats" / "p7_1" / "structures"
DOMAINS = ROOT / "résultats" / "p7_1" / "domains"
OUT = ROOT / "résultats" / "p7_1"
TMALIGN = ROOT / "tools" / "TMalign"

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",
}

# Chaines telles que nommees par les hits Foldseek de la fiche atlas.
# La classe est ce qui est teste : HTH_MerR (ADN) vs S6 (ARN).
PDB_REFS = {
    "3hh0": ("A", "MerR", "MerR family regulator, B. cereus"),
    "5d8c": ("A", "MerR", "HiNmlR MerR-family regulator, DNA-bound"),
    "5d90": ("C", "MerR", "HiNmlR MerR-family regulator, DNA-bound"),
    "7b90": ("E", "S6", "circular permutant of ribosomal protein S6"),
    "7bff": ("E", "S6", "circular permutant of ribosomal protein S6"),
    "7bfd": ("K", "S6", "circular permutant of ribosomal protein S6"),
}


def parse_pdb_chains(path: Path) -> dict[str, list[str]]:
    """Lignes ATOM par chaine, acides amines standards seulement.

    Les structures MerR sont co-cristallisees avec leur promoteur : sans ce filtre
    on alignerait de l'ADN, et TM-align rendrait un score sur des nucleotides.
    """
    chains: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[17:20].strip() not in AA3:
            continue
        chains.setdefault(line[21], []).append(line)
    return chains


def write_chain(lines: list[str], dest: Path) -> int:
    dest.write_text("\n".join(lines) + "\nEND\n")
    return len({ln[22:27] for ln in lines})


def slice_residues(src: Path, lo: int, hi: int, dest: Path) -> int:
    """Extrait [lo, hi] (numerotation du modele) -- utilise UNIQUEMENT pour
    decouper le domaine r4 des facteurs sigma de reference, jamais la requete."""
    keep = []
    for line in src.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[17:20].strip() not in AA3:
            continue
        try:
            n = int(line[22:26])
        except ValueError:
            continue
        if lo <= n <= hi:
            keep.append(line)
    if not keep:
        return 0
    return write_chain(keep, dest)


def tmalign(a: Path, b: Path) -> dict | None:
    """TM-align(a, b). Renvoie les scores et l'alignement brut."""
    try:
        r = subprocess.run([str(TMALIGN), str(a), str(b)],
                           capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if r.returncode != 0:
        return None
    out = r.stdout
    res: dict = {"raw": out}
    for line in out.splitlines():
        s = line.strip()
        if s.startswith("Aligned length="):
            parts = s.replace("=", ",").split(",")
            try:
                res["aligned"] = int(parts[1])
                res["rmsd"] = float(parts[3])
                res["seqid"] = float(parts[-1])
            except (ValueError, IndexError):
                pass
        elif s.startswith("TM-score=") and "Chain_1" in s:
            res["tm_ref"] = float(s.split("=")[1].split("(")[0])
        elif s.startswith("TM-score=") and "Chain_2" in s:
            res["tm_query"] = float(s.split("=")[1].split("(")[0])
    return res if "tm_ref" in res else None


def aligned_span_on_second(raw: str) -> tuple[int, int, int] | None:
    """Ou l'alignement tombe-t-il sur la SECONDE structure ?

    TM-align imprime trois lignes : sequence 1 avec gaps, marqueurs, sequence 2
    avec gaps. On compte les residus non-gap de la sequence 2 pour retrouver la
    numerotation, et on borne les positions effectivement appariees. C'est la
    sortie qui dit si une reference HTH se pose ou non sur la fenetre 98-121.
    """
    lines = raw.splitlines()
    for i, line in enumerate(lines):
        if set(line.strip()) <= set(":. ") and line.strip() and i + 1 < len(lines):
            marks, seq2 = line, lines[i + 1]
            pos = 0
            hits = []
            for j, ch in enumerate(seq2):
                if ch != "-":
                    pos += 1
                    if j < len(marks) and marks[j] in ":.":
                        hits.append(pos)
            if hits:
                return min(hits), max(hits), len(hits)
    return None


def main() -> None:
    print("== P7.1 : comparaison structurale ciblee de Rv2516c ==")
    DOMAINS.mkdir(parents=True, exist_ok=True)
    meta = json.load(open(OUT / "sigma_meta.json"))

    entries: list[tuple[str, str, str, Path]] = []  # (id, classe, description, pdb)

    query = STRUCT / "query_Rv2516c.pdb"
    if not query.exists():
        raise SystemExit("Modele de la requete absent : lancer phase1_p7_1_fetch_structures.sh")

    # -- domaines region 4 des facteurs sigma (decoupes aux bornes Pfam) --
    for rv, m in sorted(meta.items()):
        src = STRUCT / f"sigma_{rv}.pdb"
        if not src.exists():
            print(f"  [absent] {rv} ({m['gene']}) : modele non recupere")
            continue
        dest = DOMAINS / f"r4_{rv}.pdb"
        n = slice_residues(src, m["r4_from"], m["r4_to"], dest)
        if n < 20:
            print(f"  [court]  {rv} : {n} residus extraits, ignore")
            continue
        entries.append((f"{m['gene']}_r4", "sigma_r4", f"{rv} {m['r4_name']} {m['r4_from']}-{m['r4_to']}", dest))

    # -- references experimentales MerR et S6 --
    for pid, (chain, klass, desc) in PDB_REFS.items():
        src = STRUCT / f"pdb_{pid}.pdb"
        if not src.exists():
            print(f"  [absent] {pid} : non telecharge")
            continue
        chains = parse_pdb_chains(src)
        use = chain if chain in chains else (sorted(chains)[0] if chains else None)
        if use is None:
            print(f"  [vide]   {pid} : aucune chaine proteique")
            continue
        if use != chain:
            print(f"  [note]   {pid} : chaine {chain} absente, repli sur {use}")
        dest = DOMAINS / f"{pid}_{use}.pdb"
        n = write_chain(chains[use], dest)
        entries.append((f"{pid}_{use}", klass, f"{desc} ({n} res.)", dest))

    print(f"\n  {len(entries)} references retenues "
          f"({sum(1 for e in entries if e[1] == 'sigma_r4')} sigma_r4, "
          f"{sum(1 for e in entries if e[1] == 'MerR')} MerR, "
          f"{sum(1 for e in entries if e[1] == 'S6')} S6)")

    # -- 1. requete vs chaque reference --
    rows = []
    for name, klass, desc, path in entries:
        r = tmalign(path, query)
        if not r:
            print(f"  [echec TM-align] {name}")
            continue
        span = aligned_span_on_second(r["raw"])
        rows.append({
            "reference": name, "classe": klass, "description": desc,
            "tm_norm_ref": r["tm_ref"], "tm_norm_query": r.get("tm_query"),
            "aligned": r.get("aligned"), "rmsd": r.get("rmsd"), "seqid": r.get("seqid"),
            "span_query": f"{span[0]}-{span[1]}" if span else "?",
            "n_span": span[2] if span else 0,
        })
    rows.sort(key=lambda d: -d["tm_norm_ref"])

    with open(OUT / "tmalign_query.tsv", "w") as fh:
        fh.write("reference\tclasse\ttm_norm_ref\ttm_norm_query\taligned\trmsd\tseqid\tspan_sur_query\tn_span\tdescription\n")
        for d in rows:
            fh.write(f"{d['reference']}\t{d['classe']}\t{d['tm_norm_ref']:.4f}\t"
                     f"{d['tm_norm_query']:.4f}\t{d['aligned']}\t{d['rmsd']:.2f}\t"
                     f"{d['seqid']:.3f}\t{d['span_query']}\t{d['n_span']}\t{d['description']}\n")

    print("\n-- Rv2516c vs references (TM-score normalise par la reference) --")
    print(f"{'reference':16} {'classe':9} {'TM':>6} {'RMSD':>6} {'alig':>5}  span sur la requete")
    for d in rows:
        print(f"{d['reference']:16} {d['classe']:9} {d['tm_norm_ref']:6.3f} "
              f"{d['rmsd']:6.2f} {d['aligned']:5d}  {d['span_query']}")

    # -- 2. ligne de base : references entre elles --
    base = []
    for (n1, k1, _, p1), (n2, k2, _, p2) in combinations(entries, 2):
        r = tmalign(p1, p2) if len(open(p1).readlines()) <= len(open(p2).readlines()) else tmalign(p2, p1)
        if not r:
            continue
        base.append({"a": n1, "b": n2, "ka": k1, "kb": k2, "tm": r["tm_ref"]})

    with open(OUT / "tmalign_matrix.tsv", "w") as fh:
        fh.write("a\tb\tclasse_a\tclasse_b\tpaire\ttm_norm_court\n")
        for d in base:
            paire = "intra" if d["ka"] == d["kb"] else "inter"
            fh.write(f"{d['a']}\t{d['b']}\t{d['ka']}\t{d['kb']}\t{paire}\t{d['tm']:.4f}\n")

    def stats(vals):
        if not vals:
            return "n/a"
        vals = sorted(vals)
        return f"n={len(vals)} med={vals[len(vals)//2]:.3f} [{vals[0]:.3f}-{vals[-1]:.3f}]"

    print("\n-- Ligne de base entre references (l'echelle de lecture) --")
    for ka in ("sigma_r4", "MerR", "S6"):
        intra = [d["tm"] for d in base if d["ka"] == ka and d["kb"] == ka]
        print(f"  intra-{ka:9} : {stats(intra)}")
    for ka, kb in (("sigma_r4", "MerR"), ("sigma_r4", "S6"), ("MerR", "S6")):
        inter = [d["tm"] for d in base if {d["ka"], d["kb"]} == {ka, kb}]
        print(f"  {ka} vs {kb:9} : {stats(inter)}")

    print("\n-- Rv2516c par classe --")
    for k in ("sigma_r4", "MerR", "S6"):
        v = [d["tm_norm_ref"] for d in rows if d["classe"] == k]
        print(f"  {k:9} : {stats(v)}")

    print(f"\nEcrit {OUT/'tmalign_query.tsv'} et {OUT/'tmalign_matrix.tsv'}")


if __name__ == "__main__":
    main()
