#!/usr/bin/env python3
"""phase25_p10_2_unit_c_fusion_timing.py -- P10.2 : l'unité C est-elle une acquisition récente,
restreinte à la branche MTBAP/MTBC, absente des NTM vrais ?

POURQUOI CE SCRIPT, ET CE QU'IL NE PREND PAS POUR ACQUIS.
P3.1 (`phase13_p3_1_ntm_orthologs.py`) a déjà tourné tblastn sur les 6 génomes NTM porteurs et gardé
UNIQUEMENT le meilleur HSP par bitscore par espèce (`best_hit`, max_hsps=3 mais un seul retenu). En
relisant `résultats/p3_1_ntm_orthologs/resume.json` sous l'angle des bornes qend plutôt que de
l'identité moyenne, un motif saute aux yeux : 5 des 6 espèces s'arrêtent avant ou tout juste au seuil
de l'unité C (178-267) -- trois d'entre elles (clade MAC) à qend=175 EXACTEMENT identique -- alors que
M_shinjukuense (MTBAP, le plus proche parent du MTBC parmi les 6) couvre les 267 résidus à 100%.

CE QUE CE MOTIF PEUT VOULOIR DIRE, ET CE QUI DISTINGUE CES LECTURES.
  1. Absence biologique réelle : l'unité C a été acquise APRÈS la divergence NTM/MTBAP -- fusion plus
     récente que celle des unités A+B, restreinte à la branche menant au MTBC.
  2. Artefact d'ANNOTATION/ASSEMBLAGE : le gène NTM est en réalité aussi long, mais l'assemblage est
     fragmenté à cet endroit (le contig s'arrête), ou l'annotation automatique tronque le modèle de
     gène.
  3. Simple décrochage de DIVERGENCE : l'unité C existe chez ces NTM mais a divergé au point de ne
     plus être détectable par tblastn à evalue=1e-5 sur la requête ENTIÈRE (267 aa) -- un signal faible
     dilué par le reste de la requête, exactement le type de dilution déjà documenté ailleurs sur ce
     projet (P2.1, interroger par domaine plutôt que la protéine entière).
Ce script cherche à départager ces trois lectures, PAS à confirmer la première par défaut.

QUATRE VÉRIFICATIONS INDÉPENDANTES, dans l'ordre où elles permettent d'éliminer une lecture :
  A. Ré-exécuter tblastn en gardant TOUS les HSP (pas seulement le meilleur) : un second HSP, même
     faible, qui reprendrait après un gap dans la fenêtre 178-267 serait un signal direct de lecture 3
     (diverge mais présent), pas de lecture 1.
  B. Recherche RESTREINTE à l'unité C seule (178-267) comme requête, evalue relâchée -- teste
     spécifiquement la dilution de lecture 3 : un signal qui n'apparaît qu'une fois l'unité C isolée de
     la protéine entière serait la même leçon méthodologique que P2.1, transposée ici.
  C. Fragmentation d'assemblage (lecture 2) : la position sstart/send du meilleur HSP est-elle proche
     de la fin du contig/scaffold sur lequel il tombe ? Un hit qui s'arrête à quelques centaines de
     paires de bases de la fin d'un contig est suspect d'assemblage tronqué, pas de vraie absence.
  D. Codon stop en phase (lecture 1 la plus directe) : traduire, dans le cadre de lecture et le brin
     de l'alignement, les ~300 nt du génome immédiatement après la fin du HSP -- un codon stop précoce
     dans ce cadre serait la signature la plus directe d'une vraie fin de gène à cet endroit ; son
     absence n'exclut rien (l'ORF pourrait continuer hors cadre ou sur l'autre brin) mais son PRÉSENCE
     est un argument fort pour la lecture 1.

Entrées : mêmes génomes que P3.1 (`bdd/hors_mtbc/<espèce>/**/genome.fna`), même requête
(`data/Rv2516c.faa`), mêmes seuils de présence (pident>=30, qcov>=50) pour le hit principal.
Sorties : résultats/p10_2_unit_c_fusion/{tous_hsp.json, recherche_unite_c_seule.json,
fragmentation_contig.json, codon_stop.json, verdict.md}
Run: python analyses/phase25_p10_2_unit_c_fusion_timing.py
"""
from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BDD_HORS = ROOT.parent / "bdd" / "hors_mtbc"
OUT = ROOT / "résultats" / "p10_2_unit_c_fusion"
QUERY_FULL = ROOT / "data" / "Rv2516c.faa"

SPECIES = ["M_avium", "M_bouchedurhonense", "M_celatum", "M_shinjukuense", "M_simiae", "M_timonense"]
MTBAP_MEMBER = "M_shinjukuense"

MIN_PIDENT = 30.0
MIN_QCOV = 50.0
QLEN = 267
UNIT_C = (178, 267)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def revcomp(seq: str) -> str:
    comp = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def translate(seq: str) -> str:
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3].upper()
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def find_genome(species: str) -> Path | None:
    hits = sorted(BDD_HORS.glob(f"{species}/**/genome.fna"))
    return hits[0] if hits else None


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    buf: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def run_tblastn_all_hsps(query_fasta: Path, genome: Path, tmp: Path, evalue: str) -> list[dict]:
    db = tmp / "db"
    subprocess.run(["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db)],
                    check=True, capture_output=True)
    out = tmp / "hits.tsv"
    fmt = "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq qlen slen"
    r = subprocess.run(
        ["tblastn", "-query", str(query_fasta), "-db", str(db), "-out", str(out),
         "-evalue", evalue, "-max_target_seqs", "20", "-max_hsps", "20", "-outfmt", fmt],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        raise RuntimeError(f"tblastn a echoue : {r.stderr[:300]}")
    rows = []
    text = out.read_text()
    if not text.strip():
        return rows
    for line in text.splitlines():
        f = line.split("\t")
        qstart, qend = int(f[4]), int(f[5])
        qcov = 100.0 * (qend - qstart + 1) / QLEN
        rows.append({"sseqid": f[1], "pident": float(f[2]), "length": int(f[3]),
                     "qstart": qstart, "qend": qend, "sstart": int(f[6]), "send": int(f[7]),
                     "evalue": float(f[8]), "bitscore": float(f[9]), "qseq": f[10], "sseq": f[11],
                     "qlen": int(f[12]), "slen": int(f[13]), "qcov": qcov})
    return rows


def best_hit(rows: list[dict]) -> dict | None:
    passed = [r for r in rows if r["pident"] >= MIN_PIDENT and r["qcov"] >= MIN_QCOV]
    pool = passed or rows
    return max(pool, key=lambda r: r["bitscore"]) if pool else None


def part_a_all_hsps(genomes: dict[str, Path]) -> dict:
    """Tous les HSP (pas seulement le meilleur) contre la protéine entière -- un second HSP qui
    reprend dans 178-267 après le meilleur serait un signal direct de divergence-mais-présence."""
    print("\n[A] Tous les HSP par espèce (protéine entière comme requête) --------------------------")
    result = {}
    for sp, genome in genomes.items():
        with tempfile.TemporaryDirectory() as td:
            rows = run_tblastn_all_hsps(QUERY_FULL, genome, Path(td), evalue="1e-5")
        best = best_hit(rows)
        others = [r for r in rows if r is not best]
        into_c = [r for r in others
                  if max(r["qstart"], UNIT_C[0]) <= min(r["qend"], UNIT_C[1])]
        print(f"  {sp:20s} {len(rows)} HSP au total, meilleur span={best['qstart']}-{best['qend']}"
              f" ({best['bitscore']:.0f} bits)" if best else f"  {sp:20s} aucun HSP")
        if into_c:
            for r in into_c:
                print(f"      -> HSP secondaire dans l'unité C : span={r['qstart']}-{r['qend']}"
                      f" pident={r['pident']:.1f}% evalue={r['evalue']:.2g} bits={r['bitscore']:.0f}")
        else:
            print("      -> aucun HSP secondaire, même faible, n'entre dans la fenêtre 178-267")
        result[sp] = {
            "n_hsp_total": len(rows),
            "best": {k: best[k] for k in ("qstart", "qend", "sstart", "send", "sseqid", "pident",
                                           "evalue", "bitscore")} if best else None,
            "hsp_secondaires_dans_unite_C": [
                {k: r[k] for k in ("qstart", "qend", "sstart", "send", "pident", "evalue", "bitscore")}
                for r in into_c
            ],
        }
    return result


def part_b_unit_c_only(genomes: dict[str, Path]) -> dict:
    """Requête restreinte à l'unité C seule (178-267), evalue relâchée -- teste si le signal était
    dilué par la protéine entière (même leçon que P2.1, transposée à la conservation NTM)."""
    print("\n[B] Recherche restreinte à l'unité C seule (178-267), evalue relâchée -------------------")
    full = load_fasta(QUERY_FULL)
    seq_full = next(iter(full.values()))
    seq_c = seq_full[UNIT_C[0] - 1:UNIT_C[1]]
    assert len(seq_c) == UNIT_C[1] - UNIT_C[0] + 1, "découpage unité C incohérent avec Rv2516c.faa"
    query_c = OUT / "unite_C_seule.faa"
    query_c.write_text(f">Rv2516c_uniteC_{UNIT_C[0]}-{UNIT_C[1]}\n{seq_c}\n")

    result = {}
    for sp, genome in genomes.items():
        with tempfile.TemporaryDirectory() as td:
            rows = run_tblastn_all_hsps(query_c, genome, Path(td), evalue="10")
        rows.sort(key=lambda r: -r["bitscore"])
        top = rows[:3]
        if top:
            print(f"  {sp:20s} {len(rows)} HSP (evalue<=10) contre l'unité C seule -- meilleur : "
                  f"pident={top[0]['pident']:.1f}% evalue={top[0]['evalue']:.2g} "
                  f"bits={top[0]['bitscore']:.0f} span_query={top[0]['qstart']}-{top[0]['qend']}")
        else:
            print(f"  {sp:20s} AUCUN HSP même à evalue<=10 -- négatif franc pour l'unité C isolée")
        result[sp] = [
            {k: r[k] for k in ("qstart", "qend", "sstart", "send", "pident", "evalue", "bitscore")}
            for r in top
        ]
    return result


def part_c_contig_fragmentation(genomes: dict[str, Path], best_hits: dict) -> dict:
    """La position du meilleur HSP est-elle proche de la fin de son contig/scaffold ? Un hit qui
    s'arrête a quelques centaines de pb de la fin d'un contig est suspect d'assemblage tronque."""
    print("\n[C] Distance entre la fin du meilleur HSP et la fin de son contig/scaffold -------------")
    result = {}
    for sp, genome in genomes.items():
        best = best_hits[sp]["best"]
        if best is None:
            continue
        seqs = load_fasta(genome)
        contig_len = len(seqs.get(best["sseqid"], ""))
        plus_strand = best["send"] >= best["sstart"]
        dist_to_contig_end = (contig_len - best["send"]) if plus_strand else (best["send"] - 1)
        suspect = dist_to_contig_end < 500
        print(f"  {sp:20s} contig {best['sseqid']} ({contig_len} pb), hit brin "
              f"{'plus' if plus_strand else 'moins'} se termine à {dist_to_contig_end} pb de la fin "
              f"du contig" + ("  <-- SUSPECT, assemblage possiblement tronqué ici" if suspect else ""))
        result[sp] = {"sseqid": best["sseqid"], "contig_len": contig_len,
                      "brin": "plus" if plus_strand else "moins",
                      "distance_fin_contig_pb": dist_to_contig_end, "suspect_troncature": suspect}
    return result


def part_d_stop_codon(genomes: dict[str, Path], best_hits: dict) -> dict:
    """Traduit ~300 nt du genome immediatement apres la fin du meilleur HSP, dans le cadre et le brin
    de l'alignement -- un codon stop precoce est le signal le plus direct d'une vraie fin de gene."""
    print("\n[D] Codon stop en phase juste après la fin du meilleur HSP -------------------------------")
    result = {}
    for sp, genome in genomes.items():
        best = best_hits[sp]["best"]
        if best is None:
            continue
        seqs = load_fasta(genome)
        contig = seqs.get(best["sseqid"], "")
        plus_strand = best["send"] >= best["sstart"]
        window = 303  # multiple de 3, ~101 codons
        if plus_strand:
            start = best["send"]  # 1-based inclusive fin du HSP -> premier nt apres = send (0-based)
            downstream = contig[start:start + window]
        else:
            end = best["send"] - 1  # 1-based -> 0-based, hit va vers les positions decroissantes
            raw = contig[max(0, end - window):end]
            downstream = revcomp(raw)
        aa = translate(downstream)
        stop_pos = aa.find("*")
        n_codons_avant_stop = stop_pos if stop_pos >= 0 else len(aa)
        print(f"  {sp:20s} {n_codons_avant_stop} codons traduits avant le premier stop "
              f"(sur {len(aa)} traduits)"
              + (f" -- stop au codon {stop_pos + 1}" if stop_pos >= 0 else " -- AUCUN stop dans la fenêtre"))
        result[sp] = {"aa_traduits_downstream": aa, "position_premier_stop": stop_pos,
                      "n_codons_avant_stop": n_codons_avant_stop}
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P10.2 : chronologie de fusion de l'unité C -- artefact ou acquisition récente ? ==")

    genomes = {}
    for sp in SPECIES:
        g = find_genome(sp)
        if g is None:
            print(f"  [ABSENT] {sp} introuvable")
            continue
        genomes[sp] = g

    res_a = part_a_all_hsps(genomes)
    (OUT / "tous_hsp.json").write_text(json.dumps(res_a, indent=1, ensure_ascii=False))

    res_b = part_b_unit_c_only(genomes)
    (OUT / "recherche_unite_c_seule.json").write_text(json.dumps(res_b, indent=1, ensure_ascii=False))

    res_c = part_c_contig_fragmentation(genomes, res_a)
    (OUT / "fragmentation_contig.json").write_text(json.dumps(res_c, indent=1, ensure_ascii=False))

    res_d = part_d_stop_codon(genomes, res_a)
    (OUT / "codon_stop.json").write_text(json.dumps(res_d, indent=1, ensure_ascii=False))

    print("\n[TERMINE] Sorties dans résultats/p10_2_unit_c_fusion/")


if __name__ == "__main__":
    main()
