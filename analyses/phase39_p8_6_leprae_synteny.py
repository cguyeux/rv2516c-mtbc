#!/usr/bin/env python3
"""P8.6 -- la cassette Rv2516c-Rv2517c survit-elle a la reduction genomique de M. leprae ?

Approche par SYNTENIE (genes flanquants), pas par le seul hit tblastn de P3.1 : le PAF
genome-entier deja sur disque (global_supplementary/hors_mtbc_align/M_leprae.paf) montre un
GRAND TROU non aligne (H37Rv 2 784 061-2 839 565, 55.5 kb) qui engloutit toute la fenetre
Rv2505c-Rv2523c, y compris orn/Rv2511 pourtant tres conserve (53/53 NTM) : l'aligneur
genome-entier (minimap2 asm) ne resout pas cette region, ce qui n'implique PAS une absence
biologique -- juste que la chaine de colinearite locale est trop perturbee (indels/pseudogenes
denses) pour etre chainee automatiquement. D'ou l'approche gene par gene de cette piste.

Methode : tblastn de chaque CDS proteique de la fenetre Rv2503c-Rv2524c (dont Rv2516c/Rv2517c)
contre le genome M. leprae complet (bdd/hors_mtbc/M_leprae/ref/genome.fna), meilleur HSP par
gene, verification d'un codon stop en phase immediatement apres la fin de l'alignement (meme
methode que P10.2 pour les orthologues NTM), classification ORF intact / pseudogene / absent.

Lit  : Canettii/NC_000962.3_CDS.fasta, bdd/hors_mtbc/M_leprae/ref/genome.fna
Ecrit: résultats/p8_6_leprae_synteny/synteny.tsv, résumé.md
"""
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
CDS = ROOT / "Canettii/NC_000962.3_CDS.fasta"
GENOME = ROOT / "bdd/hors_mtbc/M_leprae/ref/genome.fna"
OUT_DIR = Path(__file__).resolve().parent.parent / "résultats" / "p8_6_leprae_synteny"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# fenetre Rv2503c (housekeeping en amont) a Rv2524c (fas, housekeeping en aval), dans l'ordre
# genomique H37Rv croissant
WINDOW = [
    "Rv2503c", "Rv2504c", "Rv2505c", "Rv2506", "Rv2507", "Rv2508c", "Rv2509", "Rv2510c",
    "Rv2511", "Rv2512c", "Rv2513", "Rv2514c", "Rv2515c", "Rv2516c", "Rv2517c", "Rv2518c",
    "Rv2519", "Rv2520c", "Rv2521", "Rv2522c", "Rv2523c", "Rv2524c",
]

GENOME_BG_RATE = 0.50  # Cole et al. 2001 : ~50 % du contenu codant est pseudogenise / perdu


def load_cds():
    """locus_tag -> (proteine traduite, gene_name, location_str, strand, start, end)."""
    entries = {}
    header = None
    buf = []

    def flush():
        if header is None:
            return
        m = re.search(r"\[locus_tag=([^\]]+)\]", header)
        if not m:
            return
        lt = m.group(1)
        gm = re.search(r"\[gene=([^\]]+)\]", header)
        gene = gm.group(1) if gm else None
        loc = re.search(r"\[location=([^\]]+)\]", header)
        locstr = loc.group(1) if loc else ""
        strand = "minus" if "complement" in locstr else "plus"
        nums = [int(x) for x in re.findall(r"\d+", locstr)]
        start, end = (min(nums), max(nums)) if nums else (None, None)
        nt = "".join(buf)
        entries[lt] = {"nt": nt, "gene": gene, "loc": locstr, "strand": strand, "start": start, "end": end}

    with open(CDS) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header, buf = line, []
            else:
                buf.append(line.strip())
    flush()
    return entries


def translate(nt: str) -> str:
    from Bio.Seq import Seq

    try:
        return str(Seq(nt).translate(table=11))
    except Exception:  # noqa: BLE001
        return ""


def run_tblastn(query_fasta: Path) -> list[dict]:
    fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen sstrand"
    cmd = [
        "tblastn", "-query", str(query_fasta), "-subject", str(GENOME),
        "-outfmt", fmt, "-evalue", "1e-3", "-max_target_seqs", "5",
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    rows = []
    for line in res.stdout.strip().splitlines():
        p = line.split("\t")
        rows.append({
            "sseqid": p[1], "pident": float(p[2]), "length": int(p[3]),
            "qstart": int(p[6]), "qend": int(p[7]),
            "sstart": int(p[8]), "send": int(p[9]),
            "evalue": float(p[10]), "bitscore": float(p[11]),
            "qlen": int(p[12]), "sstrand": p[13],
        })
    return rows


def load_genome() -> str:
    seq = []
    with open(GENOME) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def revcomp(s: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return s.translate(comp)[::-1]


def check_stop_after(genome: str, hit: dict, extra_codons: int = 40) -> dict:
    """Traduit jusqu'a `extra_codons` codons en phase apres la fin du HSP et cherche un stop."""
    if hit["sstrand"] == "plus":
        pos = hit["send"]  # 1-based, dernier nt aligne
        window = genome[pos: pos + extra_codons * 3]
    else:
        pos = hit["send"]  # sstart>send pour le brin moins en tblastn sortie sstrand=minus
        window = revcomp(genome[max(0, pos - 1 - extra_codons * 3): pos - 1])
    codons = [window[i:i + 3] for i in range(0, len(window) - 2, 3)]
    stop_codons = {"TAA", "TAG", "TGA"}
    for i, c in enumerate(codons):
        if c in stop_codons:
            return {"stop_found": True, "codon_offset": i + 1}
    return {"stop_found": False, "codon_offset": None}


def classify(qcov_pct: float, pident: float, stop_info: dict) -> str:
    if qcov_pct < 30:
        return "absent (pas de hit couvrant)"
    if stop_info["stop_found"] and stop_info["codon_offset"] is not None and stop_info["codon_offset"] <= 3:
        return "pseudogene (stop immediat en phase)"
    if qcov_pct >= 80 and pident >= 30:
        return "ORF probable intact (couverture et identite suffisantes, pas de stop immediat)"
    return "fragment divergent / ambigu"


def main() -> None:
    entries = load_cds()
    genome = load_genome()

    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for lt in WINDOW:
            e = entries.get(lt)
            if e is None:
                rows.append({"locus": lt, "verdict": "CDS H37Rv introuvable"})
                continue
            prot = translate(e["nt"])
            qfa = Path(tmp) / f"{lt}.faa"
            qfa.write_text(f">{lt}\n{prot}\n")
            hits = run_tblastn(qfa)
            if not hits:
                rows.append({
                    "locus": lt, "gene": e["gene"], "h37rv_loc": e["loc"],
                    "qlen_aa": len(prot), "verdict": "absent (aucun hit tblastn e<=1e-3)",
                })
                continue
            best = max(hits, key=lambda h: h["bitscore"])
            qcov = 100.0 * (best["qend"] - best["qstart"] + 1) / best["qlen"]
            stop_info = check_stop_after(genome, best)
            verdict = classify(qcov, best["pident"], stop_info)
            rows.append({
                "locus": lt, "gene": e["gene"], "h37rv_loc": e["loc"], "qlen_aa": len(prot),
                "leprae_sstart": best["sstart"], "leprae_send": best["send"],
                "leprae_strand": best["sstrand"], "pident": round(best["pident"], 1),
                "qcov_pct": round(qcov, 1), "evalue": best["evalue"],
                "stop_apres_alignement_codon": stop_info["codon_offset"],
                "verdict": verdict,
            })

    # --- écriture TSV ---
    cols = ["locus", "gene", "h37rv_loc", "qlen_aa", "leprae_sstart", "leprae_send",
            "leprae_strand", "pident", "qcov_pct", "evalue", "stop_apres_alignement_codon", "verdict"]
    with open(OUT_DIR / "synteny.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # --- résumé ---
    lines = ["# P8.6 -- synténie Rv2516c-Rv2517c chez M. leprae (GCF_003584725.1)", ""]
    lines.append("| locus H37Rv | gène | position M. leprae | %identité | %couverture | stop après (codon) | verdict |")
    lines.append("|---|---|---|---|---|---|---|")
    for r in rows:
        pos = f"{r.get('leprae_sstart','?')}-{r.get('leprae_send','?')} ({r.get('leprae_strand','?')})" if "leprae_sstart" in r else "—"
        lines.append(
            f"| {r['locus']} | {r.get('gene','') or ''} | {pos} | {r.get('pident','')} "
            f"| {r.get('qcov_pct','')} | {r.get('stop_apres_alignement_codon','')} | {r['verdict']} |"
        )
    lines.append("")

    r2516 = next(r for r in rows if r["locus"] == "Rv2516c")
    r2517 = next(r for r in rows if r["locus"] == "Rv2517c")
    lines.append(f"**Rv2516c** : {r2516['verdict']}")
    lines.append(f"**Rv2517c** : {r2517['verdict']}")
    lines.append("")
    n_present = sum(1 for r in rows if "hit" not in r["verdict"] and "introuvable" not in r["verdict"] and "absent" not in r["verdict"])
    lines.append(
        f"Sur {len(WINDOW)} loci de la fenêtre Rv2503c-Rv2524c, {n_present} rendent un hit "
        f"tblastn exploitable chez M. leprae (seuil E<=1e-3). Taux de fond de pseudogénisation "
        f"du génome entier (Cole et al. 2001) : ~{int(GENOME_BG_RATE*100)} %."
    )
    (OUT_DIR / "resume.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
