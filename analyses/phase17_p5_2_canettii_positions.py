#!/usr/bin/env python3
"""phase17_p5_2_canettii_positions.py -- P5.2, vérification secondaire annoncée par P5.2.1 :
où tombent les 5 substitutions canettii-consensus par rapport aux fenêtres de domaine ?

Explicitement une vérification FAIBLE (n=5), pas un test : P5.2.1 l'annonçait comme "à mentionner
sans en faire un argument". Ce script se contente de LOCALISER les 5 substitutions déjà comptées par
l'atlas (P10, data/Rv2516c_canettii_conservation.json), sans reproduire tout le calcul dN/dS.

Réutilise la machinerie de traduction déjà validée par le projet Canettii (`annotate_variants`), pour
ne pas dériver à la main l'effet synonyme/missense sur un gène BRIN MOINS -- piège classique de ce
genre de conversion (décalage de coordonnées, sens de lecture du codon).

Run: python analyses/phase17_p5_2_canettii_positions.py
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
CANETTII = MTBC / "Canettii"
sys.path.insert(0, str(CANETTII))
import analyse_dnds_categories as adc  # type: ignore[import-not-found]  # noqa: E402

GENOME = MTBC / "investigate_phylo" / "resources" / "NC_000962.3.fasta"
GFF3 = MTBC / "investigate_phylo" / "resources" / "NC_000962.3.gff3"
CONSENSUS = CANETTII / "canettii_vs_mtbc" / "canettii_consensus_variants.txt"

RV = "Rv2516c"
CDS_START, CDS_END = 2832710, 2833513   # coordonnees H37Rv verifiees sur le GFF3, brin moins

DOMAINS = {"unite_A_ferredoxin (1-87)": (1, 87), "unite_B_wHTH (88-147)": (88, 147),
           "linker (148-177)": (148, 177), "unite_C_Ig (178-267)": (178, 267)}


def domain_of(res: int) -> str:
    for name, (lo, hi) in DOMAINS.items():
        if lo <= res <= hi:
            return name
    return "hors_bornes"


def load_genome(path: Path) -> str:
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq)


def main() -> None:
    seq = load_genome(GENOME)
    genes = adc.load_gene_index(str(GFF3))          # LISTE triee par start, pas un dict
    spdis = adc.load_variant_file(str(CONSENSUS))    # liste de (spdi, freqs)

    annotated = adc.annotate_variants(spdis, genes, seq, "canettii")
    rows = [a for a in annotated if a.get("locus_tag") == RV]
    print(f"Substitutions canettii-consensus annotées sur {RV} : {len(rows)}\n")

    counts: dict[str, int] = {}
    for a in sorted(rows, key=lambda x: x["position"]):
        m = re.match(r"^([A-Z*])(\d+)([A-Z*])$", a.get("aa_change") or "")
        aa_pos = int(m.group(2)) if m else None
        dom = domain_of(aa_pos) if aa_pos else "?"
        print(f"  pos_genomique={a['position']}  {a['ref']}>{a['alt']}  "
              f"aa_change={a.get('aa_change') or '(non-coding/complex)'}  effet={a['effect']}  "
              f"résidu={aa_pos}  domaine={dom}")
        if dom != "?":
            counts[dom] = counts.get(dom, 0) + 1

    print("\nComptage par domaine :")
    for d in DOMAINS:
        print(f"  {d} : {counts.get(d, 0)}")

    # Verification de coherence avec l'agrege deja publie par l'atlas (P10).
    n_syn = sum(1 for a in rows if a["effect"] == "synonymous_variant")
    n_mis = sum(1 for a in rows if a["effect"] == "missense_variant")
    print(f"\nsyn={n_syn} missense={n_mis} total_codant={n_syn+n_mis} "
          f"(atlas P10 rapportait 3 syn + 2 nonsyn = 5)")


if __name__ == "__main__":
    main()
