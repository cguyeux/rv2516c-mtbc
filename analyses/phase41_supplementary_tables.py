#!/usr/bin/env python3
"""phase41_supplementary_tables.py -- Export four publication-ready supplementary tables from
already-computed, already-verified result files, with no new computation beyond sorting and
formatting: cross-genus conservation by non-tuberculous mycobacterial lineage, the McDonald-Kreitman
test contingency data, the IS6110 insertions found within the Rv2512c-Rv2518c genomic window, and the
genome-wide co-expression top-50/top-100 partner lists for Rv2516c in each of the two expression
compendia used in the Methods.

Inputs (relative to the project root):
  resultats/p3_1_ntm_orthologs/resume.json
  resultats/p3_1_1_gradient_par_unite/gradient.json
  resultats/p10_7_mk_test/mk_resultats.json
  resultats/p10_4_is6110_rd/is_positions.tsv
  resultats/p10_3_coexpression/gse166501_correlations.tsv
  resultats/p10_3_coexpression/gse71200_correlations.tsv

Outputs: article/supplementary_materials/Table_S1..S4_*.tsv
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "resultats" if (ROOT / "resultats").exists() else ROOT / "résultats"
OUT = ROOT / "article" / "supplementary_materials"
OUT.mkdir(parents=True, exist_ok=True)

# --- Table S1: cross-genus conservation ------------------------------------------------------

def table_s1():
    resume = json.loads((RES / "p3_1_ntm_orthologs" / "resume.json").read_text())
    gradient = json.loads((RES / "p3_1_1_gradient_par_unite" / "gradient.json").read_text())
    mac_clade = {"M_avium", "M_bouchedurhonense", "M_timonense"}

    rows = []
    for hit in resume["hits"]:
        sp = hit["species"]
        independent_lineage = (
            "MAC_clade_avium_bouchedurhonense_timonense" if sp in mac_clade else sp
        )
        rows.append({
            "species": sp,
            "independent_phylogenetic_lineage": independent_lineage,
            "is_mtbap": hit["is_mtbap"],
            "blast_percent_identity": hit["pident_blast"],
            "query_coverage_percent": hit["qcov"],
            "evalue": hit["evalue"],
            "n_aligned_positions": hit["n_aligned_positions"],
            "n_identical_positions": hit["n_identical_positions"],
        })

    with open(OUT / "Table_S1_cross_genus_conservation.tsv", "w", newline="") as f:
        f.write(
            "# Table S1. Cross-genus conservation of Rv2516c against six non-tuberculous "
            "mycobacterial orthologs.\n"
            "# M. avium, M. bouchedurhonense and M. timonense share an identical per-position "
            "identity profile (M. avium complex) and are treated as one independent "
            "phylogenetic lineage, leaving four independent lineages overall.\n"
            "# is_mtbap = member of the M. shinjukuense/MTBAP clade, the closest outgroup used "
            "in the housekeeping-gene control tree.\n"
        )
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
        f.write("\n# Per-unit identity, deduplicated to four independent lineages "
                 "(main text figure/table source).\n")
        f.write("unit\tlabel\tn_positions\tidentity_percent_4_lineages_dedup\t"
                 "identity_percent_6_species_raw\n")
        for g in gradient["gradient_par_unite_dedup_4_lignees"]:
            f.write(f"{g['unite']}\t{g['label']}\t{g['n_positions']}\t"
                     f"{g['identite_dedup_4lignees_pct']}\t{g['identite_brute_6especes_pct']}\n")
    print("Table_S1_cross_genus_conservation.tsv written")


# --- Table S2: McDonald-Kreitman test ---------------------------------------------------------

DEFINITION_LABELS = {
    "tous les sites variables": "all variable sites (Pn=243, Ps=78)",
    "les 9 SNP de la fiche d'atlas": "curated nine-polymorphic-site set (Pn=7, Ps=2)",
}


def table_s2():
    mk = json.loads((RES / "p10_7_mk_test" / "mk_resultats.json").read_text())

    with open(OUT / "Table_S2_mcdonald_kreitman_test.tsv", "w", newline="") as f:
        f.write(
            "# Table S2. McDonald-Kreitman test, Rv2516c intra-MTBC polymorphism versus "
            "divergence to the M. canettii outgroup.\n"
            "# Divergence substitutions (Dn=2 non-synonymous, Ds=3 synonymous), each mapped to "
            "its structural unit (A: ferredoxin-like, B: helix-turn-helix, C: Ig-like).\n"
        )
        f.write("substitution\ttype\tunit\n")
        for sub, kind, unit in mk["divergence"]["substitutions"]:
            f.write(f"{sub}\t{kind}\t{unit}\n")

        f.write("\n# Contingency table and test statistics, two alternative definitions of the "
                 "polymorphism denominator.\n")
        f.write("definition\tPn\tPs\tDn\tDs\tNI\tDoS\talpha\tp_fisher\n")
        for t in mk["tests"]:
            label = DEFINITION_LABELS.get(t["definition"], t["definition"])
            f.write(f"{label}\t{t['Pn']}\t{t['Ps']}\t{t['Dn']}\t{t['Ds']}\t"
                     f"{t['NI']}\t{t['DoS']}\t{t['alpha']}\t{t['p_fisher']}\n")

        f.write("\n# Exhaustive enumeration of statistical power: p-value for every possible "
                 "(Dn, Ds) outcome at Dn+Ds=5 fixed, computed before the observed value was read, "
                 "for each polymorphism definition.\n")
        f.write("definition\tDn\tDs\tNI\tp\n")
        for definition, block in mk["enumeration_issues"].items():
            label = DEFINITION_LABELS.get(definition, definition)
            for issue in block["issues"]:
                f.write(f"{label}\t{issue['Dn']}\t{issue['Ds']}\t{issue['NI']}\t{issue['p']}\n")
    print("Table_S2_mcdonald_kreitman_test.tsv written")


# --- Table S3: IS6110 insertions in the locus window ------------------------------------------

def table_s3():
    # Same constants as analyses/phase34_p10_4_is6110_rd.py (LOCUS +/- 10 kb; Rv2512c-Rv2518c
    # cassette bounds).
    locus_start, locus_end = 2_832_710, 2_833_761
    window_start, window_end = locus_start - 10_000, locus_end + 10_000
    cassette_start, cassette_end = 2_827_000, 2_836_000

    rows = []
    with open(RES / "p10_4_is6110_rd" / "is_positions.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["IS"] != "IS6110":
                continue
            pos = int(r["position"])
            if window_start <= pos <= window_end:
                rows.append({
                    "clade": r["clade"],
                    "strain": r["souche"],
                    "position_h37rv": pos,
                    "in_9kb_cassette_window": cassette_start <= pos <= cassette_end,
                })

    rows.sort(key=lambda r: r["position_h37rv"])
    with open(OUT / "Table_S3_is6110_insertions_locus_window.tsv", "w", newline="") as f:
        f.write(
            "# Table S3. Non-reference IS6110 insertions found within +/-10 kb of the "
            "Rv2516c-Rv2517c pair (H37Rv coordinates 2,827,000-2,836,000 for the 9 kb "
            "Rv2512c-Rv2518c cassette itself), out of 18,435 genome-wide non-reference IS6110 "
            "insertions scored across 5,169 MTBC genomes spanning 614 phylogenetic clades.\n"
            "# The Rv2516c-Rv2517c pair itself (2,832,710-2,833,761) carries zero insertions in "
            "this sample; the counts below fall in the flanking genes of the same cassette.\n"
        )
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        else:
            f.write("clade\tstrain\tposition_h37rv\tin_9kb_cassette_window\n")
    print(f"Table_S3_is6110_insertions_locus_window.tsv written ({len(rows)} rows)")


# --- Table S4: co-expression top partners -------------------------------------------------------

def table_s4():
    compendia = {
        "GSE166501": RES / "p10_3_coexpression" / "gse166501_correlations.tsv",
        "GSE71200": RES / "p10_3_coexpression" / "gse71200_correlations.tsv",
    }
    with open(OUT / "Table_S4_coexpression_top_partners.tsv", "w", newline="") as f:
        f.write(
            "# Table S4. Top-100 co-expression partners of Rv2516c in each of the two "
            "expression compendia used for the genome-wide, unbiased co-expression scan "
            "(3,911 genes each), ranked by Spearman correlation within each compendium. "
            "The top-50 subset used for the hypergeometric enrichment test against the four "
            "curated regulons (DosR, SigH, SigF, and the union of Rv2516c's five known "
            "upstream regulators) is flagged in the top50 column.\n"
        )
        f.write("compendium\trank\tpartner_locus\tproduct\tspearman_r\tp_value\tq_BH\ttop50\n")
        for name, path in compendia.items():
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                partners = [r for r in reader if r["focus"] == "Rv2516c"]
            partners.sort(key=lambda r: float(r["spearman"]), reverse=True)
            for rank, r in enumerate(partners[:100], start=1):
                f.write(f"{name}\t{rank}\t{r['partenaire']}\t{r['produit']}\t"
                         f"{r['spearman']}\t{r['p']}\t{r['q_BH']}\t{'yes' if rank <= 50 else 'no'}\n")
    print("Table_S4_coexpression_top_partners.tsv written")


if __name__ == "__main__":
    table_s1()
    table_s2()
    table_s3()
    table_s4()
