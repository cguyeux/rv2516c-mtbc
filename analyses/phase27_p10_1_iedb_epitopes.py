#!/usr/bin/env python3
"""P10.1 -- Confondant epitope T sur la fenetre HTH invariante 98-121 de Rv2516c.

QUESTION. La fenetre 98-121, et en son coeur les positions 111-118 (RQRVHQLR)
invariantes sur 4 lignees NTM independantes (P3.1), est lue dans ce dossier comme
une conservation FONCTIONNELLE liee a la liaison a l'ADN. Or le MTBC est connu pour
maintenir ses epitopes T hyperconserves plutot que de les faire varier (paradoxe de
Comas et al. 2010) : une hyperconservation IMMUNOLOGIQUE produirait exactement le
meme signal de conservation, pour une raison entierement differente. Tant que ce
confondant n'est pas teste, il pese sur l'argument le plus important du dossier.

METHODE. Interrogation de l'API IEDB (query-api.iedb.org) sur l'accession UniProt
I6YDM0, volets lymphocytes T ET lymphocytes B. Les positions rendues par IEDB se
referent a l'antigene source CURE de chaque enregistrement (parfois un orthologue
M. bovis, parfois un construct fragmentaire) : elles ne sont PAS fiables comme
coordonnees H37Rv. Chaque peptide est donc RELOCALISE ici par recherche exacte dans
la sequence H37Rv de 267 aa, et l'ecart avec la position annoncee par IEDB est
rapporte.

GARDE-FOU, ecrit avant le calcul (enonce de la piste) : meme si un epitope
chevauchait la fenetre, cela ne refuterait pas la fonction de liaison a l'ADN --
les deux pressions peuvent coexister sur le meme patch de surface. Le resultat
nuance l'interpretation, il ne la tranche pas seul. Symetriquement, un
enregistrement "Negative" signifie teste-et-non-reactif dans UN panel d'alleles et
de donneurs, pas l'absence d'epitope pour tout HLA.

Sortie : resultats/p10_1_iedb_epitopes/
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p10_1_iedb_epitopes"

# Sequence H37Rv de Rv2516c (267 aa, UniProt I6YDM0, verifiee contre AFDB AF-I6YDM0-F1)
SEQ = ("MTADWVVTFTFDADPSMETMDAWETQLEGFDALVSRVPGHGIDVTVYAPGDWSVFDALAKMAGEVMPVVQAKSPIAVQ"
       "IISEPEHRLRAEAFTTPELMSAAEIADELGVSRQRVHQLRSTAGFPAPLADLRGGAVWDAAAVRRFAETWERKPGRPH"
       "TGTAKFAYSWAVGPAVGRSGKAPNVRWRVENPDKIRFVLRNIGDDIAEDVEIDLSRIDAITRNVPKKTVIRPGEGLNM"
       "VLIAAWGHPLPNQLYVRWAGQDEWAAVPLHPAH")

HTH_WINDOW = (98, 121)      # fenetre HTH conservee (P3.1/P2.3)
INVARIANT_CORE = (111, 118)  # RQRVHQLR, invariant sur les 4 lignees NTM (P3.1)
UNITS = {"A_ferredoxin": (1, 87), "B_wHTH": (88, 147),
         "linker": (148, 177), "C_Ig_like": (178, 267)}


def unit_of(start, end):
    hit = [n for n, (a, b) in UNITS.items() if not (end < a or start > b)]
    return "+".join(hit)


def overlap(r1, r2):
    lo, hi = max(r1[0], r2[0]), min(r1[1], r2[1])
    return max(0, hi - lo + 1)


def main():
    assert len(SEQ) == 267, len(SEQ)
    records = []
    for kind in ("tcell", "bcell"):
        p = OUT / f"{kind}_raw.json"
        rows = json.loads(p.read_text()) if p.exists() else []
        for r in rows:
            records.append((kind, r))

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("P10.1 -- Epitopes T/B references (IEDB) sur Rv2516c / I6YDM0")
    emit("=" * 78)
    n_t = sum(1 for k, _ in records if k == "tcell")
    n_b = sum(1 for k, _ in records if k == "bcell")
    emit(f"Enregistrements IEDB : {n_t} lymphocyte T, {n_b} lymphocyte B.")
    emit("")

    # --- deduplication par peptide, relocalisation dans H37Rv ---------------
    peptides = {}
    for kind, r in records:
        pep = r.get("linear_sequence")
        if not pep:
            continue
        cs = r.get("curated_source_antigen") or {}
        e = peptides.setdefault(pep, {"kind": kind, "assays": [], "iedb_pos": set()})
        e["assays"].append({
            "mesure": r.get("qualitative_measure"),
            "mhc_class": r.get("mhc_class"),
            "mhc": r.get("mhc_restriction"),
            "hote": r.get("host_organism_name"),
            "pmid": r.get("pubmed_id"),
            "antigene_source": cs.get("name"),
            "organisme_source": cs.get("source_organism_name"),
        })
        if cs.get("starting_position"):
            e["iedb_pos"].add((cs["starting_position"], cs["ending_position"]))

    rows_out = []
    emit("PEPTIDES DISTINCTS, RELOCALISES DANS LA SEQUENCE H37Rv")
    emit("-" * 78)
    for pep, e in sorted(peptides.items(), key=lambda kv: SEQ.find(kv[0])):
        idx = SEQ.find(pep)
        if idx < 0:
            emit(f"  {pep}  ABSENT de la sequence H37Rv (orthologue divergent ?) -- exclu")
            continue
        start, end = idx + 1, idx + len(pep)
        ov_w = overlap((start, end), HTH_WINDOW)
        ov_c = overlap((start, end), INVARIANT_CORE)
        pos_msg = ""
        if e["iedb_pos"]:
            declared = sorted(e["iedb_pos"])
            if all(d[0] != start for d in declared):
                pos_msg = f"   [IEDB annoncait {', '.join(f'{a}-{b}' for a, b in declared)}]"
        pos_all = [a["mesure"] for a in e["assays"] if (a["mesure"] or "").lower() == "positive"]
        neg_all = [a["mesure"] for a in e["assays"] if (a["mesure"] or "").lower() == "negative"]
        emit(f"  {pep}")
        emit(f"    H37Rv {start}-{end}  unite {unit_of(start, end)}{pos_msg}")
        emit(f"    chevauchement fenetre HTH {HTH_WINDOW[0]}-{HTH_WINDOW[1]} : {ov_w} residus"
             f"   |  coeur invariant {INVARIANT_CORE[0]}-{INVARIANT_CORE[1]} : {ov_c} residus")
        emit(f"    dosages : {len(pos_all)} POSITIF, {len(neg_all)} negatif")
        for a in e["assays"]:
            emit(f"      - {a['mesure']:<8s} hote {a['hote']}, MHC {a['mhc_class']}/{a['mhc']},"
                 f" PMID {a['pmid']}")
        emit("")
        rows_out.append({
            "peptide": pep, "h37rv_start": start, "h37rv_end": end,
            "unite": unit_of(start, end),
            "recouvre_fenetre_hth": ov_w, "recouvre_coeur_invariant": ov_c,
            "n_positif": len(pos_all), "n_negatif": len(neg_all),
            "hotes": ";".join(sorted({a["hote"] or "?" for a in e["assays"]})),
            "pmids": ";".join(sorted({str(a["pmid"]) for a in e["assays"] if a["pmid"]})) or "aucun",
        })

    # --- verdict ------------------------------------------------------------
    emit("=" * 78)
    emit("VERDICT SUR LE CONFONDANT")
    emit("")
    tested_window = [r for r in rows_out if r["recouvre_fenetre_hth"] > 0]
    tested_core = [r for r in rows_out if r["recouvre_coeur_invariant"] > 0]
    pos_window = [r for r in tested_window if r["n_positif"] > 0]

    if not tested_window:
        emit("  La fenetre HTH n'a JAMAIS ete testee. Absence de preuve, pas preuve d'absence :")
        emit("  le confondant reste ouvert.")
    else:
        cov = set()
        for r in tested_window:
            cov |= set(range(max(r["h37rv_start"], HTH_WINDOW[0]),
                             min(r["h37rv_end"], HTH_WINDOW[1]) + 1))
        emit(f"  La fenetre HTH {HTH_WINDOW[0]}-{HTH_WINDOW[1]} A ETE TESTEE experimentalement :")
        emit(f"  {len(tested_window)} peptides distincts la chevauchent, couvrant "
             f"{len(cov)}/{HTH_WINDOW[1]-HTH_WINDOW[0]+1} de ses residus "
             f"({min(cov)}-{max(cov)}).")
        core_cov = set()
        for r in tested_core:
            core_cov |= set(range(max(r["h37rv_start"], INVARIANT_CORE[0]),
                                  min(r["h37rv_end"], INVARIANT_CORE[1]) + 1))
        emit(f"  Le coeur invariant {INVARIANT_CORE[0]}-{INVARIANT_CORE[1]} (RQRVHQLR) est couvert "
             f"a {len(core_cov)}/{INVARIANT_CORE[1]-INVARIANT_CORE[0]+1} residus.")
        emit("")
        if pos_window:
            emit("  >>> AU MOINS UN DOSAGE POSITIF chevauche la fenetre : le confondant")
            emit("      immunologique est REEL et doit etre porte dans l'interpretation.")
        else:
            emit("  >>> TOUS les dosages chevauchant la fenetre sont NEGATIFS.")
            emit("      Le confondant 'hyperconservation immunologique' n'est pas soutenu")
            emit("      POUR CETTE FENETRE : la conservation de 98-121 n'est pas explicable")
            emit("      par une pression de maintien d'epitope T documentee. La lecture")
            emit("      fonctionnelle (liaison a l'ADN) en sort renforcee -- par elimination")
            emit("      d'une alternative, pas par preuve directe.")
    emit("")
    pos_any = [r for r in rows_out if r["n_positif"] > 0]
    if pos_any:
        emit("  Epitope(s) POSITIF(s) ailleurs dans la proteine :")
        for r in pos_any:
            emit(f"    {r['peptide']}  H37Rv {r['h37rv_start']}-{r['h37rv_end']}, "
                 f"unite {r['unite']}, hote {r['hotes']}, PMID {r['pmids']}")
    emit("")
    emit("  RESERVES, a porter dans toute citation :")
    emit("   - 'Negative' = teste et non reactif dans UN panel d'alleles/donneurs, pas")
    emit("     l'absence d'epitope pour tout HLA. Effectif tres faible (quelques dosages).")
    emit("   - Les peptides testes sont des 15-meres se chevauchant : ils couvrent la")
    emit("     fenetre mais n'epuisent pas tous les cadres de presentation possibles.")

    hdr = list(rows_out[0].keys())
    (OUT / "epitopes_mappes.tsv").write_text(
        "\t".join(hdr) + "\n"
        + "\n".join("\t".join(str(r[h]) for h in hdr) for r in rows_out) + "\n")
    (OUT / "rapport.txt").write_text("\n".join(lines) + "\n")
    emit("")
    emit(f"Ecrit : {OUT/'epitopes_mappes.tsv'}")
    emit(f"Ecrit : {OUT/'rapport.txt'}")


if __name__ == "__main__":
    main()
