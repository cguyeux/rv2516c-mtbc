#!/usr/bin/env python3
"""P3.2.1 -- Les 9 paires ex-aequo du null N2b sont-elles elles-memes des modules
fonctionnels connus (systemes toxine-antitoxine, operons caracterises) ?

CONTEXTE. P3.2 a montre que la co-occurrence parfaite Rv2516c/Rv2517c a travers 53
NTM (Jaccard = 1,000) n'est PAS unique dans sa classe de reference la mieux
appariee : sur les 49 paires de genes a la fois VOISINES et de meme bande de rarete
(3-10/53), 9 sont ex-aequo a 1,000. La conclusion tiree etait deflationniste : la
concordance parfaite serait un effet generique de voisinage, pas la signature d'un
module fonctionnel.

P3.2 avait elle-meme identifie le contre-argument qui pourrait renverser cette
lecture, sans le tester : "si plusieurs de ces 9 paires sont elles-memes des
systemes TA ou des operons connus, cela reformulerait la conclusion en
'concordance parfaite frequente PARCE QUE les voisins rares sont souvent des
modules fonctionnels'". C'est exactement ce que teste ce script.

METHODE. Le pool N2b est REGENERE a l'identique (memes sources, memes bornes,
meme fenetre exclue) plutot que relu d'un fichier, car P3.2 n'avait sauvegarde que
son resume. Chaque paire des 49 est ensuite annotee depuis `db.sqlite` de l'atlas,
et le test compare les 9 ex-aequo aux 40 non-ex-aequo par un Fisher exact bilateral.

DEUX PIEGES TRAITES EXPLICITEMENT :
 1. NON-INDEPENDANCE. Les "paires" sont des fenetres glissantes sur l'ordre
    genomique : Rv1678-1679, Rv1679-1680 et Rv1680-1681 sont trois fenetres du MEME
    bloc de quatre genes contigus, pas trois observations independantes. Le script
    rapporte les deux comptages, par paire ET par bloc genomique fusionne.
 2. BIAIS D'ANNOTATION. Un "module connu" est detecte par l'annotation existante ;
    les paires de deux hypothetiques sont donc NON CLASSABLES, pas negatives. Le
    test est fait sur le sous-ensemble annote, et la fraction non classable est
    rapportee, car elle plafonne ce que ce test peut conclure.

Sortie : resultats/p3_2_1_exaequo_annotation/
"""
from __future__ import annotations

import csv
import json
import sqlite3
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT.parent / "annotation_mtbc"
NTM_TSV = ATLAS / "résultats" / "phase57_ntm" / "ntm_presence_per_gene.tsv"
DB = ATLAS / "site" / "data" / "db.sqlite"
OUT = ROOT / "résultats" / "p3_2_1_exaequo_annotation"

GENE_A, GENE_B = "Rv2516c", "Rv2517c"
PREVALENCE_BAND = (3, 10)
EXCLUDE_WINDOW = {"Rv2513", "Rv2514c", "Rv2515c", "Rv2516c", "Rv2517c", "Rv2518c"}


def load_presence():
    rows = {}
    with open(NTM_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sp = frozenset(r["species_present"].split(",")) if r["species_present"] else frozenset()
            rows[r["rv"]] = {"species": sp, "n": int(r["n_ntm_present"])}
    return rows


def jaccard(a, b):
    if not a and not b:
        return float("nan")
    return len(a & b) / len(a | b)


def fisher_exact_2x2(a, b, c, d):
    """p bilateral par somme des tables de probabilite <= celle observee."""
    n = a + b + c + d
    r1, r2 = a + b, c + d
    c1 = a + c

    def p_of(x):
        y, z, w = r1 - x, c1 - x, r2 - (c1 - x)
        if min(y, z, w) < 0:
            return 0.0
        return comb(r1, x) * comb(r2, z) / comb(n, c1)

    p_obs = p_of(a)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return min(1.0, sum(p_of(x) for x in range(lo, hi + 1) if p_of(x) <= p_obs * (1 + 1e-9)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    presence = load_presence()

    con = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    order = [r[0] for r in con.execute(
        "SELECT rv FROM genes WHERE start_mtbc0 IS NOT NULL ORDER BY start_mtbc0")]
    meta = {r["rv"]: dict(r) for r in con.execute(
        "SELECT rv, gene_name, len_aa, strand, product_h37rv FROM genes")}
    con.close()

    # --- regeneration du pool N2b, a l'identique de P3.2 ---------------------
    pairs = [(order[i], order[i + 1]) for i in range(len(order) - 1)]
    pairs = [(g1, g2) for g1, g2 in pairs
             if g1 not in EXCLUDE_WINDOW and g2 not in EXCLUDE_WINDOW
             and g1 in presence and g2 in presence]
    n2b = [(g1, g2) for g1, g2 in pairs
           if PREVALENCE_BAND[0] <= presence[g1]["n"] <= PREVALENCE_BAND[1]
           and PREVALENCE_BAND[0] <= presence[g2]["n"] <= PREVALENCE_BAND[1]]

    scored = []
    for g1, g2 in n2b:
        j = jaccard(presence[g1]["species"], presence[g2]["species"])
        if j == j:
            scored.append({"g1": g1, "g2": g2, "jaccard": j})

    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    emit("P3.2.1 -- Les paires ex-aequo de N2b sont-elles des modules fonctionnels connus ?")
    emit("=" * 84)
    emit(f"Pool N2b regenere : {len(scored)} paires "
         f"(attendu 49 d'apres le resume de P3.2 : "
         f"{'CONCORDE' if len(scored) == 49 else 'DIVERGE, verifier'})")
    exa = [p for p in scored if p["jaccard"] == 1.0]
    emit(f"Paires ex-aequo a Jaccard = 1,000 : {len(exa)} (attendu 9 : "
         f"{'CONCORDE' if len(exa) == 9 else 'DIVERGE'})")
    emit("")

    # --- annotation ----------------------------------------------------------
    def classify(p):
        m1, m2 = meta.get(p["g1"], {}), meta.get(p["g2"], {})
        prods = f"{m1.get('product_h37rv') or ''} | {m2.get('product_h37rv') or ''}"
        names = [m1.get("gene_name"), m2.get("gene_name")]
        low = prods.lower()
        is_ta = "toxin" in low  # couvre 'toxin' et 'antitoxin'
        n_named = sum(1 for x in names if x)
        n_hypo = low.count("hypothetical protein")
        if is_ta:
            cls = "TA_valide"
        elif n_hypo == 2:
            cls = "NON_CLASSABLE_2_hypothetiques"
        elif n_named >= 1:
            cls = "au_moins_1_gene_nomme"
        else:
            cls = "annote_sans_nom"
        return cls, prods, names

    for p in scored:
        p["classe"], p["produits"], p["noms"] = classify(p)
        p["ex_aequo"] = p["jaccard"] == 1.0

    emit("LES 9 PAIRES EX-AEQUO, ANNOTEES")
    emit("-" * 84)
    for p in sorted(exa, key=lambda x: x["g1"]):
        emit(f"  {p['g1']:<9s} + {p['g2']:<9s}  [{p['classe']}]")
        emit(f"      {p['produits']}")
    emit("")

    # --- piege 1 : non-independance -----------------------------------------
    idx = {rv: i for i, rv in enumerate(order)}
    exa_sorted = sorted(exa, key=lambda p: idx[p["g1"]])
    blocks, cur = [], []
    for p in exa_sorted:
        if cur and idx[p["g1"]] <= idx[cur[-1]["g2"]]:
            cur.append(p)
        else:
            if cur:
                blocks.append(cur)
            cur = [p]
    if cur:
        blocks.append(cur)
    emit("PIEGE 1 -- NON-INDEPENDANCE DES PAIRES (fenetres glissantes)")
    emit("-" * 84)
    emit(f"  Les {len(exa)} paires ex-aequo se reduisent a {len(blocks)} BLOCS genomiques disjoints :")
    for bl in blocks:
        genes = sorted({g for p in bl for g in (p["g1"], p["g2"])}, key=lambda g: idx[g])
        tag = "TA_valide" if any(p["classe"] == "TA_valide" for p in bl) else bl[0]["classe"]
        emit(f"    {'-'.join(genes):<45s} ({len(bl)} paire(s))  [{tag}]")
    n_ta_blocks = sum(1 for bl in blocks if any(p["classe"] == "TA_valide" for p in bl))
    emit(f"  Blocs qui sont des systemes TA valides : {n_ta_blocks}/{len(blocks)}")
    emit("")

    # --- test d'enrichissement -----------------------------------------------
    nex = [p for p in scored if not p["ex_aequo"]]
    a = sum(1 for p in exa if p["classe"] == "TA_valide")
    b = len(exa) - a
    c = sum(1 for p in nex if p["classe"] == "TA_valide")
    d = len(nex) - c
    pval = fisher_exact_2x2(a, b, c, d)
    emit("TEST D'ENRICHISSEMENT EN SYSTEMES TA (Fisher exact bilateral)")
    emit("-" * 84)
    emit(f"                        TA valide    autre")
    emit(f"  ex-aequo (J=1,000)    {a:^9d}    {b:^5d}")
    emit(f"  reste du pool N2b     {c:^9d}    {d:^5d}")
    emit(f"  p = {pval:.4f}")
    emit(f"  Fraction TA : {a}/{len(exa)} = {100*a/len(exa):.0f}% chez les ex-aequo "
         f"contre {c}/{len(nex)} = {100*c/len(nex):.0f}% dans le reste du pool")
    emit("")

    # --- piege 2 : plafond d'annotation --------------------------------------
    nc_exa = sum(1 for p in exa if p["classe"] == "NON_CLASSABLE_2_hypothetiques")
    nc_all = sum(1 for p in scored if p["classe"] == "NON_CLASSABLE_2_hypothetiques")
    emit("PIEGE 2 -- PLAFOND D'ANNOTATION")
    emit("-" * 84)
    emit(f"  Paires faites de DEUX hypothetiques, donc non classables : "
         f"{nc_exa}/{len(exa)} chez les ex-aequo, {nc_all}/{len(scored)} dans tout le pool.")
    emit("  Ces paires ne sont PAS des negatifs : ce sont des inconnues. Elles plafonnent")
    emit("  ce que ce test peut etablir, dans les deux sens.")
    emit("")

    # --- le BON niveau de comparaison : la classe N2b vs le fond genomique ----
    # Le test ci-dessus compare ex-aequo vs reste DE N2b. Il est negatif. Mais la
    # question de P3.2 portait sur la classe de reference elle-meme : est-ce un
    # echantillon de paires voisines quelconques, ou une classe deja particuliere ?
    n2b_keys = {(p["g1"], p["g2"]) for p in scored}
    bg = [(g1, g2) for g1, g2 in pairs if (g1, g2) not in n2b_keys]

    def is_ta_pair(g1, g2):
        prods = f"{meta.get(g1, {}).get('product_h37rv') or ''} " \
                f"{meta.get(g2, {}).get('product_h37rv') or ''}"
        return "toxin" in prods.lower()

    ta_n2b = sum(1 for p in scored if p["classe"] == "TA_valide")
    ta_bg = sum(1 for g1, g2 in bg if is_ta_pair(g1, g2))
    p_class = fisher_exact_2x2(ta_n2b, len(scored) - ta_n2b, ta_bg, len(bg) - ta_bg)
    emit("LE BON NIVEAU DE COMPARAISON : la classe N2b vs toutes les autres paires voisines")
    emit("-" * 84)
    emit(f"  systemes TA parmi les paires voisines ET rares (N2b) : "
         f"{ta_n2b}/{len(scored)} = {100*ta_n2b/len(scored):.1f} %")
    emit(f"  systemes TA parmi toutes les autres paires voisines  : "
         f"{ta_bg}/{len(bg)} = {100*ta_bg/len(bg):.1f} %")
    emit(f"  Fisher exact bilateral : p = {p_class:.3g}   "
         f"(rapport de frequences x{(ta_n2b/len(scored))/(ta_bg/len(bg)):.0f})")
    emit("")

    # --- verdict --------------------------------------------------------------
    emit("=" * 84)
    emit("VERDICT -- en deux temps, car les deux niveaux de test ne disent PAS la meme chose")
    emit("")
    emit("  1. LA QUESTION LITTERALE DE LA PISTE : parmi les 49 paires de N2b, les 9 ex-aequo")
    emit(f"     sont-elles plus souvent des modules connus que les 40 autres ? NON. {a}/{len(exa)} = "
         f"{100*a/len(exa):.0f} % de systemes TA")
    emit(f"     chez les ex-aequo contre {c}/{len(nex)} = {100*c/len(nex):.0f} % dans le reste du pool, "
         f"p = {pval:.3f} : la")
    emit("     fraction est LA MEME. Atteindre une co-occurrence PARFAITE n'est donc pas, a")
    emit("     l'interieur de cette classe, un marqueur de modularite fonctionnelle. Le")
    emit("     contre-argument que P3.2 avait laisse ouvert n'est PAS confirme sous cette forme.")
    emit("")
    emit("  2. MAIS LE TEST DEPLACE AU BON NIVEAU EST, LUI, TRES POSITIF, et c'est le")
    emit("     resultat qui compte. La classe de reference ELLE-MEME -- les paires de genes")
    emit(f"     a la fois voisines et rares -- contient {100*ta_n2b/len(scored):.0f} % de systemes")
    emit(f"     toxine-antitoxine, contre {100*ta_bg/len(bg):.1f} % pour les autres paires voisines du")
    emit(f"     genome, soit un enrichissement d'un facteur {(ta_n2b/len(scored))/(ta_bg/len(bg)):.0f} "
         f"(p = {p_class:.3g}).")
    emit("     Autrement dit : la deflation de P3.2 reposait sur l'idee que le groupe de")
    emit("     comparaison etait banal. Il ne l'est pas. Etre un couple de genes rares et")
    emit("     voisins dans H37Rv, c'est deja appartenir a une population fortement enrichie")
    emit("     en modules fonctionnels co-herites.")
    emit("")
    emit("  CE QUE CELA CHANGE POUR Rv2516c-Rv2517c : la lecture deflationniste de P3.2")
    emit("  ('9/49 ex-aequo, donc banal') doit etre requalifiee. La paire n'est pas banale")
    emit("  parce qu'elle rejoint un groupe quelconque ; elle rejoint un groupe deja")
    emit("  atypique. Cela ne vient PAS de sa co-occurrence parfaite (point 1), mais de sa")
    emit("  rarete conjointe et de son adjacence. Le voisinage immediat contient d'ailleurs")
    emit("  deja un systeme TA valide (Rv2514c-Rv2515c, Tandon 2019).")
    emit("")
    emit("  CE QUE CELA NE CHANGE PAS, et il faut le dire aussi nettement : appartenir a une")
    emit("  classe enrichie en modules ne demontre PAS que Rv2516c-Rv2517c EST un module, ni")
    emit("  qu'il serait un systeme TA -- rien dans ce dossier ne suggere une toxine ou une")
    emit(f"  antitoxine. L'effectif reste minuscule ({len(exa)} paires, {len(blocks)} blocs "
         f"independants), et {nc_exa} des {len(exa)}")
    emit("  paires ex-aequo sont deux hypothetiques dont on ne sait rien. C'est un")
    emit("  deplacement de plausibilite, pas une preuve -- exactement le statut annonce")
    emit("  par la piste.")

    hdr = ["g1", "g2", "jaccard", "ex_aequo", "classe", "produits"]
    (OUT / "n2b_paires_annotees.tsv").write_text(
        "\t".join(hdr) + "\n"
        + "\n".join("\t".join(str(p[h]) for h in hdr)
                    for p in sorted(scored, key=lambda x: -x["jaccard"])) + "\n")
    (OUT / "resume.json").write_text(json.dumps({
        "n_pool_n2b": len(scored), "n_ex_aequo": len(exa),
        "ta_dans_n2b": ta_n2b, "ta_fond_voisins": ta_bg, "n_fond_voisins": len(bg),
        "fisher_p_classe_vs_fond": float(f"{p_class:.3g}"),
        "n_blocs_independants": len(blocks),
        "n_ta_parmi_ex_aequo": a, "n_ta_blocs": n_ta_blocks,
        "n_ta_reste_pool": c, "fisher_p": round(pval, 4),
        "n_non_classables_ex_aequo": nc_exa, "n_non_classables_pool": nc_all,
    }, indent=1, ensure_ascii=False))
    (OUT / "rapport.txt").write_text("\n".join(lines) + "\n")
    emit("")
    emit(f"Ecrit : {OUT}/")


if __name__ == "__main__":
    main()
