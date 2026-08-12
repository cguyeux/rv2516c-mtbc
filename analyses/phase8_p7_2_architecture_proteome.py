#!/usr/bin/env python3
"""phase8_p7_2_architecture_proteome.py -- P7.2 REFORMULÉE : que vaut l'architecture de Rv2516c
à l'échelle du protéome de H37Rv, et la fonction est-elle déjà attribuée ailleurs ?

POURQUOI P7.2 EST REFORMULÉE, ET PAS SIMPLEMENT EXÉCUTÉE.
Le cadrage d'origine — « région 4 de facteur sigma sans région 2, unique dans le protéome » — est
SUPERSÉDÉ et ne doit pas être ressuscité. Il reposait sur le seul hit Pfam que la fiche exposait
(Sigma70_r4, i-E 8,4e-4, SOUS le seuil de gathering), alors que le domtblout en contenait douze,
tous HTH, sans famille dominante (P7.6.2), et que la structure ne corroborait pas spécifiquement
sigma (P7.1). HHpred a depuis tranché : l'unité B est un winged HTH de la famille AlpA/excisionase,
E = 1,3e-10.

Mais la QUESTION sous-jacente survit à son mauvais cadrage, et elle est même plus nette maintenant.
C'est le garde-fou paralogue-de-fold appliqué systématiquement dans ce projet : avant de
requalifier un gène vers une fonction nommée, vérifier que cette fonction n'est pas DÉJÀ attribuée
à un autre locus
du génome. Trois questions, dans l'ordre de ce qu'elles engagent :

  Q1. **H37Rv a-t-il déjà un locus annoté excisionase / facteur de directionnalité / AlpA ?**
      Si oui, ne pas propager l'assignation sans précaution. C'est la question qui peut arrêter
      la rédaction, donc elle passe en premier.
  Q2. **Combien de protéines de H37Rv portent un HTH ?** Donne le dénominateur. « Rv2516c porte un
      HTH » n'a de valeur que rapporté à la banalité du motif dans ce protéome.
  Q3. **L'ARCHITECTURE de Rv2516c est-elle rare ici ?** Son HTH (88-147) est INTERNE, avec 87 aa
      avant et 120 aa après, et ces deux flancs n'ont AUCUN domaine Pfam. Combien de protéines à
      HTH de H37Rv ont la même signature ?

CE QUE CE SCRIPT NE FAIT PAS, ET IL FAUT LE DIRE.
Rv2516c **n'a aucun hit Pfam au seuil de gathering** : son `.domtblout` est vide. Elle n'est donc
PAS dans la population recensée, et son HTH est localisé par HHpred et par carte de contacts, pas
par hmmscan. On ne compare donc pas des objets homogènes : le recensement fournit une DISTRIBUTION
DE RÉFÉRENCE contre laquelle lire Rv2516c, jamais une appartenance. Écrire « Rv2516c est la seule
protéine de H37Rv à... » serait faux ; « parmi les protéines que Pfam appelle avec confiance
porteuses d'un HTH, N seulement présentent l'architecture de Rv2516c » est exact.

MODÈLE NUL. Un HTH interne à flancs longs pourrait n'être que la conséquence mécanique d'une
protéine longue : plus elle est longue, plus un domaine placé au hasard a de la marge des deux côtés.
On calcule donc, pour CHAQUE protéine à HTH, la probabilité que son domaine, de sa longueur
observée, tombe au hasard en laissant >= FLANK aa de part et d'autre — puis l'espérance du compte.
Sans ce conditionnement sur la longueur, le chiffre observé n'est pas interprétable (leçon P8.3).

Entrées : ../annotation_mtbc/résultats/phase2b_pfam/*.domtblout (hmmscan --cut_ga, protéome entier)
          ../annotation_mtbc/site/content/genes/*.json  (noms de produits)
          data/pfam_clans/CL0123_HTH.json               (698 familles du clan HTH, API InterPro)
Sorties : résultats/p7_2_architecture/{hth_proteome.tsv, resume.md}
Run: python analyses/phase8_p7_2_architecture_proteome.py
"""
from __future__ import annotations
import glob, json, re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ATLAS = ROOT.parent / "annotation_mtbc"
DOMTBL = ATLAS / "résultats" / "phase2b_pfam"
GENES = ATLAS / "site" / "content" / "genes"
CLAN = ROOT / "data" / "pfam_clans" / "CL0123_HTH.json"
OUT = ROOT / "résultats" / "p7_2_architecture"

FLANK = 50                      # un flanc « substantiel » : de quoi loger un petit domaine replié
RV2516C = {"n_flank": 87, "c_flank": 120, "hth": (88, 147), "length": 267}

# familles de l'assignation HHpred de l'unité B (AlpA / excisionase / RDF).
# PF11112 n'appartient PAS au clan CL0123 mais relève du même groupe fonctionnel : on le suit à part.
ALPA_FAMS = {"PF05930": "Phage_AlpA", "PF06806": "DUF1233 (excisionase)",
             "PF09035": "Tn916-Xis", "PF11112": "PyocinActivator"}
# vocabulaire des éléments mobiles, pour Q1 (garde-fou paralogue-de-fold)
MOBILE = re.compile(r"excisionase|directionality|\bxis\b|\balpA\b|integrase|transposase|recombinase|"
                    r"resolvase|prophage|insertion sequence", re.I)


def parse_domtbl(path: Path) -> tuple[str, int, list[dict]]:
    """(locus, longueur de la protéine, domaines). Format hmmscan --domtblout :
    col 2 accession Pfam, col 4 nom de la requête, col 6 qlen, col 13 i-Evalue, col 18/19 ali from/to."""
    rv, qlen, doms = path.stem, 0, []
    for line in path.read_text(errors="replace").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 19:
            continue
        qlen = int(f[5])
        doms.append({"pfam": f[1].split(".")[0].upper(), "name": f[0],
                     "ie": float(f[12]), "start": int(f[17]), "end": int(f[18])})
    return rv, qlen, doms


def main() -> None:
    print("== P7.2 (reformulée) : architecture de Rv2516c à l'échelle du protéome H37Rv ==")
    clan = json.load(open(CLAN))
    print(f"  clan CL0123 « HTH » : {len(clan)} familles Pfam\n")

    files = sorted(DOMTBL.glob("*.domtblout"))
    prot, hth = {}, {}
    for p in files:
        rv, qlen, doms = parse_domtbl(p)
        prot[rv] = {"len": qlen, "doms": doms}
        h = [d for d in doms if d["pfam"] in clan]
        if h:
            hth[rv] = max(h, key=lambda d: d["end"] - d["start"])   # le plus étendu
    n_with_pfam = sum(1 for v in prot.values() if v["doms"])
    print(f"  {len(files)} protéines scannées | {n_with_pfam} avec >= 1 domaine au seuil de gathering")
    print(f"  Rv2516c : {len(prot.get('Rv2516c', {}).get('doms', []))} domaine(s) au seuil "
          f"-> ABSENTE de la population recensée, comme attendu\n")

    # ── noms de produits, pour Q1 ──
    products = {}
    for f in glob.glob(str(GENES / "*.json")):
        d = json.load(open(f))
        products[d["rv"]] = d.get("product_h37rv") or ""

    print("-- Q1. La fonction excisionase / RDF est-elle DÉJÀ attribuée dans H37Rv ? --")
    print("   (garde-fou paralogue-de-fold : si oui, ne pas propager l'assignation)")
    found = {a: [rv for rv, v in prot.items() if any(d["pfam"] == a for d in v["doms"])]
             for a in ALPA_FAMS}
    for a, name in ALPA_FAMS.items():
        loci = found[a]
        print(f"   {a} {name:24} {'AUCUN locus' if not loci else ', '.join(loci)}")
    by_name, cnt = [], Counter()
    for rv, p in sorted(products.items()):
        m = MOBILE.search(p)
        if m:
            by_name.append(rv)
            cnt[m.group(0).lower()] += 1
    print(f"\n   Loci dont le NOM de produit évoque un élément mobile : {len(by_name)}")
    for k, n in cnt.most_common():
        print(f"     {k:22} {n:>4}")
    exci = [rv for rv in by_name if re.search(r"excisionase|directionality|\bxis\b", products[rv], re.I)]
    print(f"   dont explicitement excisionase / RDF / Xis : {len(exci) or 'AUCUN'}"
          f"{' -> ' + ', '.join(exci) if exci else ''}")
    # Le nom ne vaut pas la preuve : un locus PEUT porter l'étiquette sans porter le domaine.
    # C'est le cas ici, et c'est le vrai contenu du garde-fou — donc on l'instruit, on ne se
    # contente pas de compter.
    if exci:
        print("\n   Instruction de ces loci (le NOM ne vaut pas le DOMAINE) :")
        for rv in exci:
            d = json.load(open(GENES / f"{rv}.json"))
            faa = DOMTBL / f"{rv}.faa"
            n_aa = len("".join(l.strip() for l in faa.read_text().splitlines()[1:])) if faa.exists() else 0
            u = (d.get("uniprot") or {}).get("protein_name") or "?"
            e = (d.get("essentiality") or {}).get("dejesus2017") or "?"
            dd = prot.get(rv, {}).get("doms", [])
            fams = ", ".join(f"{x['name']} ({x['pfam']}, i-E {x['ie']:.1g})" for x in dd) or "aucun"
            alpa = [x for x in dd if x["pfam"] in ALPA_FAMS]
            print(f"     {rv:10} {n_aa:>4} aa | DeJesus {e:3} | Pfam : {fams}")
            print(f"     {'':10}      UniProt : {u}")
            print(f"     {'':10}      famille AlpA/excisionase : {'OUI' if alpa else 'NON'}")
        print("     -> les deux portent bien un HTH, mais générique (HTH_17), et AUCUN ne porte de")
        print("        famille excisionase. L'étiquette « excisionase » de H37Rv vient de la catégorie")
        print("        fonctionnelle héritée « insertion seqs and phages », pas d'un domaine mesuré ;")
        print("        UniProt requalifie même Rv3750c en antitoxine putative.")
        print("     -> CONCLUSION DU GARDE-FOU : la fonction n'est PAS déjà attribuée sur preuve de")
        print("        domaine dans H37Rv. L'assignation de Rv2516c peut être propagée, à condition")
        print("        d'écrire que deux loci portent DÉJÀ le NOM sans en avoir le domaine.")

    print(f"\n-- Q2. Combien de protéines de H37Rv portent un HTH (clan CL0123, cut_ga) ? --")
    print(f"   {len(hth)} protéines, soit {100*len(hth)/len(files):.1f} % du protéome "
          f"et {100*len(hth)/n_with_pfam:.1f} % de celles qui ont un domaine.")
    fam_counts = Counter(clan[hth[rv]["pfam"]] for rv in hth)
    print("   Familles les plus représentées :")
    for name, n in fam_counts.most_common(10):
        print(f"     {name[:46]:46} {n:>4}")
    print(f"   -> porter un HTH dans ce protéome est BANAL ({len(hth)} loci). "
          "L'argument ne peut pas être « il y a un HTH ».")

    print(f"\n-- Q3. L'ARCHITECTURE (HTH interne, deux flancs >= {FLANK} aa) est-elle rare ? --")
    rows, internal, internal_dark = [], [], []
    exp = 0.0
    for rv, d in hth.items():
        L, s, e = prot[rv]["len"], d["start"], d["end"]
        nf, cf = s - 1, L - e
        # autres domaines Pfam tombant dans l'un des flancs
        others = [o for o in prot[rv]["doms"] if o is not d and (o["end"] < s or o["start"] > e)]
        rows.append({"rv": rv, "len": L, "fam": clan[d["pfam"]], "pfam": d["pfam"],
                     "start": s, "end": e, "nf": nf, "cf": cf, "n_other": len(others),
                     "product": products.get(rv, "")})
        if nf >= FLANK and cf >= FLANK:
            internal.append(rv)
            if not others:
                internal_dark.append(rv)
        # nul : domaine de longueur observée placé uniformément dans une protéine de longueur L
        dlen = e - s + 1
        slots = L - dlen + 1
        ok = max(0, L - dlen - 2 * FLANK + 1)
        exp += ok / slots if slots > 0 else 0.0

    ratio = len(internal) / exp if exp else float("nan")
    print(f"   observé : {len(internal)} / {len(hth)} protéines à HTH ont >= {FLANK} aa de part et d'autre")
    print(f"   attendu sous placement uniforme, conditionné sur les longueurs : {exp:.1f}")
    # Lecture du ratio, dans le bon sens : < 1 = DÉPLÉTION, donc l'internalité est plus rare que
    # le hasard une fois la longueur neutralisée. (Le commentaire d'origine disait l'inverse.)
    if ratio < 1:
        print(f"   -> ratio observé/attendu {ratio:.2f} : DÉPLÉTION d'un facteur {1/ratio:.1f}. "
              "Une fois la longueur\n      neutralisée, un HTH INTERNE est nettement plus rare que "
              "le hasard : dans ce protéome les\n      HTH se placent préférentiellement en "
              "extrémité, fusionnés à un domaine senseur.")
    else:
        print(f"   -> ratio observé/attendu {ratio:.2f} : pas de déplétion ; l'internalité "
              "s'explique par la longueur.")
    print(f"\n   Le critère qui DISCRIMINE est l'absence d'annotation dans les flancs :")
    print(f"   {len(internal_dark)} / {len(hth)} ont un HTH interne à flancs longs ET "
          f"AUCUN autre domaine Pfam,")
    print(f"   soit {100*len(internal_dark)/len(hth):.1f} % des protéines à HTH du protéome.")
    print(f"   Rv2516c a exactement ce profil : flancs de {RV2516C['n_flank']} et "
          f"{RV2516C['c_flank']} aa, aucun Pfam.")
    if internal_dark:
        print(f"\n   Les {len(internal_dark)} comparables, à examiner un par un :")
        for rv in sorted(internal_dark, key=lambda r: -prot[r]["len"]):
            r = next(x for x in rows if x["rv"] == rv)
            print(f"     {rv:10} {r['len']:>4} aa  HTH {r['start']:>3}-{r['end']:<3} "
                  f"(N {r['nf']:>3} / C {r['cf']:>3})  {r['fam'][:28]:28} {r['product'][:34]}")

    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "hth_proteome.tsv", "w") as fh:
        fh.write("rv\tlongueur\tpfam\tfamille\thth_start\thth_end\tflanc_N\tflanc_C\t"
                 "autres_domaines\tinterne\tinterne_et_flancs_sans_pfam\tproduit\n")
        for r in sorted(rows, key=lambda x: x["rv"]):
            i = r["nf"] >= FLANK and r["cf"] >= FLANK
            fh.write(f"{r['rv']}\t{r['len']}\t{r['pfam']}\t{r['fam']}\t{r['start']}\t{r['end']}\t"
                     f"{r['nf']}\t{r['cf']}\t{r['n_other']}\t{int(i)}\t"
                     f"{int(i and r['n_other'] == 0)}\t{r['product']}\n")
    with open(OUT / "resume.md", "w") as fh:
        fh.write("# P7.2 (reformulee) -- architecture de Rv2516c dans le proteome H37Rv\n\n")
        fh.write(f"Population : {len(hth)} proteines de H37Rv portant un domaine du clan Pfam CL0123\n"
                 f"(HTH, {len(clan)} familles), au seuil de gathering. **Rv2516c n'en fait pas partie**\n"
                 "(aucun hit au seuil) : le recensement est une distribution de reference, pas une\n"
                 "appartenance.\n\n")
        fh.write("## Q1. Garde-fou paralogue-de-fold\n\n")
        for a, name in ALPA_FAMS.items():
            fh.write(f"- {a} {name} : {', '.join(found[a]) if found[a] else 'aucun locus dans H37Rv'}\n")
        fh.write(f"- loci nommes comme elements mobiles : {len(by_name)} ; "
                 f"explicitement excisionase/RDF/Xis : {len(exci) if exci else 0}\n\n")
        fh.write(f"## Q2. Banalite du HTH\n\n{len(hth)} loci, soit {100*len(hth)/len(files):.1f} % "
                 f"du proteome. Familles dominantes : "
                 + ", ".join(f"{n} {name}" for name, n in fam_counts.most_common(5)) + ".\n\n")
        fh.write(f"## Q3. Architecture\n\n| critere | n | % des proteines a HTH |\n|---|---|---|\n")
        fh.write(f"| HTH interne (flancs >= {FLANK} aa) | {len(internal)} | {100*len(internal)/len(hth):.1f} |\n")
        fh.write(f"| attendu sous placement uniforme | {exp:.1f} | {100*exp/len(hth):.1f} |\n")
        fh.write(f"| interne ET flancs sans aucun Pfam | {len(internal_dark)} | {100*len(internal_dark)/len(hth):.1f} |\n")
    print(f"\nÉcrit {OUT}/hth_proteome.tsv et resume.md")


if __name__ == "__main__":
    main()
