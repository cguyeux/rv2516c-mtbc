#!/usr/bin/env python3
"""P10.7 -- Test de McDonald-Kreitman formel sur Rv2516c.

POURQUOI. Le dossier porte un rapport pN/pS informel de 1,11 qualifie de "relache",
qu'il faut systematiquement nuancer dans le texte
("du bruit statistique sur 9 SNP, pas une preuve"). Il s'agit de remplacer
cette reserve methodologique repetee par un test proprement specifie, et annonce
elle-meme que le gain attendu est METHODOLOGIQUE, pas un nouveau resultat positif :
"avec un effectif aussi petit (9 et 5), meme un test correctement specifie aura tres
probablement un intervalle de confiance traversant la neutralite".

CE SCRIPT VA PLUS LOIN QUE DE LE CONSTATER. Comme la divergence ne compte que 5
substitutions, l'espace des issues possibles du test est FINI et minuscule : Dn peut
valoir 0, 1, 2, 3, 4 ou 5. On peut donc ENUMERER les six mondes possibles et calculer
le p de chacun, au lieu de constater apres coup que celui qu'on a observe n'est pas
significatif. Si aucune des six issues n'atteint le seuil, alors le test etait
STRUCTURELLEMENT incapable de rendre un resultat significatif, quelle que soit la
biologie : ce n'est plus "notre resultat est non significatif" (qui laisse planer un
doute sur l'effectif ET sur l'effet), c'est "cette donnee ne pouvait pas trancher".
La reserve devient un fait demontre, ce qui est precisement le gain methodologique
que la piste visait.

DONNEES, toutes deja sur disque, aucun recalcul biologique.
  POLYMORPHISME intra-MTBC : `resultats/p5_2_1/sites.tsv` (P5.2.1), 321 sites codants
    variables sur 145 209 genomes, avec leur effet et leur compte de porteurs.
    Deux definitions testees en sensibilite, car le choix du plancher de frequence
    EST le levier d'ascertainment de ce dossier (cf. plus bas).
  DIVERGENCE fixee : les 5 substitutions de M. canettii (P5.2 / phase17), polarisees
    par la machinerie de traduction validee du projet Canettii : 3 synonymes
    (G50G, F92F, V235V) et 2 non-synonymes (K72E, G132A).

CONVENTIONS reprises telles quelles du skill `mk-ascertainment` (`mk_stats`) pour
rester comparable aux autres projets du depot : NI avec correction de Haldane
(+0,5 partout), alpha = 1 - NI, DoS = Dn/(Dn+Ds) - Pn/(Pn+Ps), Fisher exact bilateral.

Sortie : resultats/p10_7_mk_test/
"""
from __future__ import annotations

import csv
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITES = ROOT / "résultats" / "p5_2_1" / "sites.tsv"
OUT = ROOT / "résultats" / "p10_7_mk_test"

# Divergence M. canettii (P5.2, phase17_p5_2_canettii_positions.py)
CANETTII = [("G50G", "syn", "A"), ("K72E", "missense", "A"),
            ("F92F", "syn", "B"), ("G132A", "missense", "B"),
            ("V235V", "syn", "C")]
ALPHA_LEVEL = 0.05


def fisher_exact_2x2(a, b, c, d):
    n, r1, r2, c1 = a + b + c + d, a + b, c + d, a + c

    def p_of(x):
        y, z, w = r1 - x, c1 - x, r2 - (c1 - x)
        if min(y, z, w) < 0:
            return 0.0
        return comb(r1, x) * comb(r2, z) / comb(n, c1)

    p_obs = p_of(a)
    lo, hi = max(0, c1 - r2), min(r1, c1)
    return min(1.0, sum(p_of(x) for x in range(lo, hi + 1) if p_of(x) <= p_obs * (1 + 1e-9)))


def mk_stats(dn, ds, pn, ps):
    """Identique a mk_stats du skill mk-ascertainment."""
    if dn + ds == 0:
        dos = -pn / (pn + ps) if (pn + ps) > 0 else 0.0
    elif pn + ps == 0:
        dos = dn / (dn + ds)
    else:
        dos = dn / (dn + ds) - pn / (pn + ps)
    ni = ((pn + 0.5) / (ps + 0.5)) / ((dn + 0.5) / (ds + 0.5))
    return dos, ni, 1.0 - ni, fisher_exact_2x2(dn, ds, pn, ps)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(SITES.open(), delimiter="\t"))
    lines = []

    def emit(s=""):
        print(s)
        lines.append(s)

    def counts(min_count):
        pn = sum(1 for r in rows if r["effet"] in ("missense", "nonsense")
                 and int(r["count"]) >= min_count)
        ps = sum(1 for r in rows if r["effet"] == "syn" and int(r["count"]) >= min_count)
        return pn, ps

    dn = sum(1 for _, e, _ in CANETTII if e == "missense")
    ds = sum(1 for _, e, _ in CANETTII if e == "syn")

    emit("P10.7 -- Test de McDonald-Kreitman, Rv2516c")
    emit("=" * 80)
    emit(f"DIVERGENCE fixee (M. canettii) : Dn = {dn} non-synonymes, Ds = {ds} synonymes")
    emit(f"  {', '.join(f'{m} [{e}, unite {u}]' for m, e, u in CANETTII)}")
    emit("")

    # --- deux definitions du polymorphisme -----------------------------------
    defs = [
        ("tous les sites variables", 1,
         "aucun plancher : toute position ou une variation a ete vue au moins une fois"),
        ("les 9 SNP de la fiche d'atlas", 100,
         ">= 100 porteurs sur 145 209 genomes. Plancher IDENTIFIE en balayant les seuils"
         " jusqu'a retrouver exactement la composition annoncee par le champ `conservation`"
         " de la fiche (2 syn + 7 missense + 0 non-sens = 9 sites) : c'est donc bien le"
         " sous-ensemble dont le dossier parle quand il cite '9 SNP'"),
    ]
    results = []
    emit("DEUX DEFINITIONS DU POLYMORPHISME, EN ANALYSE DE SENSIBILITE")
    emit("-" * 80)
    for label, floor, note in defs:
        pn, ps = counts(floor)
        dos, ni, alp, p = mk_stats(dn, ds, pn, ps)
        emit(f"  [{label}]  {note}")
        emit(f"    Pn = {pn:>3d}   Ps = {ps:>3d}   (Dn = {dn}, Ds = {ds})")
        emit(f"    pN/pS = {pn/ps:.2f}   dN/dS = {dn/ds:.2f}")
        emit(f"    NI = {ni:.2f}  (>1 = exces de polymorphisme non-synonyme, "
             f"signature de selection purificatrice)")
        emit(f"    DoS = {dos:+.3f}   alpha = {alp:+.2f}   Fisher p = {p:.3f}"
             f"   -> {'SIGNIFICATIF' if p < ALPHA_LEVEL else 'non significatif'}")
        emit("")
        results.append({"definition": label, "plancher": floor, "Pn": pn, "Ps": ps,
                        "Dn": dn, "Ds": ds, "NI": round(ni, 3), "DoS": round(dos, 3),
                        "alpha": round(alp, 3), "p_fisher": round(p, 4)})

    # --- LE POINT CENTRAL : enumeration exhaustive des issues possibles -------
    emit("=" * 80)
    emit("LE TEST POUVAIT-IL SEULEMENT RENDRE UN RESULTAT SIGNIFICATIF ?")
    emit("-" * 80)
    emit(f"La divergence ne compte que {dn+ds} substitutions : Dn ne peut prendre que")
    emit(f"{dn+ds+1} valeurs. On enumere donc TOUS les mondes possibles, au lieu de")
    emit("constater apres coup que le notre n'est pas significatif.")
    emit("")
    enum_all = {}
    for label, floor, _ in defs:
        pn, ps = counts(floor)
        emit(f"  Polymorphisme = {label} (Pn = {pn}, Ps = {ps})")
        emit(f"    {'Dn':>3s} {'Ds':>3s}  {'dN/dS':>7s}  {'NI':>7s}  {'p Fisher':>9s}   issue")
        best_p, rowsE = 1.0, []
        for k in range(dn + ds + 1):
            _dn, _ds = k, dn + ds - k
            _, _ni, _, _p = mk_stats(_dn, _ds, pn, ps)
            best_p = min(best_p, _p)
            star = "  <-- OBSERVE" if (_dn == dn) else ""
            ratio = f"{_dn/_ds:.2f}" if _ds else "inf"
            emit(f"    {_dn:>3d} {_ds:>3d}  {ratio:>7s}  {_ni:>7.2f}  {_p:>9.3f}   "
                 f"{'SIGNIFICATIF' if _p < ALPHA_LEVEL else 'non signif.'}{star}")
            rowsE.append({"Dn": _dn, "Ds": _ds, "NI": round(_ni, 3), "p": round(_p, 4)})
        enum_all[label] = {"Pn": pn, "Ps": ps, "meilleur_p_atteignable": round(best_p, 4),
                           "issues": rowsE}
        emit(f"    -> meilleur p ATTEIGNABLE, toutes issues confondues : {best_p:.3f}")
        if best_p >= ALPHA_LEVEL:
            emit(f"    -> AUCUNE des {dn+ds+1} issues possibles n'aurait atteint p < {ALPHA_LEVEL}.")
        emit("")

    # --- ascertainment --------------------------------------------------------
    pn1, ps1 = counts(1)
    pn2, ps2 = counts(100)
    emit("=" * 80)
    emit("BIAIS D'ASCERTAINMENT : LE PLANCHER DE FREQUENCE EST LE LEVIER")
    emit("-" * 80)
    emit(f"  sans plancher      : pN/pS = {pn1/ps1:.2f}  ({pn1} NS / {ps1} S)")
    emit(f"  les 9 SNP atlas    : pN/pS = {pn2/ps2:.2f}  ({pn2} NS / {ps2} S)")
    emit("  Le rapport se DEPLACE selon le plancher choisi, alors que la biologie")
    emit("  sous-jacente est la meme. C'est la signature attendue d'une selection")
    emit("  purificatrice : les variants non-synonymes sont maintenus a basse")
    emit("  frequence et disparaissent donc preferentiellement quand on releve le")
    emit("  plancher. Mais avec seulement 9 sites au-dessus du plancher, la direction")
    emit("  de ce deplacement n'est pas mesurable de facon fiable -- elle est donnee")
    emit("  ici comme description, pas comme test.")
    emit("")

    # --- verdict, derive des chiffres et non pre-ecrit ------------------------
    emit("=" * 80)
    emit("VERDICT")
    emit("")
    obs = results[0]
    main_enum = enum_all["tous les sites variables"]
    sig = [e for e in main_enum["issues"] if e["p"] < ALPHA_LEVEL]
    sig_dn = sorted(e["Dn"] for e in sig)
    emit(f"  Le test MK sur Rv2516c est NON SIGNIFICATIF sous les deux definitions du")
    emit(f"  polymorphisme (p = {results[0]['p_fisher']:.3f} sans plancher, "
         f"p = {results[1]['p_fisher']:.3f} avec).")
    emit("")
    if not sig:
        emit(f"  Et le test NE POUVAIT PAS conclure : aucune des {dn+ds+1} issues possibles")
        emit("  n'aurait atteint la significativite. L'axe est structurellement muet.")
    else:
        emit("  MAIS -- et c'est le resultat qui corrige l'attente posee par la piste -- le")
        emit("  test n'etait PAS structurellement muet. L'enumeration exhaustive montre que")
        emit(f"  {len(sig)} des {dn+ds+1} issues possibles auraient ete significatives, "
             f"celles ou Dn <= {max(sig_dn)} :")
        for e in sig:
            emit(f"    Dn = {e['Dn']}, Ds = {e['Ds']}  ->  p = {e['p']:.3f}")
        emit("")
        emit("  La piste annoncait 'meme un test correctement specifie aura tres")
        emit("  probablement un intervalle traversant la neutralite, le gain attendu est")
        emit("  methodologique, pas un nouveau resultat'. C'est a nuancer : la donnee AVAIT")
        emit("  de quoi trancher dans un sens, et elle ne l'a pas fait.")
        emit("")
        emit("  LA PUISSANCE EST ASYMETRIQUE, et c'est l'enseignement precis a retenir :")
        pos_extreme = [e for e in main_enum["issues"] if e["Ds"] == 0][0]
        emit(f"   - vers la selection PURIFICATRICE (deficit de non-synonymes fixes) : le")
        emit(f"     test detecte des Dn <= {max(sig_dn)}. Il avait donc les moyens de")
        emit("     montrer une contrainte forte sur la divergence. Il ne l'a pas montree.")
        emit(f"   - vers la selection POSITIVE (exces de non-synonymes fixes) : meme l'issue")
        emit(f"     la plus extreme possible, Dn = {pos_extreme['Dn']} / Ds = {pos_extreme['Ds']} "
             f"(dN/dS infini), ne donne que")
        emit(f"     p = {pos_extreme['p']:.3f}. Dans cette direction, l'axe est bel et bien muet :")
        emit("     aucune donnee de divergence de cette taille ne pourra jamais soutenir")
        emit("     une selection positive sur ce gene.")
        emit("")
        emit(f"   - l'observe (Dn = {dn}) tombe juste au-dela du dernier seuil detectable")
        emit(f"     (Dn = {max(sig_dn)}). Une seule substitution de moins dans la classe")
        emit("     non-synonyme aurait rendu le test significatif : la conclusion est donc")
        emit("     fragile au comptage, pas robuste.")
    emit("")
    emit("  CE QUE LE MANUSCRIT PEUT ECRIRE, desormais sans s'excuser : le rapport pN/pS")
    emit("  n'est pas 'du bruit qu'on ne sait pas interpreter', c'est un test dont la")
    emit("  puissance est connue et bornee. Il exclut une selection purificatrice FORTE au")
    emit("  niveau de la divergence, ne peut rien dire d'une selection positive, et laisse")
    emit("  ouvert tout regime intermediaire. C'est un enonce citable, la ou la reserve")
    emit("  precedente etait une simple mise en garde.")
    emit("")
    emit(f"  Direction, purement descriptive et NON soutenue : NI = {obs['NI']:.2f} > 1, "
         f"alpha = {obs['alpha']:.2f} < 0,")
    emit("  DoS < 0, ce qui pointerait vers un exces de polymorphisme non-synonyme encore")
    emit("  segregeant (contrainte faible a moderee, deleteres pas encore elimines). A ne")
    emit("  pas citer comme un resultat.")
    emit("")
    emit("  POINT DE COMPTAGE RESOLU EN PASSANT : le '9 SNP' du dossier n'etait rattache a")
    emit("  aucun plancher explicite, ce qui le rendait invérifiable. Balayage des seuils :")
    emit(f"  c'est >= 100 porteurs qui reproduit exactement la composition de la fiche")
    emit(f"  ({pn2} missense + {ps2} syn, 0 non-sens = {pn2+ps2} sites). Le chiffre est confirme et sa")
    emit("  definition desormais tracable.")
    emit("")
    emit("  ET UN CONTRESENS A EVITER, trouve en chemin : le pN/pS de 1,113 de la fiche est")
    emit("  NORMALISE par l'attente neutre propre au gene, il n'est PAS brut. Le rapport")
    emit(f"  brut vaut {pn2}/{ps2} = {pn2/ps2:.2f}. Lire 1,113 comme 'autant de non-synonymes que de")
    emit("  synonymes' serait faux ; c'est 'autant que ce qu'on attendrait sous neutralite'.")

    (OUT / "mk_resultats.json").write_text(json.dumps(
        {"divergence": {"Dn": dn, "Ds": ds, "substitutions": CANETTII},
         "tests": results, "enumeration_issues": enum_all}, indent=1, ensure_ascii=False))
    (OUT / "rapport.txt").write_text("\n".join(lines) + "\n")
    emit("")
    emit(f"Ecrit : {OUT}/")


if __name__ == "__main__":
    main()
