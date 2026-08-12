#!/usr/bin/env python3
"""phase31b_p9_3_test_differentiel.py -- P9.3, second temps : les doublons divergents
donnent-ils les MÊMES RÉSULTATS ?

phase31 dit quels noms de fonction existent en plusieurs versions ; il ne dit pas si la divergence
est cosmétique ou substantielle. Une empreinte d'AST différente peut ne refléter qu'un style
(`for` vs compréhension), ou cacher une erreur de calcul. Seul un test différentiel tranche : on
extrait chaque implémentation, on l'exécute EN ISOLATION sur les mêmes entrées, et on compare.

Ce que le test peut conclure, et ce qu'il ne peut pas :
  - sorties identiques sur tous les cas -> divergence COSMÉTIQUE (les deux versions sont
    interchangeables ; le doublon reste un coût de maintenance, pas un risque de résultat) ;
  - sorties différentes -> il faut lire, et l'une des deux au moins est fausse OU les deux
    répondent à des questions différentes sous le même nom (le cas le plus dangereux) ;
  - non exécutable en isolation (dépend de globales du module) -> non conclu, signalé comme tel,
    jamais silencieusement compté comme « identique ».

Les jeux d'entrées incluent délibérément les cas limites où les implémentations naïves divergent :
ensembles vides (Jaccard 0/0), séquences de longueurs inégales (Hamming), lettres ambiguës et
minuscules (revcomp), codons incomplets et codons stop (translate).

Aucun fichier hors de ce projet n'est modifié : cf. la consigne de P9.3.

Entrée  : résultats/p9_3_audit_doublons/fonctions.tsv (produit par phase31)
Sortie  : résultats/p9_3_audit_doublons/test_differentiel.md
Run: python analyses/phase31b_p9_3_test_differentiel.py
"""
from __future__ import annotations

import ast
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
AUDIT = ROOT / "résultats" / "p9_3_audit_doublons"

# arbres tiers embarqués : ce n'est pas notre code, un doublon n'y a aucun sens
VENDORED = ("SPAdes-", "googletest", "site-packages", "/build/")

CASES: dict[str, list[tuple]] = {
    "revcomp": [("ATGC",), ("atgc",), ("",), ("ATGCN",), ("ACGTRYKM",), ("AAAATTTT",)],
    "translate": [("ATGGCCTAA",), ("ATGGCC",), ("ATG",), ("",), ("ATGGCCTA",), ("TTTTTAATG",)],
    "jaccard": [({1, 2, 3}, {2, 3, 4}), (set(), set()), ({1}, set()), ({1, 2}, {1, 2})],
    "jaccard_dist": [({1, 2, 3}, {2, 3, 4}), (set(), set()), ({1}, set()), ({1, 2}, {1, 2})],
    "hamming": [("ATGC", "ATGG"), ("ATGC", "ATGC"), ("ATGC", "AT"), ("", "")],
}


def load_impls(name: str) -> list[tuple[str, int, object | str]]:
    """(fichier, ligne, fonction ou message d'échec) pour chaque implémentation du nom."""
    out = []
    for line in (AUDIT / "fonctions.tsv").read_text().splitlines()[1:]:
        p = line.split("\t")
        if p[0] != name or any(v in p[2] for v in VENDORED):
            continue
        src_lines = (MTBC / p[2]).read_text(errors="replace").splitlines()
        start, n = int(p[3]) - 1, int(p[4])
        block = "\n".join(src_lines[start:start + n])
        # dédente si la fonction est imbriquée / méthode
        indent = len(block) - len(block.lstrip())
        if indent:
            block = "\n".join(l[indent:] if len(l) > indent else l.lstrip()
                              for l in block.splitlines())
        ns: dict = {}
        try:
            exec("from collections import Counter, defaultdict\nimport re, math, itertools, "
                 "statistics\n" + block, ns)
            fn = ns[name]
        except Exception as e:
            fn = f"non exécutable en isolation ({type(e).__name__}: {e})"
        out.append((p[2], int(p[3]), fn))
    return out


def main() -> None:
    report = ["# P9.3 — test différentiel des doublons à enjeu", "",
              "Chaque implémentation est extraite et exécutée en isolation sur les mêmes entrées, "
              "cas limites inclus. `=` signifie sorties identiques à l'implémentation de "
              "référence (la première exécutable) ; `≠` une divergence de RÉSULTAT.", ""]
    for name, cases in CASES.items():
        impls = load_impls(name)
        if len(impls) < 2:
            continue
        report.append(f"## `{name}` — {len(impls)} implémentations")
        results: dict[str, list] = {}
        for path, line, fn in impls:
            key = f"{path}:{line}"
            if isinstance(fn, str):
                report.append(f"  - `{key}` — {fn}")
                continue
            outs = []
            for c in cases:
                try:
                    outs.append(repr(fn(*c)))
                except Exception as e:
                    outs.append(f"<{type(e).__name__}>")
            results[key] = outs
        if not results:
            report.append("")
            continue
        ref_key = next(iter(results))
        ref = results[ref_key]
        report.append(f"  - référence : `{ref_key}`")
        report.append("")
        report.append("    | entrée | " + " | ".join(f"impl{i+1}" for i in range(len(results))) + " |")
        report.append("    |---|" + "---|" * len(results))
        for i, c in enumerate(cases):
            row = " | ".join(v[i] for v in results.values())
            report.append(f"    | `{c}` | {row} |")
        report.append("")
        # Un `<NameError>` partout signifie que la fonction dépend d'une table globale du module :
        # elle n'a PAS été testée. La confondre avec une divergence de résultat gonflerait
        # artificiellement le nombre de doublons suspects — c'est le piège que ce script doit éviter.
        def testee(v: list[str]) -> bool:
            return any(not x.startswith("<") for x in v)

        n_test = sum(1 for v in results.values() if testee(v))
        n_div = sum(1 for v in results.values() if testee(v) and v != ref)
        for k, v in results.items():
            mark = "?" if not testee(v) else ("=" if v == ref else "≠")
            label = {"?": "non testable en isolation (dépend d'une globale du module)",
                     "=": "résultats identiques", "≠": "RÉSULTATS DIVERGENTS"}[mark]
            report.append(f"  - {mark} `{k}` — {label}")
        report.append("")
        print(f"{name:14s} {n_test} réellement testées / {len(impls)} — {n_div} divergentes")
        for k, v in results.items():
            if testee(v) and v != ref:
                print(f"    DIVERGENT : {k}")
                for i, c in enumerate(cases):
                    if v[i] != ref[i]:
                        print(f"        {c} -> {v[i]}   (réf : {ref[i]})")
    (AUDIT / "test_differentiel.md").write_text("\n".join(report))
    print(f"\nrapport : {AUDIT / 'test_differentiel.md'}")


if __name__ == "__main__":
    main()
