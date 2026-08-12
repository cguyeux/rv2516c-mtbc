#!/usr/bin/env python3
"""phase31_p9_3_audit_doublons.py -- P9.3 : inventaire des routines scientifiques dupliquées
entre projets MTBC, et détection des doublons DIVERGENTS.

POURQUOI. Pendant P4.3bis-e, `dinuc_shuffle` a été reforgé ici alors que
`mtbc/SpacerEgalVirus/analyses/phase1_test_nul_permutation.py` contenait DÉJÀ l'algorithme correct
d'Altschul & Erickson 1985 (tirage d'arbre couvrant), pendant que Rv2516c dupliquait deux fois une
marche eulérienne naïve. Le défaut n'a rien coûté à trouver une fois cherché : il a coûté de ne pas
chercher. Un doublon silencieux est une occasion d'hériter de la PLUS MAUVAISE des deux versions.

CE QUE FAIT CE SCRIPT, ET CE QU'IL NE FAIT PAS.
Il INVENTORIE : toutes les fonctions Python de ~/docs/codes/mtbc/, groupées par nom, avec le projet
d'origine, la taille, et une empreinte du corps NORMALISÉ (AST re-dumpé sans docstring ni constantes
de commentaire, donc insensible au renommage de variables locales : deux implémentations réellement
identiques ont la même empreinte même si l'une est commentée et l'autre non).
Il CLASSE : un nom présent dans ≥ 2 projets est un doublon inter-projets ; si les empreintes
diffèrent, c'est un doublon DIVERGENT, le seul cas qui demande une lecture humaine.
Il PRIORISE : un lexique de racines scientifiques (modèles nuls, coordonnées, distances, parseurs de
confiance structurale, statistiques) sépare les routines à enjeu des utilitaires d'affichage.

Il NE FACTORISE RIEN et NE MODIFIE AUCUN FICHIER. Le contre-argument posé par la piste tient :
mutualiser une routine entre dépôts crée une dépendance croisée qui casse l'autonomie de chaque
projet et peut rendre irreproductible un résultat déjà publié si la version partagée évolue. Et la
piste interdit explicitement de toucher aux scripts de SpacerEgalVirus depuis Rv2516c.

Sorties : résultats/p9_3_audit_doublons/{fonctions.tsv, doublons.tsv, divergents.md, resume.md}
Run: python analyses/phase31_p9_3_audit_doublons.py
"""
from __future__ import annotations

import ast
import hashlib
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
OUT = ROOT / "résultats" / "p9_3_audit_doublons"

SKIP_DIRS = {"__pycache__", ".git", "node_modules", ".venv", "venv", "site-packages",
             "zenodo_bundle"}
# `bdd/` contient ~146 000 répertoires de souches et aucun code d'analyse : un rglob qui le traverse
# ne rend jamais la main (constaté). Il est écarté à la racine, pas par SKIP_DIRS (qui filtre APRÈS
# le parcours et ne ferait donc rien gagner).
SKIP_TOP = {"bdd", "investigate_phylo"}

# racines de noms à enjeu scientifique : une divergence y change un RÉSULTAT, pas un affichage.
SCIENTIFIC = [
    ("modèle nul", r"shuffle|permut|null|randomi|bootstrap|resampl|monte_?carlo"),
    ("coordonnées/séquence", r"coord|spdi|translate|revcomp|reverse_comp|codon|gc_content|"
                             r"anchored|offset|to_h37rv|mtbc0"),
    ("distance/similarité", r"distance|hamming|jaccard|kimura|pdist|tm_?score|identity|pident|"
                            r"align|nj_|neighbor"),
    ("confiance structurale", r"plddt|_pae|pae_|iptm|ptm|confidence|boltz|af3|alphafold"),
    ("statistique", r"fisher|fdr|benjamini|bonferroni|pvalue|p_value|chi2|mann|welch|"
                    r"spearman|pearson|enrich|odds"),
    ("phylogénie/parcimonie", r"fitch|dollo|parsimon|tree|clade|mrca|newick|lineage|barcode"),
]


def project_of(p: Path) -> str:
    rel = p.relative_to(MTBC).parts
    return rel[0] if rel else "?"


def normalise(node: ast.AST) -> str:
    """Empreinte du corps : AST re-dumpé sans docstring, insensible aux commentaires."""
    body = list(node.body)  # type: ignore[attr-defined]
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
            and isinstance(body[0].value.value, str):
        body = body[1:]
    dumped = "\n".join(ast.dump(b, annotate_fields=False) for b in body)
    return hashlib.sha1(dumped.encode()).hexdigest()[:12]


def category(name: str) -> str:
    low = name.lower()
    for label, pat in SCIENTIFIC:
        if re.search(pat, low):
            return label
    return ""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    funcs: list[dict] = []
    n_files = n_bad = 0
    tops = [d for d in sorted(MTBC.iterdir())
            if d.is_dir() and d.name not in SKIP_TOP and not d.name.startswith(".")]
    for py in (p for top in tops for p in top.rglob("*.py")):
        if any(part in SKIP_DIRS for part in py.parts):
            continue
        try:
            tree = ast.parse(py.read_text(errors="replace"))
        except Exception:
            n_bad += 1
            continue
        n_files += 1
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                funcs.append({
                    "name": node.name,
                    "project": project_of(py),
                    "file": str(py.relative_to(MTBC)),
                    "line": node.lineno,
                    "n_lines": (node.end_lineno or node.lineno) - node.lineno + 1,
                    "hash": normalise(node),
                    "cat": category(node.name),
                })
    print(f"{n_files} fichiers .py analysés ({n_bad} illisibles), {len(funcs)} fonctions")

    with (OUT / "fonctions.tsv").open("w") as fh:
        fh.write("name\tproject\tfile\tline\tn_lines\thash\tcategorie\n")
        for f in sorted(funcs, key=lambda x: (x["name"], x["project"])):
            fh.write(f"{f['name']}\t{f['project']}\t{f['file']}\t{f['line']}\t"
                     f"{f['n_lines']}\t{f['hash']}\t{f['cat']}\n")

    by_name: dict[str, list[dict]] = defaultdict(list)
    for f in funcs:
        by_name[f["name"]].append(f)

    dups = []
    for name, occ in by_name.items():
        projects = {o["project"] for o in occ}
        if len(projects) < 2:
            continue
        hashes = {o["hash"] for o in occ}
        dups.append({
            "name": name, "n_proj": len(projects), "n_occ": len(occ),
            "projects": ",".join(sorted(projects)), "n_impl": len(hashes),
            "cat": occ[0]["cat"] or next((o["cat"] for o in occ if o["cat"]), ""),
            "occ": occ,
        })

    dups.sort(key=lambda d: (d["cat"] == "", -d["n_impl"], -d["n_proj"], d["name"]))
    with (OUT / "doublons.tsv").open("w") as fh:
        fh.write("name\tcategorie\tn_projets\tn_occurrences\tn_implementations\tprojets\tstatut\n")
        for d in dups:
            statut = "IDENTIQUE" if d["n_impl"] == 1 else "DIVERGENT"
            fh.write(f"{d['name']}\t{d['cat']}\t{d['n_proj']}\t{d['n_occ']}\t{d['n_impl']}\t"
                     f"{d['projects']}\t{statut}\n")

    sci_div = [d for d in dups if d["cat"] and d["n_impl"] > 1]
    lines = ["# P9.3 — doublons inter-projets DIVERGENTS à enjeu scientifique", "",
             f"{len(dups)} noms de fonction partagés par ≥2 projets ; "
             f"{sum(1 for d in dups if d['n_impl'] > 1)} ont des implémentations divergentes ; "
             f"**{len(sci_div)}** de ceux-là portent un nom à enjeu scientifique.", ""]
    for d in sci_div:
        lines.append(f"## `{d['name']}` — {d['cat']} — {d['n_impl']} implémentations "
                     f"dans {d['n_proj']} projets")
        for o in sorted(d["occ"], key=lambda x: (x["project"], x["file"])):
            lines.append(f"  - `{o['file']}:{o['line']}` ({o['n_lines']} lignes, "
                         f"empreinte {o['hash']})")
        lines.append("")
    (OUT / "divergents.md").write_text("\n".join(lines))

    print(f"{len(dups)} noms partagés entre projets, "
          f"{sum(1 for d in dups if d['n_impl'] > 1)} divergents, "
          f"{len(sci_div)} divergents à enjeu scientifique")
    for d in sci_div[:25]:
        print(f"  {d['name']:32s} {d['cat']:22s} {d['n_impl']} impl / {d['n_proj']} proj "
              f"[{d['projects'][:70]}]")


if __name__ == "__main__":
    main()
