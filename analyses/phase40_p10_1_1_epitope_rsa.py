#!/usr/bin/env python3
"""P10.1.1 -- positionner l'épitope bovin 227-246 (seul fait expérimental sur l'unité C) dans le
modèle 3D : boucle exposée ou brin enfoui ?

Aucun nouveau calcul de structure : réutilise résultats/p2_3_charge_hth/sasa_par_residu.json
(P2.3, Shrake-Rupley déjà calculé sur tout le monomère AlphaFold, mêmes seuils RSA>=0.25 exposé).

Lit  : résultats/p2_3_charge_hth/sasa_par_residu.json
Écrit: résultats/p10_1_1_epitope_rsa/resume.md, epitope_rsa.json
"""
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "résultats" / "p2_3_charge_hth" / "sasa_par_residu.json"
OUT = ROOT / "résultats" / "p10_1_1_epitope_rsa"
OUT.mkdir(parents=True, exist_ok=True)

EPITOPE = (227, 246)   # RPGEGLNMVLIAAWGHPLPN, Farrell et al. 2016, PMID 28348866
UNIT_C = (178, 267)


def main() -> None:
    rows = json.loads(SRC.read_text())
    by_pos = {r["pos"]: r for r in rows}

    epi = [by_pos[p] for p in range(EPITOPE[0], EPITOPE[1] + 1) if p in by_pos]
    unit_c = [by_pos[p] for p in range(UNIT_C[0], UNIT_C[1] + 1) if p in by_pos]
    unit_c_rest = [r for r in unit_c if not (EPITOPE[0] <= r["pos"] <= EPITOPE[1])]

    def stats(group, label):
        rsas = [r["rsa"] for r in group]
        n_exposed = sum(1 for r in group if r["exposed"])
        return {
            "label": label, "n": len(group),
            "rsa_moyen": round(statistics.mean(rsas), 3),
            "rsa_median": round(statistics.median(rsas), 3),
            "n_exposed": n_exposed,
            "frac_exposed": round(n_exposed / len(group), 3) if group else None,
        }

    s_epi = stats(epi, "épitope 227-246")
    s_unitc = stats(unit_c, "unité C entière (178-267)")
    s_unitc_rest = stats(unit_c_rest, "unité C hors épitope")

    seq = "".join(r["aa"] for r in epi)
    detail = [
        f"    {r['aa']}{r['pos']} ({r['classe']:10s}) RSA={r['rsa']:.2f} "
        f"{'EXPOSÉ' if r['exposed'] else 'enfoui'}"
        for r in epi
    ]

    lines = [
        "# P10.1.1 -- position structurale de l'épitope bovin 227-246",
        "",
        f"Séquence (H37Rv 227-246) : `{seq}` (attendu `RPGEGLNMVLIAAWGHPLPN`, Farrell et al. "
        "2016, PMID 28348866, BoLA-DRB3, IFNγ)",
        "",
        "| groupe | n résidus | RSA moyen | RSA médian | % exposés (RSA≥0,25) |",
        "|---|---|---|---|---|",
        f"| {s_epi['label']} | {s_epi['n']} | {s_epi['rsa_moyen']} | {s_epi['rsa_median']} "
        f"| {100*s_epi['frac_exposed']:.0f} % |",
        f"| {s_unitc_rest['label']} | {s_unitc_rest['n']} | {s_unitc_rest['rsa_moyen']} "
        f"| {s_unitc_rest['rsa_median']} | {100*s_unitc_rest['frac_exposed']:.0f} % |",
        f"| {s_unitc['label']} | {s_unitc['n']} | {s_unitc['rsa_moyen']} | {s_unitc['rsa_median']} "
        f"| {100*s_unitc['frac_exposed']:.0f} % |",
        "",
        "## Détail par résidu",
        "",
    ] + detail

    verdict_expose = s_epi["frac_exposed"] >= s_unitc_rest["frac_exposed"]
    lines += [
        "",
        "## Verdict",
        "",
        (
            f"L'épitope (fraction exposée {100*s_epi['frac_exposed']:.0f} %, RSA moyen "
            f"{s_epi['rsa_moyen']}) est {'PLUS' if verdict_expose else 'MOINS'} exposé que le "
            f"reste de l'unité C (fraction exposée {100*s_unitc_rest['frac_exposed']:.0f} %, RSA "
            f"moyen {s_unitc_rest['rsa_moyen']}). "
            + (
                "Cohérent avec une boucle de surface, position candidate d'interaction plutôt "
                "que cœur structural enfoui."
                if verdict_expose else
                "L'épitope n'est pas préférentiellement exposé par rapport au reste du domaine ; "
                "sa présentation par le CMH bovin n'implique pas nécessairement une localisation "
                "de surface distinctive dans la structure prédite."
            )
        ),
        "",
        "Réserve (déjà posée par la piste) : la géométrie seule ne dit rien de la fonction, "
        "seulement de l'accessibilité ; un épitope T CD8/CD4 bovin n'a pas besoin d'être sur une "
        "surface fonctionnellement pertinente pour Rv2516c elle-même, seulement présentable par le "
        "CMH de l'hôte.",
    ]

    (OUT / "resume.md").write_text("\n".join(lines) + "\n")
    (OUT / "epitope_rsa.json").write_text(json.dumps(
        {"epitope": s_epi, "unite_C_hors_epitope": s_unitc_rest, "unite_C_entiere": s_unitc,
         "detail_par_residu": epi},
        indent=1, ensure_ascii=False))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
