#!/usr/bin/env python3
"""phase11_p4_3bis_dimere.py -- P4.3bis, volet co-repliement : le module HTH dimérise-t-il ?

POURQUOI CE TEST-LÀ, ET PAS CELUI QUE P4.3bis ANNONÇAIT.
P4.3bis proposait de co-plier Rv2516c avec la machinerie d'élément mobile (Rv2512c transposase
IS1081, intégrases Rv2646 et Rv1586c). Deux raisons de ne pas commencer par là :
  1. **Mécanisme.** AlpA, chez *E. coli*, est un activateur TRANSCRIPTIONNEL du gène d'excision du
     prophage CP4-57 : il se fixe à l'ADN, il n'est pas un partenaire protéique de la recombinase.
     La famille ne prédit donc pas d'interaction protéine-protéine avec la transposase.
  2. **Coût.** Rv2516c (267 aa) + Rv2512c (415 aa) = 682 résidus. Sur CPU, avec un criblage déjà en
     cours, cela ne termine pas dans un horizon utile. Un test qu'on ne peut pas finir n'est pas un
     test.

Le co-repliement qui INFORME vraiment est ailleurs, et il est petit. Les winged HTH de cette famille
se fixent à l'ADN **en dimère**, ce qui est précisément la prémisse de la recherche de répétition
inversée menée en `phase10_p4_3bis_operateur.py`. Les deux volets se testent donc l'un l'autre :

    si l'unité B dimérise ET qu'un opérateur palindromique ressort -> le modèle est cohérent
    si elle ne dimérise pas -> le modèle « opérateur palindromique » perd sa base structurale,
                               et une éventuelle répétition inversée trouvée devient anecdotique

DESIGN, entièrement fait de contrôles (discipline du skill `boltz` et héritage Rv1025).

    D1  TEST              unité B (88-147) x2, le module wHTH AlpA/excisionase
    D2  contrôle NÉGATIF de spécificité : unité A (1-87) x2, domaine ferredoxin-like/DUF8830
        qu'aucune analyse ne rattache à une interface de dimérisation
        -> mesure ce qu'un domaine NON attendu obtient dans le même protocole
    D3  contrôle POSITIF  BldC de *S. venezuelae* (6AMA) x2 — régulateur qui se fixe à l'ADN en
        dimère, actinobactérie, et l'un des meilleurs hits Foldseek de l'unité B (TM 0,748)
        -> calibre le plafond atteignable ; sans lui un ipTM médiocre est ininterprétable

Tailles : 120, 174 et 136 résidus. Réalisable sur CPU même chargé, contrairement au co-repliement
avec la transposase.

LECTURE DES SORTIES. Le discriminant est le **PAE inter-chaînes minimum**, pas l'ipTM absolu ; et
c'est le CONTRASTE D1 vs D2, calibré par D3, qui informe. Parseur bi-schéma AF3/Boltz :
`../Rv1025/analyses/phase3_afmultimer_parse.py`. Méthode à ÉTIQUETER comme du criblage Boltz-2, à ne
jamais comparer à un chiffre AF3 sans le dire.

LIMITE. Un homodimère prédit n'est pas un homodimère physiologique, et Boltz place volontiers deux
copies d'un petit domaine l'une contre l'autre. C'est pour cela que D2 existe : si l'unité A obtient
autant que l'unité B, le test n'a rien montré.

Run: python analyses/phase11_p4_3bis_dimere.py        # écrit les YAML + le runner
"""
from __future__ import annotations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p4_3bis_dimere"
BOLTZ = Path.home() / "venvs" / "boltz" / "bin" / "boltz"

FAA = ROOT / "data" / "Rv2516c.faa"
SEQ = "".join(l.strip() for l in FAA.read_text().splitlines()[1:])
UNIT_A = SEQ[0:87]          # 1-87   ferredoxin-like / DUF8830
UNIT_B = SEQ[87:147]        # 88-147 winged HTH AlpA/excisionase
# BldC de S. venezuelae, 6AMA chaîne A, tag GSHM de purification retiré (même construction qu'en P2.4)
BLDC = "MTARTPDAEPLLTPAEVATMFRVDPKTVTRWAKAGKLTSIRTLGGHRRYREAEVRALLAGIPQQRSEA"

JOBS = {
    "d1_unitB_dimere": ("unité B (88-147) x2 -- TEST", UNIT_B),
    "d2_unitA_dimere": ("unité A (1-87) x2 -- CONTRÔLE NÉGATIF DE SPÉCIFICITÉ", UNIT_A),
    "d3_bldc_dimere":  ("BldC (6AMA) x2 -- CONTRÔLE POSITIF", BLDC),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P4.3bis, volet co-repliement : dimérisation du module HTH ==")
    print("  Prémisse testée : un wHTH de cette famille se fixe à l'ADN en DIMÈRE.")
    print("  Ce volet et la recherche d'opérateur (phase10) se valident mutuellement.\n")
    for name, (desc, seq) in JOBS.items():
        (OUT / f"{name}.yaml").write_text(
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: [A, B]\n"          # deux chaînes identiques = homodimère
            f"      sequence: {seq}\n"
        )
        print(f"  {name:18} {2*len(seq):>3} résidus au total  {desc}")

    runner = OUT / "run_dimeres.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "# P4.3bis -- criblage d'homodimérisation. À lancer APRÈS les contrôles de P2.4 :\n"
        "# le CPU est partagé, et P2.4 décide de la branche ADN dont ce volet dépend.\n"
        "set -u\n"
        "trap '' HUP  # survit a la fermeture du shell qui l'a lance (P9.2 : vecu le 2026-08-01,\n"
        "             # un runner attache au shell est mort 4 min apres son lancement, en silence)\n"
        'cd "$(dirname "$0")"\n'
        f"B={BOLTZ}\n"
        "STATUS=dimeres_status.log\n"
        "\n"
        "if [ ! -x \"$B\" ]; then  # echouer BRUYAMMENT plutot que job par job (P9.2 : vecu le\n"
        "  # 2026-08-01, le venv Boltz avait disparu et 3 des 4 jobs de P2.4 ont echoue en silence,\n"
        "  # laissant un seul ipTM sans aucun de ses controles -- un resultat sans ses controles\n"
        "  # ressemble exactement a un resultat.)\n"
        "  echo \"ERREUR : binaire Boltz introuvable ou non executable : $B\" | tee -a \"$STATUS\" >&2\n"
        "  echo \"  Le venv a-t-il disparu ? cf. skill boltz, section Installation.\" >&2\n"
        "  exit 1\n"
        "fi\n"
        "\n"
        "for y in d3_bldc_dimere d1_unitB_dimere d2_unitA_dimere; do\n"
        "  echo \"=== $y $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
        "  \"$B\" predict $y.yaml --out_dir out_$y --accelerator cpu --devices 1 \\\n"
        "      --output_format mmcif --diffusion_samples 1 --num_workers 4 --use_msa_server \\\n"
        "      >> boltz_dimeres.log 2>&1 && echo \"  $y OK\" | tee -a \"$STATUS\" || echo \"  $y ECHEC\" | tee -a \"$STATUS\"\n"
        "done\n"
        "echo \"=== TERMINE $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
    )
    runner.chmod(0o755)
    print(f"\n  runner : {runner}")
    print("  ordre : contrôle POSITIF d'abord, c'est lui qui rend les deux autres lisibles.")
    print("  lecture : ../Rv1025/analyses/phase3_afmultimer_parse.py (bi-schéma AF3/Boltz)")


if __name__ == "__main__":
    main()
