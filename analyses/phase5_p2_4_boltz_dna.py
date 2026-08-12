#!/usr/bin/env python3
"""phase5_p2_4_boltz_dna.py -- P2.4 : le module HTH de Rv2516c peut-il loger un ADN double brin ?

Contexte. HHpred assigne l'unité B (88-147) à la famille AlpA / excisionase / facteur de
directionnalité de recombinaison, avec E = 1,3e-10 : ce sont des winged HTH de liaison à
l'ADN. P2.4 teste la conséquence structurale de cette assignation, en local avec Boltz-2
(AlphaFold Server n'étant pas automatisable, cf. skill `boltz`).

DESIGN, et il est presque entièrement fait de contrôles — c'est le CONTRASTE qui informe,
jamais la valeur absolue d'un ipTM (garde-fou du skill, et discipline anti-circularité
héritée de Rv1025).

  J1  test          unité B (88-147, 60 aa)              + ADN
  J2  contrôle NÉGATIF MUTANT : unité B avec R113A/H115A/R118A, les résidus basiques de
      l'hélice de reconnaissance (114-121 VHQLRSTA)      + ADN
      -> si le signal de J1 est réel et porté par le HTH, il doit s'effondrer ici.
  J3  contrôle de SPÉCIFICITÉ : unité A (1-87), le domaine ferredoxin-like, qu'aucune
      analyse ne rattache à l'ADN                        + ADN
      -> mesure ce qu'un domaine NON attendu obtient sur le même ADN.
  J4  contrôle POSITIF : BldC de *Streptomyces venezuelae* (6AMA), régulateur de type
      MerR co-cristallisé AVEC cet ADN                   + ADN
      -> calibre le plafond atteignable. Sans lui, un ipTM médiocre est ininterprétable.
      BldC est aussi un des meilleurs hits Foldseek de l'unité B (TM 0,748) et vient d'une
      actinobactérie, donc le comparateur taxonomiquement le plus proche disponible.

**Le même ADN pour les quatre jobs** : duplex de 22 pb tiré de l'opérateur de BldC (6AMA,
répétition ATTCGGGTA). La protéine est ainsi la SEULE variable.

LIMITE À ÉNONCER. Aucun site de liaison n'est connu pour Rv2516c. Un résultat positif dirait
donc « cette surface peut accommoder un ADN B », ce qui pour un wHTH est presque acquis et
donc peu informatif en soi. **Toute l'information est dans J2 et J3** : le mutant perd-il le
signal, et un domaine non-ADN fait-il aussi bien ?

MÉTHODE À ÉTIQUETER : Boltz-2, c'est-à-dire du CRIBLAGE. À ne pas comparer à un chiffre AF3
sans le dire, et à ne pas présenter comme méthode de référence d'un manuscrit.

Run: python analyses/phase5_p2_4_boltz_dna.py        # écrit les YAML + le runner
"""
from __future__ import annotations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p2_4_boltz"
BOLTZ = Path.home() / "venvs" / "boltz" / "bin" / "boltz"

FAA = ROOT / "data" / "Rv2516c.faa"
SEQ = "".join(l.strip() for l in FAA.read_text().splitlines()[1:])

UNIT_A = SEQ[0:87]          # 1-87   ferredoxin-like
UNIT_B = SEQ[87:147]        # 88-147 winged HTH AlpA/excisionase

# mutant du contrôle négatif : résidus basiques de l'hélice de reconnaissance -> Ala.
# Indices protéine (1-based) 113 R, 115 H, 118 R ; dans l'unité B ils tombent en 113-87=26, 28, 31.
MUT_POS = [113, 115, 118]
_b = list(UNIT_B)
for p in MUT_POS:
    i = p - 88
    assert _b[i] in "RH", f"résidu {p} attendu R/H, trouvé {_b[i]}"
    _b[i] = "A"
UNIT_B_MUT = "".join(_b)

# BldC de S. venezuelae, 6AMA chaîne A, tag GSHM de purification retiré
BLDC = ("MTARTPDAEPLLTPAEVATMFRVDPKTVTRWAKAGKLTSIRTLGGHRRYREAEVRALLAGIPQQRSEA")

# duplex de 22 pb tiré de l'opérateur de BldC (6AMA : répétition ATTCGGGTA)
DNA_F = "ATTCGGGTAATTCGGGTAATTC"
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
DNA_R = "".join(COMP[c] for c in reversed(DNA_F))

JOBS = {
    "j1_unitB_dna":      ("unité B (88-147), le module HTH -- TEST", UNIT_B),
    "j2_unitB_mut_dna":  ("unité B R113A/H115A/R118A -- CONTRÔLE NÉGATIF MUTANT", UNIT_B_MUT),
    "j3_unitA_dna":      ("unité A (1-87), ferredoxin-like -- CONTRÔLE DE SPÉCIFICITÉ", UNIT_A),
    "j4_bldc_dna":       ("BldC de S. venezuelae (6AMA) -- CONTRÔLE POSITIF", BLDC),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P2.4 : liaison à l'ADN, criblage Boltz-2 ==")
    print(f"  ADN commun aux 4 jobs : {DNA_F} / {DNA_R} ({len(DNA_F)} pb)")
    print(f"  mutant : {'/'.join(f'{p}A' for p in MUT_POS)}\n")
    for name, (desc, prot) in JOBS.items():
        y = OUT / f"{name}.yaml"
        y.write_text(
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {prot}\n"
            "  - dna:\n"
            "      id: B\n"
            f"      sequence: {DNA_F}\n"
            "  - dna:\n"
            "      id: C\n"
            f"      sequence: {DNA_R}\n"
        )
        print(f"  {name:20} {len(prot):>3} aa  {desc}")

    runner = OUT / "run_boltz.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "# P2.4 -- criblage Boltz-2. Lancer en TÂCHE DE FOND : plusieurs dizaines de minutes par job.\n"
        "# Les séquences sont publiques (génome H37Rv, PDB 6AMA), donc --use_msa_server est acceptable.\n"
        "set -u\n"
        "trap '' HUP  # survit a la fermeture du shell qui l'a lance (P9.2 : vecu le 2026-08-01,\n"
        "             # un runner attache au shell est mort 4 min apres son lancement, en silence)\n"
        f"cd {OUT}\n"
        "STATUS=boltz_status.log\n"
        "\n"
        f"if [ ! -x {BOLTZ} ]; then  # echouer BRUYAMMENT plutot que job par job (P9.2 : vecu le\n"
        "  # 2026-08-01, le venv Boltz avait disparu et 3 des 4 jobs ont echoue en silence, laissant\n"
        "  # un seul ipTM sans aucun de ses controles -- un resultat sans ses controles ressemble\n"
        "  # exactement a un resultat.)\n"
        f"  echo \"ERREUR : binaire Boltz introuvable ou non executable : {BOLTZ}\" | tee -a \"$STATUS\" >&2\n"
        "  echo \"  Le venv a-t-il disparu ? cf. skill boltz, section Installation.\" >&2\n"
        "  exit 1\n"
        "fi\n"
        "\n"
        "for y in j1_unitB_dna j2_unitB_mut_dna j3_unitA_dna j4_bldc_dna; do\n"
        "  echo \"=== $y $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
        f"  {BOLTZ} predict $y.yaml --out_dir out_$y --accelerator cpu --devices 1 \\\n"
        "      --output_format mmcif --diffusion_samples 1 --num_workers 4 --use_msa_server \\\n"
        "      >> boltz.log 2>&1 && echo \"  $y OK\" | tee -a \"$STATUS\" || echo \"  $y ECHEC\" | tee -a \"$STATUS\"\n"
        "done\n"
        "echo \"=== TERMINE $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
    )
    runner.chmod(0o755)
    print(f"\n  runner : {runner}")
    print("  lecture des sorties : ../Rv1025/analyses/phase3_afmultimer_parse.py (bi-schéma AF3/Boltz)")
    print("  RAPPEL : le discriminant est le PAE inter-chaînes minimum, pas l'ipTM absolu ;")
    print("           et c'est le CONTRASTE J1 vs J2/J3, calibré par J4, qui informe.")


if __name__ == "__main__":
    main()
