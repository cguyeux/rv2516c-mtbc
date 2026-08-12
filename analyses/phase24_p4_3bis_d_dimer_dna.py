#!/usr/bin/env python3
"""phase24_p4_3bis_d_dimer_dna.py -- P4.3bis-d : le module HTH dimérise-t-il en présence de l'ADN ?

POURQUOI CE TEST, ET CE QU'IL AJOUTE AUX DEUX DÉJÀ FAITS.
P4.3bis-b (apo, deux copies de l'unité B seule) n'a pas pu trancher : l'ordre D2 > D1 > D3 était
inversé par rapport à l'attendu, pour deux raisons indépendantes -- un biais générique
d'empilement de Boltz sur petits domaines compacts, et l'indice que BldC lui-même dimérise
probablement sous l'effet de l'ADN plutôt que de façon constitutive (BldC apo 0,564/4,48 Å/16
contacts contre BldC+ADN, déjà testé en P2.4, 0,827/0,64 Å/113 contacts -- écart massif sur la
MÊME protéine et le MÊME pipeline). P2.4 a testé UNE copie de l'unité B + ADN (le test de liaison
simple). Aucun des deux jobs déjà faits ne teste directement le mécanisme suggéré : DEUX copies de
l'unité B, EN PRÉSENCE de l'ADN.

DESIGN, allégé pour tenir compte des contraintes de ressources actuelles (cf. cahier 2026-08-04 :
swap saturé, contention mémoire concurrente confirmée avoir tué P4.1/J2).
    E1  TEST   deux copies de l'unité B (88-147) + le même duplex de 22 pb que P2.4/P4.3bis-b
               (opérateur de BldC, 6AMA) -- teste si le dimère se forme/stabilise EN PRÉSENCE de l'ADN,
               par contraste avec l'apo (P4.3bis-b/D1, ipTM 0,767) et le monomère+ADN (P2.4/J1, ipTM 0,890).
Le contrôle positif idéal (2×BldC + ADN, dans le même contexte moléculaire complet) est différé : un
troisième front sur une machine déjà en tension mémoire ajouterait un risque d'échec silencieux
(OOM) plutôt qu'une réponse plus vite. À lancer en complément si E1 aboutit et que la question se
pose encore.

LECTURE. Comparer E1 à D1 (apo, P4.3bis-b) et à J1 (monomère+ADN, P2.4) sur ipTM, PAE inter-chaînes
minimum, et contacts d'interface -- le CONTRASTE informe, jamais la valeur absolue seule (discipline
du skill boltz). Un ipTM E1 nettement supérieur à D1 apo, du même ordre que J1, serait cohérent avec
une dimérisation couplée à l'ADN. Parseur : ../Rv1025/analyses/phase3_afmultimer_parse.py.

LIMITE annoncée avant le résultat : un système à 3 chaînes (2 protéines + ADN double brin) converge
plus difficilement sur CPU en un seul échantillon de diffusion qu'un système à 2 chaînes ; un
résultat faible ici pourrait signaler un défaut de convergence plutôt qu'un vrai négatif.

Run: python analyses/phase24_p4_3bis_d_dimer_dna.py        # écrit le YAML + le runner
"""
from __future__ import annotations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p4_3bis_d_dimer_dna"
BOLTZ = Path.home() / "venvs" / "boltz" / "bin" / "boltz"

FAA = ROOT / "data" / "Rv2516c.faa"
SEQ = "".join(l.strip() for l in FAA.read_text().splitlines()[1:])
UNIT_B = SEQ[87:147]  # 88-147 winged HTH AlpA/excisionase

# Même duplex que P2.4/P4.3bis-b (opérateur de BldC, 6AMA), pour rester comparable aux deux jobs déjà faits.
DNA_F = "ATTCGGGTAATTCGGGTAATTC"
COMP = {"A": "T", "T": "A", "G": "C", "C": "G"}
DNA_R = "".join(COMP[c] for c in reversed(DNA_F))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P4.3bis-d : dimère de l'unité B COUPLÉ à l'ADN ==")
    print(f"  ADN commun aux jobs déjà faits (P2.4, P4.3bis-b) : {DNA_F} / {DNA_R} ({len(DNA_F)} pb)")
    print(f"  unité B x2 : {len(UNIT_B)*2} résidus au total + {2*len(DNA_F)} pb d'ADN\n")

    y = OUT / "e1_unitB_dimer_dna.yaml"
    y.write_text(
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: [A, B]\n"          # deux copies de l'unité B = dimere candidat
        f"      sequence: {UNIT_B}\n"
        "  - dna:\n"
        "      id: C\n"
        f"      sequence: {DNA_F}\n"
        "  - dna:\n"
        "      id: D\n"
        f"      sequence: {DNA_R}\n"
    )
    print(f"  e1_unitB_dimer_dna   {len(UNIT_B)*2:>3} aa (2x unité B) + ADN 22 pb  -- TEST")

    runner = OUT / "run_dimer_dna.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "# P4.3bis-d -- dimere de l'unite B en presence de l'ADN. Pre-vol ressources integre :\n"
        "# P4.1/J2 a ete tue silencieusement DEUX FOIS (2026-08-04 et 2026-08-05) par contention\n"
        "# memoire concurrente avec des jobs boltz d'AUTRES projets sur cette machine partagee --\n"
        "# meme categorie de risque pour ce job a 3 chaines. Le pre-vol ci-dessous SAUTE le job\n"
        "# plutot que de le lancer dans l'espoir que ca passe.\n"
        "set -u\n"
        "trap '' HUP  # survit a la fermeture du shell qui l'a lance (P9.2, vecu le 2026-08-01)\n"
        'cd "$(dirname "$0")"\n'
        f"B={BOLTZ}\n"
        "STATUS=dimer_dna_status.log\n"
        "\n"
        "if [ ! -x \"$B\" ]; then\n"
        "  echo \"ERREUR : binaire Boltz introuvable ou non executable : $B\" | tee -a \"$STATUS\" >&2\n"
        "  exit 1\n"
        "fi\n"
        "\n"
        "check_resources() {  # cf. skill boltz -- verifie AVANT chaque job, pas seulement au demarrage\n"
        "  local other free_mib swap_free_mib\n"
        "  other=$(pgrep -f \"boltz predict\" 2>/dev/null | grep -vx \"$$\" || true)\n"
        "  if [ -n \"$other\" ]; then\n"
        "    echo \"  PRE-VOL ECHEC : un autre processus 'boltz predict' tourne deja (PID $other).\" >&2\n"
        "    return 1\n"
        "  fi\n"
        "  free_mib=$(free -m | awk '/^Mem:/{print $7}')\n"
        "  swap_free_mib=$(free -m | awk '/^Swap:/{print $4}')\n"
        "  if [ \"${free_mib:-0}\" -lt 4000 ] || [ \"${swap_free_mib:-0}\" -lt 500 ]; then\n"
        "    echo \"  PRE-VOL ECHEC : RAM libre ${free_mib:-?} Mio, swap libre ${swap_free_mib:-?} Mio (seuils 4000/500).\" >&2\n"
        "    return 1\n"
        "  fi\n"
        "  return 0\n"
        "}\n"
        "\n"
        "for y in e1_unitB_dimer_dna; do\n"
        "  if ! check_resources; then\n"
        "    echo \"  $y SAUTE (pre-vol ressources echoue, $(date +%H:%M:%S)) -- a relancer manuellement plus tard\" | tee -a \"$STATUS\"\n"
        "    continue\n"
        "  fi\n"
        "  echo \"=== $y $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
        "  \"$B\" predict $y.yaml --out_dir out_$y --accelerator cpu --devices 1 \\\n"
        "      --output_format mmcif --diffusion_samples 1 --num_workers 4 --use_msa_server \\\n"
        "      >> boltz_dimer_dna.log 2>&1 && echo \"  $y OK\" | tee -a \"$STATUS\" || echo \"  $y ECHEC\" | tee -a \"$STATUS\"\n"
        "done\n"
        "echo \"=== TERMINE $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
    )
    runner.chmod(0o755)
    print(f"\n  runner : {runner}")
    print("  pré-vol ressources intégré : SAUTE le job plutôt que de le lancer sous contention.")
    print("  lecture : comparer a D1 (apo, P4.3bis-b) et J1 (monomere+ADN, P2.4).")


if __name__ == "__main__":
    main()
