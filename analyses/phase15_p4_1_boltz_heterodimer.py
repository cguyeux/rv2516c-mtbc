#!/usr/bin/env python3
"""phase15_p4_1_boltz_heterodimer.py -- P4.1 : Rv2516c forme-t-il un hétérodimère avec son partenaire
d'opéron Rv2517c ?

PRÉDICTION ÉNONCÉE AVANT LE CALCUL, à ne pas oublier une fois le résultat en main.
Rv2517c est NON essentiel (DeJesus NE) et NON vulnérable (VI CRISPRi +0,54, IC95 -1,88..+4,58,
traversant zéro), alors que Rv2516c est ESD et VI -9,91. Un hétérodimère OBLIGATOIRE est donc peu
probable : si Rv2517c était le partenaire structural indispensable de Rv2516c, sa propre extinction
devrait coûter cher, ce qui n'est pas mesuré. **Un résultat négatif ici serait cohérent, pas
décevant.** Mais un ipTM faible SANS contrôle positif ne prouve rien : un ipTM max non reproductible
n'est interprétable qu'une fois un contrôle positif ayant retrouvé une interface structuralement
connue.

DESIGN, deux jobs seulement (la piste ne demande que le contrôle positif, pas un panel complet ;
le CPU est déjà saturé par P2.4 + P4.3bis-b en file, ne pas ouvrir un troisième front sans besoin
justifié -- cf. garde-fous du skill `boltz` sur le coût réel d'un run) :

    J1  TEST     Rv2516c (267 aa) + Rv2517c (83 aa), l'opéron dark lui-même.
    J2  CONTRÔLE POSITIF   Rv2514c (153 aa) + Rv2515c (415 aa), le couple toxine-antitoxine de la
        MÊME cassette (Rv2513-Rv2518c), validé EXPÉRIMENTALEMENT (Tandon et al. 2019, J Biol Chem,
        doi 10.1074/jbc.ra118.006814) : surexpression de la toxine Rv2514c -> bactériostase,
        co-expression de l'antitoxine Rv2515c -> restauration de la croissance. Une antitoxine
        neutralise sa toxine cognate par contact DIRECT, donc ce couple a un mécanisme structural
        attendu, contrairement à Rv2516c-Rv2517c dont le lien n'est pour l'instant QUE transcriptionnel
        (opéron, chevauchement de 4 pb, co-expression r=+0,409 en P7.4). Choisi plutôt qu'un complexe
        générique hors MTBC : même organisme, même cassette génomique, donc le contrôle calibre le
        plafond ATTEIGNABLE PAR CE PIPELINE DANS CE CONTEXTE PRÉCIS, pas seulement "AF/Boltz peut
        détecter une interface quelconque".

LECTURE. Le discriminant est le PAE inter-chaînes minimum, pas l'ipTM absolu (cf. skill boltz,
garde-fous d'interprétation). Si J2 (positif) montre un signal net et J1 (test) n'en montre aucun,
le négatif de J1 devient interprétable : hétérodimère obligatoire écarté, cohérent avec la prédiction.
Si J2 lui-même ne montre rien, le pipeline est simplement insensible à ce type d'interaction dans ce
contexte (petites protéines, CPU, MSA serveur) et J1 reste ININTERPRÉTABLE quel que soit son chiffre.

MÉTHODE À ÉTIQUETER : Boltz-2 = criblage local, jamais comparé à un chiffre AF3 sans le dire (cf.
skill boltz). Les deux entrées sont indépendantes de tout MSA existant (aucun run AF3 antérieur sur
CES paires), donc `--use_msa_server` (ColabFold public) est utilisé, séquences H37Rv publiques.

FILE D'ATTENTE CPU (déjà 2 jobs Boltz en cours/en file sur cette machine) : ce runner NE DOIT PAS être
lancé avant que P2.4 (contrôles J2/J3/J4) ET P4.3bis-b (dimérisation) aient terminé. P4.1 n'a AUCUNE
dépendance LOGIQUE envers ces deux pistes (question orthogonale : interaction stable Rv2516c-Rv2517c,
indépendante de la liaison à l'ADN et de l'auto-dimérisation du module HTH), seulement une dépendance
PHYSIQUE de ressource -- lancer trois jobs Boltz de front sur 16 coeurs déjà saturés (charge mesurée
14-18 avec un seul job actif) allonge tout le monde sans rien décider plus tôt.

Run: python analyses/phase15_p4_1_boltz_heterodimer.py     # écrit les YAML + le runner, NE LANCE RIEN
"""
from __future__ import annotations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats" / "p4_1_boltz_heterodimer"
BOLTZ = Path.home() / "venvs" / "boltz" / "bin" / "boltz"

RV2516C = ("".join(l.strip() for l in (ROOT / "data" / "Rv2516c.faa").read_text().splitlines()[1:]))
RV2517C = "MNSAIIKIAKWAQSQQWTVEDDASGYTRFYNPQGVYIARFPATPSNEYRRMRDLLGALKKAGLTWPPPSKKERRAQHRKEGAQ"
RV2514C = ("MLYSFDTSAILNGRRDLFRPAVFRSLWGRVEDAISAGQIRSVDEVQRELARRDDDAKRWADGQTGLFCPLDEQIQQAARHILR"
           "LHPNMVRQGGRRSAADPFVIALAMVNNATVVTQETASGNIEKPRIPDVCDALGVPWLTLMGYIEAQGWTF")
RV2515C = ("MGIGHPMWVGWCIIIAMRSIPASVESSVLRWARESCGLTEVAAARKLGLPDDRVAAWEVGEVVPTIAQLRKAAEVYKRSLAVFF"
           "LSEPPEGFDTLRDFRRLDGAASGQWTPGLHEEFRRAHTQRDFALELADAEDREIPGAWRLPLSGDEADADIAARIRKALIEVSP"
           "LPIPVASVDPYEHLNAWVSAIETSGVLVLATRGGKVAIDEMRGMCLYFDELPVIVLNGSDHPRPRLFSLLHEFVHVVLHTEGLC"
           "DVIADAHPSTQDRSLEARCNAIAAAVLMPADVVRARPEVIVRSETPSSWDYESLRPVAAHFGVSAEAFLRRLSTLGIVPVEVYR"
           "QRRAEFIAAHEDEAERARSAGGGNWYRNTVRDLGKGYVRAVTDAHRRRVIDSNTAAIYLDAKVSQIPKLAESAELRSVV")

assert len(RV2516C) == 267 and len(RV2517C) == 83 and len(RV2514C) == 153 and len(RV2515C) == 415

JOBS = {
    "j1_rv2516c_rv2517c": ("Rv2516c + Rv2517c -- TEST, l'opéron dark", RV2516C, RV2517C),
    "j2_rv2514c_rv2515c": ("Rv2514c + Rv2515c -- CONTRÔLE POSITIF, TA validé (Tandon 2019)",
                            RV2514C, RV2515C),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P4.1 : Rv2516c forme-t-il un hétérodimère avec Rv2517c ? ==")
    print("  Prédiction : NÉGATIF probable (Rv2517c non essentiel, non vulnérable).")
    print("  Un négatif n'est interprétable QU'AVEC le contrôle positif (leçon Rv1025).\n")
    for name, (desc, seq_a, seq_b) in JOBS.items():
        (OUT / f"{name}.yaml").write_text(
            "version: 1\n"
            "sequences:\n"
            "  - protein:\n"
            "      id: A\n"
            f"      sequence: {seq_a}\n"
            "  - protein:\n"
            "      id: B\n"
            f"      sequence: {seq_b}\n"
        )
        print(f"  {name:22} {len(seq_a):>3}+{len(seq_b):<3} aa  {desc}")

    runner = OUT / "run_heterodimer.sh"
    runner.write_text(
        "#!/usr/bin/env bash\n"
        "# P4.1 -- controle positif Rv2514c/Rv2515c (TA validé). J1 (Rv2516c/Rv2517c) a DEJA reussi\n"
        "# le 2026-08-05 (ipTM 0,778) -- ne pas le relancer, seul J2 reste dans la boucle ci-dessous.\n"
        "# NE PAS LANCER si un autre boltz predict tourne deja ou si la memoire est sous tension :\n"
        "# J2 a ete tue en silence DEUX FOIS (2026-08-04 et 2026-08-05) par contention memoire\n"
        "# concurrente avec des jobs boltz d'AUTRES projets sur cette machine partagee -- d'ou le\n"
        "# pre-vol ci-dessous, qui SAUTE le job plutot que de le relancer dans l'espoir que ca passe.\n"
        "set -u\n"
        "trap '' HUP  # survit a la fermeture du shell qui l'a lance (P9.2 : vecu le 2026-08-01,\n"
        "             # un runner attache au shell est mort 4 min apres son lancement, en silence)\n"
        'cd "$(dirname "$0")"\n'
        f"B={BOLTZ}\n"
        "STATUS=heterodimer_status.log\n"
        "\n"
        "if [ ! -x \"$B\" ]; then  # echouer BRUYAMMENT plutot que job par job (P9.2 : vecu le\n"
        "  # 2026-08-01, le venv Boltz avait disparu et 3 des 4 jobs de P2.4 ont echoue en silence.)\n"
        "  echo \"ERREUR : binaire Boltz introuvable ou non executable : $B\" | tee -a \"$STATUS\" >&2\n"
        "  echo \"  Le venv a-t-il disparu ? cf. skill boltz, section Installation.\" >&2\n"
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
        "for y in j2_rv2514c_rv2515c; do\n"
        "  if ! check_resources; then\n"
        "    echo \"  $y SAUTE (pre-vol ressources echoue, $(date +%H:%M:%S)) -- a relancer manuellement plus tard\" | tee -a \"$STATUS\"\n"
        "    continue\n"
        "  fi\n"
        "  echo \"=== $y $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
        "  \"$B\" predict $y.yaml --out_dir out_$y --accelerator cpu --devices 1 \\\n"
        "      --output_format mmcif --diffusion_samples 1 --num_workers 4 --use_msa_server \\\n"
        "      >> boltz_heterodimer.log 2>&1 && echo \"  $y OK\" | tee -a \"$STATUS\" || echo \"  $y ECHEC\" | tee -a \"$STATUS\"\n"
        "done\n"
        "echo \"=== TERMINE $(date +%H:%M:%S) ===\" | tee -a \"$STATUS\"\n"
    )
    runner.chmod(0o755)
    print(f"\n  runner : {runner}")
    print("  J1 (test) a déjà réussi le 2026-08-05 -- seul J2 (contrôle positif) reste dans la boucle.")
    print("  pré-vol ressources intégré au runner : SAUTE le job plutôt que de le lancer sous contention.")
    print("  lecture : ../Rv1025/analyses/phase3_afmultimer_parse.py (bi-schéma AF3/Boltz)")


if __name__ == "__main__":
    main()
