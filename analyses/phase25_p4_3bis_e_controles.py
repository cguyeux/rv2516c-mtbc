#!/usr/bin/env python3
"""P4.3bis-e -- le triangle de contrôles qui rend E1 interprétable (ou le tue).

E1 (P4.3bis-d, 2026-08-09) a rendu ipTM 0,930 pour 2×unité B + duplex 22 pb, le meilleur
score de toute la série. Il n'est PAS lisible en l'état, pour deux raisons pré-enregistrées
avant tout calcul :

  (i)  le duplex est l'opérateur de BldC (6AMA), séquence HÉTÉROLOGUE : rien ne dit que
       le score mesure une reconnaissance de séquence plutôt que l'affinité générique
       d'un wHTH pour un sillon d'ADN B ;
  (ii) le pipeline SUR-NOTE de façon démontrée : D2 (2×unité A apo, domaine sans aucune
       fonction de liaison à l'ADN connue) a obtenu 0,873, soit mieux que le contrôle
       positif J4 (0,827).

Trois jobs, tous au même paramétrage que E1, à comparer sur ipTM / PAE inter-chaînes /
nombre de contacts d'interface :

  E2  2×unité B + duplex BROUILLÉ   -- négatif de SÉQUENCE. Le test décisif.
  E3  2×unité A + duplex intact     -- négatif de DOMAINE, le plus sévère (D2 a sur-noté).
  E4  2×BldC     + duplex intact    -- positif direct : le co-cristal 6AMA reconstitué,
                                       qui fixe enfin le plafond dans le bon contexte.

LECTURE PRÉ-ENREGISTRÉE, avant tout calcul : E1 n'est un signal
spécifique QUE si E1 > E2 nettement ET E1 > E3 nettement, E4 servant de plafond. Si E2 ou
E3 talonne E1, le verdict est « le pipeline sur-note les complexes wHTH-ADN », résultat de
méthode négatif, à écrire tel quel sans le maquiller.

Le brouillage réutilise `dinuc_shuffle` de phase19 (elle-même identique à phase10) : c'est le
MÊME modèle nul que celui des recherches de site P4.3bis-a et -c, la cohérence de méthode est
donc assurée par import, pas par copie.

RÉSERVE MÉTHODOLOGIQUE TROUVÉE ICI, 2026-08-10, et sa parade. L'implémentation de phase19 est
une marche eulérienne NAÏVE : quand elle se bloque sur un nucléotide dont tous les successeurs
sont consommés, elle SAUTE vers un autre nucléotide encore disponible. Ce saut invente un
dinucléotide absent de l'original et en perd un autre, donc la composition n'est pas
strictement préservée à chaque tirage (c'est le défaut que corrige l'algorithme d'Altschul &
Erikson 1985 par tirage préalable d'un arbre couvrant). Sur 400 permutations servant à bâtir
une distribution nulle, ce biais occasionnel est dilué et conservatif : phase10 et phase19 ne
sont PAS à refaire. Mais E2 est un contrôle négatif UNIQUE, dont un reviewer pourrait dire que
la comparaison confond composition et arrangement. On garde donc la fonction de phase19 et on
ajoute un FILTRE D'ACCEPTATION explicite : tirages successifs jusqu'au premier qui préserve
exactement la composition dinucléotidique et diffère de l'original sur au moins 6 positions.
La graine retenue est imprimée et consignée. Le score de répétition directe est RAPPORTÉ mais
jamais filtré : filtrer dessus reviendrait à fabriquer le contrôle qui arrange la conclusion.

Le brin complémentaire est recalculé après brouillage, pour que E2 reste un duplex apparié :
sinon on remplacerait un artefact par un autre (mésappariements).
"""
from __future__ import annotations

import importlib.util
import random
import sys
from pathlib import Path

ICI = Path(__file__).resolve().parent
RACINE = ICI.parent
SORTIE = RACINE / "résultats" / "p4_3bis_e_controles"

# --- séquences, reprises VERBATIM des YAML déjà lancés (aucune re-dérivation) -------------
# unité B : résultats/p4_3bis_dimere/d1_unitB_dimere.yaml et p4_3bis_d_dimer_dna/e1_*.yaml
UNITE_B = "RAEAFTTPELMSAAEIADELGVSRQRVHQLRSTAGFPAPLADLRGGAVWDAAAVRRFAET"
# unité A : résultats/p4_3bis_dimere/d2_unitA_dimere.yaml
UNITE_A = "MTADWVVTFTFDADPSMETMDAWETQLEGFDALVSRVPGHGIDVTVYAPGDWSVFDALAKMAGEVMPVVQAKSPIAVQIISEPEHRL"
# BldC : résultats/p4_3bis_dimere/d3_bldc_dimere.yaml
BLDC = "MTARTPDAEPLLTPAEVATMFRVDPKTVTRWAKAGKLTSIRTLGGHRRYREAEVRALLAGIPQQRSEA"
# duplex opérateur BldC (6AMA), 22 pb, brin sens tel qu'utilisé dans E1 et J4
DUPLEX = "ATTCGGGTAATTCGGGTAATTC"

COMPL = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(COMPL)[::-1]


def charge_phase19():
    """Importe phase19 pour réutiliser SA fonction de brouillage, pas une copie."""
    chemin = ICI / "phase19_p4_3bis_c_direct_repeats.py"
    spec = importlib.util.spec_from_file_location("phase19", chemin)
    if spec is None or spec.loader is None:
        raise ImportError(f"impossible de charger {chemin}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # guard __main__ présent : rien ne s'exécute
    return mod


def dinucs(seq: str) -> dict[str, int]:
    d: dict[str, int] = {}
    for a, b in zip(seq, seq[1:]):
        d[a + b] = d.get(a + b, 0) + 1
    return d


def yaml_complexe(prot: str, adn_sens: str | None) -> str:
    bloc = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: [A, B]",
        f"      sequence: {prot}",
    ]
    if adn_sens is not None:
        bloc += [
            "  - dna:",
            "      id: C",
            f"      sequence: {adn_sens}",
            "  - dna:",
            "      id: D",
            f"      sequence: {rc(adn_sens)}",
        ]
    return "\n".join(bloc) + "\n"


def main() -> int:
    ph19 = charge_phase19()

    # --- vérification préalable : le duplex de E1 est-il bien apparié ? ------------------
    attendu_brin_D = "GAATTACCCGAATTACCCGAAT"      # brin D du YAML de E1, vérifié ici
    if rc(DUPLEX) != attendu_brin_D:
        print("ERREUR : le brin complémentaire recalculé ne reproduit pas celui de E1.")
        print(f"  recalculé {rc(DUPLEX)}\n  E1        {attendu_brin_D}")
        return 1
    print(f"duplex E1 apparié, contrôle OK : {DUPLEX} / {attendu_brin_D}")

    # --- brouillage dinucléotidique, fonction de P4.3bis-a/-c + filtre d'acceptation ------
    # Critères d'acceptation, fixés AVANT de regarder les tirages :
    #   (1) composition dinucléotidique strictement identique à l'original ;
    #   (2) au moins 6 positions changées sur 22, sans quoi le « négatif » n'en est pas un.
    # Le score de répétition directe n'entre PAS dans les critères : on le rapporte, on ne
    # sélectionne pas dessus.
    ref_dinucs = dinucs(DUPLEX)
    brouille, graine, rejets = None, None, 0
    for g in range(1000):
        cand = ph19.dinuc_shuffle(DUPLEX, random.Random(g))
        if dinucs(cand) != ref_dinucs:
            rejets += 1
            continue
        if sum(1 for a, b in zip(DUPLEX, cand) if a != b) < 6:
            rejets += 1
            continue
        brouille, graine = cand, g
        break
    if brouille is None:
        print("ERREUR : aucun brouillage acceptable en 1000 tirages. Ne rien lancer.")
        return 1

    hamming = sum(1 for a, b in zip(DUPLEX, brouille) if a != b)
    dr_orig = ph19.best_dr(DUPLEX)
    dr_brou = ph19.best_dr(brouille)

    print("\n--- contrôle négatif de séquence (E2) ---")
    print(f"  original  {DUPLEX}   meilleure répétition directe : score {dr_orig[0]}"
          f" (bras {dr_orig[1]}, espaceur {dr_orig[2]})")
    print(f"  brouillé  {brouille}   meilleure répétition directe : score {dr_brou[0]}"
          f" (bras {dr_brou[1]}, espaceur {dr_brou[2]})")
    print(f"  graine retenue random.Random({graine}) après {rejets} tirage(s) rejeté(s)")
    print(f"  distance de Hamming : {hamming}/{len(DUPLEX)}"
          f"   composition dinucléotidique : strictement préservée (vérifiée)")

    # --- écriture des trois YAML ---------------------------------------------------------
    SORTIE.mkdir(parents=True, exist_ok=True)
    jobs = [
        ("e2_unitB_dimer_dna_brouille", UNITE_B, brouille,
         "négatif de SÉQUENCE : 2×unité B + duplex brouillé"),
        ("e3_unitA_dimer_dna", UNITE_A, DUPLEX,
         "négatif de DOMAINE : 2×unité A + duplex intact"),
        ("e4_bldc_dimer_dna", BLDC, DUPLEX,
         "positif DIRECT : 2×BldC + duplex intact (6AMA reconstitué)"),
    ]
    print("\n--- jobs écrits ---")
    for nom, prot, adn, desc in jobs:
        (SORTIE / f"{nom}.yaml").write_text(yaml_complexe(prot, adn))
        tokens = 2 * len(prot) + 2 * len(adn)
        besoin = 5367 + 0.02639 * tokens * tokens         # modèle du pré-vol, num_workers=0
        print(f"  {nom}.yaml  {tokens:>4} tokens  besoin ~{besoin:.0f} Mio"
              f"  seuil pré-vol ~{besoin * 1.2:.0f} Mio")
        print(f"      {desc}")
    print(f"\nÉcrit dans {SORTIE}")
    print("Rappel E1, pour comparaison : 164 tokens, ipTM 0,930, PAE inter-ch. 0,48 Å, 205 contacts.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
