#!/usr/bin/env python3
"""Modèles nuls de séquence, partagés par les analyses du projet.

Créé le 2026-08-10 en REFORGEANT `dinuc_shuffle`, qui vivait dupliquée dans
`phase10_p4_3bis_operateur.py` et `phase19_p4_3bis_c_direct_repeats.py` et dont le contrat
annoncé était FAUX. Les deux copies promettaient « préservant la composition en dinucléotides »
et ne la préservaient pas toujours.

LE DÉFAUT, ET POURQUOI IL A ÉCHAPPÉ SI LONGTEMPS. Ces implémentations tiraient un ordre
aléatoire des successeurs de chaque nucléotide puis parcouraient le graphe. Quand la marche
arrivait sur un nucléotide dont tous les successeurs étaient consommés AVANT d'avoir épuisé
toutes les arêtes, elles « redémarraient sur un résidu disponible » : ce saut invente un
dinucléotide absent de l'original et en perd un autre. Le docstring de phase10 le concédait
d'ailleurs à demi-mot, « marche eulérienne APPROCHÉE », mais la ligne de sortie affichée à
l'écran, elle, annonçait un nul à composition préservée. C'est exactement le problème que
résout l'algorithme d'Altschul & Erickson (1985, Mol. Biol. Evol. 2:526-538) en tirant d'abord
un arbre couvrant orienté vers le nucléotide terminal, ce qui garantit qu'un chemin eulérien
existe avant de le parcourir.

ET LA LEÇON LA PLUS COÛTEUSE DE L'AFFAIRE : cet algorithme correct était DÉJÀ IMPLÉMENTÉ dans
le dépôt, à deux répertoires d'ici, dans `mtbc/SpacerEgalVirus/analyses/phase1_test_nul_permutation.py`
(tirage d'arbre couvrant, docstring citant la référence). Rv2516c a hérité d'une version naïve
et l'a dupliquée deux fois sans jamais aller voir si le voisin avait résolu le problème. Le
défaut n'a donc rien coûté à trouver une fois cherché : il a coûté de ne pas chercher. Réflexe
à garder : `grep -rl "def <nom_de_la_fonction>"` sur `mtbc/` AVANT de recopier une routine
d'un script à l'autre.

CE QUI EST FAIT ICI. Même marche, mais l'échec est traité par REJET plutôt que par saut : si la
marche n'épuise pas toutes les arêtes, on retire un nouvel ordonnancement, jusqu'à `essais`.
Conditionnellement à l'acceptation, cela échantillonne uniformément parmi les séquences valides,
même distribution que le tirage d'arbre couvrant d'Altschul-Erickson, pour un coût négligeable
sur un alphabet de 4 lettres. Si aucun tirage n'aboutit, on lève une exception : mieux vaut un
arrêt bruyant qu'un nul silencieusement faux, qui est précisément ce qui s'est passé jusqu'ici.
Vérifié le 2026-08-10 : 0 altération sur 500 tirages du duplex 6AMA, contre 94 % avec la
version naïve. Qui préfère l'algorithme de référence sans rejet a le modèle sous les yeux dans
SpacerEgalVirus, cité ci-dessus ; les deux sont équivalents en distribution.

PORTÉE SUR LES RÉSULTATS DÉJÀ ACQUIS — MESURÉE, car la première explication venue était fausse.
Le réflexe est de dire « biais rare, donc dilué sur 400 permutations ». C'est faux : mesuré le
2026-08-10, le défaut frappe **20 à 94 % des tirages selon la séquence** (94 % sur le duplex
6AMA de 22 pb, ~70 % sur des séquences de 50 à 400 nt à 65 % GC). Il n'est pas rare du tout.

Ce qui sauve `phase10` et `phase19` est autre chose, et il fallait le mesurer plutôt que le
supposer : le rapiéçage est fréquent mais **minuscule**. Il ne déplace que ~1 % des
dinucléotides par tirage (0,23 sur 21 pour une séquence de 22 nt, 1,34 sur 149 pour 150 nt,
maximum observé 3). Conséquence sur la statistique réellement utilisée : le score de meilleure
répétition directe sous le nul vaut 6,990 (ancien) contre 6,955 (exact) sur 200 tirages d'une
séquence de 150 nt, soit un écart de -0,035 pour des écarts-types de 0,93 à 1,07. Les deux nuls
sont donc interchangeables en pratique, et l'ancien est très légèrement plus DISPERSÉ, donc
plutôt plus conservatif. **P4.3bis-a et -c ne sont pas à refaire et leurs conclusions négatives
tiennent** — pour cette raison mesurée, pas pour la raison intuitive qui aurait été erronée.

En revanche tout usage FUTUR, et en particulier tout nul en exemplaire UNIQUE (comme le contrôle
négatif E2 de P4.3bis-e), doit passer par la version corrigée. Non parce que l'effet statistique
serait grand, il ne l'est pas, mais parce qu'un contrôle unique dont on ne peut pas AFFIRMER que
la composition est identique se réfute en une phrase : « votre comparaison confond composition
et arrangement ». Sur une distribution on discute une moyenne ; sur un exemplaire unique on ne
discute rien, on vérifie.

Conséquence pratique : les tirages produits pour une graine donnée DIFFÈRENT de ceux de
l'ancienne version. Relancer phase10 ou phase19 ne redonnera donc pas au chiffre près les
p-values consignées dans le cahier au 2026-08-03. C'est le prix d'un contrat enfin honoré, et
c'est sans effet sur les conclusions, qui étaient négatives et le restent.
"""
from __future__ import annotations

import random


class NulImpossible(RuntimeError):
    """Aucun mélange dinucléotidique valide trouvé : on refuse de rendre un nul approximatif."""


def dinuc_shuffle(seq: str, rng: random.Random, essais: int = 200) -> str:
    """Mélange préservant EXACTEMENT la composition en dinucléotides.

    Préserver le seul %GC sous-estimerait le nul : la structure dinucléotidique fabrique déjà
    des palindromes et des répétitions, surtout à 65 % GC comme chez le MTBC.

    Lève `NulImpossible` plutôt que de rendre une séquence à composition altérée.
    """
    if len(seq) < 3:
        return seq

    succ: dict[str, list[str]] = {}
    for a, b in zip(seq, seq[1:]):
        succ.setdefault(a, []).append(b)
    n_aretes = len(seq) - 1

    for _ in range(essais):
        pool = {k: list(v) for k, v in succ.items()}
        for v in pool.values():
            rng.shuffle(v)
        out, cur = [seq[0]], seq[0]
        for _ in range(n_aretes):
            nxt = pool.get(cur)
            if not nxt:
                break                    # impasse : ce tirage est REJETÉ, jamais rapiécé
            c = nxt.pop()
            out.append(c)
            cur = c
        if len(out) == len(seq):         # marche complète : toutes les arêtes consommées
            return "".join(out)

    raise NulImpossible(
        f"aucune marche eulérienne complète en {essais} tirages pour une séquence de "
        f"{len(seq)} nt. Séquence trop contrainte : revoir le modèle nul, ne pas relâcher "
        f"la contrainte de composition en silence."
    )
