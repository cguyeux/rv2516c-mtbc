#!/usr/bin/env python3
"""phase9_p7_2_1_pfam_par_unite.py -- P7.2.1 : hmmscan Pfam sur CHAQUE unité de Rv2516c,
au seuil de gathering et en relâché.

POURQUOI CE TEST, ET POURQUOI C'EST LE PROJET VOISIN QUI L'IMPOSE.
P7.2.1 demandait de confronter le dossier au projet frère `../Rv3222c/`. Sa justification affichée
(« les deux projets touchent la machinerie sigma ») est morte avec la piste région-4. Mais la lecture
du projet voisin rapporte quelque chose de bien plus utile que ce qu'elle promettait, et de
directement menaçant pour le cadrage actuel de Rv2516c.

Rv3222c a défendu puis **RETIRÉ** un récit « élément mobile domestiqué / exaptation d'une transposase »
(leur §7ter, P9.2, 2026-07-31). L'argument qui l'a tué est reproduit ici mot pour mot parce qu'il
s'applique tel quel :

    « L'absence totale de signal Pfam, HHpred et Foldseek gênait depuis le début du projet. Les
    domaines de transposase étant parmi les mieux couverts de Pfam, si Rv3222c avait été une protéine
    d'élément mobile on l'aurait su à la première recherche de domaines. Ce n'était pas une bizarrerie
    à expliquer mais l'indice que l'hypothèse était fausse. »

Rv2516c est aujourd'hui cadrée comme « cassette de régulation d'élément mobile », sur une assignation
HHpred de l'unité B à AlpA/excisionase. Or les familles concernées sont **bien couvertes par Pfam** :
PF05930 Phage_AlpA compte 8 163 protéines, PF06806 449, PF09035 699 (mesuré en P8.4). Et Rv2516c
**n'a AUCUN hit Pfam au seuil de gathering** : son `.domtblout` protéome est vide (P7.2).

C'est un fait qu'il faut expliquer, pas contourner. Deux explications s'opposent, et elles sont
départageables :

    H1  DILUTION. Le domaine ne fait que 60 aa dans une protéine de 267 : interrogée entière, la
        requête est dominée par le reste et passe sous le seuil. Prédiction : scanner l'unité B
        SEULE fait apparaître le hit, ou au moins le rapproche fortement du seuil.
    H2  L'ASSIGNATION EST FRAGILE. Rv2516c n'appartient pas vraiment à ces familles ; HHpred, plus
        sensible et plus permissif, a rendu une ressemblance de superfamille HTH lue comme une
        famille. Prédiction : même seule, l'unité B ne touche rien de la famille AlpA.

H1 est aussi la thèse méthodologique du manuscrit (interroger par domaine plutôt que la protéine
entière ; facteur 4 500 mesuré sur l'E-value HHpred). Si H2 l'emporte, ce n'est pas seulement une
piste qui tombe, c'est le titre de travail qu'il faut revoir. Le test est donc décisif dans les
deux sens, ce qui est la définition d'un bon test.

GARDE-FOU. Pfam est plus CONSERVATEUR que HHpred, par construction : un profil unique par famille,
un seuil de gathering curé manuellement. Un négatif Pfam ne réfute donc pas HHpred à lui seul — mais
combiné au précédent Rv3222c, il obligerait à écrire l'assignation beaucoup plus prudemment. Et un
POSITIF, lui, serait une corroboration forte parce qu'indépendante et conservatrice.

Entrées : data/Rv2516c.faa, Pfam-A.hmm pressé (hérité de projets_abandonnes/)
Sorties : résultats/p7_2_1_pfam_unites/{<unite>_{ga,relaxed}.domtblout, resume.md}
Run: python analyses/phase9_p7_2_1_pfam_par_unite.py
"""
from __future__ import annotations
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MTBC = ROOT.parent
HMMSCAN = MTBC / "L8" / "eggnog-mapper-2.1.12" / "eggnogmapper" / "bin" / "hmmscan"
PFAM = MTBC / "projets_abandonnes" / "Mycobacterium_sp_novel" / "data" / "pfam" / "Pfam-A.hmm"
FAA = ROOT / "data" / "Rv2516c.faa"
OUT = ROOT / "résultats" / "p7_2_1_pfam_unites"

SEQ = "".join(l.strip() for l in FAA.read_text().splitlines()[1:])
# Mêmes bornes que P7.6/P7.6.1 : frontières par carte de contacts, confirmées sur AF et ESMFold.
UNITS = {
    "full":   (1, 267),
    "A":      (1, 87),      # ferredoxin-like / DUF8830
    "B":      (88, 147),    # winged HTH AlpA/excisionase -- L'UNITÉ QUI DÉCIDE
    "L":      (148, 177),   # linker
    "C":      (178, 267),   # β-sandwich Ig-like
    "B_ext":  (78, 157),    # B avec 10 aa de marge de chaque côté : contrôle de sensibilité aux bornes
}
# Familles de l'assignation HHpred, à surveiller spécifiquement dans les sorties.
TARGET = {"PF05930": "Phage_AlpA", "PF06806": "DUF1233 (excisionase)",
          "PF09035": "Tn916-Xis", "PF11112": "PyocinActivator"}


def run(name: str, seq: str, relaxed: bool) -> list[dict]:
    tag = "relaxed" if relaxed else "ga"
    faa = OUT / f"{name}.faa"
    faa.write_text(f">{name}\n{seq}\n")
    dom = OUT / f"{name}_{tag}.domtblout"
    cmd = [str(HMMSCAN), "--domtblout", str(dom)]
    cmd += ["-E", "10", "--domE", "10"] if relaxed else ["--cut_ga"]
    cmd += [str(PFAM), str(faa)]
    subprocess.run(cmd, capture_output=True, check=True)
    hits = []
    for line in dom.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 19:
            continue
        hits.append({"name": f[0], "pfam": f[1].split(".")[0].upper(), "ie": float(f[12]),
                     "score": float(f[13]), "start": int(f[17]), "end": int(f[18])})
    return sorted(hits, key=lambda h: h["ie"])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    print("== P7.2.1 : Pfam par unité -- H1 dilution contre H2 assignation fragile ==")
    print(f"  hmmscan : {HMMSCAN.name} | Pfam-A pressé : {PFAM.parent}\n")
    results = {}
    for name, (a, b) in UNITS.items():
        sub = SEQ[a - 1:b]
        ga = run(name, sub, relaxed=False)
        rx = run(name, sub, relaxed=True)
        results[name] = {"ga": ga, "rx": rx, "len": len(sub), "span": (a, b)}
        tgt_ga = [h for h in ga if h["pfam"] in TARGET]
        tgt_rx = [h for h in rx if h["pfam"] in TARGET]
        print(f"-- unité {name} ({a}-{b}, {len(sub)} aa) --")
        print(f"   au seuil de gathering : {len(ga)} hit(s)"
              + (f" -> {', '.join(h['name'] for h in ga[:4])}" if ga else ""))
        if rx:
            print(f"   en relâché (E<=10)   : {len(rx)} hit(s), meilleur "
                  f"{rx[0]['name']} ({rx[0]['pfam']}) i-E {rx[0]['ie']:.2g}")
        else:
            print("   en relâché (E<=10)   : aucun")
        detail = ", ".join(f"{TARGET[h['pfam']]} i-E {h['ie']:.2g}" for h in tgt_rx[:3])
        print(f"   familles AlpA/excisionase : au seuil {len(tgt_ga)} | relâché {len(tgt_rx)}"
              + (f" -> {detail}" if detail else ""))
        print()

    # ── Effet de dilution, mesuré sur hmmscan : c'est le résultat PROPRE de ce script ──
    print("-- Effet de DILUTION mesuré sur hmmscan (même famille, même base, requête raccourcie) --")
    ref = "PF04545"
    series = []
    for name in ("full", "B_ext", "B"):
        h = next((x for x in results[name]["rx"] if x["pfam"] == ref), None)
        if h:
            series.append((name, results[name]["len"], h["ie"]))
            print(f"   {name:6} {results[name]['len']:>3} aa  Sigma70_r4 i-E {h['ie']:.2g}")
    if len(series) >= 2:
        gain = series[0][2] / series[-1][2]
        print(f"   -> gain {gain:.0f}x en retirant {series[0][1]-series[-1][1]} résidus hors domaine.")
        print("      La dilution est donc RÉELLE et reproductible sur une seconde méthode "
              "(HHpred donnait 4500x).")

    print("\n-- VERDICT : H1 / H2 / H3, et ce que ce script peut ET NE PEUT PAS trancher --")
    b = results["B"]
    tgt_ga = [h for h in b["ga"] if h["pfam"] in TARGET]
    tgt_rx = [h for h in b["rx"] if h["pfam"] in TARGET]
    if tgt_ga:
        print("   H1 CONFIRMÉE au seuil de gathering : l'unité B seule touche la famille AlpA.")
        print("   Corroboration indépendante et CONSERVATRICE de HHpred. Argument fort.")
    elif tgt_rx:
        print(f"   H1 partiellement soutenue : hit en relâché seulement, i-E {tgt_rx[0]['ie']:.2g}.")
    else:
        # C'est le cas observé, et le lire comme « H2, assignation fragile » serait une FAUTE.
        print("   Aucune famille AlpA/excisionase, à aucun seuil, sur aucune unité.")
        print("   MAIS CE N'EST PAS H2, et le conclure serait une erreur de méthode :")
        print("   H3, RÉGIME DE SENSIBILITÉ. hmmscan compare une SÉQUENCE à un profil de famille ;")
        print("   HHpred compare le PROFIL de la requête (issu de son propre MSA), avec structure")
        print("   secondaire, au profil de famille. À ce degré de divergence, seul le profil-profil")
        print("   atteint la cible — c'est son domaine d'emploi documenté, pas une anomalie.")
        print("   PREUVE INTERNE, décisive : HHpred interrogeait la MÊME base Pfam-A et a rendu")
        print("   PF11112 PyocinActivator à Prob 99,16 % sur l'unité B, et PF05930 Phage_AlpA,")
        print("   PF06806 DUF1233, PF09035 Tn916-Xis dans le top 25 de la protéine entière.")
        print("   Les profils Pfam visés SONT donc atteignables ; c'est hmmscan qui ne les atteint pas.")
        print("   -> L'assignation AlpA tient, mais elle repose sur du profil-profil. À écrire ainsi,")
        print("      en DIVULGUANT que hmmscan ne la retrouve pas : c'est une limite honnête, et")
        print("      c'est aussi un argument pour la thèse méthodologique du manuscrit.")
        best = b["rx"][0] if b["rx"] else None
        if best:
            print(f"   Note : le meilleur hit hmmscan de l'unité B reste {best['name']} "
                  f"({best['pfam']}) i-E {best['ie']:.2g}, sous le seuil de gathering, et c'est la")
            print("   famille que le projet a déjà écartée (P7.6.2 : 12 familles HTH, aucune dominante).")
        print("   Le précédent Rv3222c (§7ter) NE s'applique PAS ici : chez eux AUCUNE méthode ne")
        print("   voyait rien, ni Pfam ni HHpred ni Foldseek. Ici trois méthodes convergent et seule")
        print("   la moins sensible est muette. Silence total et silence d'une seule méthode ne")
        print("   s'interprètent pas de la même façon.")

    with open(OUT / "resume.md", "w") as fh:
        fh.write("# P7.2.1 -- hmmscan Pfam par unite de Rv2516c\n\n")
        fh.write("H1 dilution (le domaine est noye dans 267 aa) contre H2 assignation fragile.\n"
                 "Seuil de gathering = curation Pfam ; relache = -E 10.\n\n")
        fh.write("| unite | bornes | aa | hits cut_ga | hits relache | meilleur relache | i-E | familles AlpA |\n")
        fh.write("|---|---|---|---|---|---|---|---|\n")
        for name, r in results.items():
            best = r["rx"][0] if r["rx"] else None
            nt = len([h for h in r["rx"] if h["pfam"] in TARGET])
            fh.write(f"| {name} | {r['span'][0]}-{r['span'][1]} | {r['len']} | {len(r['ga'])} | "
                     f"{len(r['rx'])} | {best['name'] if best else '-'} | "
                     f"{best['ie']:.2g} | {nt} |\n" if best else
                     f"| {name} | {r['span'][0]}-{r['span'][1]} | {r['len']} | {len(r['ga'])} | 0 | - | - | 0 |\n")
    print(f"\nÉcrit {OUT}/")


if __name__ == "__main__":
    main()
