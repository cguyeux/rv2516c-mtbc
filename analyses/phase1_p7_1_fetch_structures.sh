#!/usr/bin/env bash
# phase1_p7_1_fetch_structures.sh -- P7.1 : recuperation des structures de reference.
#
# Le reseau n'est joignable QUE par curl appele directement en Bash dans ce sandbox
# (urllib et curl-via-subprocess sont bloques, cf. KB python-patterns) : la partie
# reseau est donc isolee ici, et l'analyse (phase1_p7_1_tmalign.py) travaille hors ligne.
#
# AlphaFold DB (alphafold.ebi.ac.uk) est INJOIGNABLE depuis ce sandbox (TLS
# "unexpected eof", teste en v4 et v6, avec et sans User-Agent). On modelise donc
# la requete ET les facteurs sigma par ESMFold (api.esmatlas.com), ce qui a
# l'avantage d'etre une methode HOMOGENE pour tout le jeu compare. Les references
# experimentales (MerR, permutants S6) viennent de la RCSB, qui repond.
#
# Run: bash analyses/phase1_p7_1_fetch_structures.sh
set -u
OUT="résultats/p7_1/structures"
mkdir -p "$OUT"

echo "== 1. Requete : Rv2516c (ESMFold) =="
if [ ! -s "$OUT/query_Rv2516c.pdb" ]; then
  SEQ=$(tail -n +2 data/Rv2516c.faa | tr -d '\n')
  curl -sS -X POST --data "$SEQ" -o "$OUT/query_Rv2516c.pdb" \
    -w "  Rv2516c : HTTP %{http_code}  %{size_download} octets\n" \
    https://api.esmatlas.com/foldSequence/v1/pdb/
else
  echo "  Rv2516c : deja present"
fi

echo "== 2. Facteurs sigma de M. tuberculosis (ESMFold) =="
# ESMFold public a une limite de longueur : sigA (528 aa) et Rv0890c (882 aa)
# echoueront probablement -- c'est attendu et signale, pas une erreur du script.
python3 - <<'PY' > résultats/p7_1/_seqs.tsv
seqs = {}
name = None
for line in open("résultats/p7_1/sigma_factors.faa"):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split("|")[0]
        seqs[name] = ""
    elif name:
        seqs[name] += line
for k, v in seqs.items():
    print(f"{k}\t{v}")
PY

while IFS=$'\t' read -r RV SEQ; do
  [ -z "$RV" ] && continue
  F="$OUT/sigma_${RV}.pdb"
  if [ -s "$F" ]; then echo "  $RV : deja present"; continue; fi
  CODE=$(curl -sS -X POST --data "$SEQ" -o "$F" -w "%{http_code}" \
         --max-time 300 https://api.esmatlas.com/foldSequence/v1/pdb/ || echo "000")
  SZ=$(stat -c%s "$F" 2>/dev/null || echo 0)
  printf "  %-9s HTTP %s  %s octets  (len %s)\n" "$RV" "$CODE" "$SZ" "${#SEQ}"
  # une reponse non-PDB (erreur JSON/HTML) est jetee pour ne pas polluer l'analyse
  head -1 "$F" 2>/dev/null | grep -q "^HEADER" || { rm -f "$F"; echo "    -> reponse non-PDB, ecartee"; }
  sleep 2
done < résultats/p7_1/_seqs.tsv

echo "== 3. References experimentales (RCSB) =="
# MerR lies a leur promoteur + permutants circulaires de S6 : ce sont exactement
# les cibles des 8 hits Foldseek non significatifs de la fiche atlas.
for ID in 3hh0 5d8c 5d90 7b90 7bff 7bfd; do
  F="$OUT/pdb_${ID}.pdb"
  if [ -s "$F" ]; then echo "  $ID : deja present"; continue; fi
  curl -sS -o "$F" -w "  $ID : HTTP %{http_code}  %{size_download} octets\n" \
    "https://files.rcsb.org/download/${ID}.pdb"
  sleep 1
done

echo "== Fait. Contenu de $OUT =="
ls -1 "$OUT" | sed 's/^/  /'
