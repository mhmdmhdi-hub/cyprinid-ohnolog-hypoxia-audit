#!/usr/bin/env bash
set -euo pipefail
WORKDIR="${1:-$HOME/carassius_native_2025}"
CARP_PAIRS="${2:-$HOME/carp_hypoxia/strict_homeolog_pairs.csv}"
DIAMOND="${DIAMOND_BIN:-$HOME/miniforge3/envs/rbh_tools/bin/diamond}"
THREADS="${THREADS:-8}"
cd "$WORKDIR"
"$DIAMOND" makedb --in carp_longest_50chr.fa -d carp_db --quiet
"$DIAMOND" makedb --in ccar_v5_longest_50chr.fa -d ccar_db --quiet
"$DIAMOND" blastp -q carp_longest_50chr.fa -d ccar_db -o carp_vs_ccar.tsv --evalue 1e-10 --max-target-seqs 5 --outfmt 6 qseqid sseqid pident length evalue bitscore --threads "$THREADS" --quiet
"$DIAMOND" blastp -q ccar_v5_longest_50chr.fa -d carp_db -o ccar_vs_carp.tsv --evalue 1e-10 --max-target-seqs 5 --outfmt 6 qseqid sseqid pident length evalue bitscore --threads "$THREADS" --quiet
python3 - "$CARP_PAIRS" << 'PYEOF'
import csv, re, sys
def load_best(path):
    best = {}
    with open(path) as f:
        for ln in f:
            ps = ln.rstrip("\n").split("\t")
            q, s, pid, bs = ps[0], ps[1], float(ps[2]), float(ps[5])
            if q not in best or bs > best[q][2]: best[q] = (s, pid, bs)
    return best
k2c, c2k = load_best("carp_vs_ccar.tsv"), load_best("ccar_vs_carp.tsv")
rbh = {kq.split("|")[0]: (cs.split("|")[0], cs.split("|")[1], pid)
       for kq, (cs, pid, _) in k2c.items() if cs in c2k and c2k[cs][0] == kq}
rows_out = []
with open(sys.argv[1]) as f:
    for r in csv.DictReader(f):
        gA, gB = r["gene_id_A"], r["gene_id_B"]
        if gA in rbh and gB in rbh:
            cA, chrA, pidA = rbh[gA]; cB, chrB, pidB = rbh[gB]
            mA, mB = re.match(r"^ccar_u([ab])(\d{2})$", chrA), re.match(r"^ccar_u([ab])(\d{2})$", chrB)
            if mA and mB and mA.group(1) != mB.group(1) and mA.group(2) == mB.group(2):
                rows_out.append({
                    "gene_name": r["gene_name"], "carp_gene_A": gA, "carp_chr_A": r["chr_A"],
                    "carp_gene_B": gB, "carp_chr_B": r["chr_B"],
                    "ccar_orth_of_carpA": cA, "ccar_chr_of_carpA": chrA, "pid_A": round(pidA, 2),
                    "ccar_orth_of_carpB": cB, "ccar_chr_of_carpB": chrB, "pid_B": round(pidB, 2),
                    "ccar_subgenome_of_carpA": mA.group(1).upper(), "ccar_chr_num": int(mA.group(2)),
                })
with open("native_verified_quartets.csv", "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)
PYEOF
