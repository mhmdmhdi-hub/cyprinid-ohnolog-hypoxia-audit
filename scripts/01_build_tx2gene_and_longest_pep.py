#!/usr/bin/env python3
import argparse, gzip, os, re

def build_ccar_tx2gene(gff3_path: str, out_tsv: str) -> None:
    tx2gene = {}
    with open(gff3_path, errors="ignore") as f:
        for ln in f:
            if ln.startswith("#") or not ln.strip(): continue
            ps = ln.rstrip("\n").split("\t")
            if len(ps) < 9 or ps[2] not in ("mRNA", "transcript"): continue
            mid = re.search(r"ID=([^;]+)", ps[8])
            par = re.search(r"Parent=([^;]+)", ps[8])
            if mid and par:
                tx2gene[mid.group(1).strip()] = par.group(1).split(",")[0].strip()
    with open(out_tsv, "w") as out:
        out.write("TXNAME\tGENEID\n")
        for tx, g in sorted(tx2gene.items()):
            out.write(f"{tx}\t{g}\n")

def extract_ccar_longest_50chr(pep_path: str, out_fa: str) -> None:
    best_seq, best_hdr = {}, {}
    with open(pep_path, errors="ignore") as f:
        cur_id = cur_gene = cur_chr = None
        seq = []
        def flush():
            if cur_gene and cur_chr and re.match(r"^ccar_u[ab]\d{2}$", cur_chr):
                s = "".join(seq).rstrip("*")
                if cur_gene not in best_seq or len(s) > len(best_seq[cur_gene]):
                    best_seq[cur_gene] = s
                    best_hdr[cur_gene] = f"{cur_gene}|{cur_chr}|{cur_id}"
        for ln in f:
            ln = ln.strip()
            if not ln: continue
            if ln.startswith(">"):
                flush()
                cur_id = ln[1:].split()[0]
                gm = re.search(r"gene=([^\s]+)", ln)
                cm = re.search(r"seq_id=([^\s]+)", ln)
                cur_gene = gm.group(1) if gm else None
                cur_chr  = cm.group(1) if cm else None
                seq = []
            else:
                seq.append(ln)
        flush()
    with open(out_fa, "w") as out:
        for g in sorted(best_seq):
            out.write(f">{best_hdr[g]}\n{best_seq[g]}\n")

def extract_carp_longest_50chr(carp_pep_gz: str, out_fa: str) -> None:
    best_seq, best_chr = {}, {}
    opener = gzip.open if carp_pep_gz.endswith((".gz", ".bgz")) else open
    with opener(carp_pep_gz, "rt", errors="ignore") as f:
        cur_g = cur_c = None
        seq = []
        def flush():
            if cur_g and cur_c:
                s = "".join(seq).rstrip("*")
                if cur_g not in best_seq or len(s) > len(best_seq[cur_g]):
                    best_seq[cur_g] = s
                    best_chr[cur_g] = cur_c
        for ln in f:
            ln = ln.strip()
            if not ln: continue
            if ln.startswith(">"):
                flush()
                cm = re.search(r"ASM1834038v1:([AB]\d+):", ln)
                gm = re.search(r"gene:(ENSCCRG\d+)", ln)
                cur_c = cm.group(1) if cm else None
                cur_g = gm.group(1) if gm else None
                seq = []
            else:
                seq.append(ln)
        flush()
    with open(out_fa, "w") as out:
        for g in sorted(best_seq):
            out.write(f">{g}|{best_chr[g]}\n{best_seq[g]}\n")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--ccar-gff3", default=os.path.expanduser("~/carassius_native_2025/ccar_v5.gff3"))
    p.add_argument("--ccar-pep", default=os.path.expanduser("~/carassius_native_2025/ccar_v5_proteins.fa"))
    p.add_argument("--carp-pep", default=os.path.expanduser("~/rbh_verify/carp_pep.fa.bgz"))
    p.add_argument("--outdir", default=os.path.expanduser("~/carassius_native_2025"))
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    build_ccar_tx2gene(args.ccar_gff3, os.path.join(args.outdir, "tx2gene_native.tsv"))
    extract_ccar_longest_50chr(args.ccar_pep, os.path.join(args.outdir, "ccar_v5_longest_50chr.fa"))
    extract_carp_longest_50chr(args.carp_pep, os.path.join(args.outdir, "carp_longest_50chr.fa"))
