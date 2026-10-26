"""
Helpers for a variant calling benchmark on a targeted panel on a GIAB reference sample:
build a panel BED from a GTF, summarize coverage and uniformity, compare calls to a truth set,
and collect metrics across a depth titration. Pure Python (pandas, numpy) so they can be unit tested.
"""

import gzip
import re

import numpy as np
import pandas as pd


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# ---------------------------------------------------------------- panel BED
def merge_intervals(df):
    """Merge overlapping or adjacent intervals per chromosome (and gene, if present). Returns a sorted DataFrame."""
    keys = ["chrom", "gene"] if "gene" in df.columns else ["chrom"]
    out = []
    for key, g in df.sort_values(keys + ["start", "end"]).groupby(keys, sort=False):
        cur_s, cur_e = None, None
        for s, e in zip(g["start"], g["end"]):
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e:
                cur_e = max(cur_e, e)
            else:
                out.append((*((key,) if isinstance(key, str) else key), cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            out.append((*((key,) if isinstance(key, str) else key), cur_s, cur_e))
    res = pd.DataFrame(out, columns=keys + ["start", "end"])
    return res.sort_values(["chrom", "start"]).reset_index(drop=True)


def gtf_to_panel_bed(gtf_path, genes, feature="CDS", pad=10):
    """
    Coding exons (feature rows) of the given genes from a GENCODE/Ensembl GTF, padded by `pad` bp
    and merged. Returns chrom, start (0-based), end, gene. Genes not found are reported in .attrs['missing'].
    """
    wanted = set(genes)
    rows = []
    with _open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != feature:
                continue
            m = re.search(r'gene_name "([^"]+)"', parts[8])
            if not m or m.group(1) not in wanted:
                continue
            start = max(0, int(parts[3]) - 1 - pad)   # GTF is 1-based inclusive; BED is 0-based, end exclusive
            rows.append((parts[0], start, int(parts[4]) + pad, m.group(1)))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene"])
    bed = merge_intervals(df) if len(df) else df
    bed = bed[["chrom", "start", "end", "gene"]]   # BED column order; merge_intervals puts group keys first
    bed.attrs["missing"] = sorted(wanted - set(df["gene"]))
    return bed


# ---------------------------------------------------------------- coverage
def per_base_depth_in_panel(perbase, panel):
    """
    Depth values weighted by length restricted to panel intervals.
    perbase: DataFrame chrom,start,end,depth (mosdepth per base runs); panel: chrom,start,end[,gene].
    Returns DataFrame of (gene, length, depth) pieces covering the panel.
    """
    pieces = []
    for chrom, pg in panel.groupby("chrom"):
        runs = perbase[perbase["chrom"] == chrom].sort_values("start")
        rs, re_, rd = runs["start"].to_numpy(), runs["end"].to_numpy(), runs["depth"].to_numpy()
        for _, iv in pg.iterrows():
            lo = np.searchsorted(re_, iv["start"], side="right")
            covered = iv["start"]
            for k in range(lo, len(rs)):
                if rs[k] >= iv["end"]:
                    break
                s, e = max(rs[k], iv["start"]), min(re_[k], iv["end"])
                if s > covered:          # gap in the per base output = zero depth
                    pieces.append((iv.get("gene", chrom), s - covered, 0))
                if e > s:
                    pieces.append((iv.get("gene", chrom), e - s, rd[k]))
                covered = max(covered, e)
            if covered < iv["end"]:
                pieces.append((iv.get("gene", chrom), iv["end"] - covered, 0))
    return pd.DataFrame(pieces, columns=["gene", "length", "depth"])


def coverage_summary(pieces, thresholds=(20, 50, 100)):
    """Mean depth, % of bases at or above each threshold, and uniformity (% of bases >= 0.2x mean)."""
    total = pieces["length"].sum()
    mean = (pieces["length"] * pieces["depth"]).sum() / total if total else 0.0
    out = {"bases": int(total), "mean_depth": float(mean)}
    for t in thresholds:
        out[f"pct_ge_{t}x"] = float(100 * pieces.loc[pieces["depth"] >= t, "length"].sum() / total) if total else 0.0
    out["uniformity_pct_ge_0.2x_mean"] = (
        float(100 * pieces.loc[pieces["depth"] >= 0.2 * mean, "length"].sum() / total) if total else 0.0)
    return out


def read_mosdepth_perbase(path):
    return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "depth"])


# ---------------------------------------------------------------- calls vs truth
def read_vcf_variants(path, bed=None):
    """Variants as a set of (chrom, pos, ref, alt), split per ALT allele, optionally restricted to BED intervals."""
    out = set()
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c, pos, _, ref, alts = line.split("\t")[:5]
            for alt in alts.split(","):
                if alt not in (".", "*"):
                    out.add((c, int(pos), ref, alt))
    if bed is not None:
        out = {v for v in out if _in_bed(v, bed)}
    return out


def _in_bed(v, bed):
    sub = bed[bed["chrom"] == v[0]]
    p0 = v[1] - 1
    return bool(((sub["start"] <= p0) & (p0 < sub["end"])).any())


def variant_type(v):
    return "SNP" if len(v[2]) == 1 and len(v[3]) == 1 else "INDEL"


def compare_calls(truth, query):
    """
    Exact-match comparison by type (after both VCFs are normalized with `bcftools norm`).
    A fallback when hap.py is unavailable; hap.py handles representation differences better.
    """
    rows = []
    for t in ("SNP", "INDEL"):
        T = {v for v in truth if variant_type(v) == t}
        Q = {v for v in query if variant_type(v) == t}
        tp, fn, fp = len(T & Q), len(T - Q), len(Q - T)
        rec = tp / (tp + fn) if tp + fn else float("nan")
        prec = tp / (tp + fp) if tp + fp else float("nan")
        f1 = 2 * rec * prec / (rec + prec) if rec + prec else float("nan")
        rows.append({"type": t, "truth_total": len(T), "tp": tp, "fn": fn, "fp": fp,
                     "recall": rec, "precision": prec, "f1": f1})
    return pd.DataFrame(rows)


def parse_happy_summary(path):
    """Recall/precision/F1 for PASS SNPs and indels from a hap.py *.summary.csv."""
    df = pd.read_csv(path)
    df = df[df["Filter"] == "PASS"]
    return pd.DataFrame({
        "type": df["Type"], "truth_total": df["TRUTH.TOTAL"], "tp": df["TRUTH.TP"], "fn": df["TRUTH.FN"],
        "fp": df["QUERY.FP"], "recall": df["METRIC.Recall"], "precision": df["METRIC.Precision"],
        "f1": df["METRIC.F1_Score"]}).reset_index(drop=True)
