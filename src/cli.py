"""Command line steps for the panel benchmark (called from run_benchmark.sh)."""

import argparse
import json
import os

import pandas as pd

import panelbench as pb


def make_panel(a):
    genes = [g.strip() for g in open(a.genes) if g.strip()]
    bed = pb.gtf_to_panel_bed(a.gtf, genes, pad=a.pad)
    if bed.attrs.get("missing"):
        print("Genes not found in GTF:", ", ".join(bed.attrs["missing"]))
    bed[["chrom", "start", "end", "gene"]].to_csv(a.out, sep="\t", header=False, index=False)
    print(f"Panel: {bed['gene'].nunique()} genes, {len(bed)} intervals, {int((bed['end'] - bed['start']).sum())} bp")


def coverage(a):
    panel = pd.read_csv(a.panel, sep="\t", header=None, names=["chrom", "start", "end", "gene"])
    pieces = pb.per_base_depth_in_panel(pb.read_mosdepth_perbase(a.perbase), panel)
    summ = pb.coverage_summary(pieces)
    per_gene = (pieces.assign(w=pieces["length"] * pieces["depth"]).groupby("gene")
                .agg(bases=("length", "sum"), w=("w", "sum")))
    per_gene["mean_depth"] = per_gene["w"] / per_gene["bases"]
    per_gene.drop(columns="w").to_csv(a.out_prefix + "_per_gene.csv")
    json.dump(summ, open(a.out_prefix + "_summary.json", "w"), indent=2)
    print(json.dumps(summ, indent=2))


def compare(a):
    bed = pd.read_csv(a.bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    res = pb.compare_calls(pb.read_vcf_variants(a.truth, bed), pb.read_vcf_variants(a.query, bed))
    res.to_csv(a.out, index=False)
    print(res.to_string(index=False))


def collect(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = []
    for frac in a.fractions:
        d = os.path.join(a.results, f"frac_{frac}")
        cov = json.load(open(os.path.join(d, "coverage_summary.json")))
        happy = os.path.join(d, "happy.summary.csv")
        met = pb.parse_happy_summary(happy) if os.path.exists(happy) else pd.read_csv(os.path.join(d, "compare.csv"))
        for _, m in met.iterrows():
            rows.append({"fraction": float(frac), "mean_depth": cov["mean_depth"],
                         "uniformity": cov["uniformity_pct_ge_0.2x_mean"], **m.to_dict()})
    df = pd.DataFrame(rows).sort_values(["type", "mean_depth"])
    df.to_csv(os.path.join(a.results, "titration_metrics.csv"), index=False)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for t, g in df.groupby("type"):
        ax.plot(g["mean_depth"], g["recall"], marker="o", label=f"{t} recall")
        ax.plot(g["mean_depth"], g["precision"], marker="s", linestyle="--", label=f"{t} precision")
    ax.set_xlabel("Mean depth over panel"); ax.set_ylabel("Value"); ax.set_ylim(0, 1.02)
    ax.set_title("Variant calling performance vs. depth"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(a.results, "sensitivity_vs_depth.png"), dpi=200)

    full = os.path.join(a.results, "frac_1.0", "coverage_per_gene.csv")
    if os.path.exists(full):
        pg = pd.read_csv(full).sort_values("mean_depth")
        fig, ax = plt.subplots(figsize=(5.5, 8))
        ax.barh(pg["gene"], pg["mean_depth"]); ax.set_xlabel("Mean depth")
        ax.set_title("Coverage by gene (full depth)"); ax.tick_params(axis="y", labelsize=7)
        fig.tight_layout(); fig.savefig(os.path.join(a.results, "coverage_by_gene.png"), dpi=200)
    print(df.to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("make-panel"); p.add_argument("--gtf", required=True); p.add_argument("--genes", default="config/genes.txt")
    p.add_argument("--pad", type=int, default=10); p.add_argument("--out", default="data/panel.bed"); p.set_defaults(f=make_panel)
    p = sub.add_parser("coverage"); p.add_argument("--perbase", required=True); p.add_argument("--panel", default="data/panel.bed")
    p.add_argument("--out-prefix", required=True); p.set_defaults(f=coverage)
    p = sub.add_parser("compare"); p.add_argument("--truth", required=True); p.add_argument("--query", required=True)
    p.add_argument("--bed", required=True); p.add_argument("--out", required=True); p.set_defaults(f=compare)
    p = sub.add_parser("collect"); p.add_argument("--results", default="results")
    p.add_argument("--fractions", nargs="+", default=["1.0", "0.5", "0.25", "0.1"]); p.set_defaults(f=collect)
    a = ap.parse_args(); a.f(a)


if __name__ == "__main__":
    main()
