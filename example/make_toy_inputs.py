"""
Write tiny inputs, small enough to check by hand, for the Python steps of the benchmark
(make-panel, coverage, compare, collect), standing in for the outputs of
mosdepth and bcftools at four read depths.

Planted design:
  - GTF: 3 panel genes (GENE_A, GENE_B, GENE_C) with CDS exons, plus a UTR row
    and a gene outside the panel that must be ignored. GENE_MISSING is listed in
    genes.txt but absent from the GTF, so make-panel should report it.
  - Depth at full coverage: GENE_A 200x, GENE_B 100x, GENE_C 40x, with a 20 bp
    gap with zero coverage in GENE_B. Lower fractions scale depth down.
  - Truth: 10 SNVs + 4 indels in the panel, 1 SNV outside it (must not count).
  - Calls: all found at full depth; progressively more missed, and a few
    false positives added, as depth drops.
"""

import gzip
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "toy_inputs")
FRACTIONS = ["1.0", "0.5", "0.25", "0.1"]

# chrom, gene, CDS exons (1-based inclusive, as in a GTF)
GENES = [
    ("chr1", "GENE_A", [(1001, 1100), (1201, 1300)]),
    ("chr1", "GENE_B", [(5001, 5200)]),
    ("chr2", "GENE_C", [(2001, 2150)]),
]
FULL_DEPTH = {"GENE_A": 200, "GENE_B": 100, "GENE_C": 40}

TRUTH_SNVS = [("chr1", 1010, "A", "G"), ("chr1", 1050, "C", "T"), ("chr1", 1090, "G", "A"),
              ("chr1", 1220, "T", "C"), ("chr1", 1280, "A", "C"), ("chr1", 5020, "G", "T"),
              ("chr1", 5100, "C", "A"), ("chr1", 5180, "T", "G"), ("chr2", 2030, "A", "T"),
              ("chr2", 2120, "G", "C")]
TRUTH_INDELS = [("chr1", 1060, "A", "AT"), ("chr1", 5150, "GC", "G"),
                ("chr2", 2050, "T", "TAA"), ("chr2", 2100, "CAG", "C")]
OFF_PANEL = [("chr1", 9000, "A", "G")]

# Per fraction: truth variants missed (by index into SNVs / indels) and false positives added.
MISSED = {"1.0": ([], []), "0.5": ([9], []), "0.25": ([9, 8, 4], [3]), "0.1": ([9, 8, 4, 7, 6, 2], [3, 2])}
FALSE_POS = {"1.0": [], "0.5": [], "0.25": [("chr1", 5060, "C", "G")],
             "0.1": [("chr1", 5060, "C", "G"), ("chr2", 2010, "A", "G")]}


def write_gtf(path):
    lines = ["##description: toy GTF for the example\n"]
    for chrom, gene, exons in GENES:
        attrs = f'gene_id "{gene}.1"; gene_name "{gene}";'
        lines.append(f"{chrom}\ttoy\tgene\t{exons[0][0] - 50}\t{exons[-1][1] + 50}\t.\t+\t.\t{attrs}\n")
        lines.append(f"{chrom}\ttoy\tUTR\t{exons[0][0] - 50}\t{exons[0][0] - 1}\t.\t+\t.\t{attrs}\n")
        for s, e in exons:
            lines.append(f"{chrom}\ttoy\tCDS\t{s}\t{e}\t.\t+\t0\t{attrs}\n")
    lines.append('chr3\ttoy\tCDS\t100\t400\t.\t+\t0\tgene_id "OFF.1"; gene_name "OFF_PANEL_GENE";\n')
    with gzip.open(path, "wt") as f:
        f.writelines(lines)


def write_perbase(path, frac):
    """Per base depth runs in mosdepth format (0-based, end exclusive), covering the padded panel."""
    runs = []
    for chrom, gene, exons in GENES:
        depth = round(FULL_DEPTH[gene] * float(frac))
        start, end = exons[0][0] - 1 - 10, exons[-1][1] + 10
        if gene == "GENE_B":   # 20 bp dropout in the middle of GENE_B
            runs += [(chrom, start, 5100, depth), (chrom, 5100, 5120, 0), (chrom, 5120, end, depth)]
        else:
            runs.append((chrom, start, end, depth))
    with gzip.open(path, "wt") as f:
        for r in sorted(runs):
            f.write("\t".join(map(str, r)) + "\n")


def write_vcf(path, variants):
    with gzip.open(path, "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, r, a in sorted(variants):
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    write_gtf(os.path.join(OUT, "toy.gtf.gz"))
    with open(os.path.join(OUT, "genes.txt"), "w") as f:
        f.write("GENE_A\nGENE_B\nGENE_C\nGENE_MISSING\n")
    write_vcf(os.path.join(OUT, "truth.vcf.gz"), TRUTH_SNVS + TRUTH_INDELS + OFF_PANEL)
    for frac in FRACTIONS:
        d = os.path.join(OUT, f"frac_{frac}")
        os.makedirs(d, exist_ok=True)
        write_perbase(os.path.join(d, "cov.per-base.bed.gz"), frac)
        miss_snv, miss_indel = MISSED[frac]
        calls = ([v for i, v in enumerate(TRUTH_SNVS) if i not in miss_snv]
                 + [v for i, v in enumerate(TRUTH_INDELS) if i not in miss_indel]
                 + FALSE_POS[frac] + OFF_PANEL)
        write_vcf(os.path.join(d, "calls.vcf.gz"), calls)
    print(f"Toy inputs written to {OUT}")


if __name__ == "__main__":
    main()
