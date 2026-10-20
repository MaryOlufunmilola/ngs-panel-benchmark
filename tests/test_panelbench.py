import gzip
import math

import pandas as pd
import pytest

import panelbench as pb


def test_merge_intervals_merges_overlaps_per_gene():
    df = pd.DataFrame({"chrom": ["chr1"] * 3, "start": [0, 50, 200], "end": [100, 150, 300], "gene": ["A"] * 3})
    out = pb.merge_intervals(df)
    assert list(zip(out["start"], out["end"])) == [(0, 150), (200, 300)]


def test_gtf_to_panel_bed_pads_converts_and_reports_missing(tmp_path):
    gtf = tmp_path / "x.gtf.gz"
    lines = [
        'chr17\tHAVANA\tCDS\t101\t200\t.\t-\t0\tgene_id "g1"; gene_name "TP53";',
        'chr17\tHAVANA\tCDS\t205\t300\t.\t-\t0\tgene_id "g1"; gene_name "TP53";',
        'chr17\tHAVANA\texon\t1\t1000\t.\t-\t.\tgene_id "g1"; gene_name "TP53";',
        'chr12\tHAVANA\tCDS\t500\t600\t.\t-\t0\tgene_id "g2"; gene_name "OTHER";',
    ]
    with gzip.open(gtf, "wt") as f:
        f.write("\n".join(lines) + "\n")

    bed = pb.gtf_to_panel_bed(str(gtf), ["TP53", "KRAS"], pad=10)

    # 101-200 -> [90, 210), 205-300 -> [194, 310); overlap after padding merges them
    assert list(zip(bed["chrom"], bed["start"], bed["end"], bed["gene"])) == [("chr17", 90, 310, "TP53")]
    assert bed.attrs["missing"] == ["KRAS"]
    # BED column order: chrom, start, end, name (samtools/bedtools/mosdepth read columns by position)
    assert list(bed.columns) == ["chrom", "start", "end", "gene"]


def test_coverage_counts_gaps_as_zero_and_computes_uniformity():
    panel = pd.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200], "gene": ["A"]})
    perbase = pd.DataFrame({"chrom": ["chr1", "chr1"], "start": [0, 150], "end": [150, 175], "depth": [100, 10]})
    pieces = pb.per_base_depth_in_panel(perbase, panel)
    assert pieces["length"].sum() == 100            # 50 @100x, 25 @10x, 25 gap @0x
    s = pb.coverage_summary(pieces, thresholds=(20, 50))
    assert s["mean_depth"] == pytest.approx((50 * 100 + 25 * 10) / 100)
    assert s["pct_ge_20x"] == pytest.approx(50.0)
    # 0.2 x mean = 10.5 -> only the 100x bases qualify
    assert s["uniformity_pct_ge_0.2x_mean"] == pytest.approx(50.0)


def _vcf(path, records):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c, p, r, a in records:
            f.write(f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t.\n")


def test_compare_calls_by_type_within_bed(tmp_path):
    truth = tmp_path / "t.vcf"; query = tmp_path / "q.vcf"
    _vcf(truth, [("chr1", 101, "A", "G"), ("chr1", 120, "C", "T"), ("chr1", 150, "AT", "A"), ("chr1", 900, "G", "C")])
    _vcf(query, [("chr1", 101, "A", "G"), ("chr1", 130, "G", "A"), ("chr1", 150, "AT", "A")])
    bed = pd.DataFrame({"chrom": ["chr1"], "start": [100], "end": [200]})

    res = pb.compare_calls(pb.read_vcf_variants(truth, bed), pb.read_vcf_variants(query, bed)).set_index("type")

    assert res.loc["SNP", "tp"] == 1 and res.loc["SNP", "fn"] == 1 and res.loc["SNP", "fp"] == 1
    assert res.loc["SNP", "recall"] == pytest.approx(0.5) and res.loc["SNP", "precision"] == pytest.approx(0.5)
    assert res.loc["INDEL", "recall"] == pytest.approx(1.0)
    assert res.loc["SNP", "truth_total"] == 2      # chr1:900 is outside the BED


def test_parse_happy_summary_keeps_pass_rows(tmp_path):
    p = tmp_path / "h.summary.csv"
    pd.DataFrame({
        "Type": ["INDEL", "INDEL", "SNP", "SNP"], "Filter": ["ALL", "PASS", "ALL", "PASS"],
        "TRUTH.TOTAL": [10, 10, 100, 100], "TRUTH.TP": [8, 8, 99, 98], "TRUTH.FN": [2, 2, 1, 2],
        "QUERY.FP": [1, 1, 3, 2], "METRIC.Recall": [0.8, 0.8, 0.99, 0.98],
        "METRIC.Precision": [0.89, 0.89, 0.97, 0.98], "METRIC.F1_Score": [0.84, 0.84, 0.98, 0.98],
    }).to_csv(p, index=False)
    res = pb.parse_happy_summary(p)
    assert list(res["type"]) == ["INDEL", "SNP"]
    assert res.loc[res["type"] == "SNP", "recall"].item() == pytest.approx(0.98)


def test_toy_example_matches_expected_output():
    """End-to-end run of the Python steps on the example inputs (see example/README in the main README)."""
    import pathlib
    import subprocess
    root = pathlib.Path(__file__).resolve().parents[1]
    subprocess.run(["bash", str(root / "example" / "run_example.sh")], check=True, cwd=root)
