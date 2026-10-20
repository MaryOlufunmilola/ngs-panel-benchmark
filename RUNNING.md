# Running the NGS Panel Benchmark: Step-by-Step

A walkthrough for running `ngs-panel-benchmark` on a Mac and publishing the results on GitHub.

| Part | Time |
|---|---|
| 1. Put the updated code on GitHub | 15 minutes |
| 2. Set up the environment | 15–30 minutes |
| 3. Download the inputs | 20–40 minutes (about 4 GB) |
| 4. Run the benchmark | 30–90 minutes (mostly streaming reads) |
| 5. Check the results | 10 minutes |
| 6. Publish the results | 15 minutes |

You need about 10 GB of free disk space and a stable internet connection.

---

## 1. Put the updated code on GitHub

Do this first, so the toy example and the BED fix are public before the full run.

```bash
cd ~/Projects
git clone https://github.com/MaryOlufunmilola/ngs-panel-benchmark.git
cd ngs-panel-benchmark

# Copy the updated files over your clone (the zip's top folder is ngs-panel-benchmark/)
unzip -o ~/Downloads/ngs-panel-benchmark.zip -d /tmp/npb_update
rsync -av --exclude .git /tmp/npb_update/ngs-panel-benchmark/ ./

git status          # review what changed
git diff            # read the changes, especially src/panelbench.py (BED column order)
```

**What changed in this update:**

- `src/panelbench.py`, `src/cli.py`: the panel BED now has columns in BED order
  (chrom, start, end, gene). Before, interval merging put the gene name first, which breaks
  mosdepth, bcftools and bedtools.
- `example/`: tiny hand-made inputs with planted answers and a saved expected output.
  `bash example/run_example.sh` runs the panel, coverage, comparison, and summary steps in seconds.
- `tests/`: a test for the BED column order and an end-to-end toy-example test.
- `requirements-dev.txt`: adds matplotlib, which the summary step needs (CI would fail without it).
- `config/settings.sh`: points to a public GIAB HG002 BAM and sets the depth titration.
- `run_benchmark.sh`: runs hap.py under amd64 emulation on Apple Silicon, and raises the
  bcftools depth cap so deep samples aren't truncated.
- `.gitignore`: lets the final results table and plots be committed.
- `README.md`, `RUNNING.md`: inputs and outputs, toy example, runtime, and these instructions.

Check the code, then push:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt
pytest                       # 6 tests should pass
bash example/run_example.sh  # ends with "Output matches expected/titration_metrics.csv"
deactivate

git add -A
git commit -m "Fix panel BED column order; add toy example, tests, and run instructions"
git push
```

Open the repo on GitHub, click **Actions**, and confirm the CI run is green.

---

## 2. Set up the environment

If you don't have conda, install Miniforge first (github.com/conda-forge/miniforge).

```bash
cd ~/Projects/ngs-panel-benchmark
conda env create -f environment.yml
conda activate panelbench
samtools --version | head -1 && bcftools --version | head -1 && mosdepth --version && bedtools --version
```

**On an Apple Silicon Mac (M1/M2/M3/M4),** if `conda env create` can't find a package, build the
environment for Intel and let Rosetta run it:

```bash
CONDA_SUBDIR=osx-64 conda env create -f environment.yml
conda activate panelbench
conda config --env --set subdir osx-64
```

**hap.py (optional, recommended):** start Docker Desktop, then

```bash
docker pull --platform linux/amd64 pkrusche/hap.py
```

Without it, the pipeline uses a simpler exact-match comparison automatically.

---

## 3. Download the inputs

```bash
mkdir -p data/ref data/annotation data/truth
```

**Reference genome:** GRCh38 "no alt" analysis set, about 900 MB compressed.

```bash
cd data/ref
curl -L -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GRCh38_no_alt_analysis_set.fasta
samtools faidx GRCh38_no_alt_analysis_set.fasta
cd ../..
```

**Gene annotation:** GENCODE "basic" GTF. Open
`https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/`, go into the newest `release_XX`
folder, and download `gencode.vXX.basic.annotation.gtf.gz`. Then:

```bash
mv ~/Downloads/gencode.v*.basic.annotation.gtf.gz data/annotation/gencode.basic.annotation.gtf.gz
```

Write down the release number (for example, v49); it goes in the README.

**GIAB truth set:** HG002, v4.2.1, GRCh38.

```bash
cd data/truth
B=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38
curl -L -O $B/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
curl -L -O $B/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
curl -L -O $B/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
cd ../..
ls -lh data/truth      # three files, none of them tiny error pages
```

If a file is only a few KB, the folder name changed: open the `GRCh38` folder in a browser and
download the three files by hand.

**ClinVar (optional, about 100 MB):**

```bash
curl -L -o data/annotation/clinvar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
curl -L -o data/annotation/clinvar.vcf.gz.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

**HG002 reads: don't download them.** `config/settings.sh` already points to GIAB's indexed
~300x Illumina BAM for HG002 on GRCh38; the pipeline streams only the panel regions. Test the link
and the chromosome names:

```bash
source config/settings.sh
samtools view -H "$READS" | grep "^@SQ" | head -3
```

You should see `SN:chr1`. If you get an error, or `SN:1` without "chr", choose another HG002
Illumina GRCh38 BAM that has a `.bai` index, from GIAB's index files
(github.com/genome-in-a-bottle/giab_data_indexes, folder `AshkenazimTrio`), and put its https URL in
`READS`. If you change to a ~30–60x file, also change `FRACTIONS` to `"1.0 0.5 0.25 0.1"`.

---

## 4. Run the benchmark

Check `config/settings.sh`. With the steps above, only these matter:

```bash
REF="data/ref/GRCh38_no_alt_analysis_set.fasta"
GTF="data/annotation/gencode.basic.annotation.gtf.gz"
READS="https://ftp-trace.ncbi.nlm.nih.gov/.../HG002.GRCh38.300x.bam"
CLINVAR_VCF="data/annotation/clinvar.vcf.gz"     # or "" to skip
FRACTIONS="1.0 0.2 0.1 0.05 0.02"                # ~300x, 60x, 30x, 15x, 6x
```

Run it, keeping the Mac awake:

```bash
conda activate panelbench
mkdir -p logs
caffeinate -i bash run_benchmark.sh 2>&1 | tee logs/run.log
```

What you'll see:

1. "Panel: 50 genes, ... intervals, ... bp" (any gene missing from the GTF is listed).
2. Reads fetched over the panel (the slow step; it shows no progress, so be patient).
3. "Evaluation regions: ... bp."
4. Coverage, calls, and comparison for each depth fraction.
5. The titration table.
6. "Calls with a ClinVar record: N" (if enabled).

If it stops partway, fix the problem and rerun the same command. Fetched reads are kept in
`data/hg002_panel.bam` and aren't downloaded again.

---

## 5. Check the results

| File | What it shows |
|---|---|
| `results/titration_metrics.csv` | Recall, precision, and F1 for SNVs and indels at each depth |
| `results/sensitivity_vs_depth.png` | How detection changes as depth falls |
| `results/coverage_by_gene.png` | Mean depth per gene at full depth |
| `results/frac_1.0/coverage_summary.json` | Mean depth, % bases at 20x/50x/100x, uniformity |
| `results/clinvar_matches.tsv` | Full-depth calls with a ClinVar record (optional) |

**Sanity checks:**

- **Full-depth SNV recall and precision should be high** (typically above 0.95). Very low values
  usually mean mismatched builds or chromosome names between reads, reference and truth set.
- **Indels usually score lower than SNVs.** That's normal; say so in the README.
- **Recall should fall at the lowest depths** (15x and 6x). If nothing changes, check that
  `mean_depth` differs between fractions in `titration_metrics.csv`.
- **Some genes will have low coverage** in `coverage_by_gene.png` (often GC-rich first exons).
  Name the lowest two or three in the README.

---

## 6. Publish the results

**1. Print the results as a Markdown table:**

```bash
python - <<'EOF'
import pandas as pd
d = pd.read_csv("results/titration_metrics.csv")
cols = ["type", "fraction", "mean_depth", "tp", "fn", "fp", "recall", "precision", "f1"]
print("| " + " | ".join(cols) + " |"); print("|" + "---|" * len(cols))
for _, r in d[cols].iterrows():
    print("| " + " | ".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in r) + " |")
EOF
```

**2. Edit `README.md`:**

- Replace the line `> **Status:** code and tests complete; results will be added after the full run.`
  with `> **Status:** complete; results below are from a full run on GIAB HG002 (GENCODE vXX, run on <date>).`
- Add a `## Results` section after "Pipeline" with:
  - the table you just printed;
  - `![Sensitivity vs depth](results/sensitivity_vs_depth.png)` and
    `![Coverage by gene](results/coverage_by_gene.png)`;
  - two or three sentences of interpretation: full-depth SNV and indel recall and precision,
    the depth at which SNV recall drops below 0.95, and the lowest-coverage genes.
- In "Time and resources", replace the estimates with your actual run time and disk use.

**3. Commit and push** (`.gitignore` already allows these three result files):

```bash
git add README.md results/titration_metrics.csv results/sensitivity_vs_depth.png results/coverage_by_gene.png
git commit -m "Add full HG002 benchmark results"
git push
```

Confirm CI is still green, and that the plots appear on the README page on GitHub.

**4. Update your applications:** in the resume portfolio line, replace "In progress: ... full run
pending" with the real numbers (for example, "SNV recall X at full depth, falling to Y at 15x").
Use only numbers from `titration_metrics.csv`.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| "Missing tool" | `conda activate panelbench` |
| "Missing input" | A path in `config/settings.sh` doesn't match the file name; check with `ls data/*` |
| Remote read access fails or hangs | Rerun `samtools view -H "$READS"`; try another network (some block https streaming) or another GIAB BAM |
| `samtools view` very slow | Normal for 50 genes at ~300x over the internet; leave it running (`caffeinate` keeps the Mac awake) |
| hap.py errors | Check Docker Desktop is running; or remove the image (`docker rmi pkrusche/hap.py`) and rerun, which uses the exact-match fallback |
| Very few variants found | Chromosome names differ between inputs (see the `SN:chr1` check in step 3) |
| CI fails after pushing | Open the failed run under **Actions** and read the last error lines; usually a missing package in `requirements-dev.txt` |
