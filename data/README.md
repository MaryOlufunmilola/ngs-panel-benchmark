# Data

All inputs are public; none are committed. Paths below match `config/settings.sh`.

| File | Source |
|---|---|
| `ref/GRCh38_no_alt_analysis_set.fasta` (+ `.fai`) | GRCh38 "no alt" analysis set (as used by GIAB and NCBI); index with `samtools faidx` |
| `annotation/gencode.basic.annotation.gtf.gz` | GENCODE human release, "basic" gene annotation (GRCh38) |
| `reads/HG002.GRCh38.bam` (+ index) | A GIAB HG002 (NA24385) Illumina alignment to GRCh38, from the GIAB data portal (ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/). A remote https URL with an index also works; only panel regions are read. |
| `truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and `..._noinconsistent.bed` | GIAB HG002 NIST v4.2.1 benchmark for GRCh38 (release folder for the Ashkenazim trio son) |
| `annotation/clinvar.vcf.gz` (optional) | NCBI ClinVar VCF for GRCh38 |

Exact file names differ between releases; edit `config/settings.sh` to match what you download.
The reference used for the reads, the truth set, and variant calling must all be the same build.
