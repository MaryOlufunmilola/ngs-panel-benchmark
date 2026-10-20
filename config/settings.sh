# Paths and inputs for run_benchmark.sh. Edit to match your downloads (see README "Data").
REF="data/ref/GRCh38_no_alt_analysis_set.fasta"            # indexed (.fai)
GTF="data/annotation/gencode.basic.annotation.gtf.gz"      # GENCODE basic, GRCh38
READS="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam"
                                                           # GIAB HG002 ~300x Illumina BAM (indexed); only panel regions are streamed. A local BAM/CRAM also works.
TRUTH_VCF="data/truth/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="data/truth/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
CLINVAR_VCF=""                                             # optional: data/annotation/clinvar.vcf.gz (GRCh38)
FRACTIONS="1.0 0.2 0.1 0.05 0.02"                         # depth titration (fraction of reads kept); from ~300x: ~300x, 60x, 30x, 15x, 6x
SEED=42
THREADS=4
