# Estimating relative abundance from BAM files

> Goal: estimate genome abundance per sample after read mapping.
>
> This is a **tutorial / recipe**, not a pipeline. Copy and paste only the blocks you need.

---

## Abundance estimation with **CoverM**

This approach uses **CoverM** to compute abundance metrics directly from BAM files. It uses "competitive mapping" to better predict relative abundances. It means that concatenates all the genomes of interest in one reference before doing the mapping. 

### 1) Requirements

* `coverm`
* Metagenomic read files (fastq.gz)
* Reference genomes in FASTA format in a single directory

---

### 2) Input paths

```bash
genomes=../data/genomes/tagged
bams=../analyses/mapping/*bam
outdir=../analyses/coverm
output=$outdir/MAST4.polar.coverage.95id.80rcov.tsv
R1=$(ls -d ../data/tara_samples/GGMM.reads/*.1.clean.fastq.gz)
R2=$(ls -d ../data/tara_samples/GGMM.reads/*.2.clean.fastq.gz)

mkdir -p "$outdir"
```

---

### 3) CoverM command

```bash
coverm genome \
  -1 $R1 -2 $R2 \
  -d $genomes \
  -p bwa \
  -x fasta \
  --min-read-percent-identity 95 \
  --min-read-aligned-percent 80 \
  -t 24 \
  --discard-unmapped \
  -m count rpkm tpm covered_bases \
  --min-covered-fraction 0 \
  -o $output
```

**Output**:

* One TSV file with genome-level metrics per sample (`count`, `RPKM`, `TPM`, `covered_bases`).

