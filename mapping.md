# Recipe: mapping paired-end reads against MAST-4 genomes with **BWA** + **samtools**

> Goal: map paired-end `*.fastq.gz` reads against a `MAST4X.fasta` genome, generate a **sorted and indexed BAM**, and (optionally) apply a **95% identity filter** using `idfilter.pl`.
>
> This is written as a **tutorial / recipe**. Copy and paste only the blocks you need.

---

## 0) Requirements

* `bwa` (genome indexing and mapping)
* `samtools` (SAM/BAM manipulation)
* `perl` + `idfilter.pl` (only if identity filtering is required)

On HPC systems these are often loaded via `module load bwa samtools perl`, but here we assume they are already available in your `PATH`.

---

## 1) Minimal variables (adapt paths and names)

Define the following:

* `GENOME`: target genome FASTA
* `R1` and `R2`: paired-end reads (gzipped)
* `THREADS`: number of threads
* `OUTDIR`: output directory
* `PREFIX`: base prefix for output files

Example:

```bash
GENOME=../data/genomes/MAST4A.fasta
R1=../data/tara_samples/GGMM.reads/TARA11.SUR.GGMM.1.clean.fastq.gz
R2=../data/tara_samples/GGMM.reads/TARA11.SUR.GGMM.2.clean.fastq.gz
THREADS=24
OUTDIR=../analyses/mapping/MAST4A
PREFIX=TA11.SUR.GGMM.MAST4A

mkdir -p "$OUTDIR"
```

> Tip: for MAST-4A/B/C/E, only change the genome letter and output directory.

---

## 2) Index the genome (once per FASTA)

This step only needs to be done **once** for each `MAST4X.fasta` file.

```bash
bwa index "$GENOME"
```

This creates the `.amb`, `.ann`, `.bwt`, `.pac`, and `.sa` index files next to the FASTA.

---

## 3) Mapping with `bwa mem` and BAM generation

### 3.1) Standard pipeline: `bwa mem` → BAM → sort → index

This block:

* decompresses `fastq.gz` on the fly
* converts SAM to BAM
* sorts alignments by coordinate
* indexes the BAM

```bash
bwa mem -t "$THREADS" "$GENOME" <(zcat "$R1") <(zcat "$R2") \
  | samtools view -bh - \
  | samtools sort -@ "$THREADS" -o "$OUTDIR/${PREFIX}.sorted.bam" -

samtools index "$OUTDIR/${PREFIX}.sorted.bam"
```

**Output**:

* `.../${PREFIX}.sorted.bam`
* `.../${PREFIX}.sorted.bam.bai`

---

## 4) 95% identity filtering with `idfilter.pl`

If your pipeline uses `idfilter.pl` to retain alignments with ≥95% identity and ≥80% read coverage, a typical workflow is:

```bash
IDFILTER=./idfilter.pl
IDENTITY=0.95

FILTERED_DIR="$OUTDIR/filtered_95id"
mkdir -p "$FILTERED_DIR"

perl "$IDFILTER" \
  -i "$OUTDIR/${PREFIX}.sorted.bam" \
  -o "$FILTERED_DIR/${PREFIX}.95id.sorted.bam" \
  -d "$IDENTITY" \

samtools index "$FILTERED_DIR/${PREFIX}.95id.sorted.bam"
```

**Output**:

* `.../filtered_95id/${PREFIX}.95id.sorted.bam`
* `.../filtered_95id/${PREFIX}.95id.sorted.bam.bai`

> Note: if `idfilter.pl` supports an explicit *coverage* parameter (-l, default 0.80), include it here.
>
> This is a custom script. Nowadays, there are software available to do this step that is actually supported by the developers, such as [bam-filter](https://github.com/genomewalker/bam-filter)
---

## 5) (Optional) Cleanup of intermediate files

If you only need the filtered BAM, you can remove the unfiltered BAM:

```bash
rm -f "$OUTDIR/${PREFIX}.sorted.bam" "$OUTDIR/${PREFIX}.sorted.bam.bai"
```
