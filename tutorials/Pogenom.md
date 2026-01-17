# Population genomics with **FreeBayes** + VCF filtering + **POGENOM**

> Goal: starting from sample BAMs mapped to a given genome (e.g. MAST-4A/B/C/E), generate a **VCF** with **FreeBayes**, filter variants overlapping masked nucleotides (`N` in the reference), optionally compute quick SNV stats, and finally run **POGENOM**.
>
> This is written as a **tutorial / recipe** (not as a SLURM submission script). Copy and paste only the blocks you need.

---

## 0) Requirements

* `samtools`
* `freebayes`
* `bcftools`
* `perl` + `pogenom.pl`

---

## 1) Minimal variables (adapt paths)

```bash
# Inputs
inbam=../analyses/mapping/filtered_95id/<GENOME>          # directory containing BAMs for one genome
ref=../data/genomes/tagged/<GENOME>.fasta                 # reference FASTA

# Outputs
outroot=../analyses/population_genomics_filtered/<GENOME>
outbam=$outroot/bam
outvcf=$outroot/vcf_file

mkdir -p "$outbam" "$outvcf"

# POGENOM
pogenom=/home/apps/POGENOM/pogenom.pl
genetic_code=/home/apps/POGENOM/standard_genetic_code.txt
```

> Replace `<GENOME>` with your target genome name (e.g. `MAST4A`).

---

## 2) Merge and sort BAM files

Merge all BAMs for this genome into a single BAM and sort it.

```bash
# Merge BAMs
samtools merge --threads 24 \
  "$outbam/TARA.<GENOME>.95id.merged.bam" \
  $(ls "$inbam"/*.addRG.bam | grep -Ff <GENOME>.filt | tr '\n' ' ')

# Sort
samtools sort --threads 24 \
  -o "$outbam/TARA.<GENOME>.95id.merged.sorted.bam" \
  "$outbam/TARA.<GENOME>.95id.merged.bam"

# Cleanup
rm -f "$outbam/TARA.<GENOME>.95id.merged.bam"
```

> Note: the `grep -Ff <GENOME>.filt` part is how BAMs were selected in this study. If you already have the exact list of BAMs, you can replace the whole `$(...)` with explicit file paths.

---

## 3) Call variants with FreeBayes

```bash
# Increase stack size (helpful on some HPC systems)
ulimit -s 81920

# Variant calling
# Notes:
#   -p 1 : haploid ploidy (adjust if needed)
#   -C 4 : minimum number of observations supporting the alternate allele (adjust if needed)
freebayes \
  -p 1 \
  -f "$ref" \
  -v "$outvcf/TARA.<GENOME>.95id.fb.vcf" \
  -C 4 \
  "$outbam/TARA.<GENOME>.95id.merged.sorted.bam"
```

---

## 4) Filter VCF to remove variants overlapping masked nucleotides (N)

Remove VCF entries where the **reference allele** contains `N` (masked positions in the reference). If your reference genome does not include them, you can skip this step. 

```bash
input_vcf="$outvcf/TARA.<GENOME>.95id.fb.vcf"
output_vcf="$outvcf/TARA.<GENOME>.95id.fb.noN.vcf"

# Copy header
grep "^#" "$input_vcf" > "$output_vcf"

# Keep only variant lines where REF does not contain N
grep -v "^#" "$input_vcf" \
  | awk 'BEGIN{FS="[[:space:]]+"}{ if ($4 !~ /N/) print $0 }' \
  >> "$output_vcf"
```

---

## 5) (Optional) SNV summary statistics with bcftools

```bash
bcftools stats --samples '-' "$outvcf/TARA.<GENOME>.95id.fb.noN.vcf" \
  > "$outvcf/TARA.<GENOME>.95id.fb.stats"

# Extract only SNV-related summary lines (PSC)
egrep "PSC" "$outvcf/TARA.<GENOME>.95id.fb.stats" | egrep -v "#" \
  > "$outvcf/TARA.<GENOME>.95id.fb.stats.psc"
```

---

## 6) Run POGENOM

POGENOM is a software designed for Prokaryotes. It requires GFF withot introns to compute all the available parameters it offers. However, if it is not given as input, it only computes them at the whole-genome, which is what we needed. 

```bash
perl "$pogenom" \
  --min_count 10 \
  --min_found 4 \
  --vcf_file "$outvcf/TARA.<GENOME>.95id.fb.noN.vcf" \
  --out "$outvcf/TARA.<GENOME>.95id.c10.s4" \
  --genetic_code_file "$genetic_code" \
  --fasta_file "$ref"

echo "POGENOM finished — outputs are in: $outvcf"
```

---

## Notes (study-specific)

* BAMs were assumed to be read-group annotated (`*.addRG.bam`) so Freebayes can differentiate between samples. This can be done with samtools as follows:
`samtools addreplacerg --threads 24 -r "ID:${sample}" -r "LB:lib1" -r "PL:illumina" -r "PU:unit1" -r "SM:${sample}" -o ${filtered95}/${sample}.${name}.95id.addRG.bam ${filtered95}/${sample}.${name}.95id.
sorted.bam`
* Filtering masked nucleotides was performed by removing VCF records where the REF allele contained `N`.
* FreeBayes parameters were chosen for haploid variant calling (`-p 1`) with a minimum of 4 supporting observations (`-C 4`).
* POGENOM was run with `--min_count 10` and `--min_found 4`.
