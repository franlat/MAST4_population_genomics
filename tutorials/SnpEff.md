# Variant effect annotation with **SnpEff**

> Goal: annotate SNVs with predicted functional effects using **SnpEff**, starting from a filtered VCF file (with freebayes), and extract synonymous and non-synonymous variants.
>
> This is written as a **tutorial / recipe** (not as a SLURM submission script). Copy and paste only the blocks you need.

---

## 0) Requirements

* `java`
* **SnpEff**
* **SnpSift** (bundled with SnpEff)
* A reference genome properly configured for SnpEff (FASTA + GFF3)

---

## 1) Minimal variables (adapt paths and names)

```bash
# Genome / dataset name (must match the SnpEff database name)
mast=<GENOME>

# Paths to software
snpeff=/home/apps/snpeff/snpEff/snpEff.jar
snpsift=/home/apps/snpeff/snpEff/SnpSift.jar

# Input VCF (from FreeBayes / POGENOM preprocessing)
vcf=/home/projects/SNPs/analyses/pogenom/vcf_files/${mast}/TA.${mast}.fb.noN.vcf

# Output directory
output=/home/projects/SNPs/analyses/snpeff/${mast}

mkdir -p "$output"
```

---

## 2) Prepare the SnpEff database (FASTA + GFF3)

Before running `snpEff build`, the reference genome and its annotation **must be installed inside the SnpEff data directory** using the expected structure.

Inside the SnpEff `data/` directory, create one folder per genome:

```
data/
└── <GENOME>/
    ├── sequences.fa   # reference genome (FASTA)
    └── genes.gff      # annotation (GFF3)
```

The file names **must be exactly** `sequences.fa` and `genes.gff`, regardless of their original names.

Example:

```bash
cd /home/apps/snpeff/snpEff/data

mkdir -p <GENOME>
cp /path/to/<GENOME>.fasta  <GENOME>/sequences.fa
cp /path/to/<GENOME>.gff3   <GENOME>/genes.gff
```

In addition, `<GENOME>` must be declared in `snpEff.config`, for example:

```
<GENOME>.genome : <GENOME description>
```

Only after this preparation step can the database be built.

---

## 3) Build the SnpEff database

Build the SnpEff database from a reference genome and its annotation (GFF3).

```bash
java -Xmx8g -jar "$snpeff" build \
  -gff3 \
  -v "$mast"
```

> This step only needs to be done once per genome, unless the reference or annotation changes.

---

## 4) Annotate variants with SnpEff

Annotate the VCF file and generate an HTML summary report.

```bash
java -Xmx8g -jar "$snpeff" \
  -v \
  -stats "$output/${mast}.html" \
  "$mast" "$vcf" \
  > "$output/${mast}.snpeff.vcf"
```

**Outputs**:

* Annotated VCF: `${mast}.snpeff.vcf`
* HTML summary report: `${mast}.html`

---

## 5) (Optional) Extract non-synonymous variants

Filter variants annotated as **missense**.

```bash
java -jar "$snpsift" filter \
  "ANN[*].EFFECT has 'missense_variant'" \
  "$output/${mast}.snpeff.vcf" \
  > "$output/${mast}.filter_non_synonymous_first.vcf"
```

---

## 6) (Optional) Extract synonymous variants

Filter variants annotated as **synonymous**.

```bash
java -jar "$snpsift" filter \
  "ANN[*].EFFECT has 'synonymous_variant'" \
  "$output/${mast}.snpeff.vcf" \
  > "$output/${mast}.filter_synonymous_first.vcf"
```

---

## Notes

* The database name (`mast`) must match the entry defined in `snpEff.config`.
* The input VCF is assumed to be already filtered (e.g. masked nucleotides removed).
* Only simple effect categories are extracted here (missense vs synonymous).
