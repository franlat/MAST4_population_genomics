# Recipe: dN/dS analysis using `dNdS.py`

> Goal: compute dN/dS statistics from a SnpEff-annotated VCF file, using a custom script (`dNdS.py`) together with the corresponding genome FASTA and GFF annotation.
>
> This is written as a **tutorial / recipe** (not as a SLURM submission script). Copy and paste only the blocks you need.

---

## 0) Requirements

* Python (as used in this study: `python/3.8.5`)
* Custom script: `dNdS.py`
* Inputs:

  * Genome annotation (`.gff`)
  * Genome sequence (`.fasta`)
  * SnpEff-annotated VCF (`.snpeff.vcf`)

---

## 1) Minimal variables (adapt paths and names)

```bash
# Genome / dataset name
genome=<GENOME>

# Input files
gff=/home/projects/SNPs/data/gff_files/${genome}.gff
fasta=/home/projects/SNPs/data/genomes/${genome}.fasta
vcf=/home/projects/SNPs/analyses/snpeff/${genome}/${genome}.snpeff.vcf

# Output file
output=/home/projects/SNPs/analyses/snpeff/${genome}/${genome}.snpeff.dnds.txt
```

---

## 2) Run the dN/dS script

```bash
python dNdS.py \
  -gff "$gff" \
  -f "$fasta" \
  -vcf "$vcf" \
  -o "$output"
```

**Output**:

* A text file with dN/dS results: `${genome}.snpeff.dnds.txt`

---

## Notes

* The input VCF is expected to be the output of the SnpEff annotation step (`*.snpeff.vcf`).
* The GFF and FASTA must correspond to the same reference genome used for variant calling and SnpEff database building.
* `dNdS.py` is a custom script used in this study; include it in the repository (e.g. under `scripts/`) to ensure reproducibility.
* If you working with Prokaryotic genomes, then, POGENOM already gives pNpS estimates per gene. However, POGENOM does not provide them for Eukaryotic genomes. 

**Important notes**:

* The script is tailored to this study and assumes:

  * A fixed number of samples (82 surface metagenomes)
  * Specific GFF formatting
  * Input files generated following the workflow described in the manuscript
* The script **must be adapted** for use with other datasets or experimental designs.
