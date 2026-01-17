# MAST4 Population Genomics

This repository contains the scripts and code used to analyze population genomic variation in MAST-4 species as described in:

> Latorre F., Jaillon O., Sieracki M. E., Cruaud C., Massana R., and Logares R.
>
> **Global population structure of a unicellular marine predator**
> 
> Manuscript under review / preprint available at: https://doi.org/10.1101/2024.11.27.625609

The analyses are based on single-cell–derived reference genomes and surface ocean metagenomes from the _Tara_ Oceans expedition, and focus on genetic diversity, population structure, and signatures of selection across global MAST-4 populations.

---

## Repository scope

This repository is intended to:

* Document the main computational steps used in the study
* Provide transparency on how population genomic metrics were computed
* Allow reuse or adaptation of specific analysis scripts by other researchers

It is **not intended** to be a fully automated or plug-and-play pipeline. Some scripts are tailored to the structure of the datasets used in this study.

---

## Contents

### `tutorials/`

Contains small receips/tutorials following the Methods section in the Manuscript, including:

* [Mapping](./tutorials/Mapping.md) - Mapping of Metagenomic reads against the reference genomes
* [Relative Abundance](./tutorials/Abundance.md) - Computation of Relative Abundance values for the genomes across metagenomic samples
* [Population Genomics analyses](./tutorials/Pogenom.md) - Variant Calling and the computation of population genomics parameters
* [Variant effects](./tutorials/SnpEff.md) - Computation of the effects of the Variants predicted
* [dNdS analysis](./tutorials/dNdS.md) - Computation of the dN/dS ratios for positive selection assessment

### `scripts/`

Contains custom scripts used for population genomic analyses and data processing, including:

* [`dNdS.py`](,/scripts/dNdS.py) - Python script used to compute gene-level dN/dS ratios across stations and populations.
* [`idfilter.pl`](./scripts/idfilter.pl) - Perl script used to filter BAM files by read identity and coverage.
* [`FST_density_histogram.R`](./scripts/FST_density_histogram.R) - R script to process FST data into density histograms (Figure 1).
* [`FST_dend_map.R`](./scripts/FST_dend_map.R) - R script to generate the FST dendrograms to infer populations along with the relative abundance per sample (Figure 2).
* [`dNdS_gclust_funcpies.R`](./scripts/dNdS_gclust_funcpies.R) - R script to plot and compute gene cluster based on dN/dS ratios across samples (Figure 3) and the Pie Charts with eggNOG functional data (Figure S2). 

---

## Data availability

Due to their size, raw metagenomic and genomic datasets are **not included** in this repository, but are available at:

- DNA sequences from _Tara_ Oceans are stored at ENA with the accession numbers PRJEB6603 for the SAGs, and PRJEB4352 for the metagenomes. 
- Genome coassemblies, coding sequence predictions, and amino acid predictions are available in FigShare (DOI: https://doi.org/10.6084/m9.figshare.13072322). 

All other study data are included in the article, supporting information, or this GitHub Repository.

---

## Software and dependencies

Analyses were performed using a combination of custom scripts and external software, including (but not limited to):

* Python
* Perl
* Samtools
* BWA
* CoverM
* FreeBayes
* SnpEff
* POGENOM
* R (for statistical analyses and visualization)

Exact versions and parameters are reported in the Methods section of the manuscript.

---

## Reproducibility notes

Because the analyses rely on:

* Large metagenomic datasets
* Coverage-based filtering
* Dataset-specific preprocessing steps

Full reproduction of the results requires access to the original data and following the workflow described in the manuscript.
Nevertheless, the scripts provided here capture the **core analytical logic** used in the study.

---

## Contact

For questions regarding the code or analyses, please contact:

**Francisco Latorre**
📧 [latorre@icm.csic.es]


