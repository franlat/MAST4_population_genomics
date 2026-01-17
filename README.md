# MAST4 Population Genomics

This repository contains the scripts and code used to analyze population genomic variation in MAST-4 species as described in:

> Latorre F., Jaillon O., Sieracki M. E., Cruaud C., Massana R., and Logares R.
> **Global population structure of a unicellular marine predator**
> (manuscript under review / preprint available at: https://doi.org/10.1101/2024.11.27.625609)

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

### `scripts/`

Contains custom scripts used for population genomic analyses, including:

* `dNdS.py`
  Python script used to compute gene-level dN/dS ratios across stations and populations.

  ⚠️ **Important notes**:

  * The script is tailored to this study and assumes:

    * A fixed number of samples (82 surface metagenomes)
    * Specific GFF formatting
    * Input files generated following the workflow described in the manuscript
  * The script **must be adapted** for use with other datasets or experimental designs.

Additional scripts for data processing and figure generation may be added as the repository evolves.

---

## Data availability

Due to their size, raw metagenomic and genomic datasets are **not included** in this repository, but are available at:

- DNA sequences from _Tara_ Oceans are stored at ENA with the accession numbers PRJEB6603 for the SAGs, and PRJEB4352 for the metagenomes. 
- Genome coassemblies, coding sequence predictions, and amino acid predictions are available in FigShare (DOI: https://doi.org/10.6084/m9.figshare.13072322). 

All other study data are included in the article, supporting information, or this GitHub Repository.

---

## Software and dependencies

Analyses were performed using a combination of custom scripts and external software, including (but not limited to):

* Python (custom scripts)
* Samtools
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

---

## License

[Optional – add a license if desired]

