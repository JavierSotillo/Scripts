# Scripts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24088985.svg)](https://doi.org/10.5281/zenodo.24088985)

Analysis scripts from the Helminth Unit, Parasitology Reference and Research Laboratory, Centro Nacional de Microbiología, Instituto de Salud Carlos III (Madrid, Spain).

The repository has two parts:

- **`general/`** holds reusable scripts for common tasks. Each script has a *Settings* block at the top where you set the input file and options before running it.
- **`papers/`** holds one folder per article, with the scripts used in that article and a README explaining how to reproduce its results.

## General scripts

| Folder | Script | Description |
|---|---|---|
| `proteomics/` | `volcano_plot.R` | Volcano plot from log2 fold changes and Q-values, with optional highlighting of selected proteins |
| | `pca_plot.R` | PCA of samples from a quantification table |
| | `heatmap.R` | Clustered heatmap (row z-scores) with two sample annotations |
| `enrichment/` | `go_enrichment_mouse_uniprot.R` | GO enrichment of mouse proteins from UniProt accessions |
| | `go_kegg_enrichment_human_symbols.R` | GO and KEGG enrichment of human genes from HGNC symbols |
| | `go_enrichment_interpro2go.R` | GO enrichment for non-model organisms through InterPro domains (e.g. WormBase ParaSite annotations) |
| `sequences/` | `extract_fasta_headers.R` | Extract sequence identifiers from a FASTA file |
| | `fastq_stats.R` | Number of reads and read-length summary of a FASTQ file |
| `statistics/` | `random_lca.R` | Latent class analysis of binary diagnostic test results |

The prevalence and confidence interval calculator is at the root of the repository, [`Prevalence_CI_Calculator.py`](Prevalence_CI_Calculator.py).

The [`Julia/`](Julia) folder contains a PCA and a heatmap of *Anisakis* excretory/secretory products over eight days of culture. Both read the Excel file `Ani_ESPs_Report_Scripts_R.xlsx` from the working directory.

## Scripts by article

| Article | Folder |
|---|---|
| Mutengo et al. (2026) *PLoS Negl Trop Dis* 20(10): e0014013. [doi:10.1371/journal.pntd.0014013](https://doi.org/10.1371/journal.pntd.0014013) | [`papers/Mutengo2026_Shaematobium-Lusaka`](papers/Mutengo2026_Shaematobium-Lusaka) |

## Requirements

R scripts install any missing CRAN or Bioconductor packages when first run. Scripts that map identifiers through Ensembl BioMart or query KEGG need an internet connection. Python scripts need Python 3 only.

## Citation

If you use these scripts, please cite the corresponding article. To cite the repository itself:

> Sotillo J (2026). Scripts: analysis scripts from the Helminth Unit, Instituto de Salud Carlos III (v1.0). Zenodo. https://doi.org/10.5281/zenodo.24088985

Each release on GitHub is archived in Zenodo with its own DOI.

## Licence

Released under the MIT Licence. See [`LICENSE`](LICENSE).
