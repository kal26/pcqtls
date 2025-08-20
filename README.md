# Focus on single gene effects limits discovery and interpretation of complex trait-associated variants


## Overview

This repository contains the code associated with our work developing a multi-gene eQTL mapping framework, termed cis-principal component expression QTL (cis-pc eQTL or pcQTL). pcQTL leverage allelic expression "proxitropy" - the phenomenon by which one variant changes the expression of mulptile, nearby genes - to map QTL effects missed by a standard single-gene eQTL approach. 

<div style="text-align: center;">
    <img src="images/pcqtl_method.png" alt="pcQTL vs eQTL methods" style="max-width: 100%; height: auto;">
</div>

### Preprint

Our preprint is available [here](https://www.biorxiv.org/content/10.1101/2025.06.06.658175v1)

## Usage

This repository contains a Snakemake workflow for performing pcQTL analysis. To run the complete workflow:

0. Create a conda environment according to [Environment Setup](#environment-setup)

1. Download the necessary data as outlined in [Required Datasets](#required-datasets)

2. Edit the example configuration file `config/config_example.yaml` to specify the paths to your input data and desired output directories.

3. Run the workflow with your configuration:

```bash
snakemake --configfile config/config_example.yaml --cores 10 --use-conda
```

### Repository Structure

* `workflow/rules`: Snakemake workflow for the pcQTL mapping framework.
* `workflow/scripts`: Python and R scripts used by the workflow.
* `workflow/figures`: Jupyter notebooks for data analysis and visualization and figure files.
* `config`: Example configuration file for the workflow.
* `references`: Small file-size references.
* `Snakefile`: Main Snakemake workflow file that orchestrates the pcQTL analysis pipeline.

### Environment Setup

This workflow can be run from a specific conda environment with all dependencies installed.

```bash
# Clone this repository
git clone <repository-url>
cd <repository-name>

# Create the conda environment from the exported file
conda env create -f environment.yml

# Activate the environment
conda activate pcqtl
```

### Workflow Components

The workflow performs the following analyses:

- **Gene clustering**: Identifies co-expressed gene clusters
- **eQTL analysis**: Maps expression quantitative trait loci for individual genes
- **pcQTL analysis**: Maps QTLs for principal components derived from gene clusters
- **Functional annotation**: Annotates variants and clusters with functional information
- **Co-localization**: Performs co-localization analysis with GWAS summary statistics


## Required Datasets

The workflow requires several input datasets. Below is a comprehensive list of required files and their sources:

### Core Analysis Data

| Dataset | Description | Source | Expected Location in Config |
|---------|-------------|---------|---------------------------|
| GTEx v8 Expression Data | Normalized gene expression by tissue | [GTEx Portal](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_expression_matrices.tar) | `expression_dir` |
| GTEx v8 Genotypes | Genotype data: .bed, .bim, and .fam | Protected access on dbGaP, access available via request at this [link](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v10.p2) | `genotype_stem` |
| GTEx v8 Covariates | Technical and biological covariates | [GTEx Portal](https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_covariates.tar.gz) | `covariates_dir` |
| GTEx v8 Sample Sizes | Sample sizes for each tissue | `references/gtex_sample_sizes.txt` | `gtex_meta` |
| Tissue IDs | List of tissues to analyze | User generated, see `references/selected_tissue_ids.txt` as an example | `tissue_id_path` |
| Chromosome List | Chromosomes to analyze | `references/chrs.txt` | `chr_list_path` |
| GWAS Metadata | Metadata from GWAS studies | [Zenodo](https://zenodo.org/records/3629742#.Y9rTQOzMIUF) | `gwas_meta` |
| GWAS Summary Stats | GWAS summary statistics files | [Zenodo](https://zenodo.org/records/3629742#.Y9rTQOzMIUF) | `gwas_folder` |

### Annotation Data

| Dataset | Description | Source | Expected Location in Config |
|---------|-------------|---------|---------------------------|
| GENCODE v26 | Gene annotations | [GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz) | `gencode_path` |
| ABC Enhancer Predictions | Activity-by-Contact enhancer predictions (hg19, requires liftOver to hg38) | [Engreitz Lab](https://www.engreitzlab.org/resources) | `full_abc_path` |
| CTCF ID-tissue | Matching between tissues and CTCF experiments | `references/ctcf_matched_gtex.txt` | `ctcf_match_path` |
| CTCF Peaks | CTCF binding data from ENTEx, the experiments from `references/ctcf_matched_gtex.txt` | [ENTEx](https://www.encodeproject.org/entex-matrix/?type=Experiment&status=released&internal_tags=ENTEx) | `ctcf_dir` |
| Cross-mappability | Cross-mappability from Saha and Battle (2019) | [figshare](https://figshare.com/ndownloader/files/13514741) | `cross_map_path` |
| Paralog Relationships | Gene paralog information | [Ensembl Biomart](https://jan2020.archive.ensembl.org/biomart/martview/853dbe49995a1a2712b77655b242db21) | `paralog_path` |
| Gene Ontology Terms | GO term annotations | [Ensembl Biomart](https://jan2020.archive.ensembl.org/biomart/martview/853dbe49995a1a2712b77655b242db21) | `go_path` |
| TAD Boundaries | Topologically Associating Domain boundaries (hg19, requires liftOver to hg38) | [TADKB](http://dna.cs.miami.edu/TADKB/download/TAD_annotations.tar.gz) | `tad_path` |


### Data Preparation Notes

1. **GENCODE Processing**: The GENCODE annotation must be processed to include only "gene" level features with columns `(chr,start,end,strand,gene_id,gene_name,tss,alternative_tss)`. The `tss`column should contain a list of the transcription start site position for all "basic" tagged transcripts: use the `start` coordinate for positive-stranded genes and the `end` coordinate for negative-stranded genes. For this analysis, only protein-coding genes were considered.

2. **Genome Assembly Conversion**: The ABC enhancer predictions and TADKB boundary databases are provided in hg19 coordinates. These must be converted to hg38 using liftOver before use in the workflow.  


## Output Files

The workflow generates results in the following directories (as specified in your config):
- `clusters_dir`: Gene expression clusters
- `eqtl_output_dir`: eQTL analysis results
- `pcqtl_output_dir`: pcQTL analysis results
- `annotations_output_dir`: Functional annotations
- `coloc_output_dir`: Co-localization results

### Output File Details

For detailed descriptions of the output file formats and their contents, refer to [output.md](output.md).

### Results Availability

* Clusters of neighboring correlated genes 
* Summary stats for pcQTL mapping [ZENODO](https://doi.org/10.5281/zenodo.15605351)


## Acknowledgments
We thank the donors and their families for their generous gifts of biospecimens to the GTEx research project. The Genotype-Tissue Expression (GTEx) project was supported by the Common Fund of the Office of the Director of the National Institutes of Health (http://commonfund.nih.gov/GTEx). Additional funds were provided by the National Cancer Institute (NCI), National Human Genome Research Institute (NHGRI), National Heart, Lung, and Blood Institute (NHLBI), National Institute on Drug Abuse (NIDA), National Institute of Mental Health (NIMH), and National Institute of Neurological Disorders and Stroke (NINDS). This research was supported by National Institutes of Health grants R01MH12524, U01AG072573, U01HG012069 to S.B.M. K.L. is supported by the Stanford Genome Training Program (SGTP; NIH/NHGRI T32HG000044). T.G. is supported by the Knight-Hennessy Scholars fellowship


<div style="text-align: center;">
    <img src="images/pcqtl_proxitropy.png" alt="make allelic proxitropy work for you" style="max-width: 100%; height: auto;">
</div>

