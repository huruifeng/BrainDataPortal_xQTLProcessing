                     ______________________________

                      NOTES ON DATASET PREPARATION

                           Christopher Zhang
                     ______________________________


Table of Contents
_________________

1. Introduction
2. Data Requirements
3. Data Organization
4. Design Reasoning
5. Scripts


## 1 Introduction

  This document describes the structure of the dataset used for eQTL
  visualization, along with the reasoning behind these choices.

  The visualization has two views:
  - *Gene View*: A gene track is displayed at the bottom showing the
     selected gene in black and nearby significant genes in gray.[1]
     Above it, scatter plots show SNPs associated with each celltype.
  - *SNP View*: A SNP track is displayed at the bottom showing the
     selected SNP in black and nearby significant SNPs in gray.[2] Above
     it, scatter plots show genes associated with each celltype.

  In both views:
  - The *x-axis* represents the genomic position (gene start and end or
    SNP position).
  - The *y-axis* represents −log_{10}(p).
  - Point color reflects the beta value, with red for positive and blue
    for negative. The strongest colors correspond to the largest
    magnitude of beta value; smaller beta values are grayer.


## 2 Data Requirements
  The visualization requires two types of data:

  *QTL Data*:
  - Gene ID
  - SNP ID
  - p-value
  - beta

  *Annotation Data*:
  - Genes: chromosome, start position, end position, strand
  - SNPs: chromosome, position


## 3 Data Organization
  The dataset is organized into the following directory structure:

    dataset/
    |-- celltype_mapping.json
    |-- celltypes
    |   |-- Astrocytes.parquet
    |   |-- ...
    |   `-- Pericytes.parquet
    |-- gene_jsons
    |   |-- A1BG.json
    |   |-- ...
    |   `-- A2M.json
    |-- gene_list.json
    |-- gene_locations
    |   |-- chr1.parquet
    |   |-- ...
    |   `-- chr22.parquet
    |-- snp_jsons
    |   |-- rs12345678.json
    |   |-- ...
    |   `-- rs7288382.json
    |-- snp_list.json
    `-- snp_locations
        |-- chr1.parquet
        |-- ...
        `-- chr22.parquet

  - `celltype_mapping.json`
        Maps display names to filenames in `celltypes/'.
  `celltypes/`
        Per-celltype Parquet files for QTL data (only including gene-SNP
        pairs significant in at least one celltype[3]).
  - `gene_jsons/` and `snp_jsons/`
        Individual JSON annotation files for each gene and SNP,
        including their associated celltypes.
  - `gene_list.json` and `snp_list.json`
        Lists of all genes and SNPs in the dataset.
  - `gene_locations/` and `snp_locations/`
        Per-chromosome Parquet files containing annotation data.


## 4 Design Reasoning
  Several choices were made to prioritize speed:

  Redundant Annotations
        Each gene and SNP has its own JSON file for annotations to allow
        fast access.
  Precomputed Celltypes
        Storing celltype information directly in annotation JSON files
        speeds up loading times by querying and displaying only
        celltypes containing significant data.
  Bulk Annotations
        `gene_locations/` and `snp_locations/` store annotation data in
        bulk so nearby genes/SNPs in the chromosome can be quickly
        loaded without opening a bunch of small files.
  Filtering
        By pre-filtering to only significant entries, we reduce the
        amount of data that needs to be processed and plotted.
  Parquet Format
        Storing files in Parquet format allows for compression and
        increases speed over text formats like CSV or TSV.

  Overall, we trade some disk space (especially individual annotation
  files) for much faster loading times. In future versions, individual
  annotation files may potentially be removed if speed is sufficient.


5 Scripts
=========

  The dataset is constructed through six scripts:

  1. *Preprocessing*
     - This script should vary depending on the input data.
     - Renames QTL data columns to `gene_id, snp_id, p_value,
       beta_value'.
     - Splits QTL data by celltype into `unfiltered_celltypes/'.
     - Renames gene and SNP annotation columns to `gene_id,
       position_start, position_end, strand' and `snp_id, position'.
     - Splits gene and SNP annotation files by chromosome into
       `gene_locations/` and `snp_locations/` (deduplicated).

  2. *Filter by Significance*
     - Filters each celltype in `unfiltered_celltypes/' for entries with
       p > 0.01.
     - Outputs filtered TSVs into `filtered_celltypes/'.

  3. *Celltype Mapping*
     - Prompts user for display names that correspond to filenames in
       `celltypes/'.
     - Generates `celltype_mapping.json' accordingly.

  4. *Global Significance Filtering*
     - Identifies gene-SNP pairs significant in at least one celltype
       using `filtered_celltypes/'.
     - Filters `unfiltered_celltypes/' accordingly.
     - Outputs filtered TSVs into `celltypes/'.

  5. *Generate Annotation JSONs*
     - Creates individual JSON files for each gene and SNP in
       `gene_locations/` and `snp_locations/`.
     - Includes list of celltypes in which each appears.

  6. *Convert to Parquet*
     - Converts all necessary TSV files to Parquet format for final use.

  Intermediate folders like `unfiltered_celltypes/` and
  `filtered_celltypes/` are not automatically deleted; users can choose
  to remove them if disk space is an issue.



Footnotes
_________

[1] In the next release, this will be replaced with all nearby genes,
with appropriate error handling if not present in the dataset.

[2] In the next release, this will be replaced with GWAS data.

[3] Entries are significant if −log_{10}(p) > 2, which is equivalent
to p < 0.01.
