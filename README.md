# Adult Extramedullary Haematopoiesis


The code in this repository pertains to the combined analysis of the data associated to the following biostudies:

[S-SUBS4](https://www.ebi.ac.uk/biostudies/studies/S-SUBS4)

and

[S-SUBS10](https://www.ebi.ac.uk/biostudies/studies/S-SUBS10)


## Fastq processing

``00_process_fastq_files/``

 - Processing of fastq sequencing files, using either ``Cellranger`` (RNA data) or ``CITE-seq-Count`` (Ab data) in order to yield feature count matrices.


## Data transformations

``01_data_transformations``

- Combining and integrating datasets (Seurat3), batch correction (ComBat), QC filtering, log-normalisations, scaling and other pre-processing data transformations required for downstream analyses.


## Multi-tissue analyses

``02_analyses``

- Exploring the transformed data through methods and techniques such as pseudotime (DPT), cell scoring (including cell cycle), Self-Assembly Manifolds (SAM) and others.


## Visualisations

``03_visualisations``

- Code focused on creating visual output such as dimensionality reduction embedding scatter plots, violin plots, heatmaps, etc.
