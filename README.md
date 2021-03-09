# Adult Extramedullary Haematopoiesis


The code in this repository pertains to the combined analysis of the data associated to the following studies:

[S-SUBS4](https://www.ebi.ac.uk/biostudies/studies/S-SUBS4)

and

[S-SUBS10](https://www.ebi.ac.uk/biostudies/studies/S-SUBS10)


## Fastq processing

``00_process_fastq_files/``

 - these scripts handle the processing of fastq sequencing files, using either ``Cellranger`` (RNA data) or ``CITE-seq-Count`` (Ab data) in order to yield feature count matrices.


## Data transformations

``01_data_transformations``


# Multi-tissue analyses

``02_analyses``

# Visualisations

``04_visualisations``
