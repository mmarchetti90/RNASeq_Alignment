# RNASeq_Alignment
Nextflow pipeline for RNASeq reads alignement.

Simple Nextflow pipeline for FastQC and alignemnt of RNASeq raw reads.

Analysis is based on a Docker image, automatically built if not present.

Edit nextflow.config file for declaring experiment type.

Works for bulk RNASeq samples, both single-read and paired-end.

Support for STARsolo scRNASeq run is still being tested.
