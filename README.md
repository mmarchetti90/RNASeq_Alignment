# Bulk RNASeq Alignment
Nextflow pipeline for bulk RNASeq reads alignement.

Simple Nextflow pipeline for FastQC and alignemnt of RNASeq raw reads.

Analysis is based on a Docker image, which can be build with the shell script in the 0_dockerfiles directory.

Edit nextflow.config file for declaring experiment type.

Works only for bulk RNASeq samples, both single-read and paired-end.
