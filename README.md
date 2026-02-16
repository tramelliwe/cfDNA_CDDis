# cfDNA_CDDis

This repository contains supporting code for the manuscript:

“Nanopore Sequencing Enables Tissue-of-Origin and Pathogen Detection in ICU Plasma cfDNA”

The scripts reproduce key analyses described in the paper:
- Tissue-of-origin deconvolution of plasma cfDNA methylation (low-coverage Oxford Nanopore sequencing).
- Microbial detection from cfDNA sequencing using Kraken2.
- Cross-validation workflows for benchmarking classifier performance.

The code is organized into separate folders corresponding to analysis components:
- loocv/ – leave-one-out cross-validation (LOOCV) setup and execution.
- deconvolution/ – tissue-of-origin fraction estimation from cfDNA methylation.
- kraken/ – microbial cfDNA detection using Kraken2 and downstream summaries.
