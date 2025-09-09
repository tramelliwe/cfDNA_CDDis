#!/bin/bash

#SBATCH --job-name=test_alignement
#SBATCH --partition=broadwell_256
#SBATCH -A ap_cds_tvb
#SBATCH --ntasks=1 --cpus-per-task=18
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=5g

set -ue

module purge
module load calcua/all
module load Kraken2/2.1.2-gompi-2022a
module load SAMtools/1.18-GCC-12.3.0
module load minimap2/2.26-GCCcore-12.3.0 
module load monitor/1.1.2

cd /user/antwerpen/209/vsc20913

for file in data/nonaligned_fastq/*fastq;
do
    sampleid=$(basename $file .fastq)

    kraken2 \
        --db scratch/pluspf_db \
        --threads 18 \
        --minimum-hit-groups 3 \
        --report scratch/${sampleid}_report.kraken2 \
        --report-minimizer-data \
        --minimum-base-quality 8 \
    scratch/${sampleid}_unaligned_part.fastq > scratch/${sampleid}.kraken2

    mv scratch/${sampleid}.kraken2 data/kraken2_results
    mv scratch/${sampleid}_report.kraken2 data/kraken2_results

done
