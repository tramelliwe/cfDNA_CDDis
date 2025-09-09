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

for file in data/bam/*bam;
do
    sampleid=$(basename $file .bam)
    cp $file scratch

    # How many aligned reads and how many unaligned reads?
    touch ${sampleid}_summary.txt
    echo "> Human genome" >> ${sampleid}_summary.txt
    samtools stats /mnt/d/all_samples/bam/I-078_E0.bam -@ 18 | head -n 50 | grep ^SN | cut -f 2- >> ${sampleid}_summary.txt

    #keep unaligned reads only
    samtools fastq -f 4 -@ 18 scratch/${sampleid}.bam > scratch/${sampleid}_human_removed.fastq

    #align to lambda genome and extract unaligned reads only
    minimap2 -ax sr -t 18 scratch/lambda.fa scratch/${sampleid}_human_removed.fastq > scratch/${sampleid}_lambda_aligned.sam
    samtools fastq -f 4 -@ 18 scratch/${sampleid}_lambda_aligned.sam > scratch/${sampleid}_unaligned_part.fastq

    # How many reads aligned to the lambda genome?
    echo "> Lambda genome" >> ${sampleid}_summary.txt
    samtools stats scratch/${sampleid}_lambda_aligned.sam -@ 18 | head -n 50 | grep ^SN | cut -f 2- >> ${sampleid}_summary.txt

    #Quality filtering
    softwares/bbmap/bbduk.sh in=scratch/${sampleid}_unaligned_part.fastq out=scratch/${sampleid}_kraken_input.fastq trimq=8 qtrim=rl maq=8 entropy=0.6 entropywindow=50 entropyk=5
    #trim MAPQ<8, average quality <8, entropy <0.6

    mv scratch/${sampleid}_unaligned_part.fastq data/nonaligned_fastq

    mv $file data/bam/processed

    rm scratch/${sampleid}_lambda_aligned.sam
    rm scratch/${sampleid}_human_removed.fastq
    rm scratch/${sampleid}.bam
done
