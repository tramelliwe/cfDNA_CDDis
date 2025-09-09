#!/bin/bash
set -ue #stop on any error

nr_atlas=$1
start_id=$2
#main script for LOOCV

cd ~/workspace/004_loocv/
mkdir -p data/markers
mkdir -p data/atlases
mkdir -p data/sets
#1. Select training and test set

#Go through each group and randomly select one sample as being in the test set. 
./generate_sets.R ${start_id} ${nr_atlas} loocv_01
						   
#2. Create atlas based on training set
source ~/miniconda3/bin/activate wgbs

for i in $(seq $start_id $((start_id+nr_atlas)));
do

wgbstools find_markers --config_file params.txt --groups_file sets/training_files_${i}.csv --out_dir markers/markers_${i}
#concatenate the markers and remove duplicates
./concat_markers.R markers/markers_${i} ${i}

~/softwares/UXM_deconv/uxm build --markers markers/all_markers_${i}.tsv --output atlases/atlas_${i}.tsv -@ 12 --rlen 3 --groups sets/training_files_${i}.csv --pats /mnt/d/pats/*pat.gz --verbose

done


#3. Generate proportions, samples, and deconvolution
for i in $(seq 3 10);
do
~/loocv_01/sampling_more_U.sh $i 50 0.5 #1000 samples per atlas at coverage 0.5X
~/loocv_01/sampling_more_U.sh $i 50 2 #1000 samples per atlas at coverage 0.5X
~/loocv_01/sampling_more_U.sh $i 50 5 #1000 samples per atlas at coverage 0.5X
done


#3. Create input for the software that is the most promising in first trials
	# i.e., multiple in silico samples and the atlas
#4. 