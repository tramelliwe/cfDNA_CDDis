#!/bin/bash
set -ue # stop on any error
source ~/miniconda3/bin/activate wgbs

cd ~
nr_atlas=$1 #which atlas (fold) are we working with
nr_samples=$2 #how many samples with coverage x and enrichment y have to be generated and deconvoluted
coverage=$3 #coverage (0.5X, 1X, 5X, 10X, ...)
if  (($(echo "$coverage == 0.5" | bc -l))); then #poor handling of decimal in folder paths
  coverage_folder_name="half"
else
  coverage_folder_name=$coverage
fi

precision=0.001
#When generating random proportions, the lowest proportion is 0.1%. When generating in-silico mixtures
#we want to have the same precision. This precision variable ensures that we don't sample more than 0.1% reads at a time

#based on real samples, we know that the CpG coverage usable for UXM for a 1X sequencing is ~14,000,000. 
#This does not correspond to the ~29,000,000 CpG in the genome because UXM only considers reads covering >= 3 CpGs.
#For example, for a coverage of 1X: If proportion of tissue x is y%, then we sample reads from this tissue until 
#the CpG covered reaches y%*14,000,000 

#C are methylated CpG (pat file syntax)
#T are unmethylated CpG (pat file syntax)
U=false #unmethylated read enrichment
M=false #methylated read enrichment
for arg in "$@"
do
    if [ "$arg" == "-U" ]; then
        U=true
        output_suffix="_U"
        total_cpg_covered=$(echo "$coverage * 24759817" | bc) # 24759817 is the average CpG-coverage for a 1X sequencing with U enrichment
        line_per_sampling=$(echo "$precision * $total_cpg_covered/10" | bc) #around 10 CpG per line on average (based on real samples)
		#line_per_sampling is the target amount of line we add per round. For this, we take into account the average
		#CpG coverage per line (10 for U reads) and we take the precision into account.
        step=$(echo "$line_per_sampling * 100/6" | bc)  #6% actually taken after sampling for enrichment
		#the number of line sampled is not equal to the number of line actually added to the file.
		#lines can be discarded in not covering enough CpG, or if not fullfiling the enrichment parameter
    elif [ "$arg" == "-M" ]; then
        M=true
        output_suffix="_M"
        total_cpg_covered=$(echo "$coverage * 13755454" | bc) # 13755454 is the average CpG-coverage for a 1X sequencing with M enrichment
        line_per_sampling=$(echo "$precision * $total_cpg_covered/6" | bc) #around 6 CpG per line on average (based on real samples)
        step=$(echo "$line_per_sampling * 100/30" | bc)  #30% actually taken after sampling for enrichment
    else
        output_suffix=""
        total_cpg_covered=$(echo "$coverage * 13755454" | bc) # 13755454 is the average CpG-coverage for a 1X sequencing without enrichment
        line_per_sampling=$(echo "$precision * $total_cpg_covered/6" | bc) #around 6 CpG per line on average (based on real samples)
        step=$(echo "$line_per_sampling * 100/60" | bc)  # number of lines to add per round when sampling reference files
    fi
done


mkdir -p workspace/004_loocv/data/tmp
mkdir -p workspace/004_loocv/data/plots
mkdir -p workspace/004_loocv/data/output

# Define an array of input files
readarray -t files < <(awk -F',' 'NR>1 {print $1}' workspace/004_loocv/data/sets/test_files_${nr_atlas}.csv) #test set as output by generate_sets.R

#If moved to another atlas, we want to erase the files that are not needed bc of storage limitations
# Remove files that are not in the array
for existing_file in workspace/004_loocv/data/tmp/*.pat; do
  # Get the base name of the existing file without the directory and extension
  base_name=$(basename "$existing_file" .pat)
  # Check if the base name is in the array of files
  if [[ ! " ${files[@]} " =~ " ${base_name} " ]]; then
    # Remove the file if it's not in the array
    rm "$existing_file"
  fi
done

unzip_if_needed() {
  file="$1"
  output_file="workspace/004_loocv/data/tmp/${file}.pat"
  if [ ! -f "$output_file" ]; then
    gzip -dc "/mnt/d/pats/${file}.pat.gz" > "$output_file"
  fi
}

export -f unzip_if_needed

# Unzip files that are not yet present using GNU Parallel
parallel --no-notice unzip_if_needed ::: "${files[@]}"

regular_sampling() { #this is sampling without any enrichment.
  local step=$1
    softwares/powershuf.py --file "$file" -n "$step" |  awk 'length($3)>=3 {print $1 "\t" $2 "\t" $3 "\t" 1}' >> "$sampling_output" 
	#sample $step lines from the pat file then only takes the ones covering >=3 CpGs and add it the output
	#references have ~2 reads per line, which is not reflecting real samples so we impose 1
    added_lines=$(($(wc -l "$sampling_output" | awk '{print $1}') - $total_lines))
    total_lines=$(wc -l "$sampling_output" | awk '{print $1}')
    cpg_count=$(($cpg_count + $(awk '{sum+=length($3)*$4} END {print sum}' <(tail -n "$added_lines" "$sampling_output"))))
    #echo "$cpg_count"
    #for a step of 10,000 --> 6000 lines & 70,000 CpGs added per round
}
U_sampling() {
  local step=$1
    softwares/powershuf.py --file "$file" -n "$step" | awk 'length($3)>=3 {print $1 "\t" $2 "\t" $3 "\t" 1}' | awk 'BEGIN {FS="\t"; OFS="\t"} {
    col3=$3
    countC = gsub(/C/, "", col3) 
    countT = gsub(/T/, "", col3)
    len = length($3)
    if (countC > 0.75*len) {
        print $0, "M"}
    else if (countT > 0.75*len) {
        print $0, "U"}
    else {
        print $0, "X"}}'  | awk 'BEGIN {FS="\t"; OFS="\t"} $5=="U" {print $1, $2, $3, $4}'   >> "$sampling_output" 
	#count the number of C (=methylated CpGs) and T (=unmethylated C).
	#Following the syntax of UXM, a read is considered "U" if more than 75% unmethylated 
	#or "M" if more than 75% methylated
	#or 'X' if in between
	#Here, only keep U reads
    added_lines=$(($(wc -l "$sampling_output" | awk '{print $1}') - $total_lines))
    total_lines=$(wc -l "$sampling_output" | awk '{print $1}')
    cpg_count=$(($cpg_count + $(awk '{sum+=length($3)*$4} END {print sum}' <(tail -n "$added_lines" "$sampling_output"))))
    #echo "$cpg_count"
}
M_sampling() {
  local step=$1
    softwares/powershuf.py --file "$file" -n "$step" | awk 'length($3)>=3 {print $1 "\t" $2 "\t" $3 "\t" 1}' | awk 'BEGIN {FS="\t"; OFS="\t"} {
    col3=$3
    countC = gsub(/C/, "", col3) 
    countT = gsub(/T/, "", col3)
    len = length($3)
    if (countC > 0.75*len) {
        print $0, "M"}
    else if (countT > 0.75*len) {
        print $0, "U"}
    else {
        print $0, "X"}}'  | awk 'BEGIN {FS="\t"; OFS="\t"} $5=="M" {print $1, $2, $3, $4}'   >> "$sampling_output" #30% of $step 
    added_lines=$(($(wc -l "$sampling_output" | awk '{print $1}') - $total_lines))
    total_lines=$(wc -l "$sampling_output" | awk '{print $1}')
    cpg_count=$(($cpg_count + $(awk '{sum+=length($3)*$4} END {print sum}' <(tail -n "$added_lines" "$sampling_output"))))
    #echo "$cpg_count"
}

generate_sample() {
  local mix=$1
    # Reassemble the files array from the exported string
  IFS=',' read -r -a files <<< "$files_str"
  # generate proportions: see R script for info
  workspace/004_loocv/scripts/generate_proportions.R workspace/004_loocv/data/tmp/proportions_mix_${mix}.csv
  readarray -t fractions < <(awk -F',' 'NR>1 {print $2}' workspace/004_loocv/data/tmp/proportions_mix_${mix}.csv)
  echo "Sampling mix $mix"
  start_time=$(date +%s)
  for i in "${!files[@]}"; do #for each tissue, sample the reference file of the test set until you reach the determined proportion
      name="${files[$i]}"
      file="workspace/004_loocv/data/tmp/${name}.pat"
      fraction="${fractions[$i]}"
      if (( $(echo "$fraction == 0" | bc -l) )); then
        continue
      fi

      cpg_count=0
      sampling_output="workspace/004_loocv/data/tmp/sampled_${name}_${mix}.pat"
      touch "$sampling_output"
      total_lines=0
      
	  #for more efficiency, each sampling starts by adding ~60% of the total number of lines at once. 
	  #The rest of the lines is added per round, little by little, in a while loop.
	  
	  
      if [ "$U" = true ]; then
        U_sampling $(echo "$fraction*$total_cpg_covered/16 * 100/6" | bc) #adding fraction * 60% of the total amount of lines (total_cpg/12) taking into account the 6% of actual lines added
        while ((cpg_count < $(echo "$fraction * $total_cpg_covered / 1" | bc) )); do
            U_sampling $step
        done
      elif [ "$M" = true ]; then
        M_sampling $(echo "$fraction*$total_cpg_covered/16 * 100/30" | bc) #adding 60% of the total amount of lines (total_cpg/12) taking into account the 30% of actual lines added
        while ((cpg_count < $(echo "$fraction * $total_cpg_covered / 1" | bc) )); do
            M_sampling $step
        done
      else
        regular_sampling $(echo "$fraction*$total_cpg_covered/16 * 100/60" | bc) #adding 60% of the total amount of lines (total_cpg/12) taking into account the 60% of actual lines added
        while ((cpg_count < $(echo "$fraction * $total_cpg_covered / 1" | bc) )); do
            regular_sampling $step
        done
      fi
	  
	  #adding a 'key' column to uniquely identify each read
      awk -F'\t' '{print $0 "\t" $2 "_" $3}' "$sampling_output" > "workspace/004_loocv/data/tmp/${name}_${mix}_keys.pat"
      rm $sampling_output
  done
  end_time=$(date +%s)
  echo "It took" $((end_time - start_time)) "sec. to generate mix" $mix "at coverage" $coverage "with U=" $U "and M=" $M 

  # Combine and process all tissues into a single file
  cat workspace/004_loocv/data/tmp/*${mix}_keys.pat | sort -k5,5 > workspace/004_loocv/data/tmp/${mix}_combined_keys.pat
  #deleting any duplicated lines. Generally does not have major impact on the precision.
  awk -F'\t' '{arr[$5]+=$4} END {for (i in arr) print i, arr[i]}' workspace/004_loocv/data/tmp/${mix}_combined_keys.pat > workspace/004_loocv/data/tmp/${mix}_aggregated.pat 
  awk 'NR==FNR {a[$1]=$2; next} {print $0, a[$5]}' workspace/004_loocv/data/tmp/${mix}_aggregated.pat workspace/004_loocv/data/tmp/${mix}_combined_keys.pat > workspace/004_loocv/data/tmp/${mix}_aggregated_fin.pat
  awk '!seen[$5]++' workspace/004_loocv/data/tmp/${mix}_aggregated_fin.pat | sort -k2n | awk '{print $1 "\t" $2 "\t" $3 "\t" $6}' > workspace/004_loocv/data/tmp/${mix}_input.pat
  gzip workspace/004_loocv/data/tmp/${mix}_input.pat

  # Run the deconvolution software
  softwares/UXM_deconv/uxm deconv workspace/004_loocv/data/tmp/${mix}_input.pat.gz --rlen 3 --atlas workspace/004_loocv/atlases/atlas_${nr_atlas}.tsv --force --output workspace/004_loocv/data/tmp/mix_${mix}.csv

  #For easier downstream analysis, we combine in a single table the real proportion and the estimated proportion
  paste -d "," workspace/004_loocv/data/tmp/mix_${mix}.csv workspace/004_loocv/data/tmp/proportions_mix_${mix}.csv > workspace/004_loocv/output/atlas_${nr_atlas}_coverage_${coverage_folder_name}${output_suffix}_$(date +"%Y%m%d%N").csv

  rm workspace/004_loocv/data/tmp/*${mix}_keys.pat
  rm workspace/004_loocv/data/tmp/${mix}_input.pat.gz
  rm workspace/004_loocv/data/tmp/${mix}_input.pat.gz.csi
  rm workspace/004_loocv/data/tmp/${mix}_combined_keys.pat
  rm workspace/004_loocv/data/tmp/${mix}_aggregated.pat
  rm workspace/004_loocv/data/tmp/${mix}_aggregated_fin.pat
  rm workspace/004_loocv/data/tmp/mix_${mix}.csv
  rm workspace/004_loocv/data/tmp/proportions_mix_${mix}.csv
}



export -f generate_sample
export nr_atlas
export total_cpg_covered
export step
files_str=$(IFS=','; echo "${files[*]}")
export files_str
export coverage
export coverage_folder_name
export M
export U
export output_suffix
export -f regular_sampling
export -f U_sampling
export -f M_sampling

# Parallelize the loop
seq 1 $nr_samples  | parallel --line-buffer generate_sample 
