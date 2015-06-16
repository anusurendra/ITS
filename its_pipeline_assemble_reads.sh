#!/bin/bash

#get input arguments
args=("$@")
file_dir=${args[0]}
file_prefix=${args[1]}

FILES="$file_dir/*_R1_*.fastq"

#create a folder to hold output

output_dir="$file_dir/UPARSE_OUTPUT"

if [ ! -d "$output_dir" ]; then
  mkdir "$file_dir/UPARSE_OUTPUT"
  else  
  echo -e "output folder $file_dir/UPARSE_OUTPUT already exists \n"
  exit 1
fi

for f in $FILES
do
 
  echo -e "----Processing $f----\n"
  
  #create the file string/sample string from the filename
  result_string=${f##*/}
  result_string=$(echo $result_string | sed 's/_R1_.*//g')
  
  sampleid=$(echo $result_string | sed 's/_L001.*//g')
  
  #forward read file
  fwreads="$file_dir/$result_string"
  fwreads+="_R1_001.fastq"
  
  #reverse read file
  revreads="$file_dir/$result_string"
  revreads+="_R2_001.fastq"
  
  #assembled fastq file
  fastqout="$output_dir/$result_string"
  fastqout+="_assembled.fastq"
   
  #assembled fasta file
  fastaout="$output_dir/$result_string"
  fastaout+="_assembled.fasta"
  
  
  #assembled fastq not merged forward file
  fastqoutnotmergedfwd="$output_dir/$result_string"
  fastqoutnotmergedfwd+="_R1_001_not_merged.fastq"
  
  #assembled fastq not merged reverse file
  fastqoutnotmergedrev="$output_dir/$result_string"
  fastqoutnotmergedrev+="_R2_001_not_merged.fastq"
  
  #alignment file
  alnout="$output_dir/$result_string"
  alnout+="_alnout.txt"
    
  #assembling reads
  options_mp=( "-fastq_mergepairs")
  command_mp=( "usearch8.0.1616_i86linux32" "${options_mp[@]}" "$fwreads" -reverse "$revreads" -fastqout "$fastqout" -fastaout "$fastaout" -fastqout_notmerged_fwd "$fastqoutnotmergedfwd" -fastqout_notmerged_rev "$fastqoutnotmergedrev" -fastq_merge_maxee "1.0" -alnout "$alnout" )
  "${command_mp[@]}"
  
  
  
  #adding the sample name to the fasta and fastq files
  
  #assembled fasta file with filename
  fastaout_sample="$output_dir/$result_string"
  fastaout_sample+="_assembled_sample.fasta"
  
  #assembled fastq file with filename
  fastqout_sample="$output_dir/$result_string"
  fastqout_sample+="_assembled_sample.fastq"
  
  #M00953:102:000000000-A8HM3:1:1101:3465:11015

  #s/^>M[0-9]\+:[0-9]\+:[0-9]\+-[a-zA-Z0-9]\+:[0-9]\+:[0-9]\+:
  sed  "s/^>M[0-9]\+:[0-9]\+:[0-9]\+-[a-zA-Z0-9]\+:[0-9]\+:[0-9]\+:/>barcodelabel=$sampleid:/g" $fastaout > $fastaout_sample
  sed  "s/^@M[0-9]\+:[0-9]\+:[0-9]\+-[a-zA-Z0-9]\+:[0-9]\+:[0-9]\+:/@barcodelabel=$sampleid:/g" $fastqout > $fastqout_sample
  
  
  #quality filter
  fastaout_sample_filtered="$output_dir/$result_string"
  fastaout_sample_filtered+="_assembled_sample_filtered.fasta"
  
  options_ff=( "-fastq_filter")
  command_ff=( "usearch8.0.1616_i86linux32" "${options_ff[@]}" "$fastqout_sample" -fastq_maxee "0.5" -fastaout "$fastaout_sample_filtered" -notrunclabels )
  "${command_ff[@]}"
  
 
  
  #concatenated read files
  fastqout_concat="$output_dir/$result_string"
  fastqout_concat+="_concat_assembled.fastq"
  
  fastaout_concat="$output_dir/$result_string"
  fastaout_concat+="_concat_assembled.fasta"

  #concatenating unassembled reads
  options_fj=( "-fastq_join")
  command_fj=( "usearch8.0.1616_i86linux32" "${options_fj[@]}" "$fastqoutnotmergedfwd" -reverse "$fastqoutnotmergedrev" -fastqout "$fastqout_concat" -fastaout "$fastaout_concat" -notrunclabels )
  "${command_fj[@]}"
  
  #adding the sample name to the concatenated fasta and fastq files
  fastaout_sample_concat="$output_dir/$result_string"
  fastaout_sample_concat+="_assembled_sample_concat.fasta"
  
  fastqout_sample_concat="$output_dir/$result_string"
  fastqout_sample_concat+="_assembled_sample_concat.fastq"
  
  sed  "s/^>M[0-9]\+:[0-9]\+:[0-9]\+-[a-zA-Z0-9]\+:[0-9]\+:[0-9]\+:/>barcodelabel=$sampleid:/g" $fastaout_concat > $fastaout_sample_concat
  sed  "s/^@M[0-9]\+:[0-9]\+:[0-9]\+-[a-zA-Z0-9]\+:[0-9]\+:[0-9]\+:/@barcodelabel=$sampleid:/g" $fastqout_concat > $fastqout_sample_concat
  
  
  #quality filter
  fastaout_sample_concat_filtered="$output_dir/$result_string"
  fastaout_sample_concat_filtered+="_assembled_sample_concat_filtered.fasta"
  
  options_cff=( "-fastq_filter")
  command_cff=( "usearch8.0.1616_i86linux32" "${options_ff[@]}" "$fastqout_sample_concat" -fastq_minlen "100" -fastq_truncqual "30" -fastaout "$fastaout_sample_concat_filtered" -notrunclabels )
  "${command_cff[@]}"
  
  #replace gap
  fastaout_sample_concat_filtered_gap_1="$output_dir/$result_string"
  fastaout_sample_concat_filtered_gap_1+="_assembled_sample_concat_filtered_gap_1.fasta"
  fastaout_sample_concat_filtered_gap="$output_dir/$result_string"
  fastaout_sample_concat_filtered_gap+="_assembled_sample_concat_filtered_gap.fasta"
  sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' $fastaout_sample_concat_filtered > $fastaout_sample_concat_filtered_gap_1
  sed -r "/^>/!s/N{8}/-/g" $fastaout_sample_concat_filtered_gap_1 > $fastaout_sample_concat_filtered_gap
  
  echo -e "----Finished Processing $f----\n"
done

#merging samples into one file
merged_file="$output_dir/$file_prefix"
merged_file+="_assembled.fa"
cat $output_dir/*_assembled_sample_filtered.fasta > $merged_file


#merging samples into one file
merged_file_concat="$output_dir/$file_prefix"
merged_file_concat+="_assembled_concat.fa"
cat $output_dir/*_assembled_sample_concat_filtered_gap.fasta > $merged_file_concat
