#!/bin/bash

#get input arguments
args=("$@")
file_dir=${args[0]}
assembled_file=${args[1]}
ref_db_file_uchime=${args[2]}
blastn_db_file=${args[3]}
taxa_file=${args[4]}
evalue_cutoff=${args[5]}
hitabundance_cutoff=${args[6]}


prefix=$(echo $assembled_file | sed 's/\..*//g')

output_dir="$file_dir/UPARSE_OUTPUT"

if [ ! -d "$output_dir" ]; then
  mkdir "$file_dir/UPARSE_OUTPUT"
else  
 echo -e "output folder $file_dir/UPARSE_OUTPUT already exists \n"
 exit 1
fi



# Dereplication

echo -e "----Dereplication $assembled_file----\n"

fastaout_derep="$output_dir/$prefix"
fastaout_derep+="_dereplicated.fasta"

options_drep=( "-derep_fulllength")
command_drep=( "usearch8.0.1616_i86linux32" "${options_drep[@]}" "$file_dir/$assembled_file" -fastaout "$fastaout_derep" -sizeout -notrunclabels )
"${command_drep[@]}"


# Abundance sort and discard singletons

echo -e "----Abundance sort and discard singletons $assembled_file----\n"

fastaout_sorted="$output_dir/$prefix"
fastaout_sorted+="_dereplicated_sorted.fasta"

options_sorted=( "-sortbysize")
command_sorted=( "usearch8.0.1616_i86linux32" "${options_sorted[@]}" "$fastaout_derep" -fastaout "$fastaout_sorted" -minsize "2" -notrunclabels )
"${command_sorted[@]}"

# OTU clustering

echo -e "----OTU clustering $assembled_file----\n"

fastaout_otus="$output_dir/$prefix"
fastaout_otus+="_dereplicated_sorted_otus.fasta"

uparseout_otus="$output_dir/$prefix"
uparseout_otus+="_dereplicated_sorted_otus_uparse.txt"

options_otus=( "-cluster_otus")
command_otus=( "usearch8.0.1616_i86linux32" "${options_otus[@]}" "$fastaout_sorted" -otus "$fastaout_otus" -uparseout "$uparseout_otus" -relabel "OTU_" -sizein -sizeout -notrunclabels)
"${command_otus[@]}"


# Chimera filtering using reference database

echo -e "----Chimera filtering using reference database $assembled_file----\n"

uchimeout_txt="$output_dir/$prefix"
uchimeout_txt+="_dereplicated_sorted_otus_uchime.txt"

options_uchime=( "-uchime_ref")
command_uchime=( "usearch8.0.1616_i86linux32" "${options_uchime[@]}" "$fastaout_otus" -db "$ref_db_file_uchime" -uchimeout "$uchimeout_txt" -minh "0.20" -mindiv "0.50" -notrunclabels -strand plus)
"${command_uchime[@]}"

uchimeout_labels_txt="$output_dir/$prefix"
uchimeout_labels_txt+="_dereplicated_sorted_otus_uchime_labels.txt"

grep "N$" $uchimeout_txt | cut -f2 > $uchimeout_labels_txt

uchimeout="$output_dir/$prefix"
uchimeout+="_dereplicated_sorted_otus_chimeras_removed.fasta"

options_subset=( "-fastx_getseqs")
command_subset=( "usearch8.0.1616_i86linux32" "${options_subset[@]}" "$fastaout_otus" -labels "$uchimeout_labels_txt" -fastaout "$uchimeout" -notrunclabels )
"${command_subset[@]}"



# Map reads (including singletons) back to OTUs

echo -e "----Map reads (including singletons) back to OTUs $assembled_file----\n"

ucout_map_reads="$output_dir/$prefix"
ucout_map_reads+="_dereplicated_sorted_otus_uparse.uc"

fastacout_unmap_reads="$output_dir/$prefix"
fastacout_unmap_reads+="_dereplicated_sorted_otus_uparse_unmapped.fa"

options_map_reads=( "-usearch_global")
command_map_reads=( "usearch8.0.1616_i86linux32" "${options_map_reads[@]}" "$file_dir/$assembled_file" -db "$uchimeout" -strand "plus" -id "0.97" -uc "$ucout_map_reads" -notmatched "$fastacout_unmap_reads"  -notrunclabels )
"${command_map_reads[@]}"


# Blast OTUs against FHiTINGS database to assign taxanomy

echo -e " Blasting OTUs against FHiTINGS database to assign taxanomy and create OTU table $assembled_file----\n"
blastoutput="$output_dir/$prefix"
blastoutput+="_allseqs.blastoutput"
command_blast=( "blastn" -query "$uchimeout" -out "$blastoutput" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -max_target_seqs "10" -db "$blastn_db_file" )
"${command_blast[@]}"


blastoutput_pathremoved="$prefix"
blastoutput_pathremoved+="_allseqs.blastoutput"

ucout_map_reads_pathremoved="$prefix"
ucout_map_reads_pathremoved+="_dereplicated_sorted_otus_uparse.uc"

#
command_fhitings=( "./fhitings_Anu.R" "--vanilla --default-packages=reshape,plyr,gtools" --filedir="$output_dir" --taxafile="$taxa_file" --evalue_cutoff="$evalue_cutoff"  --hit_abundance_cutoff="$hitabundance_cutoff" --blastfiles="$blastoutput_pathremoved" --ucfiles="$ucout_map_reads_pathremoved" )
"${command_fhitings[@]}"

biom_text_file="$output_dir/output/$prefix"
biom_text_file+="_allseqs_biom.txt"

biom_text_taxa_file="$output_dir/output/$prefix"
biom_text_taxa_file+="_allseqs_biom_taxa.txt"

biom_file="$output_dir/output/$prefix"
biom_file+="_allseqs.biom"

biom_metadata_file="$output_dir/output/$prefix"
biom_metadata_file+="_allseqs_metadata.biom"

command_biom_convert=( "biom" "convert" -i "$biom_text_file" -o "$biom_file" --table-type="otu table" )
"${command_biom_convert[@]}"

command_biom_add_metadata=( "biom" "add-metadata" --sc-separated="taxonomy" --observation-header="OTUID,taxonomy" --observation-metadata-fp="$biom_text_taxa_file" -i "$biom_file" -o "$biom_metadata_file" )
"${command_biom_add_metadata[@]}"