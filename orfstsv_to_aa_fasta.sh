#!/usr/bin/env bash

## GEDI PRICE cit.orfs.tsv to aa fasta. 
## Gedi orfs.tsv is 0 based. Neg strand exons are in an order c/w BED already. 
## clean up orfs.tsv to mimic BED. it's not sorted, so need to sort 

#awk 'BEGIN{FS="\t"; OFS="\t"} FNR==1{next} {
	### here changed to location_array[FNR]=$4 (location candidate) instead of $3 (final location)
	#location_array[FNR]=$4; 
	#n = split(location_array[FNR], new_location_array, ":");
	### split new_location_array[1] further into chr number (l) and strand (m)
	#m = substr(new_location_array[1], length(new_location_array[1]), 1);
	#l = substr(new_location_array[1], 1, (length(new_location_array[1]) - 1)); 
	### split new_location_array[2] further into start, end pos, width, and each exon start pos
	#o = split(new_location_array[2], new_location_array_sub, "|");
		#for (i = 1; i <= o; i++){
			### then get the start pos and end pos
			#if ( i == 1 ){ p = split(new_location_array_sub[i], first_exon_pos, "-")};
			#if ( i == o ){ q = split(new_location_array_sub[i], last_exon_pos, "-")};
			### then get width and relative pos from the start pos
			#r = split(new_location_array_sub[i], any_exon_pos, "-");
				#exon_width[FNR]=(exon_width[FNR] "," (any_exon_pos[2] - any_exon_pos[1]));
				#distance_from_start[FNR]=(distance_from_start[FNR] "," (any_exon_pos[1] - first_exon_pos[1]));
		#}
	#print "chr" l, first_exon_pos[1], last_exon_pos[2], \
	#($1 ":" $2), 0, m, first_exon_pos[1], last_exon_pos[2], 0, o, \
	### remove "," in front of the first width and distance  
	#gensub(/,(.+)$/, "\\1", "g", exon_width[FNR]), \
	#gensub(/,(.+)$/, "\\1", "g", distance_from_start[FNR]), \
	### keep some other columns for now
	#$5, $6, $7, $8, $9, $16}' \
#DM1_5cit.orfs.tsv > DM1_5cit.orfs_int_1.bed

## add geneID based on ensembl tx. To do this, first, re-split the gene name.
#awk 'BEGIN{FS="\t"; OFS="\t"} {my_array[FNR]=$4;
#n = split(my_array[FNR], new_my_array, ":");
#### split further to get ENSMUST ID
#m = split(new_my_array[2], new_my_array_sub, "_");
#print $1,$2,$3, new_my_array[1], new_my_array_sub[1], $5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' \
#DM1_5cit.orfs_int_1.bed > DM1_5cit.orfs_int_2.bed

## prepare mart_export_ENSMUST_transcript_ID.txt which looks like this 
## awk 'BEGIN{FS="\t"; OFS="\t"} $0 ~ /Ppp1r15a-201/{print $2}' mart_export_ENSMUST_transcript_ID.txt
	## ENSMUSG00000040435.12	ENSMUST00000042105.10	Ppp1r15a-201	Ppp1r15a

### register file1 ENSMUST (biomart file) with ENSMUST as an index ($2). 
### Then search for the index in file2 using file2 $5 ENSMUST as an index.  
#awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_array[$2]=$0; next} \
#($5 in my_array){ print my_array[$5], $0 }' \
#mart_export_ENSMUST_transcript_ID.txt DM1_5cit.orfs_int_2.bed > DM1_5cit.orfs_int_3.bed

## clean up and sort 
#awk 'BEGIN{FS="\t"; OFS="\t"} \
#{print $5,$6,$7,($8 ":" $9 ":" $3 ":" $4),$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' \
#DM1_5cit.orfs_int_3.bed > DM1_5cit.orfs_int_4.bed

#sort -k 1,1 -k 2,2n -k 3,3n DM1_5cit.orfs_int_4.bed > DM1_5cit.orfs_int_5.bed

##awk '{print $1}' DM1_5cit.orfs_int_5.bed | uniq

## at this point col13=codon, col14=type, col15=Start, col16=Range, col17=p value, col18=Total
## col2: col3 pos 0 of start and pos 3 of stop (for + strand)

## select uoORF, uORF, iORF, dORF and also combine ORF type ($14) to $4
#awk 'BEGIN{FS="\t"; OFS="\t"} $14 ~ /uoORF|uORF|iORF|dORF/{print $1,$2,$3,($4 ":" $14),$5,$6,$7,$8,$9,$10,$11,$12}' \
#DM1_5cit.orfs_int_5.bed > DM1_5cit.orfs_uo_i_dORF_final.bed

### Bedtools to extract coding sequences from reference genome
##wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
##gunzip GRCm38.primary_assembly.genome.fa.gz

### transcribe
#bedtools getfasta -name -s -split -fo gencode_vM23_DM1_5cit.orfs_uo_i_dORF.fasta \
#-fi GRCm38.primary_assembly.genome.fa -bed DM1_5cit.orfs_uo_i_dORF_final.bed

### faTrans to translate from DNA to AA
###wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faTrans
###chmod u+x faTrans
#./faTrans gencode_vM23_DM1_5cit.orfs_uo_i_dORF.fasta \
#gencode_vM23_DM1_5cit.orfs_uo_i_dORF_translated.fasta

## flatten the fasta 
#awk 'FNR==1{print $0; next} $0 !~ /^>/{printf "%s", $0} $0 ~ /^>/{printf "\n%s\n", $0} END{printf "\n"}' \
#gencode_vM23_DM1_5cit.orfs_uo_i_dORF_translated.fasta > gencode_vM23_DM1_5cit.orfs_uo_i_dORF_translated_flat.fasta

## cut very short peptides (arbitrarily <8 (i.e., keep M...Z >= 8)
#awk 'BEGIN{FS="\t";OFS="\t"} FNR % 2 == 1{my_array[FNR]=$0; next} FNR % 2 == 0{ \
#	if ( length($0) >= 8 ) { printf "%s\n%s\n", my_array[(FNR - 1)], $0 }}' \
#	gencode_vM23_DM1_5cit.orfs_uo_i_dORF_translated_flat.fasta > gencode_vM23_DM1_5cit.orfs_uo_i_dORF_translated_flat_above_8aa.fasta

## this is uo_i_dORF amino acid ref (without conventional CDS) for mass spec: gencode_vM23_DM1_5cit.orfs_uo_i_dORF_translated_flat_above_8aa.fasta
