!/usr/bin/env bash

## GTF to BED conversion 
## this script contains extra intermediary and inefficient steps as it was written for some other purposes. 

## gencode GTF GRCm38: https://www.gencodegenes.org/mouse/

## remove comment lines  
awk 'BEGIN{FS="\t"; OFS="\t"} $0 !~ /^#/{print $0}' fooM23.gtf > fooM23_1.gtf

## remove $2 and keep only transcript and exon records
awk 'BEGIN{FS="\t"; OFS="\t"} $3 ~ /transcript|exon/{print $1,$3,$4,$5,$7,$9}' fooM23_1.gtf > fooM23_2.gtf

## split original col9 (current col6) and get relevant info
awk 'BEGIN{FS="\t"; OFS="\t"} {my_array0[FNR]=($1 "\t" $2 "\t" $3 "\t" $4 "\t" $5); my_array9[FNR]=$6; \
n = split(my_array9[FNR], new_array, ";"); \
print my_array0[FNR], new_array[1], new_array[2], new_array[6], new_array[7]}' fooM23_2.gtf > fooM23_3.gtf

## concatenate transcripts together (with FS=" " removing "transcript_id" etc)
awk 'BEGIN{FS=" "; OFS="\t"} \
$2 ~ /transcript/{my_array[FNR]=($1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $9 "\t" $11); print my_array[FNR]}' \
fooM23_3.gtf > fooM23_tx.gtf

## concatenate exons together per transcript and also calculate width at this time
## keep gene start and end pos (in addition to exon start and end)
awk 'BEGIN{FS=" "; OFS="\t"} \
$2 ~ /exon/{my_array[FNR]=($1 "\t" $2 "\t" $3 "\t" $4 "\t" ($4 - $3 +1) "\t" $5 "\t" $7 "\t" $9 "\t" $11); print my_array[FNR]}' \
fooM23_3.gtf > fooM23_exons.gtf

## collapse exons file per tx, which should then match up with lines with tx
## calculating exon pos with respect to gene pos is hard when exon and gene info are separate apart. 
## so will first collapse exons into per gene, then combine, then later deal with position issues.
## just to be sure, sort based on fooM23_exons.gtf col 9 first. 
sort -k 1,1 -k 9,9 -k 3,3n -k 4,4n fooM23_exons.gtf > fooM23_exons_sorted.gtf

awk 'BEGIN{FS="\t"; OFS="\t"} \
FNR == 1{my_array_exon_start_end_pos[FNR]=($3 "-" $4); my_tx[FNR]=$9; next} \
FNR != 1{my_array_exon_start_end_pos[FNR]=($3 "-" $4); my_tx[FNR]=$9} \
my_tx[(FNR - 1)] != my_tx[FNR]{ \
	my_array_exon_start_end_pos[FNR]=(my_array_exon_start_end_pos[FNR]); print my_array_exon_start_end_pos[(FNR - 1)], my_tx[(FNR - 1)]} \
my_tx[(FNR - 1)] == my_tx[FNR]{ \
	my_array_exon_start_end_pos[FNR]=(my_array_exon_start_end_pos[(FNR - 1)] "," my_array_exon_start_end_pos[FNR])} \
END{my_array_exon_start_end_pos[FNR]=(my_array_exon_start_end_pos[FNR]); print my_array_exon_start_end_pos[FNR], my_tx[FNR]}' \
fooM23_exons_sorted.gtf > fooM23_exon_start_end_pos.gtf

## combine fooM23_tx.gtf and fooM23_exon_start_end_pos.gtf. assign index to my_array based on the actual value of $8 in file1 ("Xkr4-201" etc). 
## This allows us to do index search (indx in array) during file2 NR. 
awk 'BEGIN{FS="\t"; OFS="\t"} \
FNR==NR{ my_array[$8]=$0; next } ($2 in my_array) {print my_array[$2], $0}' fooM23_tx.gtf fooM23_exon_start_end_pos.gtf > fooM23_tx_exon_combined.gtf

## figure out exon position (keep GTF position rules)
awk 'BEGIN{FS="\t"; OFS="\t"} \
## think about pos strand first
$5 == "+"{my_array[FNR]=$0; tx_start[FNR]=$3; tx_end[FNR]=$4; exon_array[FNR]=$9; \
## split $9 exon pos column
n = split(exon_array[FNR], new_array, ","); \
	for (i = 1; i <= n; i++) { \
	## to do subtraction, need one more splitting
	k = split(new_array[i], new_array_sub, "-"); \
		exon_array_width[FNR]=(exon_array_width[FNR] "," (new_array_sub[2] - new_array_sub[1] +1));\
		exon_array_dist[FNR]=(exon_array_dist[FNR] "," (new_array_sub[1] - tx_start[FNR]));\
} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, tx_start[FNR], tx_end[FNR], n, \
	## just removing "," in front of the first width and exon pos numbers. 
	gensub(/^,(.+)$/, "\\1", "g", exon_array_width[FNR]), gensub(/^,(.+)$/, "\\1", "g", exon_array_dist[FNR])}}\
$5 == "-"{my_array[FNR]=$0; tx_start[FNR]=$3; tx_end[FNR]=$4; exon_array[FNR]=$9; \
## split $9 exon pos column
n = split(exon_array[FNR], new_array, ","); \
	for (i = 1; i <= n; i++) { \
	## to do subtraction, need one more splitting
	k = split(new_array[i], new_array_sub, "-"); \
		## for - strand, reverseing the order of exons at this point 
		exon_array_width[FNR]=((new_array_sub[2] - new_array_sub[1] +1) "," exon_array_width[FNR]);\
		exon_array_dist[FNR]=((new_array_sub[1] - tx_start[FNR]) "," exon_array_dist[FNR]);\
} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, tx_start[FNR], tx_end[FNR], n, \
	## removing "," from the end of the width and exon pos nums (not from the front b/c of the reversed order). 
	gensub(/^(.+),$/, "\\1", "g", exon_array_width[FNR]), gensub(/^(.+),$/, "\\1", "g", exon_array_dist[FNR])}}' \
fooM23_tx_exon_combined.gtf > fooM23_tx_exon.fakebed


## clean up and format for BED
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $3, $4, ($6 ":" $7 ":" $8), 0, $5, $11, $12, 0, $13, $14, $15}' \
fooM23_tx_exon.fakebed > gencode_vM23_tx_gtf_to_bed_int.fakebed

awk 'BEGIN{FS="\t"; OFS="\t"} { clean=gensub(/\"/, "", "g", $4); 
print $1, $2, $3, clean, $5, $6, $7, $8, $9, $10, $11, $12}' gencode_vM23_tx_gtf_to_bed_int.fakebed > gencode_vM23_tx_gtf_to_bed_final.fakebed


## if we just need a bed for usual transcripts, convert to bed now (0-based now) 
## need to take care of - strand $11 $12 where orders need to be reversed
## convert from GTF to BED (remove 1 from col2 for both pos and neg strands)
	awk 'BEGIN{FS="\t"; OFS="\t"}\
	$6 == "+"{my_array[FNR]=$0; print $1,($2 -1),$3,$4,$5,$6,($7 - 1),$8,$9,$10,$11,$12}
	#### for neg strand, tx_start and cds_start are for BED nomenclature, they are really tx end and cds end, so be mindful
	$6 == "-"{my_array[FNR]=$0; tx_start[FNR]=$2; exon_width[FNR]=$11; exon_pos[FNR]=$12; 
	### split $12 exon pos column
	n = split(exon_pos[FNR], new_array_exon_pos, ",");
	m = split(exon_width[FNR], new_array_exon_width, ",");
	for (i = 1; i <= n; i++){
	## reversing the order of exon pos and exon width
	reversed_exon_pos[FNR]=(new_array_exon_pos[i] "," reversed_exon_pos[FNR]);
	reversed_exon_width[FNR]=(new_array_exon_width[i] "," reversed_exon_width[FNR])} { \
	my_array[FNR]=($1 "\t" ($2 - 1) "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" ($7 - 1) "\t" $8 "\t" $9 "\t" $10 "\t" \
	gensub(/(.+),$/, "\\1", "g", reversed_exon_width[FNR]) "\t" gensub(/(.+),$/, "\\1", "g", reversed_exon_pos[FNR])); 
	print my_array[FNR]}}' \
	gencode_vM23_tx_gtf_to_bed_final.fakebed > gencode_vM23_tx_gtf_to_bed_regular_tx_annot_final.bed

## this is conventional vM23 mouse gtf to bed (TSS start to end, not CDS start to end)  
## gencode_vM23_tx_gtf_to_bed_regular_tx_annot_final.bed
#sort -k 1,1 -k 2,2n -k 3,3n gencode_vM23_tx_gtf_to_bed_regular_tx_annot_final.bed > gencode_vM23_tx_gtf_to_bed_regular_tx_annot_sorted_final.bed

## add start and stop position
## remove $2 and keep only start_codon and stop_codon
awk 'BEGIN{FS="\t"; OFS="\t"} $3 ~ /start_codon|stop_codon/{print $1,$3,$4,$5,$7,$9}' fooM23_1.gtf > fooM23_4.gtf

# concatenate transcripts together (with FS=" " removing "transcript_id" etc)
awk 'BEGIN{FS=" "; OFS="\t"} \
$2 ~ /start_codon|stop_codon/{my_array[FNR]=($1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $9 "\t" $17); print my_array[FNR]}' \
fooM23_4.gtf > fooM23_5.gtf

# clean up 
awk 'BEGIN{FS="\t"; OFS="\t"} { clean6=gensub(/\;|\"/, "", "g", $6); clean7=gensub(/\;|\"/, "", "g", $7); clean8=gensub(/\;|\"/, "", "g", $8); 
print $1, $2, $3, $4, $5, clean6, clean7, clean8}' fooM23_5.gtf > fooM23_6.gtf

## the original GTF file doesn't have perfect start stop codon pairs.
## therefore need to make pairs. Split a file into start_codon and stop_codon files
awk 'BEGIN{FS="\t"; OFS="\t"} $2 ~ /start_codon/{print $0}' fooM23_6.gtf > fooM23_7_start.gtf
awk 'BEGIN{FS="\t"; OFS="\t"} $2 ~ /stop_codon/{print $0}' fooM23_6.gtf > fooM23_7_stop.gtf

## remove redundant start and stop entries then sort them again.. (must be a better way of doing this...)
sort -k 8,8 fooM23_7_start.gtf > fooM23_7_start_sorted.gtf
awk 'BEGIN{FS="\t"; OFS="\t"} {print NR, $8}' fooM23_7_start_sorted.gtf | uniq -f1 -d | awk '{print $1}' > foo_delete_lines.txt
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_delete[$1]=$1; next} FNR != my_delete[FNR]{print $0}' \
foo_delete_lines.txt fooM23_7_start_sorted.gtf > fooM23_7_start_non_redundant.gtf
sort -k 1,1 -k 3,3n -k 4,4n -k 8,8 fooM23_7_start_non_redundant.gtf > fooM23_7_start_non_redundant_sorted.gtf

sort -k 8,8 fooM23_7_stop.gtf > fooM23_7_stop_sorted.gtf
awk 'BEGIN{FS="\t"; OFS="\t"} {print NR, $8}' fooM23_7_stop_sorted.gtf | uniq -f1 -d | awk '{print $1}' > foo_delete_lines2.txt
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_delete[$1]=$1; next} FNR != my_delete[FNR]{print $0}' \
foo_delete_lines2.txt fooM23_7_stop_sorted.gtf > fooM23_7_stop_non_redundant.gtf
sort -k 1,1 -k 3,3n -k 4,4n -k 8,8 fooM23_7_stop_non_redundant.gtf > fooM23_7_stop_non_redundant_sorted.gtf

## combine start and stop: if file2-derived genename exists in my_start file (which is indexed based on file1 genename), print.
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_start[$8]=$0; next} ($8 in my_start){print my_start[$8], $0}' \
fooM23_7_start_non_redundant_sorted.gtf fooM23_7_stop_non_redundant_sorted.gtf > fooM23_8_start_stop.gtf

## clean up. And to have identical column as gencode_vM23_tx_gtf_to_bed_final.bed, concatenate names.
awk '{FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$10,$11,$12,$5,($6 ":" $7 ":" $8),$16}' fooM23_8_start_stop.gtf > fooM23_8_start_stop_clean.gtf 

## add fooM23_8_start_stop_clean.gtf to gencode_vM23_tx_gtf_to_bed_final.bed
## during file1 NR, my_array is indexed based on file1 $9 values. During fil2 NR, checking if the specified index (not the actual value) exists in file1 my_array.
## The index as written is file2 $4 but the original index was defined from file1 $8. Therefore the net effect is a search for file2 $4 element in the entire file1 fields for all entries (rows) 
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_start[$9]=$0; next} ($4 in my_start){print my_start[$4], $0}' \
fooM23_8_start_stop_clean.gtf gencode_vM23_tx_gtf_to_bed_final.fakebed > fooM23_9_start_stop_bed_combo.gtf

## clean up (make it look closer to bed)
awk 'BEGIN{FS="\t"; OFS="\t"}{print $11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$2,$3,$4,$5,$6,$7}' \
fooM23_9_start_stop_bed_combo.gtf > fooM23_9_start_stop_bed_combo_clean.gtf

## calculate positions. gencode_vM23_tx_gtf_to_bed_final.fakebed is NOT 0-based. Will keep 1 based here.
## Change col2 col3 tx start and end pos to CDS start and end pos so that later one can add nt to the ends.  
## think about pos strand first. For ease of downstream analysis, will define 3rd pos of stop as the end of CDS. (and the 1st pos of start as start of CDS)

awk 'BEGIN{FS="\t"; OFS="\t"}\
$6 == "+"{my_array[FNR]=$0; cds_start[FNR]=$14; cds_stop[FNR]=$18; tx_start[FNR]=$2; exon_pos[FNR]=$12; exon_width[FNR]=$11;
n = split(exon_pos[FNR], new_array_exon_pos, ",");
m = split(exon_width[FNR], new_array_exon_width, ",");
for (i = 1; i <= n; i++){
	each_exon_start_pos[i]=(tx_start[FNR] + new_array_exon_pos[i] );
	each_exon_end_pos[i]=(tx_start[FNR] + new_array_exon_pos[i] + new_array_exon_width[i] - 1 );
	## remove exons that come before cds_start
	if ( cds_start[FNR] > each_exon_end_pos[i] ){
		#print "this exon is before cds start"
		each_exon_start_pos[i]=NULL; each_exon_end_pos[i]=NULL;}
	if ( each_exon_start_pos[i] <= cds_start[FNR] && cds_start[FNR] < each_exon_end_pos[i] ){
		#print "this is where start belongs"
		each_exon_start_pos[i]=cds_start[FNR];}
	## remove exons that come after cds_stop
	if ( each_exon_start_pos[i] <= cds_stop[FNR] && cds_stop[FNR] < each_exon_end_pos[i] ){
		#print "this is where stop belongs"
		each_exon_end_pos[i]=cds_stop[FNR];}
	if ( each_exon_start_pos[i] > cds_stop[FNR] ){
		#print "this exon is after cds ends"
		each_exon_start_pos[i]=NULL; each_exon_end_pos[i]=NULL;}

	## revese the number back with respect to cds_start CDS start not TSS even though it reads relative to TSS here
		if (length(each_exon_start_pos[i]) > 0 && length(each_exon_end_pos[i]) > 0 ){
		revised_exon_starts[FNR]=(revised_exon_starts[FNR] "," each_exon_start_pos[i]);
		revised_exon_ends[FNR]=(revised_exon_ends[FNR] "," each_exon_end_pos[i]);
		revised_exon_starts_relative_to_TSS[FNR]=(revised_exon_starts_relative_to_TSS[FNR] "," (each_exon_start_pos[i] - cds_start[FNR])  );
		revised_exon_ends_relative_to_TSS[FNR]=(revised_exon_ends_relative_to_TSS[FNR] "," (each_exon_end_pos[i] - cds_start[FNR])  );
		revised_exon_width[FNR]=(revised_exon_width[FNR] "," (each_exon_end_pos[i] - each_exon_start_pos[i] + 1 ));
		}}
	 }
#### for neg strand, tx_start and cds_start are for BED nomenclature, they are really tx end and cds end, so be mindful
$6 == "-"{my_array[FNR]=$0; cds_start[FNR]=$15; cds_stop[FNR]=$17; tx_start[FNR]=$2; exon_pos[FNR]=$12; exon_width[FNR]=$11;
### split $12 exon pos column
n = split(exon_pos[FNR], new_array_exon_pos, ",");
m = split(exon_width[FNR], new_array_exon_width, ",");
for (i = 1; i <= n; i++){
	## cds_stop for reverse strand!
	each_exon_start_pos[i]=(tx_start[FNR] + new_array_exon_pos[i]  );
	each_exon_end_pos[i]=(tx_start[FNR] + new_array_exon_pos[i] + new_array_exon_width[i] - 1 );

	## remove exons that come before cds_start (really is after cds pos as neg strand)
	if ( cds_start[FNR] < each_exon_start_pos[i] ){
		#print "this exon is before cds start";
		each_exon_start_pos[i]=NULL; each_exon_end_pos[i]=NULL;}
	if ( each_exon_start_pos[i] < cds_start[FNR] && cds_start[FNR] <= each_exon_end_pos[i] ){
		#print "this is where start belongs"
		# note each_exon_end_pos represents cds_start
		each_exon_end_pos[i]=cds_start[FNR];}
	## remove exons that come after cds_stop
	if ( each_exon_start_pos[i] < cds_stop[FNR] && cds_stop[FNR] <= each_exon_end_pos[i] ){
		#print "this is where stop belongs"
		# note each_exon_start_pos represents cds_stop
		each_exon_start_pos[i]=cds_stop[FNR];}
	if ( each_exon_end_pos[i] < cds_stop[FNR] ){
		#print "this exon is after cds ends"
		each_exon_start_pos[i]=NULL; each_exon_end_pos[i]=NULL;}
		
		# revese the number back with respect to cds_stop. (not TSS even though it reads relative to TSS here!) 
		# it is still 1-based now. but the same for 0 based. 
		# this relative-to-TSS is really relative to the end-of-cds for neg strand just in keeping with BED nomenclature here. 
		if (length(each_exon_start_pos[i]) > 0 && length(each_exon_end_pos[i]) > 0 ){
		revised_exon_starts[FNR]=(revised_exon_starts[FNR] "," each_exon_start_pos[i]);
		revised_exon_ends[FNR]=(revised_exon_ends[FNR] "," each_exon_end_pos[i]);
		revised_exon_starts_relative_to_TSS[FNR]=((each_exon_start_pos[i] - cds_stop[FNR]) "," revised_exon_starts_relative_to_TSS[FNR]);
		revised_exon_ends_relative_to_TSS[FNR]=((each_exon_end_pos[i] - cds_stop[FNR]  ) "," revised_exon_ends_relative_to_TSS[FNR]);
		revised_exon_width[FNR]=((each_exon_end_pos[i] - each_exon_start_pos[i]  + 1 ) "," revised_exon_width[FNR] )
		}}
} {print $0, revised_exon_width[FNR], revised_exon_starts_relative_to_TSS[FNR];}' fooM23_9_start_stop_bed_combo_clean.gtf > fooM23_10_CDS.gtf


## sanity check 
awk 'BEGIN{FS="\t"; OFS="\t"} $4 ~ /Chrng-201/{print FNR}' fooM23_9_start_stop_bed_combo_clean.gtf | \
##xargs -I '{}' sed -n '{}{N;N;N;N;p}' fooM23_9_start_stop_bed_combo_clean.gtf

awk 'BEGIN{FS="\t"; OFS="\t"} $4 ~ /Chrng-201/{print FNR}' fooM23_10_CDS.gtf | \
xargs -I '{}' sed -n '{}{N;N;N;N;p}' fooM23_10_CDS.gtf

# clean up column 19 and column 20, then replace columns 11 and 12 with the current column 19, 20
# this is for pos strand
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,\
gensub(/^,(.+)$/, "\\1", "g", $19), gensub(/^,(.+)$/, "\\1", "g", $20), $13,$14,$15,$16,$17,$18}' fooM23_10_CDS.gtf > fooM23_11_CDS_int.fakebed

# this is for neg strand
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,\
gensub(/(.+),$/, "\\1", "g", $11), gensub(/(.+),$/, "\\1", "g", $12), $13,$14,$15,$16,$17,$18}' fooM23_11_CDS_int.fakebed > fooM23_11_CDS_int2.fakebed

# also fix column 10 (counts that are based on the current exon counts in 11 and 12
awk 'BEGIN{FS="\t"; OFS="\t"}{my_array11[FNR]=$11; my_array12[FNR]=$12; \
n = split(my_array11[FNR], new_my_array11, ","); m = split(my_array12[FNR], new_my_array12, ",");
print $1,$2,$3,$4,$5,$6,$7,$8,$9,n,$11,$12, $13,$14,$15,$16,$17,$18}' fooM23_11_CDS_int2.fakebed > gencode_vM23_CDS_gtf_to_bed_final.fakebed

## fix start and end pos for CDS
## if + strand, col2 and col 7 = start pos1 (col 14) and col3 col8 = end pos3 (col18)
## if - strand, col2 and col7 = end pos 1 (col17) and col3 col8 = start pos 3 (col15)
awk 'BEGIN{FS="\t"; OFS="\t"} \
$6 == "+"{cds_start[FNR]=$14; cds_stop[FNR]=$18} \
$6 == "-"{cds_start[FNR]=$17; cds_stop[FNR]=$15} \
{print $1,cds_start[FNR],cds_stop[FNR],$4,$5,$6,cds_start[FNR],cds_stop[FNR],$9,$10,$11,$12}' \
gencode_vM23_CDS_gtf_to_bed_final.fakebed > gencode_vM23_CDS_gtf_to_bed_final_ends_fixed.fakebed

# sanity check 
awk 'BEGIN{FS="\t"; OFS="\t"} $4 ~ /Chrng-201/{print FNR}' gencode_vM23_CDS_gtf_to_bed_final_ends_fixed.fakebed | \
xargs -I '{}' sed -n '{}{N;N;N;N;p}' gencode_vM23_CDS_gtf_to_bed_final_ends_fixed.fakebed

## next convert from GTF to BED (just remove 1 from col2 for both pos and neg strands!
awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,($2 -1),$3,$4,$5,$6,($7 - 1),$8,$9,$10,$11,$12}' \
gencode_vM23_CDS_gtf_to_bed_final_ends_fixed.fakebed > gencode_vM23_CDS_gtf_to_bed_final_minus_1.fakebed

## this is the final BED for GTF derived BED for annotated CDS (i.e., start to stop)
## gencode_vM23_CDS_gtf_to_bed_final_minus_1.fakebed (it's not fake, it's real BED

## Bedtools to extract coding sequences from reference genome
## wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/GRCm38.primary_assembly.genome.fa.gz
## gunzip GRCm38.primary_assembly.genome.fa.gz

## transcribe
bedtools getfasta -name -s -split -fo gencode_vM23_CDS_GRCm38fa_fakebed.fasta \
-fi GRCm38.primary_assembly.genome.fa -bed gencode_vM23_CDS_gtf_to_bed_final_minus_1.fakebed

## faTrans to translate from DNA to AA
##wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faTrans
##chmod u+x faTrans
./faTrans gencode_vM23_CDS_GRCm38fa_fakebed.fasta \
gencode_vM23_CDS_GRCm38fa_fakebed_translated.fasta

# sanity check
awk 'BEGIN{FS="\t"; OFS="\t"} $4 ~ /Chrng-201/{print FNR}' gencode_vM23_CDS_gtf_to_bed_final_minus_1.fakebed | \
xargs -I '{}' sed -n '{}{N;N;N;N;p}' gencode_vM23_CDS_gtf_to_bed_final_minus_1.fakebed

awk 'BEGIN{FS="\t"; OFS="\t"} $0 ~ /Chrng-201/{print FNR}' gencode_vM23_CDS_GRCm38fa_fakebed.fasta | \
xargs -I '{}' sed -n '{}{N;N;N;N;p}' gencode_vM23_CDS_GRCm38fa_fakebed.fasta

awk 'BEGIN{FS="\t"; OFS="\t"} $0 ~ /Chrng-201/{print FNR}' gencode_vM23_CDS_GRCm38fa_fakebed_translated.fasta | \
xargs -I '{}' sed -n '{}{N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;N;p}' gencode_vM23_CDS_GRCm38fa_fakebed_translated.fasta

## TGA stop codon
## 3nt to the end is already added (stop codon). First, characterize start and stop codon usage as a sanity check
awk 'FNR % 2 == 1{print $0} FNR % 2 == 0{n=substr($0, (length($0)-2), 3); print n}' \
gencode_vM23_CDS_GRCm38fa_fakebed.fasta > stop_codon_only.fasta
awk 'FNR % 2 == 1{print $0} FNR % 2 == 0{n=substr($0, 1, 3); print n}' \
gencode_vM23_CDS_GRCm38fa_fakebed.fasta > start_codon_only.fasta

wc -l stop_codon_only.fasta 
awk 'BEGIN{FS="\t"; OFS="\t"} FNR % 2 == 0{if ($1 == "TGA"){print $1}}' stop_codon_only.fasta | wc -l
awk 'BEGIN{FS="\t"; OFS="\t"} FNR % 2 == 0{if ($1 == "TAA"){print $1}}' stop_codon_only.fasta | wc -l
awk 'BEGIN{FS="\t"; OFS="\t"} FNR % 2 == 0{if ($1 == "TAG"){print $1}}' stop_codon_only.fasta | wc -l
awk 'BEGIN{FS="\t"; OFS="\t"} FNR % 2 == 0{if ($1 == "ATG"){print $1}}' start_codon_only.fasta | wc -l

## select TGA stop entries only. then based on that TGA list, select the corresponding BED lines. 
awk 'BEGIN{FS="\t"; OFS="\t"} {my_array[FNR]=$0} \
FNR % 2 == 0{if ($1 == "TGA"){print my_array[(FNR - 1)]; print $1;}}' stop_codon_only.fasta > TGA_stop_codon_only.fasta

## find TGA stop codon only BED entries
## first, for clarity, clean up TGA_stop_codon_only.fasta file for the ease of searching. 
awk 'BEGIN{FS="\t"; OFS="\t"} FNR % 2 == 1{n=split($0, subName, ":");
print (gensub(/^>(.+)/, "\\1", "g", subName[1]) "\t" subName[2] "\t" subName[3])}' \
TGA_stop_codon_only.fasta > TGA_ID_query

## for clarity, splitting gencode_vM23_CDS_gtf_to_bed_final.bed $4 into 3 parts also (though really excessive)
awk 'BEGIN{FS="\t"; OFS="\t"}{n=split($4, subName, ":");
print $1,$2,$3,subName[1],subName[2],subName[3],$5,$6,$7,$8,$9,$10,$11,$12}' gencode_vM23_CDS_gtf_to_bed_final_minus_1.fakebed > TGA_subject.bed

# searching for file2 $6 element in the entire file1 fields for all entries (rows) 
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_array[$3]=$0; next} ($6 in my_array) {print $0}' TGA_ID_query TGA_subject.bed > TGA_only.bed

# clean up and put it back to real BED format.
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$2,$3,($4 ":" $5 ":" $6), $7,$8,$9,$10,$11,$12,$13,$14}' \
TGA_only.bed > gencode_vM23_CDS_gtf_to_bed_TGA_only_final.bed

## next, add 45 nt to the TGA only bed.
## use bedtools flank with -s option (stranded) and -r (to add nt to the end coordinate)
## but first prepare GENOME for bedtools (and download GRCm38.primary_assembly.genome.fa)
#module load samtools
#samtools faidx GRCm38.primary_assembly.genome.fa
#cut -f 1,2 GRCm38.primary_assembly.genome.fa.fai > chrom.sizes

## adding 45nt to the end pos (-r 45) accounting for strandedness 
bedtools flank -i gencode_vM23_CDS_gtf_to_bed_TGA_only_final.bed \
-g chrom.sizes -l 0 -r 45 -s > gencode_vM23_CDS_gtf_to_bed_TGA_only_final_plus_45nt.bed

## then re-extract coding sequences
## first need to reformat the bedfile (the columns related to exons need to be simplified to 1 exon to make it work properly)
## also fix the $11 (width) and $12 (pos) columns
awk 'BEGIN{FS="\t"; OFS="\t"} {$10=1; $11=45; $12=0; print}' gencode_vM23_CDS_gtf_to_bed_TGA_only_final_plus_45nt.bed > plus_45nt_int.bed

## getfasta to transcribe the 45nt
bedtools getfasta -name -s -split -fo gencode_vM23_CDS_gtf_to_bed_TGA_only_final_plus_45nt.fasta \
-fi GRCm38.primary_assembly.genome.fa -bed plus_45nt_int.bed

## do also getfasta to transcribe the main CDS for TGA only bed. 
bedtools getfasta -name -s -split -fo gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS.fasta \
-fi GRCm38.primary_assembly.genome.fa -bed gencode_vM23_CDS_gtf_to_bed_TGA_only_final.bed

## add the 45nt to the end of each fasta
awk 'BEGIN{FS="\t"; OFS="\t"} FNR==NR{my_array[FNR]=$0; next} FNR % 2 == 1{n=split($0, new_array, ":"); \
print (my_array[FNR] new_array[6])} FNR %2 == 0{print (my_array[FNR] "---" $0)}' \
gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS.fasta gencode_vM23_CDS_gtf_to_bed_TGA_only_final_plus_45nt.fasta > \
gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL.fasta

## sanity check
awk 'BEGIN{FS="\t"; OFS="\t"} $0 ~ /Chrng-201/{print FNR}' gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL.fasta | \
xargs -I '{}' sed -n '{}{N;N;N;N;p}' gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL.fasta

## faTrans to translate into amino acid 
./faTrans gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL.fasta \
gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL_translated.fasta

## flattening the fasta 
awk 'FNR==1{print $0; next} $0 !~ /^>/{printf "%s", $0} $0 ~ /^>/{printf "\n%s\n", $0} END{printf "\n"}' \
gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL_translated.fasta > gencode_vM23_CDS_gtf_to_bed_TGA_only_final_CDS_plus_45nt_FINAL_translated_flat.fasta










