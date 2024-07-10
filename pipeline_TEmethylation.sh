#!/usr/bin/env bash
set -e
#set -o xtrace
#Check parameters
if [[ $# -ne 3 ]]; then
	echo -e "\n\tMissing parameters:"
	echo -e "\te.g. ./pipeline_TEmethylation.sh repeatmasker.out genes.gff CG.sorted.bed \n"
	exit 2
fi
#Load parameters
rm_out=$1
genes=$2
methylation_CG=$3
prop_cutoff=0.6
prefix=$(basename $rm_out | cut -f1 -d'.')
bold_start="\033[1m"
bold_end="\033[0m"

echo -e "\n\t${bold_start}Starting pipeline:${bold_end}"

#File with TE lengths from library is needed
#Can be created with: grep "^ID"  ./RepeatMasker/Libraries/RMRBMeta.embl | awk 'OFS="\t" {print $2,$6}' > RMRBMeta.lengths

### 1) TEs: Clean the RMasker output and filter the results.
echo -e "\t${bold_start}1) TEs:${bold_end}"

#Check which type of repeats there are with: grep "\b[0-9]" $rm_out |  awk '{print $11}' | sort | uniq -c

#1.1 Transform the output of RMasker into a more readable table and adding more info
out=$prefix".table"
echo -e "\t- 1.1 Ordering columns, more readable than RM output -> "$out
grep "\bChr[0-9]*" $rm_out | egrep -v "Satellite|Simple_repeat|Low_complexity|ARTEFACT" | sed 's/[()]//g' \
	| awk -v OFS="\t" 'BEGIN{print "#chromo","startQ","endQ","strand","startS","endS", \
	"map_len","id","class","repeat"}; {if ($9=="+") print $5,$6,$7,$9,$12,$13,$7-$6+1, \
	$15,$11,$10; else print $5,$6,$7,"-",$14,$13,$7-$6+1, \
	$15,$11,$10}' | awk 'OFS="\t" {sub("DNA","TIR",$9); sub("RC","Helitron",$9);print}' > $out

#1.2 Merge complete LTRs and remove incomplete ones, also merge splitted TEs
out=$prefix".table.merge"
echo -e "\t- 1.2 Merge LTRs and splitted TEs and remove incomplete LTRs -> "$out
echo -e "#chromo\tstartQ\tendQ\tstrand\tstartS\tendS\tmap_len\t\
rpt_len\tp_aln\tid\tclass\trepeat" > $out

python3 rm_merger.py $prefix".table" >> $out 

#1.3 Filter RMasker output removing unwanted TE classes and non.mergeable splitted TEs
out=$prefix".table.merge.filter"
echo -e "\t- 1.3 Filter RM output, removes nonTEs, splitted TEs and some unwanted TEs -> "$out
awk 'NR==FNR{cnt[$10]++; next} cnt[$10]==1' $prefix".table.merge" $prefix".table.merge" \
 | egrep -v "RNA|SINE|Unknown|Cassandra|Retroposon|Stowaway|ATLINE1_3A" | sed 's/LTR\/Caulimovirus/Caulimovirus/' > $out

#1.4 Filter out those that are more than certain proportion of complete in length with the reference and make a gff
out_TEs=$prefix".table.merge.filter.gt06.gff"
echo -e "\t- 1.4 Filtering by $prop_cutoff of length covered -> "$out_TEs
grep "^Chr" $prefix".table.merge.filter" | awk	'$9>=0.6 {print $1, "RepeatMasker", "dispersed_repeat", $2, $3, ".", $4,".","ID="$10";name="$12";class="$11}' | sed 's/ /\t/g' > $out_TEs 
sed -i 's/\bCaulimovirus /LTR\/Caulimovirus/' $out_TEs

### 2) TEs with methylation
echo -e "\n\t${bold_start}2) TEs and methylation:${bold_end}"

#2.1 Intersect CG coordinates with all TEs (with length gt06)
out_CG=$prefix".TEs.intersectCG.bed"
echo -e "\t- 2.1 Intersect methylation files with TEs"
echo -e "\t      - Intersecting CG with TEs -> "$out_CG
bedtools intersect -wb -a $methylation_CG -b $out_TEs \
	| awk 'OFS="\t" {print $1,$2,$3,$20,$11,$6,$18}' | awk '$6==$7' > $out_CG

#2.2 Calculating average CG methylation for each TE
out_CG_avg=$prefix".TEs.intersectCG.averageCG"
echo -e "\t- 2.2 Calculate the average methylation score for TEs (slow step)"
echo -e "\t     - Calculating average CG methylation score -> "$out_CG_avg 
for i in $(cut -f4 $out_CG | sort -k1V | uniq ); do
	grep "${i}" $out_CG | awk -v OFS="\t" \
	'{sum+=$5}END{print $1, NR, sprintf("%.1f",sum/NR), $4}'
done > $out_CG_avg

#2.3 Select only TEs with 5 or more CG
out_CG_gt4=$prefix".TEs.intersectCG.averageCG.gt4"
echo -e "\t- 2.3 Extract only TEs with 5 or more CG"
echo -e "\t      - Extract from CG average table -> "$out_CG_gt4
awk '$2>4' $prefix".TEs.intersectCG.averageCG" > $out_CG_gt4

out_TEs_gt4=$prefix".table.merge.filter.gt06.gt4.gff"
echo -e "\t      - Extract from TE annotation -> "$out_TEs_gt4
awk '$2>4' $prefix".TEs.intersectCG.averageCG" | cut -f4 > temp_list
grep -f temp_list $out_TEs > $out_TEs_gt4


#2.3 Merging tables
out_avg_all=$prefix".intersectCG.averageCG.gt4.table"
echo -e "\t- 2.3 Merging CG methylation scores and TEs annotation -> "$out_avg_all
paste $out_TEs_gt4 $out_CG_gt4 | awk '{print $1,$4,$5,$7,$11,$12,$13}'\
	| sed 's/;/ /g;s/name=//;s/class=//;s/ /\t/g' > $out_avg_all

### 3) TEs and genes. Get the TEs and genes that overlap 
echo -e "\n\t${bold_start}3) TEs and genes:${bold_end} "

#3.1 Extend coordinates 2K in 5' direction
out=$prefix".genes.extend2K5p.gff"
echo -e "\t- 3.1 Extend 2 Kbp upstream the genes coordinates -> "$out
awk 'OFS="\t" {if ($7=="+") print $1,$2,$3,$4-2000,$4,$6,$7,$8,$9; \
	else print $1,$2,$3,$5,$5+2000,$6,$7,$8,$9}' $genes \
	| awk 'OFS="\t" {if($4<1){$4=1} print}' > $out

#3.3 Intersect genes+2K5p and selected TEs
out_intersect=$prefix".TEs.intersect2K5p_genes.gff"
echo -e "\t- 3.2 Intersect the +2K5p extended genes with the TEs -> "$out_intersect
bedtools intersect -b $prefix".genes.extend2K5p.gff" -a $out_TEs_gt4 -wa \
	| uniq > $out_intersect

#3.4 Get average methylation score for selected TEs
out_avg_intersect=$prefix".TEs.intersect2K5p.averageCG"
echo -e "\t- 3.3 Get methylation scores for TEs that intersect with +2K5p genes ->" $out_avg_intersect
cut -f9 $out_intersect | cut -f1 -d';' > temp_list
grep -wf temp_list $out_avg_all > $out_avg_intersect

echo -e "\n\t${bold_start}Finishing pipeline :)${bold_end}\n"


