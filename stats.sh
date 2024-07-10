#!/usr/bin/env bash

#The prefix for the data eg. 'mitt'
if [ -z "$1" ]; then
    echo "Genome prefix must be passed as an argument. e.g ./stats.sh mitt.TEs.intersectALL"
    exit
fi
prefix=$1

########################
# TEs divided by order #
########################
outdir=$(echo $prefix | rev | cut -f1 -d'.' | rev)
mkdir -p ./stats.$outdir

#Counts of all TEs
counts_all(){
    cut -f9 $1".gff"| cut -d';' -f3| cut -f2 -d'='| cut -f1 -d'/' \
        | sed 's/?//' | sort | uniq -c | sort -k2 > aux0
}    

#Function to avoid having null/empty counts 
check_counts(){
    aux=$1
    for i in $(awk '{print $2}' aux0); do 
        a=$(grep -w $i $aux| awk '{print $1}')
        if [ -z ${a} ]; then 
       	    echo "0" $i
        else 
            echo $a $i 
        fi
    done > $aux"a"
}

#Counts of methylation for all families
counts(){
    grep -v "#" $1".averageCG" | awk '$6>80' \
        | cut -f9 | cut -f1 -d'/' | sed 's/?//' | sort | uniq -c | sort -k2 > aux1
    counts_all "$1"
    check_counts aux1
    grep -v "#" $1".averageCG" | awk '$6>20 && $6<=80' \
        | cut -f9 | cut -f1 -d'/' | sed 's/?//' | sort | uniq -c | sort -k2 > aux2
    check_counts aux2
    grep -v "#" $1".averageCG" | awk '$6<=20' \
        | cut -f9 | cut -f1 -d'/' | sed 's/?//' | sort | uniq -c | sort -k2 > aux3
    check_counts aux3
    paste aux1a aux2a aux3a | awk '{print $1,$3,$5,$2}' > "./stats."$outdir"/"$1".counts_CG"
}

out_prefix=$prefix
pre=$(echo $prefix | cut -f1 -d'.')
ln -sf $pre".table.merge.filter.gt06.gt4.gff" $out_prefix".gff"
counts "$out_prefix"
rm aux*

#########################
# TEs divided by family #
########################
mkdir -p ./stats.$outdir/fam
## Files for the count plot
#Counts of all TEs
counts_all(){
    cut -f9 $1.gff| cut -d';' -f3| cut -f2 -d'=' | sed 's/?//' \
        | grep "/" | sed 's/-.*//' | sort |uniq -c | sort -k2 | sed 's/\//./' > aux0
}

#Counts of methylation for all families
counts_fam(){
    grep -v "#" $1".averageCG" | awk '$6>80' \
        | cut -f9 | grep  "/" | sed 's/-.*//' | sort | uniq -c | sort -k2 | sed 's/\//./' > aux1
    counts_all "$1"
    check_counts aux1
    grep -v "#" $1".averageCG" | awk '$6>20 && $6<=80' \
        | cut -f9 | grep  "/" | sed 's/-.*//' | sort | uniq -c | sort -k2 | sed 's/\//./' > aux2
    check_counts aux2
    grep -v "#" $1".averageCG" | awk '$6<=20' \
        | cut -f9 | grep  "/" | sed 's/-.*//' | sort | uniq -c | sort -k2 | sed 's/\//./' > aux3
    check_counts aux3
    paste aux1a aux2a aux3a | awk '{print $1,$3,$5,$2}' > "./stats."$outdir"/fam/"$1".counts_fam_CG"
}

out_prefix=$prefix
pre=$(echo $prefix | cut -f1 -d'.')
counts_fam "$out_prefix"
rm aux*

### Files for methylation histogram
#CG
values(){
    fams=($2)

    n=0
    for i in ${fams[@]}; do
        n=$((n + 1))
        touch "aux"$n
        grep $i $1.averageCG | cut -f6 > "aux"$n
    done
    if $3; then
        echo $2|sed 's/\//\./g;s/ /\t/g' > "./stats."$outdir"/"$1".values_CG"
        list_aux=$(ls aux* | sort -V)
        paste $list_aux >> "./stats."$outdir"/"$1".values_CG"
    else
        echo $2|sed 's/\//\./g;s/ /\t/g' > "./stats."$outdir"/fam/"$1".values_fam_CG"
        list_aux=$(ls aux* | sort -V)
        paste $list_aux >> "./stats."$outdir"/fam/"$1".values_fam_CG"
    fi
    rm aux*
}

order=$(awk '{print $9}' $prefix.averageCG | grep -v class | grep "/" \
	| cut -f1 -d'/' | sort | uniq -c | sort -k1nr |awk '{print $2}')
out_prefix=$prefix
values "$out_prefix" "$order" true

fami=$(awk '{print $9}' $prefix.averageCG | grep -v class | grep "/" \
	| cut -f2 -d'/' | sed 's/-.*//' | sort | uniq -c | sort -k1nr |awk '{print $2}')
out_prefix=$prefix
values "$out_prefix" "$fami" false

