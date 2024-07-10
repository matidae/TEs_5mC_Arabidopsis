#!/usr/bin/env bash
./pipeline_TEmethylation.sh ts/ait/data/AIT.cs.fasta.out ts/ait/data/ait.c.new.genes_only.gff ts/ait/data/CG.sorted.bed
./stats.sh AIT.TEs.intersect2K5p
python plots_stats.py AIT stats.intersect2K5p
python plots_stats.py AIT stats.intersect2K5p/fam
mkdir ait
mv AIT.* ait

./pipeline_TEmethylation.sh ts/brr/data/BRR.cs.fasta.out ts/brr/data/brr.c.new.genes_only.gff ts/brr/data/CG.sorted.bed
./stats.sh BRR.TEs.intersect2K5p
python plots_stats.py BRR stats.intersect2K5p
python plots_stats.py BRR stats.intersect2K5p/fam
mkdir brr
mv BRR.* brr

./pipeline_TEmethylation.sh ts/dol/data/DOL.fs.fasta.out ts/dol/data/dol.f.new.genes_only.gff ts/dol/data/CG.sorted.bed
./stats.sh DOL.TEs.intersect2K5p
python plots_stats.py DOL stats.intersect2K5p/
python plots_stats.py DOL stats.intersect2K5p/fam
mkdir dol
mv DOL.* dol

./pipeline_TEmethylation.sh ts/oy0/data/oy0.f.fasta.out ts/oy0/data/oy0.f.new.genes_only.gff ts/oy0/data/CG.sorted.bed
./stats.sh oy0.TEs.intersect2K5p
python plots_stats.py oy0 stats.intersect2K5p/
python plots_stats.py oy0 stats.intersect2K5p/fam
mkdir oy0
mv oy0.* oy0

./pipeline_TEmethylation.sh ts/pun/data/PUN.cs.fasta.out ts/pun/data/pun.c.new.genes_only.gff ts/pun/data/CG.sorted.bed
./stats.sh PUN.TEs.intersect2K5p
python plots_stats.py PUN stats.intersect2K5p/
python plots_stats.py PUN stats.intersect2K5p/fam
mkdir pun
mv PUN.* pun

./pipeline_TEmethylation.sh ts/rsh/data/RSH.fs.fasta.out ts/rsh/data/rsh.f.new.genes_only.gff ts/rsh/data/CG.sorted.bed
./stats.sh RSH.TEs.intersect2K5p
python plots_stats.py RSH stats.intersect2K5p/
python plots_stats.py RSH stats.intersect2K5p/fam
mkdir rsh
mv RSH.* rsh

./pipeline_TEmethylation.sh ts/tad/data/tad.c.fasta.out ts/tad/data/tad.c.new.genes_only.gff ts/tad/data/CG.sorted.bed
./stats.sh tad.TEs.intersect2K5p
python plots_stats.py tad stats.intersect2K5p/
python plots_stats.py tad stats.intersect2K5p/fam
mkdir tad
mv tad.* tad

./pipeline_TEmethylation.sh ts/mit/data/MIT.fs.fasta.out ts/mit/data/mit.f.new.genes_only.gff ts/mit/data/mit.flye.call_mods.CG.bed
./stats.sh MIT.TEs.intersect2K5p
python plots_stats.py MIT stats.intersect2K5p/
python plots_stats.py MIT stats.intersect2K5p/fam
mkdir mit
mv MIT.* mit


