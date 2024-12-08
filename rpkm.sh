#!/usr/bin/sh

annot=$1
gene_lengths=$2

grep -vE "^\d*dist:|\d*antisense|split-in" $1 | cut -f1,6,7 > tmp.reads

# split up multi gene lines
python3 format_rpkm.py tmp.reads > split.reads
# sort into new file
sort -k2,2 split.reads > all_$1.reads
rm tmp.reads
rm split.reads
# create uniq file accounting for pcr index
grep -E "pcrindex: 0" all_$1.reads > pcr_$1.reads
total_reads_all=$(cat all_$1.reads | wc -l)
total_reads_pcr=$(cat pcr_$1.reads | wc -l)

python3 calc_rpkm.py all_$1.reads $2 $total_reads_all > all_$1.rpkm
python3 calc_rpkm.py pcr_$1.reads $2 $total_reads_pcr > pcr_$1.rpkm

# clean up
rm all_$1.reads
rm pcr_$1.reads
