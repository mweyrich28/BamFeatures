#!/usr/bin/sh

# bash rpkm.sh <sample>.annot <sample>.gene_lengths.txt

annot=$1
gene_lengths=$2

grep -vE "^\d*dist:|\d*antisense|split-in" $1 | cut -f1,6,7 > tmp.reads
echo "Filtered annot"
# split up multi gene lines
python3 format_rpkm.py tmp.reads > split.reads
echo "Formatted annot"
# sort into new file
sort -k2,2 split.reads > all_$1.reads
echo "Sorted annot"
rm tmp.reads
rm split.reads
# create uniq file accounting for pcr index
grep -E "pcrindex: 0" all_$1.reads > pcr_$1.reads
echo "Extracted PCR reads"
total_reads_all=$(cat all_$1.reads | wc -l)
total_reads_pcr=$(cat pcr_$1.reads | wc -l)
echo "Extracted Total read counts"

python3 calc_rpkm.py all_$1.reads $2 $total_reads_all > all_$1.rpkm
python3 calc_rpkm.py pcr_$1.reads $2 $total_reads_pcr > pcr_$1.rpkm
echo "DONE"

# clean up
rm -r all_$1.reads
rm -r pcr_$1.reads
