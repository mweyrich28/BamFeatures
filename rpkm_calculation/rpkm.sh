#!/usr/bin/sh
# bash rpkm.sh <sample>.annot <sample>.gene_lengths.txt

annot=$1
gene_lengths=$2

base_name=$(basename "$annot" | cut -d'.' -f1)

mkdir -p "$base_name/Plots"
mkdir -p "$base_name/RPKM"

grep -vE "^\d*dist:|\d*antisense|split-in" "$annot" | cut -f1,6,7 > tmp.reads
echo "LOG: Filtered annot"

python3 format_rpkm.py tmp.reads > split.reads
echo "LOG: Formatted annot"

sort -k2,2 split.reads > all_$annot.reads
echo "LOG: Sorted annot"

# rm tmp.reads
# rm split.reads

grep -E "pcrindex: 0" all_$annot.reads > pcr_$annot.reads
echo "LOG: Filtered PCR uniq reads into separate file"

total_reads_all=$(cut -f1 all_${annot}.reads | sort | wc -l)
total_reads_pcr=$(cut -f1 pcr_${annot}.reads | sort | wc -l)
echo "LOG: Extracted total read counts"

python3 calc_rpkm.py all_$annot.reads "$gene_lengths" "$total_reads_all" > "$base_name/RPKM/all_$base_name.rpkm"
python3 calc_rpkm.py pcr_$annot.reads "$gene_lengths" "$total_reads_pcr" > "$base_name/RPKM/pcr_$base_name.rpkm"

# rm all_$annot.reads
# rm pcr_$annot.reads

echo "LOG: Calculation done"
               
Rscript plot.R "$base_name/RPKM/all_$base_name.rpkm" "$base_name/RPKM/pcr_$base_name.rpkm" "$base_name/Plots"
