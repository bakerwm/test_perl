#!/bin/bash

### test the update_table.pl script:

# 1.add position (col-7 to 12) to tail of bed
echo "1.add position (col-7 to 12) to tail of bed"
echo "origin:"
head -n 2  srna_vs_known.bed

echo "output:"
perl update_table.pl -add -m 4 -d srna.txt -A 7 -B 12 srna_vs_known.bed |head -n 2
echo ""

# 2.add position (col-7 to 12) after col-6 of bed
echo "2.add position (col-7 to 12) after col-6 of bed"
echo "origin:"
head -n 2  srna_vs_known.bed

echo "output:"
perl update_table.pl -add -m 4 -a 6 -d srna.txt -A 7 -B 12 srna_vs_known.bed |head -n 2
echo ""

# 3.replace gene locus by gene name
echo "3.replace gene locus by gene name"
echo "origin:"
head -n 2 known.txt

echo "output:"
perl update_table.pl -replace -m 7 -a 7 -b 7 -d dict.ptt -n 6 -A 5 -B 5 known.txt | head -n 2
perl update_table.pl -replace -m 7 -a 9 -b 9 -d dict.ptt -n 6 -A 5 -B 5 known.txt | head -n 2
echo ""

# 4. add gene names after col-7 of sort file
echo "4.add gene names after col-7 of sort file"
echo "origin:"
head -n 2 known.txt 

echo "output:"
perl update_table.pl -add -m 7 -a 7 -d dict.ptt -n 6 -A 5 -B 5 known.txt |head -n 2
echo ""

# 5. add gene names after col-7 and col-9 of sort file (be caution about the new width)
echo "5.add gene names after col-7 and col-9 of sort file (be caution about the new width)"
echo "origin:"
head -n 2 known.txt 

echo "output:"
perl update_table.pl -add -m 7 -a 7 -d dict.ptt -n 6 -A 5 -B 5 known.txt | perl tmp.update_table.pl -add -m 10 -a 10 -d dict.ptt -n 6 -A 5 -B 5 | head -n3
echo ""

# 6. extract records in db
echo "6.extract records in db"
echo "origin:"
head -n 2 srna_vs_known.bed 

echo "output:"
perl update_table.pl -line -m 4 -d srna.txt srna_vs_known.bed | head -n 2
echo ""

echo "END"
