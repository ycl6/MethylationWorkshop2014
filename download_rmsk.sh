#!/bin/bash
# Download RepMask 3.2.7 track for hg18

echo "Fetching RepMask 3.2.7 files."
for i in {1..22}
do
	echo "chr${i}"
	wget -q -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/chr${i}_rmskRM327.txt.gz | gunzip | awk -F $'\t' 'BEGIN { OFS=FS } { print $6,$7,$8,$12,0,$10,$13,$11 }' > chr${i}_rmskRM327.bed
done

echo "chrX"
wget -q -O - http://hgdownload.soe.ucsc.edu/goldenPath/hg18/database/chrX_rmskRM327.txt.gz | gunzip | awk -F $'\t' 'BEGIN { OFS=FS } { print $6,$7,$8,$12,0,$10,$13,$11 }'  > chrX_rmskRM327.bed

cat *_rmskRM327.bed | sortBed > hg18.rmskRM327.bed
rm *_rmskRM327.bed

echo "Complete."
