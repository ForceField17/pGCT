cat HM450.hg38.manifest.gencode.v36.tsv | cut -f 1,2,3,5 | awk '($1~/^chr/)' > temp
sort-bed temp > anno_hg38.txt
bedtools map -a anno_hg38.txt -b /mnt/lab/home1/dsongad/analysis/HiFi_CpG/Element/methylation/methy_bed/200_1kb/L1_promoter_200_1kb_first500bp.bed -c 1,2,3,5 -o collapse | awk '($5!=".")' | cut -f 4,5,6,7,8 | sort > 450k_L1pro_list_2026Aug.bed

echo "A_CpG_ID" > L1pro_CpG_ID_2026Aug
cut -f 1 450k_L1pro_list_2026Aug.bed | sort >> L1pro_CpG_ID_2026Aug

