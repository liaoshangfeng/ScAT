path=/data/work/previous/ScAT_github
scat=${path}/core/bin/ScAT.py
echo $scat

# gem_file=${path}/data/FP200000495BL_A4.bin1.Lasso.gem.gz
gem_file=/data/work/previous/ScAT_github/data/FP200000553BL_F2.bin1.Lasso.gem.gz
prefix=FP200000495BL_A4

###################00 QC Control#############################
# step 1 Preview the data and obtain raw_seurat_object.rds
python3 $scat QcControl -h 

# python3 $scat QcControl \
#     --Rscript $gem_file \
    --input $gem_file \
    --out ${path}/00.QC/step1 \
    --binsize 50 \
    --prefix $prefix \
    --mtPattern '^Mt' \
    --rbPattern '^Rp[sl]' \
    --hbPattern '^Hb[^(p)]' \
    --platPattern 'Pecam1|Pf4' \
    --ptsize 50

###################00 QC Control#############################


# # step 2 Perform filteraion and normalization 
# # Note! You need to manually set the filter parameters based on the QC metrics generated in the previous step.
# # ~2.2G
# python3 $scat QcControl \
#     -i ${path}/00.QC/step1/${prefix}_raw_seuratObject.rds \
#     -o ${path}/00.QC/step2 \
#     --prefix $prefix \
#     --mFilter \
#     --mtPattern '^Mt' \
#     --rbPattern '^Rp[sl]' \
#     --hbPattern '^Hb[^(p)]' \
#     --platPattern 'Pecam1|Pf4' \
#     --minCell 3 \
#     --ptsize 50

# ###################00 QC Control#############################
