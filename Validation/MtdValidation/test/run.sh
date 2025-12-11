rm tree_vali_dev_mergetest.root
cmsRun mtd_mergedcluster_vali_dev_cfg.py \
       outputFile=tree_vali_dev_mergetest.root \
       inputFile="file:/eos/home-n/npalmeri/MTD/MTD_supercluster/association_maps_merge/src/RecoLocalFastTime/FTLClusterizer/test/mtdMergedClusters_mergetest_numEvent1000.root" &> vali_dev_mergetest.log