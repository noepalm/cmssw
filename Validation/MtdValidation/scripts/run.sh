# cd ../../; scram b -j 16; cd -

# if --plot argument passed, only run plotting scripts
if [ "$1" == "--plot" ]; then
    echo "Running only plotting scripts..."
fi

# 1: topological + historical clustering
if [ "$1" != "--plot" ]; then
    rm mtd_mergedcluster_validation_topo+history.root;
    cmsRun ../test/mtd_simmergedcluster_vali_cfg.py \
           inputFile=file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_mergedclusterDev/src/SimFastTiming/MtdSimMergedClusterProducers/test/mtdSimMergedClusters_topo+history_numEvent1000.root \
           outputFile=mtd_mergedcluster_validation_topo+history.root &> log_topo+history
fi

python3 plot_script_from_tree.py mtd_mergedcluster_validation_topo+history.root \
        --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests \
        --config configs/plot_config_topo+history.py

# 2: historical clustering only
if [ "$1" != "--plot" ]; then
    rm mtd_mergedcluster_validation_history.root;
    cmsRun ../test/mtd_simmergedcluster_vali_cfg.py \
           inputFile=file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_mergedclusterDev/src/SimFastTiming/MtdSimMergedClusterProducers/test/mtdSimMergedClusters_history_numEvent1000.root \
           outputFile=mtd_mergedcluster_validation_history.root &> log_history
fi

python3 plot_script_from_tree.py mtd_mergedcluster_validation_history.root \
        --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests/history_only \
        --config configs/plot_config_history.py

# # 3: Create cluster maps for events with topological merging
# python3 plot_map.py \
#         --history-file mtd_mergedcluster_validation_history.root \
#         --topo-file mtd_mergedcluster_validation_topo+history.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests \
#         --max-events 200