python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_291125_bugfix.root \
        --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_291125_bugfix \
        --config configs/plot_config_forDev271125.py
python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_291125_bugfix.root \
        --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_291125_bugfix \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_281125_bugfix.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_281125_bugfix \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_281125_bugfix.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_281125_bugfix \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_271125_replay.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125_replay \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_271125_replay.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125_replay \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_271125_bugfix.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125_bugfix \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_271125_bugfix.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125_bugfix \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_291125_bugfix.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_291125_bugfix \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_291125_bugfix.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_291125_bugfix \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_291125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_291125 \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_291125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_291125 \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_281125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_281125 \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_281125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_281125 \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_271125_bis.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125_bis \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_271125_bis.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125_bis \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_271125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125 \
#         --config configs/plot_config_forDev271125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_271125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_271125 \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_261125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_261125 \
#         --config configs/plot_config_forDev261125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_261125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_261125 \

# python3 plot_script_from_tree.py ../test/tree_vali_dev_forDev_251125.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_251125 \
#         --config configs/plot_config_forDev251125.py
# python3 plot_complicated_plots.py ../test/tree_vali_dev_forDev_251125.root \
        # --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forDev_251125 \

# python3 plot_script_from_tree.py ../test/tree_vali_dev.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/forPresentation_251024/photonGun \
#         --config configs/plot_config_forPresentation251024.py

# ===================================================================


# # cd ../../; scram b -j 16; cd -

# # if --plot argument passed, only run plotting scripts
# if [ "$1" == "--plot" ]; then
#     echo "Running only plotting scripts..."
# fi

# # 1: topological + historical clustering
# if [ "$1" != "--plot" ]; then
#     rm mtd_mergedcluster_validation_topo+history.root;
#     cmsRun ../test/mtd_simmergedcluster_vali_cfg.py \
#            inputFile=file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_mergedclusterDev/src/SimFastTiming/MtdSimMergedClusterProducers/test/mtdSimMergedClusters_topo+history_numEvent1000.root \
#            outputFile=mtd_mergedcluster_validation_topo+history.root &> log_topo+history
# fi

# python3 plot_script_from_tree.py mtd_mergedcluster_validation_topo+history.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests \
#         --config configs/plot_config_topo+history.py

# # 2: historical clustering only
# if [ "$1" != "--plot" ]; then
#     rm mtd_mergedcluster_validation_history.root;
#     cmsRun ../test/mtd_simmergedcluster_vali_cfg.py \
#            inputFile=file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_mergedclusterDev/src/SimFastTiming/MtdSimMergedClusterProducers/test/mtdSimMergedClusters_history_numEvent1000.root \
#            outputFile=mtd_mergedcluster_validation_history.root &> log_history
# fi

# python3 plot_script_from_tree.py mtd_mergedcluster_validation_history.root \
#         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests/history_only \
#         --config configs/plot_config_history.py

# # # 3: Create cluster maps for events with topological merging
# # python3 plot_map.py \
# #         --history-file mtd_mergedcluster_validation_history.root \
# #         --topo-file mtd_mergedcluster_validation_topo+history.root \
# #         --output-dir /eos/home-n/npalmeri/www/MTD/MergedCluster/sim_tests \
# #         --max-events 200
