# 1. topological + historical clustering
rm mtdSimMergedClusters_topo+history_numEvent1000.root
cmsRun runMtdSimMergedClusterProducer.py \
       outputFile=mtdSimMergedClusters_topo+history.root \
       useTopologicalClustering=True \
       maxEvents=1000 &> log_topo+history

# 2. historical clustering only
rm mtdSimMergedClusters_history_numEvent1000.root
cmsRun runMtdSimMergedClusterProducer.py \
       outputFile=mtdSimMergedClusters_history.root \
       useTopologicalClustering=False \
       maxEvents=1000 &> log_historyOnly