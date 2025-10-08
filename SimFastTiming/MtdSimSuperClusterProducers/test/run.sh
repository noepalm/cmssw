# 1. topological + historical clustering
rm mtdSimSuperClusters_topo+history_numEvent1000.root
cmsRun runMtdSimSuperClusterProducer.py \
       outputFile=mtdSimSuperClusters_topo+history.root \
       useTopologicalClustering=True \
       maxEvents=1000 &> log_topo+history

# 2. historical clustering only
rm mtdSimSuperClusters_history_numEvent1000.root
cmsRun runMtdSimSuperClusterProducer.py \
       outputFile=mtdSimSuperClusters_history.root \
       useTopologicalClustering=False \
       maxEvents=1000 &> log_historyOnly