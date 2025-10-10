#include "Validation/MtdValidation/plugins/SuperClusterValidationExample.h"
#include <iostream>

SuperClusterValidationExample::SuperClusterValidationExample(const edm::ParameterSet& iConfig) {
    usesResource("TFileService");
    
    superClustersToken_ = consumes<FTLSuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters"));
    clustersToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("clusters"));
    mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>(); 
    mtdgeoToken_ = esConsumes<MTDGeometry, MTDDigiGeometryRecord>();

    totalAdjacentPairs_ = 0;
    totalMergedPairs_ = 0;
}

void SuperClusterValidationExample::beginJob() {
    edm::Service<TFileService> fs;
    
    // SuperCluster histograms
    h_sc_energy_ = fs->make<TH1F>("h_sc_energy", "SuperCluster Energy;Energy [MeV];Count", 100, 0, 50);
    h_sc_time_ = fs->make<TH1F>("h_sc_time", "SuperCluster Time;Time [ns];Count", 100, -5, 20);
    h_sc_timeError_ = fs->make<TH1F>("h_sc_timeError", "SuperCluster Time Error;Time Error [ns];Count", 100, 0, 1);
    h_sc_x_ = fs->make<TH1F>("h_sc_x", "SuperCluster X;X [mm];Count", 100, -200, 200);
    h_sc_y_ = fs->make<TH1F>("h_sc_y", "SuperCluster Y;Y [mm];Count", 100, -200, 200);
    h_sc_nClusters_ = fs->make<TH1F>("h_sc_nClusters", "Number of Clusters;N_{clusters};Count", 10, 0.5, 10.5);
    
    // Cluster histograms for comparison
    h_cluster_energy_ = fs->make<TH1F>("h_cluster_energy", "Cluster Energy;Energy [MeV];Count", 100, 0, 50);
    h_cluster_time_ = fs->make<TH1F>("h_cluster_time", "Cluster Time;Time [ns];Count", 100, -5, 20);
    
    // 2D histograms
    h_sc_energy_vs_time_ = fs->make<TH2F>("h_sc_energy_vs_time", "SuperCluster Energy vs Time;Time [ns];Energy [MeV]", 100, -5, 20, 100, 0, 50);
    h_sc_xy_ = fs->make<TH2F>("h_sc_xy", "SuperCluster Position;X [mm];Y [mm]", 100, -200, 200, 100, -200, 200);
    h_sc_energy_vs_nClusters_ = fs->make<TH2F>("h_sc_energy_vs_nClusters_", "Energy vs N Clusters;N_{clusters};Energy [MeV]", 10, 0.5, 10.5, 100, 0, 50);
    h_merging_efficiency_ = fs->make<TH2F>("h_merging_efficiency", "Merging Efficiency;Cluster Energy [MeV];Merged?", 50, 0, 50, 2, -0.5, 1.5);
    
    // efficiency
    h_eta_adjacent_pairs_ = fs->make<TH1F>("h_eta_adjacent_pairs", "Adjacent Pairs vs Eta;#eta;Number of Adjacent Pairs", 100, -1.5, 1.5);
    h_eta_merged_pairs_ = fs->make<TH1F>("h_eta_merged_pairs", "Merged Pairs vs Eta;#eta;Number of Merged Pairs", 100, -1.5, 1.5);
    h_eta_merging_efficiency_ = fs->make<TH1F>("h_eta_merging_efficiency", "Merging Efficiency vs Eta;#eta;Efficiency", 100, -1.5, 1.5);
    
    // Analysis tree
    tree_ = fs->make<TTree>("MTDSuperClusters", "MTD SuperCluster Analysis Tree");
    
    tree_->Branch("evt_run", &evt_run_);
    tree_->Branch("evt_event", &evt_event_);
    
    tree_->Branch("sc_n", &sc_n_);
    tree_->Branch("sc_energy", &sc_energy_);
    tree_->Branch("sc_time", &sc_time_);
    tree_->Branch("sc_timeError", &sc_timeError_);
    tree_->Branch("sc_x", &sc_x_);
    tree_->Branch("sc_y", &sc_y_);
    tree_->Branch("sc_nClusters", &sc_nClusters_);
    tree_->Branch("sc_seedId", &sc_seedId_);
    tree_->Branch("sc_clusterIds", &sc_clusterIds_);

    tree_->Branch("cluster_n", &cluster_n_);
    tree_->Branch("cluster_energy", &cluster_energy_);
    tree_->Branch("cluster_time", &cluster_time_);
    tree_->Branch("cluster_x", &cluster_x_);
    tree_->Branch("cluster_y", &cluster_y_);
    tree_->Branch("cluster_detId", &cluster_detId_);
    tree_->Branch("cluster_inSuperCluster", &cluster_inSuperCluster_);
}

void SuperClusterValidationExample::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<FTLSuperClusterCollection> superClustersHandle;
    iEvent.getByToken(superClustersToken_, superClustersHandle);
    
    edm::Handle<FTLClusterCollection> clustersHandle;
    iEvent.getByToken(clustersToken_, clustersHandle);
    
    if (!superClustersHandle.isValid() || !clustersHandle.isValid()) {
        std::cout << "Invalid handles!" << std::endl;
        return;
    }
    
    auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
    const MTDTopology* topology = topologyHandle.product();

    auto geomHandle = iSetup.getTransientHandle(mtdgeoToken_);
    const MTDGeometry* geom = geomHandle.product();

    evt_run_ = iEvent.id().run();
    evt_event_ = iEvent.id().event();
    
    sc_energy_.clear(); sc_time_.clear(); sc_timeError_.clear();
    sc_x_.clear(); sc_y_.clear(); sc_nClusters_.clear();
    sc_seedId_.clear(); sc_clusterIds_.clear();
    cluster_energy_.clear(); cluster_time_.clear();
    cluster_x_.clear(); cluster_y_.clear(); cluster_detId_.clear();
    cluster_inSuperCluster_.clear();
    
    sc_n_ = superClustersHandle->size();
    
    std::set<uint32_t> clusterDetIds;

    std::map<uint32_t, int> detIdToSCIndex;
    int scIndex = 0;
    
    for (const auto& sc : *superClustersHandle) {

        h_sc_energy_->Fill(sc.energy());
        h_sc_time_->Fill(sc.time());
        h_sc_timeError_->Fill(sc.timeError());
        h_sc_x_->Fill(sc.x());
        h_sc_y_->Fill(sc.y());
        h_sc_nClusters_->Fill(sc.nClusters());

        h_sc_energy_vs_time_->Fill(sc.time(), sc.energy());
        h_sc_xy_->Fill(sc.x(), sc.y());
        h_sc_energy_vs_nClusters_->Fill(sc.nClusters(), sc.energy());

        sc_energy_.push_back(sc.energy());
        sc_time_.push_back(sc.time());
        sc_timeError_.push_back(sc.timeError());
        sc_x_.push_back(sc.x());
        sc_y_.push_back(sc.y());
        sc_nClusters_.push_back(sc.nClusters());
        sc_seedId_.push_back(sc.id().rawId());
        
        std::vector<uint32_t> clusterIds;
        for (const auto& detId : sc.clusterIds()) {
            clusterIds.push_back(detId.rawId());
            clusterDetIds.insert(detId.rawId());
            detIdToSCIndex[detId.rawId()] = scIndex;
        }
        sc_clusterIds_.push_back(clusterIds);
        scIndex++;

        if (evt_event_ <= 3) {
            std::cout << "SuperCluster: E=" << sc.energy() << " MeV, t=" << sc.time() 
                      << " ns, pos=(" << sc.x() << "," << sc.y() << ")" << std::endl;
            std::cout << "  Made from " << sc.nClusters() << " clusters:" << std::endl;
            for (const auto& detId : sc.clusterIds()) {
                std::cout << "    DetId: " << detId.rawId() << std::endl;
            }
        }
    }
    
    std::map<BTLDetId, const FTLCluster*> clusterMap;
    std::vector<FTLCluster> validClusters;
    cluster_n_ = 0;

    for (const auto& detSet : *clustersHandle) {
        for (const auto& cluster : detSet) {
            if (cluster.energy() < 1.0) continue; 
            
            cluster_n_++;
            validClusters.push_back(cluster);
            clusterMap[cluster.id()] = &cluster;

            h_cluster_energy_->Fill(cluster.energy());
            h_cluster_time_->Fill(cluster.time());
            
            // Check if cluster was merged
            bool inSuperCluster = clusterDetIds.count(cluster.id().rawId()) > 0;
            //h_merging_efficiency_->Fill(cluster.energy(), inSuperCluster ? 1 : 0);
            
            cluster_energy_.push_back(cluster.energy());
            cluster_time_.push_back(cluster.time());
            cluster_x_.push_back(cluster.x());
            cluster_y_.push_back(cluster.y());
            cluster_detId_.push_back(cluster.id().rawId());
            cluster_inSuperCluster_.push_back(inSuperCluster);
        }
    }
    
    std::set<std::pair<uint32_t, uint32_t>> processedPairs;

    int adjacentPairs = 0;
    int mergedPairs = 0;
    
    if (evt_event_ <= 3) {
        std::cout << "Checking " << validClusters.size() << " clusters for adjacency:" << std::endl;
    }

    //WIP
    for (const auto& cluster : validClusters) {
        BTLDetId cluId = cluster.id();
        
        // get eta of the cluster
        const MTDGeomDet* thedet = geom->idToDet(cluId);
        if (!thedet) continue;
        GlobalPoint global_point = thedet->surface().toGlobal(LocalPoint(cluster.x(), cluster.y(), 0));
        double cluster_eta = global_point.eta();
        
        // topology indices
        std::pair<uint32_t, uint32_t> indices = topology->btlIndex(cluId.rawId());
        uint32_t iphi = indices.first;
        uint32_t ieta = indices.second;
        
        if (iphi == std::numeric_limits<uint32_t>::max() || ieta == std::numeric_limits<uint32_t>::max()) {
            continue;
        }
        
        // Check edge hits for eta direction
        bool edgeHitIn0 = false;
        bool edgeHitIn15 = false;
        for (int i = 0; i < cluster.size(); ++i) {
            auto hit = cluster.hit(i);
            int hit_col = hit.y();
            if (hit_col == 0) edgeHitIn0 = true;
            else if (hit_col == 15) edgeHitIn15 = true;
        }
        bool hasEdgeHit = edgeHitIn0 || edgeHitIn15;

        //eta
        if (hasEdgeHit) {
            std::vector<int> etaOffsets = {1, -1};
            for (int etaOffset : etaOffsets) {
                uint32_t adjDetIdRaw = topology->btlidFromIndex(iphi, ieta + etaOffset);
                if (adjDetIdRaw == 0) continue;
                
                BTLDetId adjDetId(adjDetIdRaw);
                auto it = clusterMap.find(adjDetId);
                if (it == clusterMap.end()) continue;
                
                const FTLCluster* adjCluster = it->second;
                
                // Check for opposite edge hit in neighbor
                bool hasOppositeEdgeHit = false;
                for (int j = 0; j < adjCluster->size(); ++j) {
                    auto hit = adjCluster->hit(j);
                    int hit_col = hit.y();
                    if ((edgeHitIn0 && hit_col == 15) || (edgeHitIn15 && hit_col == 0)) {
                        hasOppositeEdgeHit = true;
                        break;
                    }
                }

                if (hasOppositeEdgeHit) {

                    std::pair<uint32_t, uint32_t> pairKey = 
                    (cluId.rawId() < adjDetId.rawId()) ? 
                    std::make_pair(cluId.rawId(), adjDetId.rawId()) :
                    std::make_pair(adjDetId.rawId(), cluId.rawId());
                
                    // Skip if already processed
                    if (processedPairs.find(pairKey) != processedPairs.end()) {
                        continue;
                    }
                    processedPairs.insert(pairKey);
                    
                    adjacentPairs++;
                    h_eta_adjacent_pairs_->Fill(cluster_eta);

                    // Check if both clusters are in the same SuperCluster
                    uint32_t detId1 = cluId.rawId();
                    uint32_t detId2 = adjDetId.rawId();
                    
                    auto it1 = detIdToSCIndex.find(detId1);
                    auto it2 = detIdToSCIndex.find(detId2);
                    
                    bool bothInSC = (it1 != detIdToSCIndex.end()) && (it2 != detIdToSCIndex.end());
                    bool sameIndexMerged = bothInSC && (it1->second == it2->second);
                    
                    if (sameIndexMerged) {
                        mergedPairs++;

                        h_eta_merged_pairs_->Fill(cluster_eta);
                    }

                    if (evt_event_ <= 3 && adjacentPairs <= 10) {
                        std::cout << "  ETA adjacent pair " << adjacentPairs << ": " 
                                  << detId1 << " <-> " << detId2 
                                  << ", bothInSC=" << bothInSC
                                  << ", sameIndex=" << sameIndexMerged;
                        if (bothInSC) {
                            std::cout << " (SC indices: " << it1->second << "," << it2->second << ")";
                        }
                        std::cout << std::endl;
                    }
                }
            }
        }

        //phi
        /*std::vector<int> phiOffsets = {1, -1};
        for (int phiOffset : phiOffsets) {
            uint32_t adjDetIdRaw = topology->btlidFromIndex(iphi + phiOffset, ieta);
            if (adjDetIdRaw == 0) continue;
            
            BTLDetId adjDetId(adjDetIdRaw);
            auto it = clusterMap.find(adjDetId);
            if (it == clusterMap.end()) continue;

            std::pair<uint32_t, uint32_t> pairKey = 
                (cluId.rawId() < adjDetId.rawId()) ? 
                std::make_pair(cluId.rawId(), adjDetId.rawId()) :
                std::make_pair(adjDetId.rawId(), cluId.rawId());
            
            // Skip if already processed
            if (processedPairs.find(pairKey) != processedPairs.end()) {
                continue;
            }
            processedPairs.insert(pairKey);
            
            adjacentPairs++;

            h_eta_adjacent_pairs_->Fill(cluster_eta);

            // Check if both clusters are in the same SuperCluster
            uint32_t detId1 = cluId.rawId();
            uint32_t detId2 = adjDetId.rawId();
            
            auto it1 = detIdToSCIndex.find(detId1);
            auto it2 = detIdToSCIndex.find(detId2);
            
            bool bothInSC = (it1 != detIdToSCIndex.end()) && (it2 != detIdToSCIndex.end());
            bool sameIndexMerged = bothInSC && (it1->second == it2->second);
            
            if (sameIndexMerged) {
                mergedPairs++;

                h_eta_merged_pairs_->Fill(cluster_eta);
            }
            
            if (evt_event_ <= 3 && adjacentPairs <= 10) {
                std::cout << "  PHI adjacent pair " << adjacentPairs << ": " 
                          << detId1 << " <-> " << detId2 
                          << ", bothInSC=" << bothInSC
                          << ", sameIndex=" << sameIndexMerged;
                if (bothInSC) {
                    std::cout << " (SC indices: " << it1->second << "," << it2->second << ")";
                }
                std::cout << std::endl;
            }
        }*/
    }

    totalAdjacentPairs_ += adjacentPairs;
    totalMergedPairs_ += mergedPairs;

    tree_->Fill();

    if (evt_event_ <= 3) {
        std::cout << "Event " << evt_event_ << ": " << sc_n_ << " SuperClusters, " 
                  << cluster_n_ << " individual clusters" << std::endl;
        std::cout << "Adjacent pairs: " << adjacentPairs 
                  << ", Merged pairs: " << mergedPairs << std::endl;
    }
}

void SuperClusterValidationExample::endJob() {

    h_eta_merging_efficiency_->Divide(h_eta_merged_pairs_, h_eta_adjacent_pairs_, 1, 1, "B");

    std::cout << "\n=== MTD SuperCluster Validation Summary ===" << std::endl;
    std::cout << "Total events processed: " << tree_->GetEntries() << std::endl;
    std::cout << "Average SuperClusters per event: " << h_sc_energy_->GetEntries() / tree_->GetEntries() << std::endl;
    std::cout << "Total adjacent cluster pairs found: " << totalAdjacentPairs_ << std::endl;
    std::cout << "Total merged pairs: " << totalMergedPairs_ << std::endl;
    if (totalAdjacentPairs_ > 0) {
        std::cout << "Pair merging efficiency: " << (double)totalMergedPairs_ / totalAdjacentPairs_ * 100 << "%" << std::endl;
    }
    //std::cout << "Merging efficiency: " << h_merging_efficiency_->GetBinContent(2, 2) / h_merging_efficiency_->GetEntries() * 100 << "%" << std::endl;
}

DEFINE_FWK_MODULE(SuperClusterValidationExample);