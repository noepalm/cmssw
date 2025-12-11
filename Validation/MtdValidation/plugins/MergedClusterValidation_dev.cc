#include "Validation/MtdValidation/plugins/MergedClusterValidation_dev.h"
#include <iostream>
#include <CLHEP/Units/SystemOfUnits.h>
#include "DataFormats/Math/interface/GeantUnits.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#define DEBUG 0

MergedClusterValidation_dev::MergedClusterValidation_dev(const edm::ParameterSet& iConfig):
    mtdgeoToken_(esConsumes<MTDGeometry, MTDDigiGeometryRecord>()),
    mtdtopoToken_(esConsumes<MTDTopology, MTDTopologyRcd>()) {
    usesResource("TFileService");
    
    mergedClustersToken_ = consumes<FTLMergedClusterCollection>(iConfig.getParameter<edm::InputTag>("mergedClusters"));
    clustersToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("clusters"));

    simMergedClustersToken_ = consumes<MtdSimMergedClusterCollection>(iConfig.getParameter<edm::InputTag>("simMergedClusters"));
    simClustersToken_ = consumes<MtdSimLayerClusterCollection>(iConfig.getParameter<edm::InputTag>("simLayerClusters"));

    totalAdjacentPairs_ = 0;
    totalMergedPairs_ = 0;
}

void MergedClusterValidation_dev::beginJob() {
    edm::Service<TFileService> fs;
    
    // -------------------- //
    // ------- RECO ------- //
    // -------------------- //

    // MergedCluster histograms
    h_mc_energy_ = fs->make<TH1F>("h_mc_energy", "MergedCluster Energy;Energy [MeV];Count", 100, 0, 50);
    h_mc_time_ = fs->make<TH1F>("h_mc_time", "MergedCluster Time;Time [ns];Count", 100, -5, 20);
    h_mc_timeError_ = fs->make<TH1F>("h_mc_timeError", "MergedCluster Time Error;Time Error [ns];Count", 100, 0, 1);
    h_mc_x_ = fs->make<TH1F>("h_mc_x", "MergedCluster X;X [mm];Count", 100, -200, 200);
    h_mc_y_ = fs->make<TH1F>("h_mc_y", "MergedCluster Y;Y [mm];Count", 100, -200, 200);
    h_mc_nClusters_ = fs->make<TH1F>("h_mc_nClusters", "Number of Clusters;N_{clusters};Count", 10, 0.5, 10.5);
    
    // Cluster histograms for comparison
    h_cluster_energy_ = fs->make<TH1F>("h_cluster_energy", "Cluster Energy;Energy [MeV];Count", 100, 0, 50);
    h_cluster_time_ = fs->make<TH1F>("h_cluster_time", "Cluster Time;Time [ns];Count", 100, -5, 20);
    
    // 2D histograms
    h_mc_energy_vs_time_ = fs->make<TH2F>("h_mc_energy_vs_time", "MergedCluster Energy vs Time;Time [ns];Energy [MeV]", 100, -5, 20, 100, 0, 50);
    h_mc_xy_ = fs->make<TH2F>("h_mc_xy", "MergedCluster Position;X [mm];Y [mm]", 100, -200, 200, 100, -200, 200);
    h_mc_energy_vs_nClusters_ = fs->make<TH2F>("h_mc_energy_vs_nClusters_", "Energy vs N Clusters;N_{clusters};Energy [MeV]", 10, 0.5, 10.5, 100, 0, 50);
    h_merging_efficiency_ = fs->make<TH2F>("h_merging_efficiency", "Merging Efficiency;Cluster Energy [MeV];Merged?", 50, 0, 50, 2, -0.5, 1.5);
    
    // efficiency
    h_eta_adjacent_pairs_ = fs->make<TH1F>("h_eta_adjacent_pairs", "Adjacent Pairs vs Eta;#eta;Number of Adjacent Pairs", 100, -1.5, 1.5);
    h_eta_merged_pairs_ = fs->make<TH1F>("h_eta_merged_pairs", "Merged Pairs vs Eta;#eta;Number of Merged Pairs", 100, -1.5, 1.5);
    h_eta_merging_efficiency_ = fs->make<TH1F>("h_eta_merging_efficiency", "Merging Efficiency vs Eta;#eta;Efficiency", 100, -1.5, 1.5);
    
    // ------------------- //
    // ------- SIM ------- //
    // ------------------- //

    // 1D histograms
    h_simmc_energy_ = fs->make<TH1F>("h_simmc_energy", "MergedCluster Energy;Energy [MeV];Count", 50, 0, 50);
    h_simmc_logEnergy_ = fs->make<TH1F>("h_simmc_logEnergy", "MergedCluster Log(Energy);Log(Energy [MeV]);Count", 100, -3, 3);
    h_simmc_time_ = fs->make<TH1F>("h_simmc_time", "MergedCluster Time;Time [ns];Count", 50, 0, 30);
    h_simmc_x_ = fs->make<TH1F>("h_simmc_x", "MergedCluster X;X [mm];Count", 50, -120, 120);
    h_simmc_y_ = fs->make<TH1F>("h_simmc_y", "MergedCluster Y;Y [mm];Count", 50, -120, 120);
    h_simmc_eta_ = fs->make<TH1F>("h_simmc_eta", "MergedCluster Eta;#eta;Count", 50, -5, 5);
    h_simmc_nClusters_ = fs->make<TH1F>("h_simmc_nClusters", "Number of Clusters;N_{clusters};Count", 11, -0.5, 10.5);
    h_simmc_n_ = fs->make<TH1F>("h_simmc_n", "Number of MergedClusters;N_{MergedClusters};Count", 51, -0.5, 50.5);

    // 1D histograms -- per simLC 
    h_simmc_logEnergy_perCluster_ = fs->make<TH1F>("h_simmc_logEnergy_perCluster", "Cluster Log(Energy);Log(Energy [MeV]);Count", 100, -3, 3);
    h_simmc_time_perCluster_ = fs->make<TH1F>("h_simmc_time_perCluster", "Cluster Time;Time [ns];Count", 50, 0, 30);
    h_simmc_clusterType_ = fs->make<TH1F>("h_simmc_clusterType", "Cluster Type;Type;Count", 4, -0.5, 3.5);
    
    // 2D histograms
    h_simmc_xy_ = fs->make<TH2F>("h_simmc_xy", "MergedCluster Position;X [mm];Y [mm]", 100, -120, 120, 100, -120, 120);
    h_simmc_energy_vs_time_ = fs->make<TH2F>("h_simmc_energy_vs_time", "MergedCluster Energy vs Time;Time [ns];Energy [MeV]", 30, 4, 27, 20, 0, 1);
    h_simmc_energy_vs_nClusters_ = fs->make<TH2F>("h_simmc_energy_vs_nClusters", "Energy vs N Clusters;N_{clusters};Energy [MeV]", 4, -0.5, 3.5, 20, 0, 1);
    h_simmc_primaryPt_vs_nClusters_ = fs->make<TH2F>("h_simmc_primaryPt_vs_nClusters", "Primary Particle pT vs N Clusters;N_{clusters};p_{T} [GeV]", 4, -0.5, 3.5, 20, 0, 20);
    h_simmc_primaryPt_vs_energy_ = fs->make<TH2F>("h_simmc_primaryPt_vs_energy", "Primary Particle pT vs MergedCluster Energy;MergedCluster Energy [MeV];p_{T} [GeV]", 20, 0, 1, 20, 0, 11);
    h_simmc_primaryEnergy_vs_energy_ = fs->make<TH2F>("h_simmc_primaryEnergy_vs_energy", "Primary Particle Energy vs MergedCluster Energy;MergedCluster Energy [MeV];Primary Particle Energy [GeV]", 20, 0, 1, 20, 0, 40);
    h_simmc_primaryEnergy_vs_nClusters_ = fs->make<TH2F>("h_simmc_primaryEnergy_vs_nClusters", "Primary Particle Energy vs N Clusters;N_{clusters};Primary Particle Energy [GeV]", 4, -0.5, 3.5, 20, 0, 40);

    // Analysis tree
    tree_ = fs->make<TTree>("MTDMergedClusters", "MTD MergedCluster Analysis Tree");
    
    tree_->Branch("evt_run", &evt_run_);
    tree_->Branch("evt_event", &evt_event_);
    
    // RECO
    tree_->Branch("mc_n", &mc_n_);
    tree_->Branch("mc_energy", &mc_energy_);
    tree_->Branch("mc_time", &mc_time_);
    tree_->Branch("mc_timeError", &mc_timeError_);
    tree_->Branch("mc_x", &mc_x_);
    tree_->Branch("mc_y", &mc_y_);
    tree_->Branch("mc_nClusters", &mc_nClusters_);
    tree_->Branch("mc_seedId", &mc_seedId_);
    tree_->Branch("mc_clusterIds", &mc_clusterIds_);

    tree_->Branch("cluster_n", &cluster_n_);
    tree_->Branch("cluster_energy", &cluster_energy_);
    tree_->Branch("cluster_time", &cluster_time_);
    tree_->Branch("cluster_x", &cluster_x_);
    tree_->Branch("cluster_y", &cluster_y_);
    tree_->Branch("cluster_detId", &cluster_detId_);
    tree_->Branch("cluster_inMergedCluster", &cluster_inMergedCluster_);

    // SIM
    tree_->Branch("simmc_n", &simmc_n_);
    tree_->Branch("simmc_energy", &simmc_energy_);
    tree_->Branch("simmc_time", &simmc_time_);
    // tree_->Branch("simmc_timeError", &simmc_timeError_);
    tree_->Branch("simmc_x", &simmc_x_);
    tree_->Branch("simmc_y", &simmc_y_);
    tree_->Branch("simmc_eta", &simmc_eta_);
    tree_->Branch("simmc_nClusters", &simmc_nClusters_);
    tree_->Branch("simmc_nModules", &simmc_nModules_);
    tree_->Branch("simmc_iphi_perCluster", &simmc_iphi_perCluster_);
    tree_->Branch("simmc_ieta_perCluster", &simmc_ieta_perCluster_);
    tree_->Branch("simmc_energy_perCluster", &simmc_energy_perCluster_);
    tree_->Branch("simmc_time_perCluster", &simmc_time_perCluster_);
    tree_->Branch("simmc_earliestHitTime_perCluster", &simmc_earliestHitTime_perCluster_);
    tree_->Branch("simmc_hitCols_perCluster", &simmc_hitCols_perCluster_);
    tree_->Branch("simmc_clusterType", &simmc_clusterType_);
    tree_->Branch("simmc_primary_energy", &simmc_primary_energy_);
    tree_->Branch("simmc_primary_et", &simmc_primary_et_);
    tree_->Branch("simmc_primary_phi", &simmc_primary_phi_);
    tree_->Branch("simmc_primary_eta", &simmc_primary_eta_);
    tree_->Branch("simmc_primary_pdgId", &simmc_primary_pdgId_);    
}

void MergedClusterValidation_dev::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace geant_units::operators; // for energy conversion

    edm::Handle<FTLMergedClusterCollection> mergedClustersHandle;
    iEvent.getByToken(mergedClustersToken_, mergedClustersHandle);
    
    edm::Handle<FTLClusterCollection> clustersHandle;
    iEvent.getByToken(clustersToken_, clustersHandle);

    edm::Handle<MtdSimMergedClusterCollection> simMergedClustersHandle;
    iEvent.getByToken(simMergedClustersToken_, simMergedClustersHandle);

    edm::Handle<MtdSimLayerClusterCollection> mtdSimLCHandle;
    iEvent.getByToken(simClustersToken_, mtdSimLCHandle);

    std::cout << "Analyzing event " << iEvent.id().event() << " in run " << iEvent.id().run() << std::endl;
    std::cout << "Number of MTD MergedClusters: " << mergedClustersHandle->size() << std::endl;
    std::cout << "Number of MTD Clusters: " << clustersHandle->size() << std::endl;
    std::cout << "Number of MTD SimMergedClusters: " << simMergedClustersHandle->size() << std::endl;
    std::cout << "Number of MTD SimLayerClusters: " << mtdSimLCHandle->size() << std::endl;

    // count BTL SimLayerClusters
    int nBTLSimLCs = 0;
    for (const auto& simLC : *mtdSimLCHandle) {
        BTLDetId detId(simLC.detIds_and_rows()[0].first);
        if (detId.mtdSubDetector() == MTDDetId::BTL) {
            nBTLSimLCs++;
        }
    }
    std::cout << "Number of BTL SimLayerClusters: " << nBTLSimLCs << std::endl;

    
    if (!mergedClustersHandle.isValid() || !clustersHandle.isValid() || !simMergedClustersHandle.isValid() || !mtdSimLCHandle.isValid()) {
        std::cout << "Invalid handles!" << std::endl;
        std::cout << "  mergedClustersHandle: " << mergedClustersHandle.isValid() << std::endl;
        std::cout << "  clustersHandle: " << clustersHandle.isValid() << std::endl;
        std::cout << "  simMergedClustersHandle: " << simMergedClustersHandle.isValid() << std::endl;
        std::cout << "  mtdSimLCHandle: " << mtdSimLCHandle.isValid() << std::endl;
        return;
    }
    
    auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
    const MTDTopology* topology = topologyHandle.product();

    auto geomHandle = iSetup.getTransientHandle(mtdgeoToken_);
    const MTDGeometry* geom = geomHandle.product();

    evt_run_ = iEvent.id().run();
    evt_event_ = iEvent.id().event();
    
    mc_energy_.clear(); mc_time_.clear(); mc_timeError_.clear();
    mc_x_.clear(); mc_y_.clear(); mc_nClusters_.clear();
    mc_seedId_.clear(); mc_clusterIds_.clear();
    cluster_energy_.clear(); cluster_time_.clear();
    cluster_x_.clear(); cluster_y_.clear(); cluster_detId_.clear();
    cluster_inMergedCluster_.clear();

    simmc_energy_.clear(); simmc_time_.clear();
    simmc_x_.clear(); simmc_y_.clear(); simmc_eta_.clear(); simmc_nClusters_.clear(); simmc_nModules_.clear();
    simmc_iphi_perCluster_.clear(); simmc_ieta_perCluster_.clear(); 
    simmc_energy_perCluster_.clear(); simmc_time_perCluster_.clear(); 
    simmc_hitCols_perCluster_.clear(); simmc_earliestHitTime_perCluster_.clear();
    simmc_clusterType_.clear(); 
    simmc_primary_energy_.clear(); 
    simmc_primary_et_.clear(); simmc_primary_phi_.clear(); 
    simmc_primary_eta_.clear(); simmc_primary_pdgId_.clear();

    // flatten detsetvector
    std::vector<const FTLMergedCluster*> mergedVec;
    mergedVec.reserve(mergedClustersHandle->size());
    for (const auto& detSet : *mergedClustersHandle) {
        for (const auto& mc : detSet) {
            mergedVec.push_back(&mc);
        }
    }

    mc_n_ = mergedVec.size();
    simmc_n_ = simMergedClustersHandle->size();
    
    std::set<uint32_t> clusterDetIds;

    std::map<uint32_t, int> detIdToSCIndex;
    int scIndex = 0;
    
    // ---------------------- //
    // -------- RECO -------- //
    // ---------------------- //

    for (const auto* mcPtr : mergedVec) {
        const auto& mc = *mcPtr;
        h_mc_energy_->Fill(mc.energy());
        h_mc_time_->Fill(mc.time());
        h_mc_timeError_->Fill(mc.timeError());
        h_mc_x_->Fill(mc.x());
        h_mc_y_->Fill(mc.y());
        h_mc_nClusters_->Fill(mc.nClusters());

        h_mc_energy_vs_time_->Fill(mc.time(), mc.energy());
        h_mc_xy_->Fill(mc.x(), mc.y());
        h_mc_energy_vs_nClusters_->Fill(mc.nClusters(), mc.energy());

        mc_energy_.push_back(mc.energy());
        mc_time_.push_back(mc.time());
        mc_timeError_.push_back(mc.timeError());
        mc_x_.push_back(mc.x());
        mc_y_.push_back(mc.y());
        mc_nClusters_.push_back(mc.nClusters());
        mc_seedId_.push_back(mc.id().rawId());
        
        std::vector<uint32_t> clusterIds;
        for (const auto& detId : mc.clusterIds()) {
            clusterIds.push_back(detId.rawId());
            clusterDetIds.insert(detId.rawId());
            detIdToSCIndex[detId.rawId()] = scIndex;
        }
        mc_clusterIds_.push_back(clusterIds);
        scIndex++;

        if (evt_event_ <= 3) {
            std::cout << "MergedCluster: E=" << mc.energy() << " MeV, t=" << mc.time() 
                      << " ns, pos=(" << mc.x() << "," << mc.y() << ")" << std::endl;
            std::cout << "  Made from " << mc.nClusters() << " clusters:" << std::endl;
            for (const auto& detId : mc.clusterIds()) {
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
            bool inMergedCluster = clusterDetIds.count(cluster.id().rawId()) > 0;
            //h_merging_efficiency_->Fill(cluster.energy(), inMergedCluster ? 1 : 0);
            
            cluster_energy_.push_back(cluster.energy());
            cluster_time_.push_back(cluster.time());
            cluster_x_.push_back(cluster.x());
            cluster_y_.push_back(cluster.y());
            cluster_detId_.push_back(cluster.id().rawId());
            cluster_inMergedCluster_.push_back(inMergedCluster);
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

                    // Check if both clusters are in the same MergedCluster
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

            // Check if both clusters are in the same MergedCluster
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


    // ---------------------- //
    // -------- SIM --------- //
    // ---------------------- //

    for (const auto& simmc : *simMergedClustersHandle) {
        #if DEBUG>0
        std::cout << "\tSC: time = " << simmc.simTime() << std::endl;
        std::cout << "\tSC: position: " << simmc.simPos() << std::endl;
        std::cout << "\tSC: energy = " << simmc.simEnergy() << std::endl;
        std::cout << "\tSC: #clusters = " << simmc.clusters().size() << std::endl;
        #endif
        
        auto energy = convertUnitsTo(0.001_MeV, simmc.simEnergy()); // convert energy from GeV to MeV
        auto time = simmc.simTime();
        auto nClusters = simmc.clusters().size();

        h_simmc_energy_->Fill(energy);
        h_simmc_logEnergy_->Fill(log10(energy));
        h_simmc_time_->Fill(time);
        h_simmc_nClusters_->Fill(nClusters);

        h_simmc_energy_vs_nClusters_->Fill(nClusters, energy);
        h_simmc_energy_vs_time_->Fill(time, energy);

        // retrieve detId (-> geometry) of earliest cluster
        BTLDetId detId = simmc.simDetId();
        DetId geoId = detId.geographicalId(BTLDetId::CrysLayout::v3);

        const MTDGeomDet* thedet = geom->idToDet(geoId);
        // Convert simLC local position to global position
        GlobalPoint global_point(0,0,0);
        if (thedet != nullptr) {
            // get global position of cluster
            global_point = thedet->toGlobal(simmc.simPos());
        } else {
            #if DEBUG>0
            std::cout << "\tSC: WARNING - no geometry for detId " << detId.rawId() << std::endl;
            #endif
            global_point = GlobalPoint(-999, -999, -999);
        }

        h_simmc_x_->Fill(global_point.x());
        h_simmc_y_->Fill(global_point.y());
        h_simmc_eta_->Fill(global_point.eta());

        h_simmc_xy_->Fill(global_point.x(), global_point.y());

        // Tree filling
        simmc_energy_.push_back(energy);
        simmc_time_.push_back(time);
        // simmc_timeError_.push_back(simmc.timeError());
        simmc_x_.push_back(global_point.x());
        simmc_y_.push_back(global_point.y());
        simmc_eta_.push_back(global_point.eta());
        simmc_nClusters_.push_back(nClusters);

        // Collect per-cluster information for this mergedcluster
        std::vector<uint32_t> iphi_perCluster;
        std::vector<uint32_t> ieta_perCluster;
        std::vector<float> energy_perCluster;
        std::vector<float> time_perCluster;
        std::vector<float> earliestHitTime;
        std::vector<uint32_t> clusterType_perCluster;
        std::vector<std::vector<int>> hitCols;
        // std::vector<std::vector<float>> hitTimes;
        // std::vector<std::vector<float>> hitEnergies;

        // Access individual clusters from the mergedcluster
        for (const auto& cluster_ref : simmc.clusters()) {
            const auto& cluster = *cluster_ref;
            
            // Get detId from first hit (following MtdSimMergedClusterProducer pattern)
            if (cluster.detIds_and_rows().empty()) continue;
            BTLDetId clusterDetId(cluster.detIds_and_rows()[0].first);
            
            // Get topology indices for this cluster
            std::pair<uint32_t, uint32_t> cluster_indices = topology->btlIndex(clusterDetId.geographicalId(BTLDetId::CrysLayout::v3).rawId());
            iphi_perCluster.push_back(cluster_indices.first);
            ieta_perCluster.push_back(cluster_indices.second);
            
            // Store cluster energy and time using MtdSimLayerCluster methods
            energy_perCluster.push_back(convertUnitsTo(0.001_MeV, cluster.simLCEnergy())); // convert GeV to MeV
            time_perCluster.push_back(cluster.simLCTime());
            
            // Store cluster type (trackIdOffset)
            clusterType_perCluster.push_back(cluster.trackIdOffset());

            std::vector<int> hitCols_singleCluster;
            // Store list of hit cols for this cluster
            for (const auto& detId_row : cluster.detIds_and_rows()) {
                hitCols_singleCluster.push_back(detId_row.second.second);
            }
            hitCols.push_back(hitCols_singleCluster);

            // // take time of earliest hit in cluster
            // std::vector<float> hitTimes_singleCluster;
            // float earliestTime_cluster = 999;
            // for (const auto& hit_and_time : cluster.hits_and_times()) {
            //     hitTimes_singleCluster.push_back(hit_and_time.second);
            //     if (hit_and_time.second < earliestTime_cluster) {
            //         earliestTime_cluster = hit_and_time.second;
            //     }
            // }
            // hitTimes.push_back(hitTimes_singleCluster);
            // earliestHitTime.push_back(earliestTime_cluster);

            // std::vector<float> hitEnergies_singleCluster;
            // for (const auto& hit_and_energy : cluster.hits_and_energies()) {
            //     hitEnergies_singleCluster.push_back(convertUnitsTo(0.001_MeV, hit_and_energy.second)); // convert GeV to MeV
            // }
            // hitEnergies.push_back(hitEnergies_singleCluster);

            // Fill histograms
            h_simmc_logEnergy_perCluster_->Fill(log10(convertUnitsTo(0.001_MeV, cluster.simLCEnergy())));
            h_simmc_time_perCluster_->Fill(cluster.simLCTime());
            h_simmc_clusterType_->Fill(cluster.trackIdOffset());            
            
            #if DEBUG>0
            std::cout << "\t\t cluster detId: " << clusterDetId.rawId() 
                        << " iPhi: " << cluster_indices.first 
                        << " iEta: " << cluster_indices.second
                        << " E: " << convertUnitsTo(0.001_MeV, cluster.simLCEnergy()) << " MeV"
                        << " t: " << cluster.simLCTime() << " ns"
                        << " trackIdOffset: " << cluster.trackIdOffset() << std::endl;
            #endif
        }
        
        simmc_iphi_perCluster_.push_back(iphi_perCluster);
        simmc_ieta_perCluster_.push_back(ieta_perCluster);
        simmc_energy_perCluster_.push_back(energy_perCluster);
        simmc_time_perCluster_.push_back(time_perCluster);
        simmc_clusterType_.push_back(clusterType_perCluster);
        simmc_hitCols_perCluster_.push_back(hitCols);
        simmc_earliestHitTime_perCluster_.push_back(earliestHitTime);

        // count the number of unique ieta entries: that's your number of modules
        std::set<uint32_t> unique_modules(ieta_perCluster.begin(), ieta_perCluster.end());
        simmc_nModules_.push_back(int(unique_modules.size()));
        
        // Find primary tracking particle for this mergedcluster
        float primary_energy = -999.0;
        float primary_et = -999.0;
        float primary_phi = -999.0;
        float primary_eta = -999.0;
        int primary_pdgId = 0;
        
        // Take any tracking particle and climb back to original ancestor
        if (!simmc.trackingParticles().empty()) {
            auto tp_ref = simmc.trackingParticles()[0]; // Take first TP
            
            #if DEBUG>0
            std::cout << "\t\t Starting TP: pdgId = " << tp_ref->pdgId() 
                      << ", status = " << tp_ref->status()
                      << ", E = " << tp_ref->energy() << " GeV"
                      << ", #g4tracks = " << tp_ref->g4Tracks().size()
                      << std::endl;
            #endif
            
            // Climb back to original ancestor
            TrackingParticleRef current = tp_ref;
            
            while (true) {
                // Check if this particle has parent vertices
                const auto& parentVertices = current->parentVertex();
                if (parentVertices.isNull() || !parentVertices.isAvailable()) {
                    // No parent vertex - this is the root ancestor
                    #if DEBUG>0
                    std::cout << "\t\t Found root ancestor (no parent vertex)" << std::endl;
                    #endif
                    break;
                }
                
                // Get the parent tracks from the parent vertex
                const auto& parentTracks = parentVertices->sourceTracks();
                if (parentTracks.empty()) {
                    // No parent tracks - this is the root ancestor
                    #if DEBUG>0
                    std::cout << "\t\t Found root ancestor (no parent tracks)" << std::endl;
                    #endif
                    break;
                }
                
                // Move to the first parent
                current = parentTracks[0];
                #if DEBUG>0
                std::cout << "\t\t Moving to parent: pdgId = " << current->pdgId() 
                          << ", status = " << current->status() << std::endl;
                #endif
            }
            
            // Now check if the ancestor passes primary particle criteria
            const auto& ancestor = *current;
            
            #if DEBUG>0
            std::cout << "\t\t Final ancestor: pdgId = " << ancestor.pdgId() 
                      << ", status = " << ancestor.status()
                      << ", E = " << ancestor.energy() << " GeV"
                      << ", #g4tracks = " << ancestor.g4Tracks().size()
                      << std::endl;
            #endif
            
            // Apply primary particle criteria to ancestor
            bool isPrimary = true;
            if (ancestor.status() != 1) isPrimary = false;
            if (ancestor.g4Tracks().size() > 0) {
                if (ancestor.g4Tracks().front().vertIndex() != 0) isPrimary = false;
            } else {
                isPrimary = false;
            }
            
            if (isPrimary) {
                // This ancestor passes the criteria - extract its properties
                primary_energy = ancestor.energy(); // GeV
                primary_et = ancestor.pt(); // GeV
                primary_phi = ancestor.phi();
                primary_eta = ancestor.eta();
                primary_pdgId = ancestor.pdgId();
                
                #if DEBUG>0
                std::cout << "\t\t Found valid primary ancestor: E=" << primary_energy << " GeV, ET=" << primary_et 
                          << " GeV, phi=" << primary_phi << ", eta=" << primary_eta << ", pdgId=" << primary_pdgId << std::endl;
                #endif
            } else {
                #if DEBUG>0
                std::cout << "\t\t Ancestor does not pass primary criteria" << std::endl;
                #endif
            }
        }
        
        h_simmc_primaryPt_vs_nClusters_->Fill(nClusters, primary_et);
        h_simmc_primaryPt_vs_energy_->Fill(energy, primary_et);
        h_simmc_primaryEnergy_vs_energy_->Fill(energy, primary_energy);
        h_simmc_primaryEnergy_vs_nClusters_->Fill(nClusters, primary_energy);

        simmc_primary_energy_.push_back(primary_energy);
        simmc_primary_et_.push_back(primary_et);
        simmc_primary_phi_.push_back(primary_phi);
        simmc_primary_eta_.push_back(primary_eta);
        simmc_primary_pdgId_.push_back(primary_pdgId);
        
        // Print detailed info for first few events
        if (evt_event_ <= 3) {
            std::cout << "MergedCluster: MeV, t=" << time << " ns;"
                      << " E = " << energy << " MeV;"
                      << " pos = (" << simmc.simPos().x() << "," << simmc.simPos().y() << ")" << std::endl;
            std::cout << "  Made from " << nClusters << " clusters:" << std::endl;
            for (const auto& detId : simmc.detIds()) {
                std::cout << "    DetId: " << detId.rawId() << std::endl;
            }
        }
    }

    h_simmc_n_->Fill(simmc_n_);

    tree_->Fill();

    if (evt_event_ <= 3) {
        std::cout << "Event " << evt_event_ << ": " << mc_n_ << " MergedClusters, " 
                  << cluster_n_ << " individual clusters" << std::endl;
        std::cout << "Adjacent pairs: " << adjacentPairs 
                  << ", Merged pairs: " << mergedPairs << std::endl;
    }
}

void MergedClusterValidation_dev::endJob() {

    h_eta_merging_efficiency_->Divide(h_eta_merged_pairs_, h_eta_adjacent_pairs_, 1, 1, "B");

    std::cout << "\n=== MTD MergedCluster Validation Summary ===" << std::endl;
    std::cout << "Total events processed: " << tree_->GetEntries() << std::endl;
    std::cout << "Average MergedClusters per event: " << h_mc_energy_->GetEntries() / tree_->GetEntries() << std::endl;
    std::cout << "Total adjacent cluster pairs found: " << totalAdjacentPairs_ << std::endl;
    std::cout << "Total merged pairs: " << totalMergedPairs_ << std::endl;
    if (totalAdjacentPairs_ > 0) {
        std::cout << "Pair merging efficiency: " << (double)totalMergedPairs_ / totalAdjacentPairs_ * 100 << "%" << std::endl;
    }
    //std::cout << "Merging efficiency: " << h_merging_efficiency_->GetBinContent(2, 2) / h_merging_efficiency_->GetEntries() * 100 << "%" << std::endl;
}

DEFINE_FWK_MODULE(MergedClusterValidation_dev);