#include "Validation/MtdValidation/plugins/MergedClusterValidation_withAssociationMaps.h"
#include <iostream>
#include <CLHEP/Units/SystemOfUnits.h>
#include "DataFormats/Math/interface/GeantUnits.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"

#define DEBUG 0

MergedClusterValidation_withAssociationMaps::MergedClusterValidation_withAssociationMaps(const edm::ParameterSet& iConfig):
    mtdgeoToken_(esConsumes<MTDGeometry, MTDDigiGeometryRecord>()),
    mtdtopoToken_(esConsumes<MTDTopology, MTDTopologyRcd>()) {
    usesResource("TFileService");
    
    mergedClustersToken_ = consumes<FTLMergedClusterCollection>(iConfig.getParameter<edm::InputTag>("mergedClusters"));
    clustersToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("clusters"));

    simMergedClustersToken_ = consumes<MtdSimMergedClusterCollection>(iConfig.getParameter<edm::InputTag>("simMergedClusters"));
    simClustersToken_ = consumes<MtdSimLayerClusterCollection>(iConfig.getParameter<edm::InputTag>("simLayerClusters"));

    mergedRecoToSimMap_ = consumes<reco::MergedRecoToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("mergedRecoToSimMap"));
    mergedSimToRecoMap_ = consumes<reco::MergedSimToRecoCollectionMtd>(iConfig.getParameter<edm::InputTag>("mergedSimToRecoMap"));

    totalAdjacentPairs_ = 0;
    totalMergedPairs_ = 0;
}

void MergedClusterValidation_withAssociationMaps::beginJob() {
    edm::Service<TFileService> fs;
    
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

    // histos -- resolution plots using association maps between reco <--> sim
    h_deltaTime_ = fs->make<TH1F>("h_deltaTime", "Time Resolution (Reco - Sim);#Delta t [ns]; Entries", 50, -0.5, .5);
    h_deltaEnergy_ = fs->make<TH1F>("h_deltaEnergy", "Energy Resolution (Reco - Sim); #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaX_ = fs->make<TH1F>("h_deltaX", "Local X Resolution (Reco - Sim); #Delta X [cm]; Entries", 40, -10, 10);
    h_deltaY_ = fs->make<TH1F>("h_deltaY", "Local Y Resolution (Reco - Sim); #Delta Y [cm]; Entries", 50, -3, 3);

    h_deltaNclu_ = fs->make<TH1F>("h_deltaNclu", "Difference in Cluster Multiplicity (Reco - Sim); #Delta N_{clusters}; Entries", 20, -10, 10);
    h_deltaNhits_ = fs->make<TH1F>("h_deltaNhits", "Difference in Hit Multiplicity (Reco - Sim); #Delta N_{hits}; Entries", 40, -20, 20);
    h_nSimPerReco_ = fs->make<TH1F>("h_nSimPerReco", "Number of SimMergedClusters associated to each RecoMergedCluster; N_{SimMergedClusters} per Reco MergedCluster; Entries", 6, -0.5, 5.5);

    h_deltaTime_multiClu_ = fs->make<TH1F>("h_deltaTime_multiClu", "Time Resolution (Reco - Sim) Multi-Cluster;#Delta t [ns]; Entries", 50, -0.5, .5);
    h_deltaEnergy_multiClu_ = fs->make<TH1F>("h_deltaEnergy_multiClu", "Energy Resolution (Reco - Sim) Multi-Cluster; #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaX_multiClu_ = fs->make<TH1F>("h_deltaX_multiClu", "Local X Resolution (Reco - Sim) Multi-Cluster; #Delta X [cm]; Entries", 40, -10, 10);
    h_deltaY_multiClu_ = fs->make<TH1F>("h_deltaY_multiClu", "Local Y Resolution (Reco - Sim) Multi-Cluster; #Delta Y [cm]; Entries", 50, -3, 3);

    h_deltaNclu_multiClu_ = fs->make<TH1F>("h_deltaNclu_multiClu", "Difference in Cluster Multiplicity (Reco - Sim) Multi-Cluster; #Delta N_{clusters}; Entries", 20, -10, 10);
    h_deltaNhits_multiClu_ = fs->make<TH1F>("h_deltaNhits_multiClu", "Difference in Hit Multiplicity (Reco - Sim) Multi-Cluster; #Delta N_{hits}; Entries", 40, -20, 20);
    h_nSimPerReco_multiClu_ = fs->make<TH1F>("h_nSimPerReco_multiClu", "Number of SimMergedClusters associated to each RecoMergedCluster Multi-Cluster; N_{SimMergedClusters} per Reco MergedCluster; Entries", 6, -0.5, 5.5);

    h_deltaTime_singleClu_ = fs->make<TH1F>("h_deltaTime_singleClu", "Time Resolution (Reco - Sim) Single-Cluster;#Delta t [ns]; Entries", 50, -0.5, 0.5);
    h_deltaEnergy_singleClu_ = fs->make<TH1F>("h_deltaEnergy_singleClu", "Energy Resolution (Reco - Sim) Single-Cluster; #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaX_singleClu_ = fs->make<TH1F>("h_deltaX_singleClu", "Local X Resolution (Reco - Sim) Single-Cluster; #Delta X [cm]; Entries", 40, -10, 10);
    h_deltaY_singleClu_ = fs->make<TH1F>("h_deltaY_singleClu", "Local Y Resolution (Reco - Sim) Single-Cluster; #Delta Y [cm]; Entries", 50, -3, 3);

    // deltaEnergy histograms by trackIdOffset (only for single sim MC matches)
    h_deltaEnergy_trackIdOffset_0_ = fs->make<TH1F>("h_deltaEnergy_trackIdOffset_0", "Energy Resolution (Primary); #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaEnergy_trackIdOffset_1_ = fs->make<TH1F>("h_deltaEnergy_trackIdOffset_1", "Energy Resolution (Secondary); #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaEnergy_trackIdOffset_2_ = fs->make<TH1F>("h_deltaEnergy_trackIdOffset_2", "Energy Resolution (Loopers); #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaEnergy_trackIdOffset_3_ = fs->make<TH1F>("h_deltaEnergy_trackIdOffset_3", "Energy Resolution (Backscatter); #Delta E [MeV]; Entries", 80, -5, 5);

    // deltaTime histograms by trackIdOffset (only for single sim MC matches)
    h_deltaTime_trackIdOffset_0_ = fs->make<TH1F>("h_deltaTime_trackIdOffset_0", "Time Resolution (Primary); #Delta t [ns]; Entries", 50, -0.5, 0.5);
    h_deltaTime_trackIdOffset_1_ = fs->make<TH1F>("h_deltaTime_trackIdOffset_1", "Time Resolution (Secondary); #Delta t [ns]; Entries", 50, -0.5, 0.5);
    h_deltaTime_trackIdOffset_2_ = fs->make<TH1F>("h_deltaTime_trackIdOffset_2", "Time Resolution (Loopers); #Delta t [ns]; Entries", 50, -0.5, 0.5);
    h_deltaTime_trackIdOffset_3_ = fs->make<TH1F>("h_deltaTime_trackIdOffset_3", "Time Resolution (Backscatter); #Delta t [ns]; Entries", 50, -0.5, 0.5);

    h_nSimPerReco_trackIdOffset_0_ = fs->make<TH1F>("h_nSimPerReco_trackIdOffset_0", "Number of Primary SimMergedClusters per Reco; N_{SimMC} (Primary); Entries", 6, -0.5, 5.5);
    h_nSimPerReco_trackIdOffset_1_ = fs->make<TH1F>("h_nSimPerReco_trackIdOffset_1", "Number of Secondary SimMergedClusters per Reco; N_{SimMC} (Secondary); Entries", 6, -0.5, 5.5);
    h_nSimPerReco_trackIdOffset_2_ = fs->make<TH1F>("h_nSimPerReco_trackIdOffset_2", "Number of Looper SimMergedClusters per Reco; N_{SimMC} (Loopers); Entries", 6, -0.5, 5.5);
    h_nSimPerReco_trackIdOffset_3_ = fs->make<TH1F>("h_nSimPerReco_trackIdOffset_3", "Number of Backscatter SimMergedClusters per Reco; N_{SimMC} (Backscatter); Entries", 6, -0.5, 5.5);

    h_deltaTime_vs_Eta_ = fs->make<TH2F>("h_deltaTime_vs_Eta", "Time Resolution vs Eta;#eta; #Delta t [ns]", 50, -1.5, 1.5, 50, -0.5, 0.5);
    h_deltaEnergy_vs_Eta_ = fs->make<TH2F>("h_deltaEnergy_vs_Eta", "Energy Resolution vs Eta;#eta; #Delta E [MeV]", 50, -1.5, 1.5, 80, -5, 5.);
}

void MergedClusterValidation_withAssociationMaps::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace geant_units::operators; // for energy conversion

    edm::Handle<FTLMergedClusterCollection> mergedClustersHandle;
    iEvent.getByToken(mergedClustersToken_, mergedClustersHandle);
    
    edm::Handle<FTLClusterCollection> clustersHandle;
    iEvent.getByToken(clustersToken_, clustersHandle);

    edm::Handle<MtdSimMergedClusterCollection> simMergedClustersHandle;
    iEvent.getByToken(simMergedClustersToken_, simMergedClustersHandle);

    edm::Handle<MtdSimLayerClusterCollection> mtdSimLCHandle;
    iEvent.getByToken(simClustersToken_, mtdSimLCHandle);

    const auto& mergedRecoToSimMap = iEvent.get(mergedRecoToSimMap_);
    const auto& mergedSimToRecoMap = iEvent.get(mergedSimToRecoMap_);
    
    // edm::Handle<reco::MergedRecoToSimCollectionMtd> mergedRecoToSimMapHandle;
    // iEvent.getByToken(mergedRecoToSimMap_, mergedRecoToSimMapHandle);
    // edm::Handle<reco::MergedSimToRecoCollectionMtd> mergedSimToRecoMapHandle;
    // iEvent.getByToken(mergedSimToRecoMap_, mergedSimToRecoMapHandle);

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

    mc_n_ = mergedClustersHandle->size();
    simmc_n_ = simMergedClustersHandle->size();
    
    std::set<uint32_t> clusterDetIds;

    std::map<uint32_t, int> detIdToSCIndex;
    int scIndex = 0;

    // ----------------------------- //
    // -------- RECO VS SIM -------- //
    // ----------------------------- //

    for (const auto& detSet : *mergedClustersHandle) {
        for (const auto& mc : detSet) {
            // look for match in association map
            FTLMergedClusterRef recoMergedClusRef = edmNew::makeRefTo(mergedClustersHandle, &mc);
            
            // check if reference is valid
            if (!recoMergedClusRef.isNonnull()) {
                std::cout << "ERROR: Invalid recoMergedClusRef!" << std::endl;
                continue;
            }

            auto itp = mergedRecoToSimMap.equal_range(recoMergedClusRef);
            if (itp.first == itp.second) {
                std::cout << "No matching SimMergedCluster found for Reco MergedCluster" << std::endl;
                continue;
            }
            
            const std::vector<MtdSimMergedClusterRef>& simMergedRefs = itp.first->second;
            
            h_nSimPerReco_->Fill(simMergedRefs.size());

            // Count sim matches by trackIdOffset type
            int nType0 = 0, nType1 = 0, nType2 = 0, nType3 = 0;
            std::cout << "#sim matches for this reco merged cluster: " << simMergedRefs.size() << std::endl;
            for (const auto& simRef : simMergedRefs) {
                if (simRef.isNonnull() && !simRef->clusters().empty()) {
                    int offset = simRef->clusters()[0]->trackIdOffset();
                    if (offset == 0) nType0++;
                    else if (offset == 1) nType1++;
                    else if (offset == 2) nType2++;
                    else if (offset == 3) nType3++;
                } else {
                    std::cout << "WARNING: invalid simMergedClusRef when counting trackIdOffset!" << std::endl;
                    std::cout << "WARNING: out of " << simMergedRefs.size() << " sim matches for this reco merged cluster." << std::endl;
                }
            }
            std::cout << "of which: " << nType0 << " primary, " << nType1 << " secondary, " << nType2 << " loopers, " << nType3 << " backscatter." << std::endl;
            if (nType0 + nType1 + nType2 + nType3 != simMergedRefs.size()) {
                std::cout << "WARNING: sum of type counts does not equal total number of sim matches!" << std::endl;
            }
            h_nSimPerReco_trackIdOffset_0_->Fill(nType0);
            h_nSimPerReco_trackIdOffset_1_->Fill(nType1);
            h_nSimPerReco_trackIdOffset_2_->Fill(nType2);
            h_nSimPerReco_trackIdOffset_3_->Fill(nType3);

            std::cout << "RECO MergedCluster with " << mc.clusterRefs().size() << " reco clusters has " << simMergedRefs.size() << " matches in SIM." << std::endl;
            std::cout << "E = " << mc.energy() << " MeV, t = " << mc.time() << " ns" << std::endl;
        
            // For single sim MC matches, we'll extract trackIdOffset for classification
            int trackIdOffset = -1;
            if (simMergedRefs.size() == 1) {
                const auto& simRef = simMergedRefs[0];
                if (simRef.isNonnull() && !simRef->clusters().empty()) {
                    trackIdOffset = simRef->clusters()[0]->trackIdOffset();
                }
            }

            // Calculate total sim energy for energy resolution (sum over all matched sim MCs)
            float totalSimEnergy = 0.0;
            int totalSimHits = 0;
            int totalSimClusters = 0;
            for (const auto& simRef : simMergedRefs) {
                if (!simRef.isNonnull()) continue;
                float simEnergy = convertUnitsTo(0.001_MeV, simRef->simEnergy());
                totalSimEnergy += simEnergy;
                totalSimClusters += simRef->clusters().size();
                for (const auto& simCluRef : simRef->clusters()) {
                    totalSimHits += simCluRef->hits_and_fractions().size();
                }
            }

            // Calculate reco hits
            int nClusters = mc.clusterIds().size();
            int recoHits = 0;
            for (const auto& cluRef : mc.clusterRefs()) {
                recoHits += cluRef->size();  
            }

            // Calculate deltas using total sim energy
            float deltaEnergy = mc.energy() - totalSimEnergy;
            int deltaNclu = nClusters - totalSimClusters;
            int deltaNhits = recoHits - totalSimHits;

            // Fill general histograms once per reco MC (using sum of sim energies)
            h_deltaEnergy_->Fill(deltaEnergy);
            h_deltaNclu_->Fill(deltaNclu);
            h_deltaNhits_->Fill(deltaNhits);

            if (nClusters > 1) {
                h_deltaEnergy_multiClu_->Fill(deltaEnergy);
                h_deltaNclu_multiClu_->Fill(deltaNclu);
                h_deltaNhits_multiClu_->Fill(deltaNhits);
                h_nSimPerReco_multiClu_->Fill(simMergedRefs.size());
            } else {
                h_deltaEnergy_singleClu_->Fill(deltaEnergy);
            }

            // Print warning for large total energy differences
            if (deltaEnergy > 1) {
                std::cout << "WARNING: Large TOTAL energy difference: DeltaE = " << deltaEnergy << " MeV" << std::endl;
                std::cout << "         Reco E = " << mc.energy() << " MeV, Total Sim E = " << totalSimEnergy << " MeV" << std::endl;
                std::cout << "         Reco nClusters = " << nClusters << ", Total Sim nClusters = " << totalSimClusters << std::endl;
                std::cout << "         Reco nHits = " << recoHits << ", Total Sim nHits = " << totalSimHits << std::endl;
                std::cout << "         Number of sim matches = " << simMergedRefs.size() << std::endl;
            }

            // iterate over matches and plot 
            for (const auto& simRef : simMergedRefs){
                if (!simRef.isNonnull()) {
                    std::cout << "ERROR: Invalid simMergedClusRef!" << std::endl;
                    continue;
                }

                float simEnergy = convertUnitsTo(0.001_MeV, simRef->simEnergy());
                if (simEnergy < 1.0) continue;

                std::cout << "    Sim MergedCluster: E=" << convertUnitsTo(0.001_MeV, simRef->simEnergy()) << " MeV, t=" << simRef->simTime() << " ns" 
                        << ", trackIdOffset = " << simRef->clusters()[0]->trackIdOffset() << " (nClusters = " << simRef->clusters().size() << ")" << std::endl;
                std::cout << "                        particles =";
                if(simRef->trackingParticles().size() > 0){
                    for(auto const& p : simRef->trackingParticles()){
                        std::cout << " " << p->pdgId();
                    }
                } else {
                    std::cout << "(no valid TP refs)";
                }
                std::cout << std::endl;

                float deltaTime = mc.time() - simRef->simTime();
                float deltaEnergyThisMatch = mc.energy() - simEnergy;
                LocalPoint simCluPos = (*simRef->clusters().begin())->simLCPos();
                float deltaX = mc.x() - simCluPos.x();
                float deltaY = mc.y() - simCluPos.y();

                int deltaNcluThisMatch = nClusters - simRef->clusters().size();
                
                int simHits = 0;
                for (const auto& simCluRef : simRef->clusters()) {
                    simHits += simCluRef->hits_and_fractions().size();
                }
                
                int deltaNhitsThisMatch = recoHits - simHits;
                
                std::cout << "    Delta t = " << deltaTime << " ns, Delta E (this match) = " << deltaEnergyThisMatch << " MeV" << "; Delta Nclu = " << deltaNcluThisMatch << ", Delta Nhits = " << deltaNhitsThisMatch << std::endl;

                h_deltaTime_->Fill(deltaTime);
                h_deltaX_->Fill(deltaX);
                h_deltaY_->Fill(deltaY);

                // Fill trackIdOffset-specific histograms for single sim MC matches
                if (simMergedRefs.size() == 1 && trackIdOffset >= 0 && trackIdOffset <= 3) {
                    if (trackIdOffset == 0) {
                        h_deltaEnergy_trackIdOffset_0_->Fill(deltaEnergy);
                        h_deltaTime_trackIdOffset_0_->Fill(deltaTime);
                    } else if (trackIdOffset == 1) {
                        h_deltaEnergy_trackIdOffset_1_->Fill(deltaEnergy);
                        h_deltaTime_trackIdOffset_1_->Fill(deltaTime);
                    } else if (trackIdOffset == 2) {
                        h_deltaEnergy_trackIdOffset_2_->Fill(deltaEnergy);
                        h_deltaTime_trackIdOffset_2_->Fill(deltaTime);
                    } else if (trackIdOffset == 3) {
                        h_deltaEnergy_trackIdOffset_3_->Fill(deltaEnergy);
                        h_deltaTime_trackIdOffset_3_->Fill(deltaTime);
                    }
                }

                //std::cout << "  deltaX = " << deltaX << " mm, deltaY = " << deltaY << " mm" << std::endl;
                //std::cout << "  reco: x=" << mc.x() << " y=" << mc.y() << std::endl;
                //std::cout << "  sim:  x=" << simRef->simPos().x() << " y=" << simRef->simPos().y() << std::endl;
                
                if (nClusters > 1) {
                    h_deltaTime_multiClu_->Fill(deltaTime);
                    h_deltaX_multiClu_->Fill(deltaX);
                    h_deltaY_multiClu_->Fill(deltaY);
                } else {
                    h_deltaTime_singleClu_->Fill(deltaTime);
                    h_deltaX_singleClu_->Fill(deltaX);
                    h_deltaY_singleClu_->Fill(deltaY);
                }

                DetId simDetId = simRef->simDetId();
                if (simDetId.rawId() != 0) {
                    BTLDetId btlId(simDetId);
                    DetId geoId = btlId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(topology->getMTDTopologyMode()));
                    const MTDGeomDet* thedet = geom->idToDet(geoId);

                    if (thedet != nullptr) {
                        GlobalPoint simGlobalPos = thedet->toGlobal(simRef->simPos());
                        float simEta = simGlobalPos.eta();

                        h_deltaTime_vs_Eta_->Fill(simEta, deltaTime);
                        h_deltaEnergy_vs_Eta_->Fill(simEta, deltaEnergy);
                    } else {
                        std::cout << "WARNING - could not get thedet for geoId " << geoId.rawId() << std::endl;
                    }
                } else {
                    std::cout << "WARNING - simDetId is zero!" << std::endl;
                }
            }
        }
    }

}

void MergedClusterValidation_withAssociationMaps::endJob() {

    // h_eta_merging_efficiency_->Divide(h_eta_merged_pairs_, h_eta_adjacent_pairs_, 1, 1, "B");

    // std::cout << "\n=== MTD MergedCluster Validation Summary ===" << std::endl;
    // std::cout << "Total events processed: " << tree_->GetEntries() << std::endl;
    // std::cout << "Average MergedClusters per event: " << h_mc_energy_->GetEntries() / tree_->GetEntries() << std::endl;
    // std::cout << "Total adjacent cluster pairs found: " << totalAdjacentPairs_ << std::endl;
    // std::cout << "Total merged pairs: " << totalMergedPairs_ << std::endl;
    // if (totalAdjacentPairs_ > 0) {
    //     std::cout << "Pair merging efficiency: " << (double)totalMergedPairs_ / totalAdjacentPairs_ * 100 << "%" << std::endl;
    // }
    //std::cout << "Merging efficiency: " << h_merging_efficiency_->GetBinContent(2, 2) / h_merging_efficiency_->GetEntries() * 100 << "%" << std::endl;
}

DEFINE_FWK_MODULE(MergedClusterValidation_withAssociationMaps);