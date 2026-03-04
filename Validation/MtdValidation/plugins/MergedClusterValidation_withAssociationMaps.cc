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
    trackingParticlesToken_ = consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"));    

    mergedRecoToSimMap_ = consumes<reco::MergedRecoToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("mergedRecoToSimMap"));
    mergedSimToRecoMap_ = consumes<reco::MergedSimToRecoCollectionMtd>(iConfig.getParameter<edm::InputTag>("mergedSimToRecoMap"));

    mergedSimToTPMap_ = consumes<reco::MergedSimToTPCollectionMtd>(iConfig.getParameter<edm::InputTag>("mergedSimToTPMap"));
    mergedTPToSimMap_ = consumes<reco::TPToMergedSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("mergedTPToSimMap"));

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
    h_deltaTime_ = fs->make<TH1F>("h_deltaTime", "Time Resolution (Reco - Sim);#Delta t [ns]; Entries", 50, -1, 1);
    h_deltaEnergy_ = fs->make<TH1F>("h_deltaEnergy", "Energy Resolution (Reco - Sim); #Delta E [MeV]; Entries", 80, -20, 20);
    h_deltaX_ = fs->make<TH1F>("h_deltaX", "Local X Resolution (Reco - Sim); #Delta X [cm]; Entries", 40, -10, 10);
    h_deltaY_ = fs->make<TH1F>("h_deltaY", "Local Y Resolution (Reco - Sim); #Delta Y [cm]; Entries", 50, -3, 3);

    h_deltaNclu_ = fs->make<TH1F>("h_deltaNclu", "Difference in Cluster Multiplicity (Reco - Sim); #Delta N_{clusters}; Entries", 21, -10, 10);
    h_deltaNhits_ = fs->make<TH1F>("h_deltaNhits", "Difference in Hit Multiplicity (Reco - Sim); #Delta N_{hits}; Entries", 41, -20, 20);
    h_nSimPerReco_ = fs->make<TH1F>("h_nSimPerReco", "Number of SimMergedClusters associated to each RecoMergedCluster; N_{SimMergedClusters} per Reco MergedCluster; Entries", 10, 0, 10);

    h_deltaTime_multiClu_ = fs->make<TH1F>("h_deltaTime_multiClu", "Time Resolution (Reco - Sim) Multi-Cluster;#Delta t [ns]; Entries", 50, -1, 1);
    h_deltaEnergy_multiClu_ = fs->make<TH1F>("h_deltaEnergy_multiClu", "Energy Resolution (Reco - Sim) Multi-Cluster; #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaX_multiClu_ = fs->make<TH1F>("h_deltaX_multiClu", "Local X Resolution (Reco - Sim) Multi-Cluster; #Delta X [cm]; Entries", 40, -10, 10);
    h_deltaY_multiClu_ = fs->make<TH1F>("h_deltaY_multiClu", "Local Y Resolution (Reco - Sim) Multi-Cluster; #Delta Y [cm]; Entries", 50, -3, 3);

    h_deltaNclu_multiClu_ = fs->make<TH1F>("h_deltaNclu_multiClu", "Difference in Cluster Multiplicity (Reco - Sim) Multi-Cluster; #Delta N_{clusters}; Entries", 21, -10, 10);
    h_deltaNhits_multiClu_ = fs->make<TH1F>("h_deltaNhits_multiClu", "Difference in Hit Multiplicity (Reco - Sim) Multi-Cluster; #Delta N_{hits}; Entries", 41, -20, 20);
    h_nSimPerReco_multiClu_ = fs->make<TH1F>("h_nSimPerReco_multiClu", "Number of SimMergedClusters associated to each RecoMergedCluster Multi-Cluster; N_{SimMergedClusters} per Reco MergedCluster; Entries", 10, 0, 10);

    h_deltaTime_singleClu_ = fs->make<TH1F>("h_deltaTime_singleClu", "Time Resolution (Reco - Sim) Single-Cluster;#Delta t [ns]; Entries", 50, -1, 1);
    h_deltaEnergy_singleClu_ = fs->make<TH1F>("h_deltaEnergy_singleClu", "Energy Resolution (Reco - Sim) Single-Cluster; #Delta E [MeV]; Entries", 80, -5, 5);
    h_deltaX_singleClu_ = fs->make<TH1F>("h_deltaX_singleClu", "Local X Resolution (Reco - Sim) Single-Cluster; #Delta X [cm]; Entries", 40, -10, 10);
    h_deltaY_singleClu_ = fs->make<TH1F>("h_deltaY_singleClu", "Local Y Resolution (Reco - Sim) Single-Cluster; #Delta Y [cm]; Entries", 50, -3, 3);

    h_deltaTime_vs_Eta_ = fs->make<TH2F>("h_deltaTime_vs_Eta", "Time Resolution vs Eta;#eta; #Delta t [ns]", 50, -1.5, 1.5, 50, -0.5, 0.5);
    h_deltaEnergy_vs_Eta_ = fs->make<TH2F>("h_deltaEnergy_vs_Eta", "Energy Resolution vs Eta;#eta; #Delta E [MeV]", 50, -1.5, 1.5, 80, -20., 20.);
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

    edm::Handle<TrackingParticleCollection> trackingParticlesHandle;
    iEvent.getByToken(trackingParticlesToken_, trackingParticlesHandle);

    const auto& mergedRecoToSimMap = iEvent.get(mergedRecoToSimMap_);
    const auto& mergedSimToRecoMap = iEvent.get(mergedSimToRecoMap_);
    
    const auto& mergedSimToTPMap = iEvent.get(mergedSimToTPMap_);
    const auto& mergedTPToSimMap = iEvent.get(mergedTPToSimMap_);

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

    /*for (const auto& detSet : *mergedClustersHandle) {
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
            
            const std::vector<MtdSimMergedClusterRef>& simMergedRefs = itp. first->second;
            
            h_nSimPerReco_->Fill(simMergedRefs.size());

            std::cout << "RECO MergedCluster with " << mc.clusterRefs().size() << " reco clusters has " << simMergedRefs.size() << " matches in SIM." << std::endl;
            std::cout << "E = " << mc.energy() << " MeV, t = " << mc.time() << " ns" << std::endl;
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
                float deltaEnergy = mc.energy() - simEnergy;
                LocalPoint simCluPos = (*simRef->clusters().begin())->simLCPos();
                float deltaX = mc.x() - simCluPos.x();
                float deltaY = mc.y() - simCluPos.y();

                int nClusters = mc.clusterIds().size();
                int deltaNclu = nClusters - simRef->clusters().size();
                
                int recoHits = 0;
                for (const auto& cluRef : mc.clusterRefs()) {
                    recoHits += cluRef->size();  
                }
                int simHits = 0;
                for (const auto& simCluRef : simRef->clusters()) {
                    simHits += simCluRef->hits_and_fractions().size();
                }
                
                int deltaNhits = recoHits - simHits;

                h_deltaTime_->Fill(deltaTime);
                h_deltaEnergy_->Fill(deltaEnergy);
                h_deltaX_->Fill(deltaX);
                h_deltaY_->Fill(deltaY);
                h_deltaNclu_->Fill(deltaNclu);
                h_deltaNhits_->Fill(deltaNhits);

                //std::cout << "  deltaX = " << deltaX << " mm, deltaY = " << deltaY << " mm" << std::endl;
                //std::cout << "  reco: x=" << mc.x() << " y=" << mc.y() << std::endl;
                //std::cout << "  sim:  x=" << simRef->simPos().x() << " y=" << simRef->simPos().y() << std::endl;
                
                if (nClusters > 1) {
                    h_deltaTime_multiClu_->Fill(deltaTime);
                    h_deltaEnergy_multiClu_->Fill(deltaEnergy);
                    h_deltaX_multiClu_->Fill(deltaX);
                    h_deltaY_multiClu_->Fill(deltaY);

                    h_deltaNclu_multiClu_->Fill(deltaNclu);
                    h_deltaNhits_multiClu_->Fill(deltaNhits);
                    h_nSimPerReco_multiClu_->Fill(simMergedRefs.size());
                } else {
                    h_deltaTime_singleClu_->Fill(deltaTime);
                    h_deltaEnergy_singleClu_->Fill(deltaEnergy);
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
    }*/

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

            std::cout << "RECO MergedCluster with " << mc.clusterRefs().size() << " reco clusters has " << simMergedRefs.size() << " matches in SIM." << std::endl;
            std::cout << "E = " << mc.energy() << " MeV, t = " << mc.time() << " ns" << std::endl;
        
            // Calculate energy-weighted averages if multiple sim matches
            float totalSimEnergy = 0.0;
            float weightedSimTime = 0.0;
            float weightedSimX = 0.0;
            float weightedSimY = 0.0;
            int totalSimClusters = 0;
            int totalSimHits = 0;
            float simEta = 0.0;
            bool hasValidEta = false;
            
            // First pass: collect valid sim objects and calculate total energy
            std::vector<std::pair<MtdSimMergedClusterRef, float>> validSimRefs;
            for (const auto& simRef : simMergedRefs) {
                if (!simRef.isNonnull()) {
                    std::cout << "ERROR: Invalid simMergedClusRef!" << std::endl;
                    continue;
                }
                
                float simEnergy = convertUnitsTo(0.001_MeV, simRef->simEnergy());
                if (simEnergy < 1.0) continue;
                
                validSimRefs.push_back(std::make_pair(simRef, simEnergy));
                totalSimEnergy += simEnergy;
            }
            
            if (validSimRefs.empty()) continue;
            
            // Second pass: calculate weighted averages
            for (const auto& [simRef, simEnergy] : validSimRefs) {
                std::cout << "    Sim MergedCluster: E=" << simEnergy << " MeV, t=" << simRef->simTime() << " ns" 
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
                
                float weight = simEnergy / totalSimEnergy;
                
                weightedSimTime += simRef->simTime() * weight;
                
                LocalPoint simCluPos = (*simRef->clusters().begin())->simLCPos();
                weightedSimX += simCluPos.x() * weight;
                weightedSimY += simCluPos.y() * weight;
                
                totalSimClusters += simRef->clusters().size();
                
                for (const auto& simCluRef : simRef->clusters()) {
                    totalSimHits += simCluRef->hits_and_fractions().size();
                }
                
                // Get eta from first valid sim (or could also do weighted average)
                if (!hasValidEta) {
                    DetId simDetId = simRef->simDetId();
                    if (simDetId.rawId() != 0) {
                        BTLDetId btlId(simDetId);
                        DetId geoId = btlId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(topology->getMTDTopologyMode()));
                        const MTDGeomDet* thedet = geom->idToDet(geoId);
                        
                        if (thedet != nullptr) {
                            GlobalPoint simGlobalPos = thedet->toGlobal(simRef->simPos());
                            simEta = simGlobalPos.eta();
                            hasValidEta = true;
                        }
                    }
                }
            }
            
            // Calculate resolutions using energy-weighted averages
            float deltaTime = mc.time() - weightedSimTime;
            float deltaEnergy = mc.energy() - totalSimEnergy;
            float deltaX = mc.x() - weightedSimX;
            float deltaY = mc.y() - weightedSimY;
            
            int nClusters = mc.clusterIds().size();
            int deltaNclu = nClusters - totalSimClusters;
            
            int recoHits = 0;
            for (const auto& cluRef : mc.clusterRefs()) {
                recoHits += cluRef->size();  
            }
            int deltaNhits = recoHits - totalSimHits;
            
            std::cout << "  Energy-weighted averages: simTime=" << weightedSimTime << " ns, simEnergy=" << totalSimEnergy 
                      << " MeV, simX=" << weightedSimX << " cm, simY=" << weightedSimY << " cm" << std::endl;
            std::cout << "  Resolutions: deltaT=" << deltaTime << " ns, deltaE=" << deltaEnergy 
                      << " MeV, deltaX=" << deltaX << " cm, deltaY=" << deltaY << " cm" << std::endl;

            h_deltaTime_->Fill(deltaTime);
            h_deltaEnergy_->Fill(deltaEnergy);
            h_deltaX_->Fill(deltaX);
            h_deltaY_->Fill(deltaY);
            h_deltaNclu_->Fill(deltaNclu);
            h_deltaNhits_->Fill(deltaNhits);
            
            if (nClusters > 1) {
                h_deltaTime_multiClu_->Fill(deltaTime);
                h_deltaEnergy_multiClu_->Fill(deltaEnergy);
                h_deltaX_multiClu_->Fill(deltaX);
                h_deltaY_multiClu_->Fill(deltaY);

                h_deltaNclu_multiClu_->Fill(deltaNclu);
                h_deltaNhits_multiClu_->Fill(deltaNhits);
                h_nSimPerReco_multiClu_->Fill(simMergedRefs.size());
            } else {
                h_deltaTime_singleClu_->Fill(deltaTime);
                h_deltaEnergy_singleClu_->Fill(deltaEnergy);
                h_deltaX_singleClu_->Fill(deltaX);
                h_deltaY_singleClu_->Fill(deltaY);
            }

            if (hasValidEta) {
                h_deltaTime_vs_Eta_->Fill(simEta, deltaTime);
                h_deltaEnergy_vs_Eta_->Fill(simEta, deltaEnergy);
            }
        }
    }

    // ----------------------------------- //
    // -------- SIM <-> TP testing --------- //
    // ----------------------------------- //
    
    // let's just check multiplicity of association map: for each sim merged cluster, how many TP matches do we have?
    for (const auto& simMergedCluster : *simMergedClustersHandle) {
        // compute ref for sim merged cluster
        MtdSimMergedClusterRef simMergedClusterRef = edm::Ref<MtdSimMergedClusterCollection>(simMergedClustersHandle, &simMergedCluster - &(*simMergedClustersHandle->begin()));
        auto itp = mergedSimToTPMap.find(simMergedClusterRef);

        if (itp != mergedSimToTPMap.end()){
            int nTPMatches = std::distance(itp->val.begin(), itp->val.end());
            std::cout << "Found " << nTPMatches << " TP matches to this SimMergedCluster." << std::endl;
            for (const auto& tpMatch : itp->val){
                if (tpMatch.isNonnull()) {
                    std::cout << "  TP MATCH : pdgId=" << tpMatch->pdgId() << ", pt=" << tpMatch->pt() << " GeV, eta=" << tpMatch->eta() << ", phi=" << tpMatch->phi() << std::endl;
                } else {
                    std::cout << "  WARNING: invalid TP ref in mergedSimToTPMap!" << std::endl;
                }
            }
        } else {
            std::cout << "No TP matches found for this SimMergedCluster." << std::endl;
        }
    }

    // and now do the opposite: for each TP, how many sim merged cluster matches do we have?
    std::cout << "TrackingParticles collection size: " << trackingParticlesHandle->size() << std::endl;
    for (const auto& tp : *trackingParticlesHandle) {
        TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticlesHandle, &tp - &(*trackingParticlesHandle->begin()));
        auto itp = mergedTPToSimMap.find(tpRef);
        if (itp != mergedTPToSimMap.end()) {
            int nSimMatches = std::distance(itp->val.begin(), itp->val.end());
            std::cout << "Found " << nSimMatches << " SimMergedCluster matches to this TP." << std::endl;
            for (const auto& simMatch : itp->val) {
                if (simMatch.isNonnull()) {
                    std::cout << "  Matched SimMergedCluster: E=" << convertUnitsTo(0.001_MeV, simMatch->simEnergy()) << " MeV, t=" << simMatch->simTime() << " ns" << std::endl;
                } else {
                    std::cout << "  WARNING: invalid SimMergedCluster ref in mergedTPToSimMap!" << std::endl;
                }
            }
        } else {
            std::cout << "No SimMergedCluster matches found for this TP." << std::endl;
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