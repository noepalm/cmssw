#include "Validation/MtdValidation/plugins/MergedClusterValidation_withAssociationMaps.h"
#include <iostream>
#include <CLHEP/Units/SystemOfUnits.h>
#include "DataFormats/Math/interface/GeantUnits.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

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
            // FTLMergedClusterRef recoMergedClusRef(mergedClustersHandle, &mc - &(*mergedClustersHandle->begin()));
            
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
            
            std::cout << "RECO MergedCluster with " << mc.clusterRefs().size() << " reco clusters has " << simMergedRefs.size() << " matches in SIM." << std::endl;
            std::cout << "E = " << mc.energy() << " MeV, t = " << mc.time() << " ns" << std::endl;
        
            // iterate over matches and print all properties
            for (const auto& simRef : simMergedRefs){
                if (!simRef.isNonnull()) {
                    std::cout << "ERROR: Invalid simMergedClusRef!" << std::endl;
                    continue;
                }

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