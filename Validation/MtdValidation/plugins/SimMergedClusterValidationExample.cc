#include "Validation/MtdValidation/plugins/SimMergedClusterValidationExample.h"
#include <iostream>
#include <CLHEP/Units/SystemOfUnits.h>
#include "DataFormats/Math/interface/GeantUnits.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

SimMergedClusterValidationExample::SimMergedClusterValidationExample(const edm::ParameterSet& iConfig):
    mtdgeoToken_(esConsumes<MTDGeometry, MTDDigiGeometryRecord>()),
    mtdtopoToken_(esConsumes<MTDTopology, MTDTopologyRcd>()) {
    usesResource("TFileService");

    // mergedClustersToken_ = consumes<FTLMergedClusterCollection>(iConfig.getParameter<edm::InputTag>("mergedClusters"));
    // clustersToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("clusters"));

    mergedClustersToken_ = consumes<MtdSimMergedClusterCollection>(iConfig.getParameter<edm::InputTag>("simMergedClusters"));
    clustersToken_ = consumes<MtdSimLayerClusterCollection>(iConfig.getParameter<edm::InputTag>("simLayerClusters"));
    
    // GenParticles are optional - only consume if specified in config
    if (iConfig.exists("genParticles")) {
        genParticlesToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
    }
}

#define DEBUG 0

void SimMergedClusterValidationExample::beginJob() {
    // Initialize event counter
    processed_event_counter_ = 0;
    
    edm::Service<TFileService> fs;
    
    // MergedCluster histograms
    h_sc_energy_ = fs->make<TH1F>("h_sc_energy", "MergedCluster Energy;Energy [MeV];Count", 50, 0, 2);
    h_sc_time_ = fs->make<TH1F>("h_sc_time", "MergedCluster Time;Time [ns];Count", 50, 0, 30);
    // h_sc_timeError_ = fs->make<TH1F>("h_sc_timeError", "MergedCluster Time Error;Time Error [ns];Count", 100, 0, 1);
    h_sc_x_ = fs->make<TH1F>("h_sc_x", "MergedCluster X;X [mm];Count", 50, -120, 120);
    h_sc_y_ = fs->make<TH1F>("h_sc_y", "MergedCluster Y;Y [mm];Count", 50, -120, 120);
    h_sc_nClusters_ = fs->make<TH1F>("h_sc_nClusters", "Number of Clusters;N_{clusters};Count", 11, -0.5, 10.5);
    
    // // Cluster histograms for comparison
    // h_cluster_energy_ = fs->make<TH1F>("h_cluster_energy", "Cluster Energy;Energy [MeV];Count", 100, 0, 50);
    // h_cluster_time_ = fs->make<TH1F>("h_cluster_time", "Cluster Time;Time [ns];Count", 100, -5, 20);
    
    // // 2D histograms
    h_sc_energy_vs_time_ = fs->make<TH2F>("h_sc_energy_vs_time", "MergedCluster Energy vs Time;Time [ns];Energy [MeV]", 30, 4, 27, 20, 0, 1);
    h_sc_xy_ = fs->make<TH2F>("h_sc_xy", "MergedCluster Position;X [mm];Y [mm]", 100, -120, 120, 100, -120, 120);
    h_sc_energy_vs_nClusters_ = fs->make<TH2F>("h_sc_energy_vs_nClusters", "Energy vs N Clusters;N_{clusters};Energy [MeV]", 4, -0.5, 3.5, 20, 0, 1);
    // h_merging_efficiency_ = fs->make<TH2F>("h_merging_efficiency", "Merging Efficiency;Cluster Energy [MeV];Merged?", 50, 0, 50, 2, -0.5, 1.5);
    
    // Analysis tree
    tree_ = fs->make<TTree>("MTDMergedClusters", "MTD MergedCluster Analysis Tree");
    
    tree_->Branch("evt_run", &evt_run_);
    tree_->Branch("evt_event", &evt_event_);
    tree_->Branch("evt_number", &evt_number_);
    
    tree_->Branch("sc_n", &sc_n_);
    tree_->Branch("sc_energy", &sc_energy_);
    tree_->Branch("sc_time", &sc_time_);
    // tree_->Branch("sc_timeError", &sc_timeError_);
    tree_->Branch("sc_x", &sc_x_);
    tree_->Branch("sc_y", &sc_y_);
    tree_->Branch("sc_nClusters", &sc_nClusters_);
    tree_->Branch("sc_iphi_perCluster", &sc_iphi_perCluster_);
    tree_->Branch("sc_ieta_perCluster", &sc_ieta_perCluster_);
    tree_->Branch("sc_energy_perCluster", &sc_energy_perCluster_);
    tree_->Branch("sc_time_perCluster", &sc_time_perCluster_);
    tree_->Branch("sc_clusterType", &sc_clusterType_);
    tree_->Branch("sc_primary_energy", &sc_primary_energy_);
    tree_->Branch("sc_primary_et", &sc_primary_et_);
    tree_->Branch("sc_primary_phi", &sc_primary_phi_);
    tree_->Branch("sc_primary_eta", &sc_primary_eta_);
    tree_->Branch("sc_primary_pdgId", &sc_primary_pdgId_);

    // tree_->Branch("cluster_n", &cluster_n_);
    // tree_->Branch("cluster_energy", &cluster_energy_);
    // tree_->Branch("cluster_time", &cluster_time_);
    // tree_->Branch("cluster_x", &cluster_x_);
    // tree_->Branch("cluster_y", &cluster_y_);
    // tree_->Branch("cluster_detId", &cluster_detId_);
    // tree_->Branch("cluster_inMergedCluster", &cluster_inMergedCluster_);
}

void SimMergedClusterValidationExample::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace geant_units::operators; // for energy conversion

    // edm::Handle<FTLMergedClusterCollection> mergedClustersHandle;
    // iEvent.getByToken(mergedClustersToken_, mergedClustersHandle);
    
    // edm::Handle<FTLClusterCollection> clustersHandle;
    // iEvent.getByToken(clustersToken_, clustersHandle);

    edm::Handle<MtdSimMergedClusterCollection> mtdSimSCHandle;
    iEvent.getByToken(mergedClustersToken_, mtdSimSCHandle);

    edm::Handle<MtdSimLayerClusterCollection> mtdSimLCHandle;
    iEvent.getByToken(clustersToken_, mtdSimLCHandle);

    // GenParticles are optional
    edm::Handle<reco::GenParticleCollection> genParticlesHandle;
    if (!genParticlesToken_.isUninitialized()) {
        iEvent.getByToken(genParticlesToken_, genParticlesHandle);
    }

    // retrieve geometry for global positions
    auto geometryHandle = iSetup.getTransientHandle(mtdgeoToken_);
    const MTDGeometry* geom = geometryHandle.product();

    // retrieve topology for geographical ids
    auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
    const MTDTopology* topology = topologyHandle.product();

    if (!mtdSimSCHandle.isValid()) {
        std::cout << "Invalid handles!" << std::endl;
        return;
    }
    
    evt_run_ = iEvent.id().run();
    evt_event_ = iEvent.id().event();
    evt_number_ = processed_event_counter_++;
    
    sc_energy_.clear();
    // sc_timeError_.clear();
    sc_time_.clear();
    sc_x_.clear(); sc_y_.clear(); sc_nClusters_.clear();
    sc_iphi_perCluster_.clear(); sc_ieta_perCluster_.clear();
    sc_energy_perCluster_.clear(); sc_time_perCluster_.clear();
    sc_clusterType_.clear();
    sc_primary_energy_.clear(); sc_primary_et_.clear(); 
    sc_primary_phi_.clear(); sc_primary_eta_.clear();
    sc_primary_pdgId_.clear();
    // cluster_energy_.clear();
    cluster_time_.clear();
    cluster_x_.clear(); cluster_y_.clear(); cluster_detId_.clear();
    cluster_inMergedCluster_.clear();
    
    sc_n_ = mtdSimSCHandle->size();
    
    // Create set of constituent DetIds for efficiency tracking
    std::set<uint32_t> clusterDetIds;

    #if DEBUG>0
    std::cout << "ACCESSING MTD SIM SUPERC" << std::endl;
    #endif

    for (const auto& sc : *mtdSimSCHandle) {
        #if DEBUG>0
        std::cout << "\tSC: time = " << sc.simTime() << std::endl;
        std::cout << "\tSC: position: " << sc.simPos() << std::endl;
        std::cout << "\tSC: energy = " << sc.simEnergy() << std::endl;
        std::cout << "\tSC: #clusters = " << sc.clusters().size() << std::endl;
        #endif
        
        auto energy = convertUnitsTo(0.001_MeV, sc.simEnergy()); // convert energy from GeV to MeV
        auto time = sc.simTime();
        auto nClusters = sc.clusters().size();

        h_sc_energy_->Fill(energy);
        h_sc_time_->Fill(time);
        h_sc_nClusters_->Fill(nClusters);

        // retrieve detId (-> geometry) of earliest cluster
        BTLDetId detId = sc.simDetId();
        // DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(topology->getMTDTopologyMode()));
        DetId geoId = detId.geographicalId(BTLDetId::CrysLayout::v3);

        const MTDGeomDet* thedet = geom->idToDet(geoId);
        // Convert simLC local position to global position
        GlobalPoint global_point(0,0,0);
        if (thedet != nullptr) {
            // get global position of cluster
            global_point = thedet->toGlobal(sc.simPos());
        } else {
            #if DEBUG>0
            std::cout << "\tSC: WARNING - no geometry for detId " << detId.rawId() << std::endl;
            #endif
            global_point = GlobalPoint(-999, -999, -999);
        }

        h_sc_x_->Fill(global_point.x());
        h_sc_y_->Fill(global_point.y());

        h_sc_xy_->Fill(global_point.x(), global_point.y());

        h_sc_energy_vs_nClusters_->Fill(nClusters, energy);
        h_sc_energy_vs_time_->Fill(time, energy);

        // // Get topology indices from geographicalId (module-level) with crystal layout
        // std::pair<uint32_t, uint32_t> indices = topology->btlIndex(detId.geographicalId(BTLDetId::CrysLayout::v3).rawId());
        // uint32_t iphi = indices.first;
        // uint32_t ieta = indices.second;

        sc_energy_.push_back(energy);
        sc_time_.push_back(time);
        // sc_timeError_.push_back(sc.timeError());
        sc_x_.push_back(global_point.x());
        sc_y_.push_back(global_point.y());
        sc_nClusters_.push_back(nClusters);

        // Collect per-cluster information for this mergedcluster
        std::vector<uint32_t> iphi_perCluster;
        std::vector<uint32_t> ieta_perCluster;
        std::vector<float> energy_perCluster;
        std::vector<float> time_perCluster;
        std::vector<uint32_t> clusterType_perCluster;
        
        // Access individual clusters from the mergedcluster
        for (const auto& cluster_ref : sc.clusters()) {
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
            
            #if DEBUG>0
            std::cout << "\t\t cluster detId: " << clusterDetId.rawId() 
                        << " iPhi: " << cluster_indices.first 
                        << " iEta: " << cluster_indices.second
                        << " E: " << convertUnitsTo(0.001_MeV, cluster.simLCEnergy()) << " MeV"
                        << " t: " << cluster.simLCTime() << " ns"
                        << " trackIdOffset: " << cluster.trackIdOffset() << std::endl;
            #endif
        }
        
        sc_iphi_perCluster_.push_back(iphi_perCluster);
        sc_ieta_perCluster_.push_back(ieta_perCluster);
        sc_energy_perCluster_.push_back(energy_perCluster);
        sc_time_perCluster_.push_back(time_perCluster);
        sc_clusterType_.push_back(clusterType_perCluster);
        
        // Find primary tracking particle for this mergedcluster
        float primary_energy = -999.0;
        float primary_et = -999.0;
        float primary_phi = -999.0;
        float primary_eta = -999.0;
        int primary_pdgId = 0;
        
        // Take any tracking particle and climb back to original ancestor
        if (!sc.trackingParticles().empty()) {
            auto tp_ref = sc.trackingParticles()[0]; // Take first TP
            
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
        
        sc_primary_energy_.push_back(primary_energy);
        sc_primary_et_.push_back(primary_et);
        sc_primary_phi_.push_back(primary_phi);
        sc_primary_eta_.push_back(primary_eta);
        sc_primary_pdgId_.push_back(primary_pdgId);
        
        // Print detailed info for first few events
        if (evt_event_ <= 3) {
            std::cout << "MergedCluster: MeV, t=" << time << " ns;"
                      << " E = " << energy << " MeV;"
                      << " pos = (" << sc.simPos().x() << "," << sc.simPos().y() << ")" << std::endl;
            std::cout << "  Made from " << nClusters << " clusters:" << std::endl;
            for (const auto& detId : sc.detIds()) {
                std::cout << "    DetId: " << detId.rawId() << std::endl;
            }
        }
    }
    
    // // Process individual clusters
    // cluster_n_ = 0;
    // for (const auto& detSet : *mtdSimLCHandle) {
    //     for (const auto& cluster : detSet) {
    //         if (cluster.energy() < 1.0) continue; 
            
    //         cluster_n_++;
            
    //         // h_cluster_energy_->Fill(cluster.energy());
    //         h_cluster_time_->Fill(cluster.time());
            
    //         // Check if cluster was merged
    //         bool inMergedCluster = clusterDetIds.count(cluster.id().rawId()) > 0;
    //         h_merging_efficiency_->Fill(cluster.energy(), inMergedCluster ? 1 : 0);
            
    //         cluster_energy_.push_back(cluster.energy());
    //         cluster_time_.push_back(cluster.time());
    //         cluster_x_.push_back(cluster.x());
    //         cluster_y_.push_back(cluster.y());
    //         cluster_detId_.push_back(cluster.id().rawId());
    //         cluster_inMergedCluster_.push_back(inMergedCluster);
    //     }
    // }
    
    tree_->Fill();
    
    if (evt_event_ <= 3) {
        std::cout << "Event " << evt_event_ << ": " << sc_n_ << " MergedClusters, " 
                  << cluster_n_ << " individual clusters" << std::endl;
    }
}

void SimMergedClusterValidationExample::endJob() {
    std::cout << "\n=== MTD MergedCluster Validation Summary ===" << std::endl;
    std::cout << "Total events processed: " << tree_->GetEntries() << std::endl;
}

DEFINE_FWK_MODULE(SimMergedClusterValidationExample);