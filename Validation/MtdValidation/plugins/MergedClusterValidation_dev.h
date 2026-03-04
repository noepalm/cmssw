#ifndef Validation_MtdValidation_MergedClusterValidation_dev_h
#define Validation_MtdValidation_MergedClusterValidation_dev_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDTopology.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "DataFormats/FTLRecHit/interface/FTLMergedClusterCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/MtdSimMergedCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimMergedClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class MergedClusterValidation_dev : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit MergedClusterValidation_dev(const edm::ParameterSet&);
    ~MergedClusterValidation_dev() = default;

private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void beginJob() override;
    void endJob() override;

    
    edm::EDGetTokenT<FTLMergedClusterCollection> mergedClustersToken_;
    edm::EDGetTokenT<FTLClusterCollection> clustersToken_;

    edm::EDGetTokenT<MtdSimMergedClusterCollection> simMergedClustersToken_;
    edm::EDGetTokenT<MtdSimLayerClusterCollection> simClustersToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
    edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;

    // RECO
    TH1F* h_mc_energy_;
    TH1F* h_mc_time_;
    TH1F* h_mc_timeError_;
    TH1F* h_mc_x_;
    TH1F* h_mc_y_;
    TH1F* h_mc_nClusters_;
    TH1F* h_cluster_energy_;
    TH1F* h_cluster_time_;
    
    TH1F* h_eta_adjacent_pairs_;
    TH1F* h_eta_merged_pairs_;
    TH1F* h_eta_merging_efficiency_;

    TH2F* h_mc_energy_vs_time_;
    TH2F* h_mc_xy_;
    TH2F* h_mc_energy_vs_nClusters_;
    TH2F* h_merging_efficiency_;

    // SIM
    TH1F* h_simmc_energy_;
    TH1F* h_simmc_logEnergy_;
    TH1F* h_simmc_time_;
    TH1F* h_simmc_x_;
    TH1F* h_simmc_y_;
    TH1F* h_simmc_eta_;
    TH1F* h_simmc_nClusters_;
    TH1F* h_simmc_n_;

    TH1F* h_simmc_logEnergy_perCluster_;
    TH1F* h_simmc_time_perCluster_;
    TH1F* h_simmc_clusterType_;

    TH2F* h_simmc_xy_;
    TH2F* h_simmc_energy_vs_time_;
    TH2F* h_simmc_energy_vs_nClusters_;
    TH2F* h_simmc_primaryPt_vs_nClusters_;
    TH2F* h_simmc_primaryPt_vs_energy_;
    TH2F* h_simmc_primaryEnergy_vs_energy_;
    TH2F* h_simmc_primaryEnergy_vs_nClusters_;

    TTree* tree_;
    
    int evt_run_, evt_event_;

    // RECO
    int mc_n_;
    std::vector<float> mc_energy_, mc_time_, mc_timeError_, mc_x_, mc_y_;
    std::vector<int> mc_nClusters_;
    std::vector<uint32_t> mc_seedId_;
    std::vector<std::vector<uint32_t>> mc_clusterIds_;

    int cluster_n_;
    std::vector<float> cluster_energy_, cluster_time_, cluster_x_, cluster_y_;
    std::vector<uint32_t> cluster_detId_;
    std::vector<bool> cluster_inMergedCluster_;

    // SIM
    int simmc_n_;
    std::vector<float> simmc_energy_, simmc_time_, simmc_timeError_, simmc_x_, simmc_y_, simmc_eta_;
    std::vector<int> simmc_nClusters_, simmc_nModules_;
    std::vector<std::vector<uint32_t>> simmc_iphi_perCluster_, simmc_ieta_perCluster_;
    std::vector<std::vector<float>> simmc_energy_perCluster_, simmc_time_perCluster_;
    std::vector<std::vector<float>> simmc_earliestHitTime_perCluster_;
    std::vector<std::vector<uint32_t>> simmc_clusterType_;
    std::vector<std::vector<std::vector<int>>> simmc_hitCols_perCluster_;
    
    // Primary particle information per mergedcluster
    std::vector<float> simmc_primary_energy_, simmc_primary_et_, simmc_primary_phi_, simmc_primary_eta_;
    std::vector<int> simmc_primary_pdgId_;



    int totalAdjacentPairs_;
    int totalMergedPairs_;
};

#endif