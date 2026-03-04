#ifndef Validation_MtdValidation_MergedClusterValidation_withAssociationMaps_h
#define Validation_MtdValidation_MergedClusterValidation_withAssociationMaps_h

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

#include "SimDataFormats/Associations/interface/MtdRecoMergedClusterToSimMergedClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimMergedClusterToRecoMergedClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdRecoMergedClusterToSimMergedClusterAssociatorBaseImpl.h"
#include "SimDataFormats/Associations/interface/MtdSimMergedClusterToTPAssociator.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class MergedClusterValidation_withAssociationMaps : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit MergedClusterValidation_withAssociationMaps(const edm::ParameterSet&);
    ~MergedClusterValidation_withAssociationMaps() = default;

private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void beginJob() override;
    void endJob() override;

    
    edm::EDGetTokenT<FTLMergedClusterCollection> mergedClustersToken_;
    edm::EDGetTokenT<FTLClusterCollection> clustersToken_;

    edm::EDGetTokenT<MtdSimMergedClusterCollection> simMergedClustersToken_;
    edm::EDGetTokenT<MtdSimLayerClusterCollection> simClustersToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
    edm::EDGetTokenT<reco::MergedRecoToSimCollectionMtd> mergedRecoToSimMap_;
    edm::EDGetTokenT<reco::MergedSimToRecoCollectionMtd> mergedSimToRecoMap_;

    edm::EDGetTokenT<reco::MergedSimToTPCollectionMtd> mergedSimToTPMap_;
    edm::EDGetTokenT<reco::TPToMergedSimCollectionMtd> mergedTPToSimMap_;
    
    edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
    edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;

    // reco - sim
    TH1F* h_deltaTime_;
    TH1F* h_deltaEnergy_;
    TH1F* h_deltaX_;
    TH1F* h_deltaY_;

    TH1F* h_deltaNclu_;
    TH1F* h_deltaNhits_;
    TH1F* h_nSimPerReco_;
    
    TH1F* h_deltaTime_multiClu_;
    TH1F* h_deltaEnergy_multiClu_;
    TH1F* h_deltaX_multiClu_;
    TH1F* h_deltaY_multiClu_;

    TH1F* h_deltaNclu_multiClu_;
    TH1F* h_deltaNhits_multiClu_;
    TH1F* h_nSimPerReco_multiClu_;

    TH1F* h_deltaTime_singleClu_;
    TH1F* h_deltaEnergy_singleClu_;
    TH1F* h_deltaX_singleClu_;
    TH1F* h_deltaY_singleClu_;

    // 2D
    TH2F* h_deltaTime_vs_Eta_;
    TH2F* h_deltaEnergy_vs_Eta_;

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