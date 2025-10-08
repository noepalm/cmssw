#ifndef Validation_MtdValidation_SimSuperClusterValidationExample_h
#define Validation_MtdValidation_SimSuperClusterValidationExample_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// #include "DataFormats/FTLRecHit/interface/FTLSuperCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimSuperCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimSuperClusterFwd.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"

class SimSuperClusterValidationExample : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit SimSuperClusterValidationExample(const edm::ParameterSet&);
    ~SimSuperClusterValidationExample() = default;

private:
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void beginJob() override;
    void endJob() override;

    edm::EDGetTokenT<MtdSimSuperClusterCollection> superClustersToken_;
    edm::EDGetTokenT<MtdSimLayerClusterCollection> clustersToken_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    const edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
    const edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;

    TH1F* h_sc_energy_;
    TH1F* h_sc_time_;
    TH1F* h_sc_timeError_;
    TH1F* h_sc_x_;
    TH1F* h_sc_y_;
    TH1F* h_sc_nClusters_;
    TH1F* h_cluster_energy_;
    TH1F* h_cluster_time_;
    
    TH2F* h_sc_energy_vs_time_;
    TH2F* h_sc_xy_;
    TH2F* h_sc_energy_vs_nClusters_;
    TH2F* h_merging_efficiency_;

    TTree* tree_;
    
    int evt_run_, evt_event_, evt_number_;
    int sc_n_;
    std::vector<float> sc_energy_, sc_time_, sc_timeError_, sc_x_, sc_y_;
    std::vector<int> sc_nClusters_;
    std::vector<std::vector<uint32_t>> sc_iphi_perCluster_, sc_ieta_perCluster_;
    std::vector<std::vector<float>> sc_energy_perCluster_, sc_time_perCluster_;
    std::vector<std::vector<uint32_t>> sc_clusterType_;
    
    // Primary particle information per supercluster
    std::vector<float> sc_primary_energy_, sc_primary_et_, sc_primary_phi_, sc_primary_eta_;
    std::vector<int> sc_primary_pdgId_;
    
    int cluster_n_;
    std::vector<float> cluster_energy_, cluster_time_, cluster_x_, cluster_y_;
    std::vector<uint32_t> cluster_detId_;
    std::vector<bool> cluster_inSuperCluster_;
    
    // Event counter for actual processed event number
    int processed_event_counter_;
};

#endif