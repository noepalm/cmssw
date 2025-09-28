#include <FWCore/Framework/interface/one/EDProducer.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerClusterFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/MtdSimSuperCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimSuperClusterFwd.h"

// MTD truth association maps
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToRecoClusterAssociationMap.h"

// Geometry and topology
#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDTopology.h"

// DetId
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"

#include <memory>
#include <set>

void traverseDecayTree(const edm::Ref<TrackingParticleCollection>& tpRef,
                       std::set<edm::Ref<TrackingParticleCollection>>& visited,
                       const std::function<void(const edm::Ref<TrackingParticleCollection>&)>& action) {

    // stop condition
    if (visited.count(tpRef))
        return;
    
    // keep track of visited particles
    visited.insert(tpRef);

    // Perform the user-defined action (e.g. fill supercluster)
    action(tpRef);

    const auto& decayVtxs = tpRef->decayVertices();
    if (decayVtxs.size() > 0) {
        // iterate using begin, end explicitly
        for (auto it = decayVtxs.begin(); it != decayVtxs.end(); ++it) {
            const auto& decayVtx = *it;
            // std::cout << "\tDecay vertex: r = " << decayVtx->position().rho() << ", z = " << decayVtx->position().z() << std::endl;
            for (const auto& daughterRef : decayVtx->daughterTracks()){
                // std::cout << "\tDaughter track: pdgId = " << daughterRef->pdgId() << std::endl;
                traverseDecayTree(daughterRef, visited, action);       
            }
        }
    }
}


class MtdSimSuperClusterProducer : public edm::one::EDProducer<edm::one::SharedResources> {
public:
    explicit MtdSimSuperClusterProducer(const edm::ParameterSet&);
    ~MtdSimSuperClusterProducer() override = default;

    void produce(edm::Event&, const edm::EventSetup&) override;

private:
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
    edm::EDGetTokenT<reco::TPToSimCollectionMtd> tpToSimClusMapToken_;
    edm::EDGetTokenT<MtdSimLayerClusterCollection> mtdSimLayerClustersToken_;    
    double minEnergy_;

    const edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
};

MtdSimSuperClusterProducer::MtdSimSuperClusterProducer(const edm::ParameterSet& iConfig)
    : trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
      tpToSimClusMapToken_(consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMap"))),
      mtdSimLayerClustersToken_(consumes<MtdSimLayerClusterCollection>(iConfig.getParameter<edm::InputTag>("mtdSimLayerClusters"))),
      minEnergy_(iConfig.getParameter<double>("minClusterEnergy")),
      mtdtopoToken_(esConsumes<MTDTopology, MTDTopologyRcd>()) {
    produces<MtdSimSuperClusterCollection>();
}

void MtdSimSuperClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // std::cout << "MtdSimSuperClusterProducer::produce() called" << std::endl;

    // Get topology for navigation
    auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
    const MTDTopology* topology = topologyHandle.product();

    // Create output collection (MtdSimSuperCluster)
    auto outputClusters = std::make_unique<MtdSimSuperClusterCollection>();

    // Retrieve collections
    edm::Handle<TrackingParticleCollection> trackingParticles;
    iEvent.getByToken(trackingParticlesToken_, trackingParticles);

    edm::Handle<reco::TPToSimCollectionMtd> tpToSimClusMap;
    iEvent.getByToken(tpToSimClusMapToken_, tpToSimClusMap);

    edm::Handle<MtdSimLayerClusterCollection> simLClusters;
    iEvent.getByToken(mtdSimLayerClustersToken_, simLClusters);

    // reserve memory for output collection
    // (worst case scenario: all TrackingParticles are primary)
    outputClusters->reserve(trackingParticles->size());

    // Create cluster map for fast lookup
    // FIXME: not BTL only
    std::map<BTLDetId, const MtdSimLayerCluster*> clusterMap;
    for (const auto& cluster : *simLClusters) {
        if (cluster.energy() >= minEnergy_) {
            // retrieve detId from first hit
            BTLDetId detId = cluster.detIds_and_rows()[0].first;
            clusterMap[detId] = &cluster;
        }
    }

    std::cout << "Found " << clusterMap.size() << " BTL clusters above energy threshold" << std::endl;

    for(size_t i = 0; i < trackingParticles->size(); ++i) {
        TrackingParticleRef tp(trackingParticles, i);

        // check if particle is primary and comes from first interaction vertex
        if (tp->status() != 1) continue;
        if (tp->genParticles().size() < 1) continue;
        if (tp->g4Tracks().size() > 0){
            if (tp->g4Tracks().front().vertIndex() != 0) continue;
        } else {
            continue;
        }

        // create supercluster
        MtdSimSuperCluster simSuperCluster(tp);

        std::set<edm::Ref<TrackingParticleCollection>> visited;

        traverseDecayTree(tp, visited, 
            [&](const TrackingParticleRef& ref) {                
                const auto& simLayerClusters = tpToSimClusMap->find(ref);
                if(simLayerClusters != tpToSimClusMap->end()){
                    for (const auto& simLayerCluster : simLayerClusters->val) {
                        if (simLayerCluster->simLCEnergy() > minEnergy_)
                            simSuperCluster.addCluster(simLayerCluster, ref);
                    }
                }
            }
        );

        outputClusters->push_back(simSuperCluster);
    }

    // print all MtdSimSuperClusters
    std::cout << "Found " << outputClusters->size() << " MtdSimSuperClusters" << std::endl;
    for (const auto& superCluster : *outputClusters) {
        std::cout << superCluster << std::endl;
        auto detids = superCluster.detIds();
    }

    iEvent.put(std::move(outputClusters));
}

DEFINE_FWK_MODULE(MtdSimSuperClusterProducer);
