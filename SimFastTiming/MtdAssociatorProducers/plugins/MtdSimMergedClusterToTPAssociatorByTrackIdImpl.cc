//
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "MtdSimMergedClusterToTPAssociatorByTrackIdImpl.h"

using namespace reco;
using namespace std;

/* Constructor */
MtdSimMergedClusterToTPAssociatorByTrackIdImpl::MtdSimMergedClusterToTPAssociatorByTrackIdImpl(
    edm::EDProductGetter const& productGetter, reco::SimToTPCollectionMtd simToTPMap, reco::TPToSimCollectionMtd tpToSimMap)
    : productGetter_(&productGetter),
      simToTPMap_(simToTPMap),
      tpToSimMap_(tpToSimMap){}

//
//---member functions
//

reco::MergedSimToTPCollectionMtd MtdSimMergedClusterToTPAssociatorByTrackIdImpl::associateSimToTP(
    const edm::Handle<MtdSimMergedClusterCollection>& simMergedClusH,
    const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {
  MergedSimToTPCollectionMtd outputCollection(productGetter_);

  const auto& simMergedClusters = *simMergedClusH.product();
  // const auto& trackingParticles = *trackingParticleH.product();

  // // -- Loop over tracking particles and build a temporary map of trackId, eventId  --> tpRef
  // // FIXME: do we still need this?
  // std::map<std::pair<unsigned int, uint32_t>, TrackingParticleRef> tpIdMap;
  // for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
  //   const auto& tp = *tpIt;
  //   EncodedEventId tpEventId = tp.eventId();
  //   for (unsigned int igt = 0; igt < tp.g4Tracks().size(); igt++) {
  //     unsigned int tpTrackId = tp.g4Tracks()[igt].trackId();
  //     TrackingParticleRef tpRef =
  //         edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIt - trackingParticles.begin());
  //     tpIdMap[std::make_pair(tpTrackId, tpEventId.rawId())] = tpRef;
  //   }
  // }

  std::vector<TrackingParticleRef> tpRefs;

  for (auto simMergedClusIt = simMergedClusters.begin(); simMergedClusIt != simMergedClusters.end(); simMergedClusIt++) {
    const auto& simMergedClus = *simMergedClusIt;
    MtdSimMergedClusterRef simMergedClusterRef = edm::Ref<MtdSimMergedClusterCollection>(simMergedClusH, &simMergedClus - &(*simMergedClusH->begin()));

    std::vector<TrackingParticleRef> associatedTPs;

    // std::cout << "Iterating on SimMergedCluster #" << simMergedClusterRef.key() << " with E=" << simMergedClus.simEnergy() << " MeV" << ", containing #" << simMergedClus.clusters().size() << "sim clusters" << std::endl;
    // for each sim merged cluster, loop over component sim clusters
    // for each component sim cluster, find associated TP and fill output collection
    // iterate over component clusters
    for (const auto& simClusterRef : simMergedClus.clusters()) {
      // query the simCluster -> TP map to find associated tracking particles
      // std::cout << "  Iterating on component simCluster #" << simClusterRef.key() << " with E=" << simClusterRef->simEnergy() << " MeV" << std::endl;
      auto simToTPIt = simToTPMap_.find(simClusterRef);
      if (simToTPIt != simToTPMap_.end()) {
        // std::cout << "    Found " << std::distance(simToTPIt->val.begin(), simToTPIt->val.end()) << " TP matches to simCluster with E=" << simClusterRef->simEnergy() << " MeV" << std::endl;
        for (const auto& tpRef : simToTPIt->val) {
          // std::cout << "    Found associated TP with pdgId=" << tpRef->pdgId() << ", pt=" << tpRef->pt() << " GeV, eta=" << tpRef->eta() << ", phi=" << tpRef->phi() << std::endl;
          associatedTPs.push_back(tpRef);
          // outputCollection.insert(simMergedClusterRef, tpRef);
        }
      } else {
        // std::cout << "    No TP matches found for simCluster with E=" << simClusterRef->simEnergy() << " MeV" << std::endl;
      }
    }

    // remove duplicates from associatedTPs
    std::sort(associatedTPs.begin(), associatedTPs.end());
    associatedTPs.erase(std::unique(associatedTPs.begin(), associatedTPs.end()), associatedTPs.end());
    for (auto const& tpRef : associatedTPs) outputCollection.insert(simMergedClusterRef, tpRef);

  }

  return outputCollection;
}

reco::TPToMergedSimCollectionMtd MtdSimMergedClusterToTPAssociatorByTrackIdImpl::associateTPToSim(
    const edm::Handle<MtdSimMergedClusterCollection>& simMergedClusH,
    const edm::Handle<TrackingParticleCollection>& trackingParticleH) const {
  TPToMergedSimCollectionMtd outputCollection(productGetter_);

  // -- get the collections
  const auto& simMergedClusters = *simMergedClusH.product();
  const auto& trackingParticles = *trackingParticleH.product();

  // make preliminary map: simCluster => simMergedCluster
  std::map<MtdSimLayerClusterRef, std::vector<MtdSimMergedClusterRef>> simClusToMergedMap;
  for (const auto& simMergedClus : simMergedClusters) {
    for (const auto& simClusRef : simMergedClus.clusters()) {
      simClusToMergedMap[simClusRef].push_back(MtdSimMergedClusterRef(simMergedClusH, &simMergedClus - &(*simMergedClusH->begin())));
    }
  }

  for (auto tpIt = trackingParticles.begin(); tpIt != trackingParticles.end(); tpIt++) {
    size_t tpIndex = tpIt - trackingParticles.begin();
    TrackingParticleRef tpRef = edm::Ref<TrackingParticleCollection>(trackingParticleH, tpIndex);

    // std::cout << "Looking for matches to TP #" << tpRef.key() << " with pt=" << tpIt->pt() << " GeV, eta=" << tpIt->eta() << ", phi=" << tpIt->phi() << std::endl;

    std::vector<MtdSimMergedClusterRef> simMergedClusRefs;

    // find SimClus associated to TPs
    auto TPtoSimIt = tpToSimMap_.find(tpRef);
    if (TPtoSimIt != tpToSimMap_.end()) {
      // std::cout << "  Found " << std::distance(TPtoSimIt->val.begin(), TPtoSimIt->val.end()) << " simCluster matches to this TP." << std::endl;
      for (auto simClusRef : TPtoSimIt->val) {
        // find simMergedCluster querying map
        auto const& simMergedClusters = simClusToMergedMap.find(simClusRef);
        if (simMergedClusters == simClusToMergedMap.end()) {
          // std::cout << "WARNING: no simMergedCluster found for simClusterRef with key " << simClusRef.key() << std::endl;
          continue;
        }
        for (const auto& simMergedClusterRef : simMergedClusters->second){
          // std::cout << "  Found associated SimMergedCluster #" << simMergedClusterRef.key() << " with E=" << simMergedClusterRef->simEnergy() << " MeV" << std::endl;
          simMergedClusRefs.push_back(simMergedClusterRef);
          // outputCollection.insert(tpRef, simMergedClusterRef);    
        }
      }
    }

    // remove duplicates from simMergedClusRefs
    std::sort(simMergedClusRefs.begin(), simMergedClusRefs.end());
    simMergedClusRefs.erase(std::unique(simMergedClusRefs.begin(), simMergedClusRefs.end()), simMergedClusRefs.end());
    for (auto const& simMergedClusterRef : simMergedClusRefs) outputCollection.insert(tpRef, simMergedClusterRef);

  }

  // save final output
  return outputCollection;
}