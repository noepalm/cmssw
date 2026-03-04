#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

using namespace reco;
using namespace std;

/* Constructor */

MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl::MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl(
    edm::EDProductGetter const& productGetter, double energyCut, double timeCut, mtd::MTDGeomUtil& geomTools, 
    reco::SimToRecoCollectionMtd simToRecoMap, reco::RecoToSimCollectionMtd recoToSimMap)
    : productGetter_(&productGetter), energyCut_(energyCut), timeCut_(timeCut), geomTools_(geomTools), simToRecoMap_(simToRecoMap),
    recoToSimMap_(recoToSimMap) {}

//
//---member functions
//

reco::MergedRecoToSimCollectionMtd MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl::associateRecoToSim(
    const edm::Handle<FTLMergedClusterCollection>& btlRecoClusH,
    const edm::Handle<FTLMergedClusterCollection>& etlRecoClusH,
    const edm::Handle<MtdSimMergedClusterCollection>& simMergedClusH) const {

  MergedRecoToSimCollectionMtd outputCollection;

  // -- get the collections
  // FIXME: reintroduce ETL when RecoMergedCluster available
  std::array<edm::Handle<FTLMergedClusterCollection>, 1> inputRecoMergedClusH{{btlRecoClusH}};

  const auto& simMergedClusters = *simMergedClusH.product();

  // make preliminary map: simCluster => simMergedCluster
  std::map<MtdSimLayerClusterRef, std::vector<MtdSimMergedClusterRef>> simClusToMergedMap;
  for (const auto& simMergedClus : simMergedClusters) {
    for (const auto& simClusRef : simMergedClus.clusters()) {
      simClusToMergedMap[simClusRef].push_back(MtdSimMergedClusterRef(simMergedClusH, &simMergedClus - &(*simMergedClusH->begin())));
    }
  }      
      
  for (auto const& recoMergedClusH : inputRecoMergedClusH) {
    for (const auto& detSet : *recoMergedClusH) {
      for (const auto& recoMergedClus : detSet) {

        // FIXME: this print
        // LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl") << "Reco cluster : " << clusId;

        FTLMergedClusterRef recoMergedClusterRef = edmNew::makeRefTo(recoMergedClusH, &recoMergedClus);
        std::vector<MtdSimMergedClusterRef> simClusterRefs;

        edm::LogWarning("MtdR2SAssoc") << "RecoMergedClus key=" << recoMergedClusterRef.key()
            << " nComponents=" << recoMergedClus.clusterRefs().size();

        // iterate over component clusters
        for (const auto& recoClusRef : recoMergedClus.clusterRefs()) {
          edm::LogWarning("MtdR2SAssoc") << "  Component id=" << recoClusRef.id() 
              << " key=" << recoClusRef.key();
          auto recoToSimIt = recoToSimMap_.equal_range(recoClusRef);
          if (recoToSimIt.first == recoToSimIt.second) {
            edm::LogWarning("MtdR2SAssoc") << "  -> NOT FOUND in recoToSimMap (size=" << recoToSimMap_.size() << ")";
            if (!recoToSimMap_.empty()) {
              edm::LogWarning("MtdR2SAssoc") << "  -> First entry id=" << recoToSimMap_.begin()->first.id()
                  << " key=" << recoToSimMap_.begin()->first.key();
            }
            // LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl") << "  No sim clusters associated to this reco cluster";
            continue;
          }
          const auto& simClusterCandidates = (*recoToSimIt.first).second;
          for (const auto& simClusterRef : simClusterCandidates){
            // retrieve simMergedClusters associated to this simLayerCluster
            auto const& simMergedClusters = simClusToMergedMap.find(simClusterRef);
            if (simMergedClusters == simClusToMergedMap.end()) {
              // Note: SHOULD NEVER HAPPEN.
              // LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl") << "  No sim merged clusters associated to this sim layer cluster";
              continue;
            }


            for (const auto& simMergedClusterRef : simMergedClusters->second){
              // LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl") << "  Found associated sim merged cluster: " << simMergedClusterRef.key();
              float dE = recoClusRef->energy() * 0.001 / simMergedClusterRef->simEnergy();  // reco cluster energy is in MeV!
              float dtSig = std::abs((recoClusRef->time() - simMergedClusterRef->simTime()) / recoClusRef->timeError());
      
              LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl")
                  << "E_recoClus = " << recoClusRef->energy() << "   E_simClus = " << simMergedClusterRef->simEnergy()
                  << "   E_recoClus/E_simClus = " << dE;
              LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl")
                  << "(t_recoClus-t_simClus)/sigma_t = " << dtSig;
      
              // FIXME: when reintroducing ETL, only consider dtSig cut for those
              if (dE < energyCut_ && dtSig < timeCut_) {  
                edm::LogWarning("MtdR2SAssoc") << "  -> MATCH PASSES dE=" << dE << " dtSig=" << dtSig;
                simClusterRefs.push_back(simMergedClusterRef);
              } else {
                edm::LogWarning("MtdR2SAssoc") << "  -> MATCH REJECTED dE=" << dE << " dtSig=" << dtSig
                    << " energyCut=" << energyCut_ << " timeCut=" << timeCut_;
              }
            }

          }
        }
        
        // -- Now fill the output collection
        // first remove duplicates from simClusterRefs
        std::sort(simClusterRefs.begin(), simClusterRefs.end());
        simClusterRefs.erase(std::unique(simClusterRefs.begin(), simClusterRefs.end()), simClusterRefs.end());
        outputCollection.emplace_back(recoMergedClusterRef, simClusterRefs);
      }
    }
  }

  outputCollection.post_insert();
  return outputCollection;
}

reco::MergedSimToRecoCollectionMtd MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl::associateSimToReco(
    const edm::Handle<FTLMergedClusterCollection>& btlRecoClusH,
    const edm::Handle<FTLMergedClusterCollection>& etlRecoClusH,
    const edm::Handle<MtdSimMergedClusterCollection>& simMergedClusH) const {
  MergedSimToRecoCollectionMtd outputCollection;

  // -- get the collections
  const auto& simMergedClusters = *simMergedClusH.product();

  // TODO: reintroduce ETL when RecoMergedCluster available
  std::array<edm::Handle<FTLMergedClusterCollection>, 1> inputH{{btlRecoClusH}};

  // make preliminary map: recoCluster => recoMergedCluster
  std::map<FTLClusterRef, std::vector<FTLMergedClusterRef>> recoClusToMergedMap;


  for (const auto& recoMergedClusH : inputH) {
    for (const auto& detSet : *recoMergedClusH) {
      for (const auto& recoMergedClus : detSet) {
        for (const auto& recoClusRef : recoMergedClus.clusterRefs()) {
          recoClusToMergedMap[recoClusRef].push_back(edmNew::makeRefTo(recoMergedClusH, &recoMergedClus));
        }
      }
    }
  }

  // -- loop over MtdSimMergedClusters
  for (auto simMergedClusIt = simMergedClusters.begin(); simMergedClusIt != simMergedClusters.end(); simMergedClusIt++) {
    const auto& simMergedClus = *simMergedClusIt;

    // query the simToReco map and retrieve list of reco clusters
    auto const& simMergedClusterRef = edm::Ref<MtdSimMergedClusterCollection>(simMergedClusH, &simMergedClus - &(*simMergedClusters.begin()));
    std::vector<FTLMergedClusterRef> recoMergedClusterRefs;

    // iterate over component clusters
    for (const auto& simClusterRef : simMergedClus.clusters()) {
      auto simToRecoIt = simToRecoMap_.equal_range(simClusterRef);
      if (simToRecoIt.first != simToRecoIt.second) {
        const auto& recoRefs = (*simToRecoIt.first).second;
        for (const auto& recoRef : recoRefs) {
          // QUESTION: no constraint on compatibility here?
          // retrieve recoMergedClusters associated to this recoCluster
          auto const& recoMergedClusters = recoClusToMergedMap.find(recoRef);
          
          if (recoMergedClusters == recoClusToMergedMap.end()) {
            // Note: SHOULD NEVER HAPPEN.
            // LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl") << "  No reco merged clusters associated to this reco cluster";
            continue;
          }

          for (const auto& recoMergedClusterRef : recoMergedClusters->second){
            LogDebug("MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsImpl") << "  Found associated reco merged cluster: " << recoMergedClusterRef.key();
            recoMergedClusterRefs.push_back(recoMergedClusterRef);
          }
          
        }
      }
    }

    std::sort(recoMergedClusterRefs.begin(), recoMergedClusterRefs.end());
    recoMergedClusterRefs.erase(std::unique(recoMergedClusterRefs.begin(), recoMergedClusterRefs.end()), recoMergedClusterRefs.end());

    outputCollection.emplace_back(simMergedClusterRef, recoMergedClusterRefs);
    
  }  // -- end loop over sim clusters

  outputCollection.post_insert();
  return outputCollection;
}
