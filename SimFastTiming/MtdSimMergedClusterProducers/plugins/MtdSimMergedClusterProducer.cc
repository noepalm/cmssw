// Enable debug logging
#define EDM_ML_DEBUG

#include <FWCore/Framework/interface/one/EDProducer.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

// Logging
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "SimDataFormats/Associations/interface/TrackAssociation.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerClusterFwd.h"

#include "SimDataFormats/CaloAnalysis/interface/MtdSimMergedCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimMergedClusterFwd.h"

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
#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
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

    // Perform the user-defined action (e.g. fill mergedcluster)
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

class MtdSimMergedClusterProducer : public edm::one::EDProducer<edm::one::SharedResources> {
public:
    explicit MtdSimMergedClusterProducer(const edm::ParameterSet&);
    ~MtdSimMergedClusterProducer() override = default;

    void produce(edm::Event&, const edm::EventSetup&) override;

private:
    edm::EDGetTokenT<TrackingParticleCollection> trackingParticlesToken_;
    edm::EDGetTokenT<reco::TPToSimCollectionMtd> tpToSimClusMapToken_;
    edm::EDGetTokenT<reco::SimToTPCollectionMtd> simClusToTPMapToken_;
    edm::EDGetTokenT<MtdSimLayerClusterCollection> mtdSimLayerClustersToken_;    
    double minEnergy_;

    bool useTopologicalClustering_ = true;

    const edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
};

MtdSimMergedClusterProducer::MtdSimMergedClusterProducer(const edm::ParameterSet& iConfig)
    : trackingParticlesToken_(consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("trackingParticles"))),
      tpToSimClusMapToken_(consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMap"))),
      simClusToTPMapToken_(consumes<reco::SimToTPCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMap"))),
      mtdSimLayerClustersToken_(consumes<MtdSimLayerClusterCollection>(iConfig.getParameter<edm::InputTag>("mtdSimLayerClusters"))),
      minEnergy_(iConfig.getParameter<double>("minClusterEnergy")),
      useTopologicalClustering_(iConfig.getParameter<bool>("useTopologicalClustering")),
      mtdtopoToken_(esConsumes<MTDTopology, MTDTopologyRcd>()) {
    produces<MtdSimMergedClusterCollection>();
}

void MtdSimMergedClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    // std::cout << "MtdSimMergedClusterProducer::produce() called" << std::endl;

    // Get topology for navigation
    auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
    const MTDTopology* topology = topologyHandle.product();

    // Create output collection (MtdSimMergedCluster)
    auto outputClusters = std::make_unique<MtdSimMergedClusterCollection>();

    // Retrieve collections
    edm::Handle<TrackingParticleCollection> trackingParticles;
    iEvent.getByToken(trackingParticlesToken_, trackingParticles);

    edm::Handle<reco::TPToSimCollectionMtd> tpToSimClusMap;
    iEvent.getByToken(tpToSimClusMapToken_, tpToSimClusMap);

    edm::Handle<reco::SimToTPCollectionMtd> simClusToTPMap;
    iEvent.getByToken(simClusToTPMapToken_, simClusToTPMap);

    edm::Handle<MtdSimLayerClusterCollection> simLClusters;
    iEvent.getByToken(mtdSimLayerClustersToken_, simLClusters);

    // BUILD ANCESTOR MAP: to check if 2 TrackingParticles have the same root ancestor
    std::map<TrackingParticleRef, TrackingParticleRef> ancestorMap;
    for(size_t i = 0; i < trackingParticles->size(); ++i) {
        TrackingParticleRef tp(trackingParticles, i);
        
        // Find earliest ancestor for this TP
        TrackingParticleRef current = tp;
        std::set<TrackingParticleRef> visited;
        
        while (true) {
            // Prevent infinite loops
            if (visited.count(current)) {
                break;
            }
            visited.insert(current);
            
            // Check if this particle has parent vertices
            const auto& parentVertices = current->parentVertex();
            if (parentVertices.isNull() || !parentVertices.isAvailable()) {
                // No parent vertex - this is the root ancestor
                break;
            }
            // Get the parent tracks from the parent vertex
            const auto& parentTracks = parentVertices->sourceTracks();
            if (parentTracks.empty()) {
                // No parent tracks - this is the root ancestor
                break;
            }
            // Move to the first parent
            current = parentTracks[0];
        }
        // Store the mapping: TP -> earliest ancestor
        ancestorMap[tp] = current;
    }

    // ANCESTOR MAP SANITY CHECK: how many unique ancestors are there?
    std::set<TrackingParticleRef> uniqueAncestors;
    for (const auto& pair : ancestorMap) {
        uniqueAncestors.insert(pair.second);
    }
    edm::LogInfo("MtdSimMergedClusterProducer") << "Total TrackingParticles: " << trackingParticles->size() 
                                          << ", Unique primary ancestors: " << uniqueAncestors.size();

    // reserve memory for output collection
    // (worst case scenario: all TrackingParticles are primary)
    outputClusters->reserve(trackingParticles->size());

    // Create cluster map for fast lookup (can have multiple clusters per DetId)
    std::vector<const MtdSimLayerCluster*> allsimLClusters;
    std::vector<const MtdSimLayerCluster*> allsimETLLClusters;
    

    for (const auto& cluster : *simLClusters) {
        if (!cluster.detIds_and_rows().empty() && MTDDetId(cluster.detIds_and_rows()[0].first).mtdSubDetector() == MTDDetId::ETL) {
            allsimETLLClusters.push_back(&cluster);
        }
        else{
            allsimLClusters.push_back(&cluster);
        }

        // if (cluster.energy() >= minEnergy_) {
        //     // retrieve detId from first hit
        //     BTLDetId detId = cluster.detIds_and_rows()[0].first;
        //     // retrieve GEOGRAPHICAL id
        //     BTLDetId geoDetId = detId.geographicalId(BTLDetId::CrysLayout::v3);

        //     clusterMap[geoDetId].push_back(&cluster);
        // }

        
    }

    // Sort simLClusters collection
    std::sort(allsimLClusters.begin(), allsimLClusters.end(), [&topology](const MtdSimLayerCluster* a, const MtdSimLayerCluster* b) {
        const auto& detIdsA = a->detIds_and_rows();
        const auto& detIdsB = b->detIds_and_rows();

        BTLDetId idA(detIdsA[0].first);
        BTLDetId idB(detIdsB[0].first);

        // side -> rod -> module
        if (idA.zside() != idB.zside()) return idA.zside() < idB.zside();
        if (idA.mtdRR() != idB.mtdRR()) return idA.mtdRR() < idB.mtdRR();

        auto [iphiA, ietaA] = topology->btlIndex(idA.geographicalId(BTLDetId::CrysLayout::v3).rawId()); //uint32_t
        auto [iphiB, ietaB] = topology->btlIndex(idB.geographicalId(BTLDetId::CrysLayout::v3).rawId());

        if (iphiA != iphiB) return iphiA < iphiB;
        if (ietaA != ietaB) return ietaA < ietaB;

        // if (idA.module() != idB.module()) return idA.module() < idB.module();
        // if same module, sort by column range

        // Get min and max columns for both clusters
        auto [minA, maxA] = std::minmax_element(detIdsA.begin(), detIdsA.end(), [](auto const& x, auto const& y){
            return x.second.second < y.second.second;
        });
        auto [minB, maxB] = std::minmax_element(detIdsB.begin(), detIdsB.end(), [](auto const& x, auto const& y){
            return x.second.second < y.second.second;
        });

        int lowest_icol_A = minA->second.second;
        int highest_icol_A = maxA->second.second;
        int lowest_icol_B = minB->second.second;
        int highest_icol_B = maxB->second.second;

        if (lowest_icol_A != lowest_icol_B) return lowest_icol_A < lowest_icol_B;
        return highest_icol_A > highest_icol_B; // larger cluster first
    });

    // DEBUG: check ordering
    // std::cout << "DEBUG: Checking cluster ordering" << std::endl;
    for(const auto* cluster : allsimLClusters ){
        // retrieve iphi, ieta
        BTLDetId cluId(cluster->detIds_and_rows()[0].first);
        std::pair<uint32_t, uint32_t> indices = topology->btlIndex(cluId.geographicalId(BTLDetId::CrysLayout::v3).rawId());
        uint32_t iphi = indices.first;
        uint32_t ieta = indices.second;
        // std::cout << "CLUSTER #" << cluster - &(*simLClusters->begin()) << " zside = " << cluId.zside() << ", RR = " << cluId.mtdRR() << ", module = " << cluId.module() << "; iphi = " << iphi << ", ieta = " << ieta << " with hits at columns: ";
        for (auto const& detId_row_col : cluster->detIds_and_rows()) {
            int col = detId_row_col.second.second;
            // std::cout << col << " ";
        }
        // std::cout << std::endl;
    }

    // Now, construct cluster map from ordered collection
    std::map<BTLDetId, std::vector<const MtdSimLayerCluster*>> clusterMap;

    for (const auto* cluster : allsimLClusters) {
        if (cluster->energy() >= minEnergy_) {
            // retrieve detId from first hit
            BTLDetId detId = cluster->detIds_and_rows()[0].first;
            // retrieve GEOGRAPHICAL id
            BTLDetId geoDetId = detId.geographicalId(BTLDetId::CrysLayout::v3);

            clusterMap[geoDetId].push_back(cluster);
        }
    }


    edm::LogInfo("MtdSimMergedClusterProducer") << "Found " << clusterMap.size() << " MTD SimLayerClusters above energy threshold";

    if (useTopologicalClustering_) {
        edm::LogInfo("MtdSimMergedClusterProducer") << "Using TOPOLOGICAL + HISTORICAL clustering algorithm";

        // ------------------------------------------------------
        // NEW: HISTORY+TOPOLOGY IMPLEMENTATION

        std::set<const MtdSimLayerCluster*> processedClusters;

        // Process clusters with full merging logic
        for (const auto* clusterPointer : allsimLClusters) {
            const auto& cluster = *clusterPointer;

            if (cluster.energy() < minEnergy_ || processedClusters.count(&cluster)) continue;

            // // TEMPORARY: forget about ETL hits
            // std::cout << "DEBUG: cluster detId: " << std::endl;
            // std::cout << BTLDetId(cluster.detIds_and_rows()[0].first) << std::endl;
            // std::cout << "Subdetector: " << BTLDetId(cluster.detIds_and_rows()[0].first).mtdSubDetector() << " (BTL = " << MTDDetId::BTL << ", ETL = " << MTDDetId::ETL << ")" << std::endl;
            if (!cluster.detIds_and_rows().empty() && MTDDetId(cluster.detIds_and_rows()[0].first).mtdSubDetector() == MTDDetId::ETL) continue;
            
            BTLDetId cluId(cluster.detIds_and_rows()[0].first);

            // LogDebug("MtdSimMergedClusterProducer") << "Processing cluster DetId " << cluId.rawId() 
            //                                     << " with energy " << cluster.simLCEnergy() << " MeV";
            // std::cout << "------------------------------------------------" << std::endl;
            // std::cout << "Processing cluster DetId " << cluId.rawId() 
            //                                    << " with energy " << cluster.simLCEnergy() << " MeV" << std::endl;
            
            // Start with current cluster
            std::vector<const MtdSimLayerCluster*> mergedClusterClusters = {&cluster};
            processedClusters.insert(&cluster);
            
            // Check for edge hits in current cluster
            bool edgeHitIn0 = false;
            bool edgeHitIn15 = false;
            int leftmost_col = 20;
            int rightmost_col = -1;

            // iterate over detIds_and_rows:
            // LogDebug("MtdSimMergedClusterProducer") << "Iterating over " << cluster.detIds_and_rows().size() << " cluster hits";
            // std::cout << "Iterating over " << cluster.detIds_and_rows().size() << " cluster hits" << std::endl;
            for (auto const& detId_row_col : cluster.detIds_and_rows()) {
                int row = detId_row_col.second.first;
                int col = detId_row_col.second.second;

                // LogDebug("MtdSimMergedClusterProducer") << "  Cluster hit at row " << row << ", col " << col;
                // std::cout << "  Cluster hit at row " << row << ", col " << col << std::endl;

                if (col == 0) {
                    edgeHitIn0 = true;
                } else if (col == 15) {
                    edgeHitIn15 = true;
                }

                if (col < leftmost_col) leftmost_col = col;
                if (col > rightmost_col) rightmost_col = col;

            }
            
            bool hasEdgeHitCurrent = edgeHitIn0 || edgeHitIn15;
            // LogDebug("MtdSimMergedClusterProducer") << "  Edge hits: col0=" << edgeHitIn0 << ", col15=" << edgeHitIn15;
            // std::cout << "  Edge hits: col0=" << edgeHitIn0 << ", col15=" << edgeHitIn15 << std::endl;
            // LogDebug("MtdSimMergedClusterProducer") << "  hasEdgeHitCurrent = " << hasEdgeHitCurrent;
            // std::cout << "  hasEdgeHitCurrent = " << hasEdgeHitCurrent << std::endl;

            // Get topology indices - use geographicalId (module-level) with crystal layout
            // std::cout << "  DEBUG: getting BTL indices from geographicalId " << cluId.geographicalId(BTLDetId::CrysLayout::v3) << std::endl;
            std::pair<uint32_t, uint32_t> indices = topology->btlIndex(cluId.geographicalId(BTLDetId::CrysLayout::v3).rawId());
            uint32_t iphi = indices.first;
            uint32_t ieta = indices.second;
            // LogDebug("MtdSimMergedClusterProducer") << "  BTL indices: iphi=" << iphi << ", ieta=" << ieta;
            // std::cout << "  BTL indices: iphi=" << iphi << ", ieta=" << ieta << std::endl;

            // ETA DIRECTION MERGING
            std::vector<int> etaOffsets = {1, 0, -1};
            for (int etaOffset : etaOffsets) {
                // std::cout << "  Eta offset = " << etaOffset << std::endl;
                if ((hasEdgeHitCurrent && iphi != std::numeric_limits<uint32_t>::max() && ieta != std::numeric_limits<uint32_t>::max()) || etaOffset == 0) {
                    // LogDebug("MtdSimMergedClusterProducer") << "  Attempting eta-direction merging...";
                    // std::cout << "  Attempting eta-direction merging..." << std::endl;
                    uint32_t adjDetIdRaw = topology->btlidFromIndex(iphi, ieta + etaOffset);
                    // LogDebug("MtdSimMergedClusterProducer") << "    Checking adjacent detId at index (" << iphi << ", " << ieta + etaOffset << "): " << adjDetIdRaw;
                    // std::cout << "    Checking adjacent detId at index (" << iphi << ", " << ieta + etaOffset << "): " << adjDetIdRaw << std::endl;
                    if (adjDetIdRaw == 0) continue;
                    
                    BTLDetId adjDetId(adjDetIdRaw);
                    auto it = clusterMap.find(adjDetId);
                    if (it == clusterMap.end()){
                        // LogDebug("MtdSimMergedClusterProducer") << "    No cluster found at this detId";
                        // std::cout << "    No cluster found at this detId" << std::endl;
                        continue;
                    }
                    
                    // Iterate over all clusters at this DetId
                    for (const MtdSimLayerCluster* adjCluster : it->second) {
                        if (processedClusters.count(adjCluster)){
                            // LogDebug("MtdSimMergedClusterProducer") << "    Cluster found but already processed";
                            // std::cout << "    Cluster found but already processed" << std::endl;
                            continue;
                        }
                        
                        // LogDebug("MtdSimMergedClusterProducer") << "    Found eta neighbor at " << adjDetId.rawId();
                        // std::cout << "    Found eta neighbor at " << adjDetId.rawId() << std::endl;

                        // Check for opposite edge hit
                        bool hasOppositeEdgeHit = false;
                        int adj_leftmost_col = 20;
                        int adj_rightmost_col = -1;

                        for (auto const& detId_row_col : adjCluster->detIds_and_rows()) {
                            int row = detId_row_col.second.first;
                            int col = detId_row_col.second.second;

                            // LogDebug("MtdSimMergedClusterProducer") << "      adjacent cluster hit at row " << row << ", col " << col;
                            // std::cout << "      adjacent cluster hit at row " << row << ", col " << col << std::endl;

                            if ((edgeHitIn0 && col == 15 && etaOffset == -1) || (edgeHitIn15 && col == 0 && etaOffset == 1)) {
                                hasOppositeEdgeHit = true;
                            }

                            if (col < adj_leftmost_col) adj_leftmost_col = col;
                            if (col > adj_rightmost_col) adj_rightmost_col = col;
                        }

                        // bool hasAdjacentCols = (abs(leftmost_col - adj_rightmost_col) == 1) || (abs(rightmost_col - adj_leftmost_col) == 1);
                        int clu_len = rightmost_col - leftmost_col + 1;

                        // bool isLeftEdgeOverlapping = (abs(adj_leftmost_col  - leftmost_col) <= clu_len) && (abs(rightmost_col - adj_rightmost_col) <= clu_len); // equality includes adjacent clusters
                        // bool isRightEdgeOverlapping = (abs(adj_rightmost_col - rightmost_col) <= clu_len) && (abs(leftmost_col - adj_leftmost_col) <= clu_len);

                        bool isLeftEdgeOverlapping = (abs(adj_leftmost_col  - leftmost_col) <= clu_len) && (abs(rightmost_col - adj_leftmost_col) <= clu_len); // equality includes adjacent clusters
                        bool isRightEdgeOverlapping = (abs(adj_rightmost_col - rightmost_col) <= clu_len) && (abs(leftmost_col - adj_rightmost_col) <= clu_len);

                        bool areClustersOverlapping = isLeftEdgeOverlapping || isRightEdgeOverlapping;

                        // if (areClustersOverlapping && etaOffset == 0) {
                            // std::cout << "      SAME MODULE MERGING: Found overlapping clusters in same module!" << std::endl;
                            // std::cout << "      SAME MODULE MERGING: leftmost_col = " << leftmost_col << ", rightmost_col = " << rightmost_col
                            //          << "; adj_leftmost_col = " << adj_leftmost_col << ", adj_rightmost_col = " << adj_rightmost_col << std::endl;
                        // }

                        if (hasOppositeEdgeHit || (areClustersOverlapping && etaOffset == 0)){ // && areTimingCompatible(&cluster, adjCluster)) { // FORGET ABOUT TIME COMPATIBILITY FOR NOW
                            // check for common ancestor
                            bool hasCommonAncestor = false;
                            const auto& simLayerClusters1 = simClusToTPMap->find(MtdSimLayerClusterRef(simLClusters, &cluster - &(*simLClusters->begin())));
                            const auto& simLayerClusters2 = simClusToTPMap->find(MtdSimLayerClusterRef(simLClusters, adjCluster - &(*simLClusters->begin())));
                            if (simLayerClusters1 != simClusToTPMap->end() && simLayerClusters2 != simClusToTPMap->end()) {
                                for (const auto& tpRef1 : simLayerClusters1->val) {
                                    for (const auto& tpRef2 : simLayerClusters2->val) {
                                        if (ancestorMap[tpRef1] == ancestorMap[tpRef2]) {
                                            hasCommonAncestor = true;
                                            break;
                                        }
                                    }
                                    if (hasCommonAncestor) break;
                                }
                            }

                            if (hasCommonAncestor) {
                                // LogDebug("MtdSimMergedClusterProducer") << "  -> MERGING ETA neighbor: " << cluId.rawId() 
                                //         << " with " << adjDetId.rawId();
                                // std::cout << "  -> MERGING ETA neighbor: " << cluId.rawId() 
                                //        << " with " << adjDetId.rawId() << std::endl;
                                bool isBackscatterMergedcluster = mergedClusterClusters[0]->trackIdOffset() == 3;
                                bool areBothBackscatter = isBackscatterMergedcluster && (adjCluster->trackIdOffset() == 3);
                                bool areBothNotBackscatter = !isBackscatterMergedcluster && (adjCluster->trackIdOffset() != 3);
                                if (areBothBackscatter || areBothNotBackscatter) {
                                    // std::cout << "  CLUSTER MERGING: offset of first = " << mergedClusterClusters[0]->trackIdOffset()
                                    //            << "( isBackscatterMerged Cluster " << isBackscatterMergedcluster << ")"
                                    //            << ", adjCluster->trackIdOffset() = " << adjCluster->trackIdOffset() 
                                    //            << ", merged cluster size = " << mergedClusterClusters.size()
                                    //            << std::endl;
                                    mergedClusterClusters.push_back(adjCluster);
                                    processedClusters.insert(adjCluster);
                                    
                                    // update lowest/highest icol
                                    if (etaOffset == 0) {
                                        // std::cout << "  SAME MODULE MERGING: current col range = ( " << leftmost_col << ", " << rightmost_col << ")" << std::endl;
                                        if (adj_leftmost_col < leftmost_col) leftmost_col = adj_leftmost_col;
                                        if (adj_rightmost_col > rightmost_col) rightmost_col = adj_rightmost_col;
                                        // std::cout << "  SAME MODULE MERGING: updated col range = ( " << leftmost_col << ", " << rightmost_col << ")" << std::endl;
                                    }

                                } 
                                // else {
                                //     // std::cout << "  NOT MERGING CLUSTER THOUGH WE SHOULD HAVE!: offset of first = " << mergedClusterClusters[0]->trackIdOffset()
                                //                 << "( isBackscatterMerged Cluster " << isBackscatterMergedcluster << ")"
                                //                 << ", adjCluster->trackIdOffset() = " << adjCluster->trackIdOffset()
                                //                 << ", merged cluster size = " << mergedClusterClusters.size() << std::endl;
                                // }
                            } else {
                                LogDebug("MtdSimMergedClusterProducer") << "    Not merging: no common ancestor found";
                            }

                        }
                    }
                } else {
                    LogDebug("MtdSimMergedClusterProducer") << "  No edge hit, invalid indices or same module check: not attempting eta merging";
                    LogDebug("MtdSimMergedClusterProducer") << "    hasEdgeHitCurrent = " << hasEdgeHitCurrent 
                                << ", iphi = " << iphi << " (is max? " << (iphi == std::numeric_limits<uint32_t>::max()) << "), ieta = " << ieta << " (is max? " << (ieta == std::numeric_limits<uint32_t>::max()) << ")";
                }
            }
            
            // // PHI DIRECTION MERGING
            // if (iphi != std::numeric_limits<uint32_t>::max() && ieta != std::numeric_limits<uint32_t>::max()) {
                
            //     std::vector<int> phiOffsets = {1, -1};
            //     for (int phiOffset : phiOffsets) {
            //         uint32_t adjDetIdRaw = topology->btlidFromIndex(iphi + phiOffset, ieta);
            //         if (adjDetIdRaw == 0) continue;
                    
            //         BTLDetId adjDetId(adjDetIdRaw);
            //         auto it = clusterMap.find(adjDetId);
            //         if (it == clusterMap.end()) continue;
            //         if (processedClusters.count(it->second)) continue;
                    
            //         const MtdSimLayerCluster* adjCluster = it->second;
                    // std::cout << "  Found phi neighbor at " << adjDetId.rawId() << std::endl;
                    
            //         // if (areTimingCompatible(&cluster, adjCluster)) {
                        // std::cout << "  -> MERGING PHI neighbor: " << cluId.rawId() 
            //         //                 << " with " << adjDetId.rawId() << std::endl;
            //         //     mergedClusterClusters.push_back(adjCluster);
            //         //     processedClusters.insert(adjCluster);
            //         // }

                    // std::cout << "  -> MERGING PHI neighbor: " << cluId.rawId() 
            //                     << " with " << adjDetId.rawId() << std::endl;
            //         mergedClusterClusters.push_back(adjCluster);
            //         processedClusters.insert(adjCluster);

            //     }
            // }
            
            // Create mergedcluster from merged clusters
            MtdSimMergedCluster simMergedCluster;

            // for each cluster, find associated TPs and add to mergedcluster using Sim to TP map
            for (const auto& simLayerCluster : mergedClusterClusters) {
                // create simLC reference by finding the index in the original collection
                size_t clusterIndex = simLayerCluster - &(*simLClusters->begin());
                MtdSimLayerClusterRef simLayerClusterRef(simLClusters, clusterIndex);
                const auto& TPs = simClusToTPMap->find(simLayerClusterRef);
                if (TPs != simClusToTPMap->end()) {
                    for (const auto& tpRef : TPs->val) {
                        simMergedCluster.addCluster(simLayerClusterRef, tpRef);
                    }
                } else {
                    LogDebug("MtdSimMergedClusterProducer") << "No TP associated to cluster index " << clusterIndex;
                    // add null reference
                    simMergedCluster.addCluster(simLayerClusterRef, TrackingParticleRef());
                }
            }

            outputClusters->push_back(simMergedCluster);
            LogDebug("MtdSimMergedClusterProducer") << "Created MergedCluster from " << mergedClusterClusters.size() 
                        << " clusters: E=" << simMergedCluster.simEnergy() 
                        << " MeV, t=" << simMergedCluster.simTime() << " ns";
        }
        for (const auto* clusterPointer : allsimETLLClusters) {
            const auto& cluster = *clusterPointer;

            // // TEMPORARY: forget about ETL hits
            // std::cout << "DEBUG: cluster detId: " << std::endl;
            // std::cout << BTLDetId(cluster.detIds_and_rows()[0].first) << std::endl;
            // std::cout << "Subdetector: " << BTLDetId(cluster.detIds_and_rows()[0].first).mtdSubDetector() << " (BTL = " << MTDDetId::BTL << ", ETL = " << MTDDetId::ETL << ")" << std::endl;
            MtdSimMergedCluster simMergedCluster;
            // create simLC reference by finding the index in the original collection
            size_t clusterIndex = clusterPointer - &(*simLClusters->begin());
            MtdSimLayerClusterRef simLayerClusterRef(simLClusters, clusterIndex);
            const auto& TPs = simClusToTPMap->find(simLayerClusterRef);
            if (TPs != simClusToTPMap->end()) {
                for (const auto& tpRef : TPs->val) {
                    simMergedCluster.addCluster(simLayerClusterRef, tpRef);
                }
            } else {
                LogDebug("MtdSimMergedClusterProducer") << "No TP associated to cluster index " << clusterIndex;
                // add null reference
                simMergedCluster.addCluster(simLayerClusterRef, TrackingParticleRef());
            }
            outputClusters->push_back(simMergedCluster);
            LogDebug("MtdSimMergedClusterProducer") << "Created ETL MergedCluster from an ETL cluster" << " E=" << simMergedCluster.simEnergy() 
                        << " MeV, t=" << simMergedCluster.simTime() << " ns";

        }
    } else {
        edm::LogInfo("MtdSimMergedClusterProducer") << "Using HISTORY-ONLY clustering algorithm";

        // OLD: HISTORY-ONLY IMPLEMENTATION
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

            // create 2 mergedclusters: one for primary+secondary+looper clusters, one for backscatter hits 
            MtdSimMergedCluster simMergedCluster(tp);
            MtdSimMergedCluster backscatterCluster(tp);

            std::set<edm::Ref<TrackingParticleCollection>> visited;

            traverseDecayTree(tp, visited, 
                [&](const TrackingParticleRef& ref) {                
                    const auto& simLayerClusters = tpToSimClusMap->find(ref);
                    if(simLayerClusters != tpToSimClusMap->end()){
                        for (const auto& simLayerCluster : simLayerClusters->val) {
                            if (simLayerCluster->simLCEnergy() > minEnergy_) {
                                //TEMPORARY: also check if cluster is in BTL
                                if (!simLayerCluster->detIds_and_rows().empty() && MTDDetId(simLayerCluster->detIds_and_rows()[0].first).mtdSubDetector() == MTDDetId::BTL){
                                    if (simLayerCluster->trackIdOffset() == 3)
                                        backscatterCluster.addCluster(simLayerCluster, ref);
                                    else
                                        simMergedCluster.addCluster(simLayerCluster, ref);
                                }
                            }
                        }
                    }
                }
            );

            outputClusters->push_back(simMergedCluster);
            outputClusters->push_back(backscatterCluster);
        }
    }

    // print all MtdSimMergedClusters
    edm::LogInfo("MtdSimMergedClusterProducer") << "Found " << outputClusters->size() << " MtdSimMergedClusters";
    for (const auto& mergedCluster : *outputClusters) {
        LogDebug("MtdSimMergedClusterProducer") << mergedCluster;
        auto detids = mergedCluster.detIds();
    }

    iEvent.put(std::move(outputClusters));
}

DEFINE_FWK_MODULE(MtdSimMergedClusterProducer);
