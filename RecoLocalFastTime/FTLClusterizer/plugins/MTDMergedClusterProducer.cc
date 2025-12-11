#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/FTLRecHit/interface/FTLMergedClusterCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDTopology.h"

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>

class MTDMergedClusterProducer : public edm::stream::EDProducer<> {
public:
    explicit MTDMergedClusterProducer(const edm::ParameterSet& conf);
    ~MTDMergedClusterProducer() override = default;
    
    void produce(edm::Event& e, const edm::EventSetup& es) override;
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions); 

private:
    edm::EDGetTokenT<FTLClusterCollection> btlClustersToken_;
    std::string btlMergedClusterInstance_;

    double timeThreshold_;  // 10 sigma timing window?
    double energyThreshold_;   
    
    edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
    edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;    

    bool areTimingCompatible(const FTLCluster* c1, const FTLCluster* c2);
    FTLMergedCluster mergeClusters(const std::vector<const FTLCluster*>& clusters, const BTLDetId& seedId, edm::Handle<FTLClusterCollection> btlClustersHandle);
};

MTDMergedClusterProducer::MTDMergedClusterProducer(const edm::ParameterSet& conf) 
    : btlClustersToken_(consumes<FTLClusterCollection>(conf.getParameter<edm::InputTag>("btlBarrel"))),
      btlMergedClusterInstance_(conf.getParameter<std::string>("btlMergedClusterInstance")),
      timeThreshold_(conf.getParameter<double>("timeThreshold")),
      energyThreshold_(conf.getParameter<double>("energyThreshold")),
      mtdgeoToken_(esConsumes<MTDGeometry, MTDDigiGeometryRecord>()),
      mtdtopoToken_(esConsumes<MTDTopology, MTDTopologyRcd>()) {

    produces<FTLMergedClusterCollection>(btlMergedClusterInstance_);

    std::cout << "MTDMergedClusterProducer: Time threshold = " << timeThreshold_ << " sigma" << std::endl;
    std::cout << "MTDMergedClusterProducer: Energy threshold = " << energyThreshold_ << " MeV" << std::endl;
    std::cout << "=== MTDMergedClusterProducer: Constructor finished ===" << std::endl;
}

bool MTDMergedClusterProducer::areTimingCompatible(const FTLCluster* c1, const FTLCluster* c2) {
    double timeDiff = std::abs(c1->time() - c2->time());
    double timeError1 = c1->timeError();
    double timeError2 = c2->timeError();
    double combinedError = std::sqrt(timeError1*timeError1 + timeError2*timeError2);
    
    bool compatible = timeDiff < (timeThreshold_ * combinedError);

    std::cout << "    Timing check: dt=" << timeDiff << " ns, combined_error=" 
              << combinedError << " ns, limit=" << (timeThreshold_ * combinedError) 
              << " ns -> " << (compatible ? "COMPATIBLE" : "INCOMPATIBLE") << std::endl;
    
    return compatible;
}

FTLMergedCluster MTDMergedClusterProducer::mergeClusters(const std::vector<const FTLCluster*>& clusters, const BTLDetId& seedId, edm::Handle<FTLClusterCollection> btlClustersHandle){
    float totalEnergy = 0;
    float weightedTime = 0;
    float weightedTimeError2 = 0;
    float weightedX = 0;
    float weightedY = 0;
    std::vector<DetId> clusterIds;
    std::vector<FTLClusterRef> clusterRefs;

    for (const auto* cluster : clusters) {
        float energy = cluster->energy();
        
        totalEnergy += energy;
        weightedTime += energy * cluster->time();
        weightedTimeError2 += energy * energy * cluster->timeError() * cluster->timeError();
        weightedX += energy * cluster->x();
        weightedY += energy * cluster->y();
        clusterIds.push_back(cluster->id());
        clusterRefs.push_back(edmNew::makeRefTo(btlClustersHandle, cluster)); // if needed
    }

    float avgTime = weightedTime / totalEnergy;
    float avgTimeError = std::sqrt(weightedTimeError2) / totalEnergy;
    float avgX = weightedX / totalEnergy;
    float avgY = weightedY / totalEnergy;

    return FTLMergedCluster(seedId, totalEnergy, avgTime, avgTimeError, avgX, avgY, clusterIds, clusterRefs);
}

void MTDMergedClusterProducer::produce(edm::Event& e, const edm::EventSetup& es) {    
    // Get topology for navigation
    auto topologyHandle = es.getTransientHandle(mtdtopoToken_);
    const MTDTopology* topology = topologyHandle.product();
    
    edm::Handle<FTLClusterCollection> btlClustersHandle;
    e.getByToken(btlClustersToken_, btlClustersHandle);
    
    std::cout << "MTDMergedClusterProducer: Processing event " << e.id() << std::endl;
    std::cout << "Time threshold: " << timeThreshold_ << " sigma, Energy threshold: " << energyThreshold_ << " MeV" << std::endl;
    std::cout << "Input BTL clusters: " << btlClustersHandle->size() << std::endl;

    auto btlOutput = std::make_unique<FTLMergedClusterCollection>();
    
    if (!btlClustersHandle.isValid() || btlClustersHandle->size() == 0) {
        std::cout << "No valid BTL clusters found in event" << std::endl;
        e.put(std::move(btlOutput), btlMergedClusterInstance_);
        return;
    }
    
    // collect clusters and sort by module ID 
    std::vector<const FTLCluster*> allClusters;
    std::map<BTLDetId, const FTLCluster*> clusterMap; // fast lookup

    for (const auto& detSet : *btlClustersHandle) {
        for (const auto& cluster : detSet) {
            if (cluster.energy() < energyThreshold_) continue;
            allClusters.push_back(&cluster);
            clusterMap[cluster.id()] = &cluster;
        }
    }

    // sorting...
    std::sort(allClusters.begin(), allClusters.end(), [](const FTLCluster* a, const FTLCluster* b) {
        BTLDetId idA = a->id();
        BTLDetId idB = b->id();

        // side -> rod -> module
        if (idA.zside() != idB.zside()) return idA.zside() < idB.zside();
        if (idA.mtdRR() != idB.mtdRR()) return idA.mtdRR() < idB.mtdRR();
        return idA.module() < idB.module();
    });

    std::cout << "Found " << allClusters.size() << " BTL clusters above energy threshold" << std::endl;

    std::set<const FTLCluster*> processedClusters;
    
    // Process clusters 
    for (const auto* cluster : allClusters) {
        if (processedClusters.count(cluster)) continue; // already part of a mergedcluster
            
        BTLDetId cluId = cluster->id();
        /*std::cout << "Cluster 1: E=" << cluster.energy() 
            << " t=" << cluster.time() 
            << " x=" << cluster.x() 
            << " y=" << cluster.y();*/
            
        std::vector<const FTLCluster*> mergedClusterClusters = {cluster};
        processedClusters.insert(cluster);
            
        // Check for edge hits in current cluster
        bool edgeHitIn0 = false;
        bool edgeHitIn15 = false;
            
        for (int i = 0; i < cluster->size(); ++i) {
            auto hit = cluster->hit(i);
            int hit_col = hit.y();
            if (hit_col == 0) {
                edgeHitIn0 = true;
            } else if (hit_col == 15) {
                edgeHitIn15 = true;
            }
        }
            
        bool hasEdgeHitCurrent = edgeHitIn0 || edgeHitIn15;
        std::cout << "  Edge hits: col0=" << edgeHitIn0 << ", col15=" << edgeHitIn15 << std::endl;
            
        // Get topology indices
        std::pair<uint32_t, uint32_t> indices = topology->btlIndex(cluId.rawId());
        uint32_t iphi = indices.first;
        uint32_t ieta = indices.second;
        std::cout << "  BTL indices: iphi=" << iphi << ", ieta=" << ieta << std::endl;
            
        // ETA DIRECTION MERGING
        if (hasEdgeHitCurrent && iphi != std::numeric_limits<uint32_t>::max() && ieta != std::numeric_limits<uint32_t>::max()) {
                
            std::vector<int> etaOffsets = {1, -1};
            for (int etaOffset : etaOffsets) {
                uint32_t adjDetIdRaw = topology->btlidFromIndex(iphi, ieta + etaOffset);
                if (adjDetIdRaw == 0) continue;
                    
                BTLDetId adjDetId(adjDetIdRaw);
                auto it = clusterMap.find(adjDetId);
                if (it == clusterMap.end()) continue;
                if (processedClusters.count(it->second)) continue;
                    
                const FTLCluster* adjCluster = it->second;
                std::cout << "  Found eta neighbor at " << adjDetId.rawId() << std::endl;
                std::cout << "  Neighbor cluster: E=" << adjCluster->energy() 
                            << " t=" << adjCluster->time() 
                            << " x=" << adjCluster->x() 
                            << " y=" << adjCluster->y();
                    
                // Check for opposite edge hit
                bool hasOppositeEdgeHit = false;
                for (int j = 0; j < adjCluster->size(); ++j) {
                    auto hit = adjCluster->hit(j);
                    int hit_col = hit.y();
                    if ((edgeHitIn0 && hit_col == 15) || (edgeHitIn15 && hit_col == 0)) {
                        hasOppositeEdgeHit = true;
                        break;
                    }
                }
                    
                if (hasOppositeEdgeHit && areTimingCompatible(cluster, adjCluster)) {
                    /*std::cout << "  -> MERGING ETA neighbor: " << cluId.rawId() 
                              << " with " << adjDetId.rawId() << std::endl;
                    std::cout << "  -> Cluster 1 module: " << cluId.module() 
                              << ", Cluster 2 module: " << adjDetId.module() << std::endl;*/
                    mergedClusterClusters.push_back(adjCluster);
                    processedClusters.insert(adjCluster);
                }
            }
        }
            
        // PHI DIRECTION MERGING
        /*if (iphi != std::numeric_limits<uint32_t>::max() && ieta != std::numeric_limits<uint32_t>::max()) {
                
            std::vector<int> phiOffsets = {1, -1};
            for (int phiOffset : phiOffsets) {
                uint32_t adjDetIdRaw = topology->btlidFromIndex(iphi + phiOffset, ieta);
                if (adjDetIdRaw == 0) continue;
                    
                BTLDetId adjDetId(adjDetIdRaw);
                auto it = clusterMap.find(adjDetId);
                if (it == clusterMap.end()) continue;
                if (processedClusters.count(it->second)) continue;
                    
                const FTLCluster* adjCluster = it->second;
                std::cout << "  Found phi neighbor at " << adjDetId.rawId() << std::endl;
                    
                if (areTimingCompatible(&cluster, adjCluster)) {
                    std::cout << "  -> MERGING PHI neighbor: " << cluId.rawId() 
                              << " with " << adjDetId.rawId() << std::endl;
                    mergedClusterClusters.push_back(adjCluster);
                    processedClusters.insert(adjCluster);
                }
            }
        }*/
            
        // Create mergedcluster from merged clusters
        if (mergedClusterClusters.size() > 0) {
            FTLMergedCluster mergedCluster = mergeClusters(mergedClusterClusters, cluId, btlClustersHandle);
            btlOutput->push_back(mergedCluster);
            
            std::cout << "Created MergedCluster from " << mergedClusterClusters.size() 
                    << " clusters: E=" << mergedCluster.energy() 
                    << " MeV, t=" << mergedCluster.time() << " ns" << std::endl;
        }        
    }
    
    std::cout << "About to put " << btlOutput->size() << " MergedClusters into event..." << std::endl;
    e.put(std::move(btlOutput), btlMergedClusterInstance_);
    std::cout << "=== Successfully put MergedClusters into event ===" << std::endl;
}
void MTDMergedClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("btlBarrel", edm::InputTag("mtdClusters", "FTLBarrel"));
    desc.add<std::string>("btlMergedClusterInstance", "FTLBarrel");    
    desc.add<double>("timeThreshold", 2.0);
    desc.add<double>("energyThreshold", 0.0);  // MeV
    descriptions.add("MTDMergedClusterProducer", desc);
}

DEFINE_FWK_MODULE(MTDMergedClusterProducer);

