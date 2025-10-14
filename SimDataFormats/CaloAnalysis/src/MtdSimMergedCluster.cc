#include "SimDataFormats/CaloAnalysis/interface/MtdSimMergedCluster.h"
#include <algorithm>
#include <iostream>

MtdSimMergedCluster::MtdSimMergedCluster(const MtdSimLayerClusterRef& clusterRef, const TrackingParticleRef& tpRef) {
    addCluster(clusterRef, tpRef);
}

void MtdSimMergedCluster::addCluster(const MtdSimLayerClusterRef& clusterRef, const TrackingParticleRef& tpRef) {
    clusters_.push_back(clusterRef);

    // check if tpRef is already in trackingParticles_
    auto it = std::find(trackingParticles_.begin(), trackingParticles_.end(), tpRef);
    if (it == trackingParticles_.end()) {
      trackingParticles_.push_back(tpRef);
    }
    // trackingParticles_.push_back(tpRef);

    // Sort clusters by time (earliest first)
    std::vector<MtdSimLayerClusterRef> sortedRefs(clusters_.begin(), clusters_.end());
    std::sort(sortedRefs.begin(), sortedRefs.end(), [](const MtdSimLayerClusterRef& a, const MtdSimLayerClusterRef& b) { return a->simLCTime() < b->simLCTime(); });
    clusters_.clear();
    for(auto& ref : sortedRefs)
        clusters_.push_back(ref);

    // std::sort(clusters_.begin(), clusters_.end(),
    //           [](MtdSimLayerClusterRef& a, MtdSimLayerClusterRef& b) {
    //             return a->simLCTime() < b->simLCTime();
    //           });
}

float MtdSimMergedCluster::simTime() const {
    if(clusters_.empty()) {
        return -999;
    } else {
        // FIRST IMPLEMENTATION: take time of earliest cluster
        return (*clusters_.begin())->simLCTime();
    }
}

LocalPoint MtdSimMergedCluster::simPos() const {
    if(clusters_.empty()) {
        return LocalPoint(-999, -999, -999);
    } else {
        // FIRST IMPLEMENTATION: take position of earliest cluster
        // (TO BE CHANGED: take energy-weighted position?)
        return (*clusters_.begin())->simLCPos();
    }
}

float MtdSimMergedCluster::simEnergy() const {
    float totalEnergy = 0;
    for (const auto& clu : clusters_) {
        // FIRST IMPLEMENTATION: sum energies of *all* clusters
        totalEnergy += clu->simLCEnergy();
    }
    return totalEnergy;
}

DetId MtdSimMergedCluster::simDetId() const {
    if(clusters_.empty()) {
        return -999;
    } else {
        // FIRST IMPLEMENTATION: take detId of earliest sim LC
        auto detIds_and_rows = (*clusters_.begin())->detIds_and_rows();
        return detIds_and_rows[0].first;
    }
}

std::vector<DetId> MtdSimMergedCluster::detIds() const {
    std::vector<DetId> ids;
    for (const auto& clu : clusters_) {
        const auto& clusterDetIds = clu->hits_and_fractions();
        for (const auto& hitFrac : clusterDetIds) {
            // check if hitFrac.first is already in ids
            if (std::find(ids.begin(), ids.end(), hitFrac.first) == ids.end()) {
                // if not, add it to the list
                ids.push_back(hitFrac.first);
            }
        }
    }

    return ids;
}

std::ostream& operator<<(std::ostream& s, const MtdSimMergedCluster& sc) {
    s << "MtdSimMergedCluster with " << sc.clusters_.size() << " clusters and TrackingParticles: = " << sc.trackingParticles_.size() << "\n";
    // for (const auto& clu : sc.clusters_) {
    //   s << "  Cluster time = " << clu->simLCTime() << "\n";
    // }
    // print also time and position using methods
    s << "Earliest time = " << sc.simTime() << "\n";
    s << "Earliest position = " << sc.simPos() << "\n";

    return s;
}