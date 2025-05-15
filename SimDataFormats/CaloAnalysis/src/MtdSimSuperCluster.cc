#include "SimDataFormats/CaloAnalysis/interface/MtdSimSuperCluster.h"
#include <algorithm>
#include <iostream>

MtdSimSuperCluster::MtdSimSuperCluster(const MtdSimLayerClusterRef& clusterRef, const TrackingParticleRef& tpRef) {
  addCluster(clusterRef, tpRef);
}

void MtdSimSuperCluster::addCluster(const MtdSimLayerClusterRef& clusterRef, const TrackingParticleRef& tpRef) {
  clusters_.push_back(clusterRef);
  trackingParticles_.push_back(tpRef);

  // Sort clusters by time (earliest first)
  std::vector<MtdSimLayerClusterRef> sortedRefs(clusters_.begin(), clusters_.end());
  std::sort(sortedRefs.begin(), sortedRefs.end(), [](const MtdSimLayerClusterRef& a, const MtdSimLayerClusterRef& b) { return a->simLCTime() > b->simLCTime(); });
  clusters_.clear();
  for(auto& ref : sortedRefs)
      clusters_.push_back(ref);

  // std::sort(clusters_.begin(), clusters_.end(),
  //           [](MtdSimLayerClusterRef& a, MtdSimLayerClusterRef& b) {
  //             return a->simLCTime() < b->simLCTime();
  //           });
}

float MtdSimSuperCluster::simTime() const {
  if(clusters_.empty()) {
    return -999;
  } else{
    return (*clusters_.begin())->simLCTime();
  }
}

LocalPoint MtdSimSuperCluster::simPos() const {
  if(clusters_.empty()) {
    return LocalPoint(-999, -999, -999);
  } else {
    return (*clusters_.begin())->simLCPos();
  }
}

std::vector<DetId> MtdSimSuperCluster::detIds() const {
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

std::ostream& operator<<(std::ostream& s, const MtdSimSuperCluster& sc) {
  s << "MtdSimSuperCluster with " << sc.clusters_.size() << " clusters and TrackingParticles: = " << sc.trackingParticles_.size() << "\n";
  // for (const auto& clu : sc.clusters_) {
  //   s << "  Cluster time = " << clu->simLCTime() << "\n";
  // }
  // print also time and position using methods
  s << "Earliest time = " << sc.simTime() << "\n";
  s << "Earliest position = " << sc.simPos() << "\n";

  return s;
}