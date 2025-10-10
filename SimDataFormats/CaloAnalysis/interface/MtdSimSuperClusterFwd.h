#ifndef CaloAnalysis_MtdSimSuperClusterFwd_h
#define CaloAnalysis_MtdSimSuperClusterFwd_h
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/RefVector.h"
#include <vector>

#include "SimDataFormats/CaloAnalysis/interface/MtdSimSuperCluster.h"
//class MtdSimSuperCluster;
typedef std::vector<MtdSimSuperCluster> MtdSimSuperClusterCollection;
typedef edm::Ref<MtdSimSuperClusterCollection> MtdSimSuperClusterRef;
typedef edm::RefVector<MtdSimSuperClusterCollection> MtdSimSuperClusterRefVector;
typedef edm::RefProd<MtdSimSuperClusterCollection> MtdSimSuperClusterRefProd;
typedef edm::RefVector<MtdSimSuperClusterCollection> MtdSimSuperClusterContainer;

std::ostream &operator<<(std::ostream &s, MtdSimSuperCluster const &tp);

#endif
