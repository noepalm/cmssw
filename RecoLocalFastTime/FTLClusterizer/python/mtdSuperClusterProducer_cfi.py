import FWCore.ParameterSet.Config as cms

mtdSuperClusters = cms.EDProducer("MTDSuperClusterProducer",
    btlBarrel = cms.InputTag("mtdClusters", "FTLBarrel"),
    timeThreshold = cms.double(5.0),
    energyThreshold = cms.double(0.0),
    btlSuperClusterInstance = cms.string("FTLBarrel")
)