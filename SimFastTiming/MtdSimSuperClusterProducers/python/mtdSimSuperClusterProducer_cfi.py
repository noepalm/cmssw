import FWCore.ParameterSet.Config as cms

mtdSimSuperClusterProducer = cms.EDProducer(
    "MtdSimSuperClusterProducer",
    trackingParticles = cms.InputTag("mix", "MergedTrackTruth"),
    tp2SimAssociationMap = cms.InputTag("mtdSimLayerClusterToTPAssociation"),
    minClusterEnergy = cms.double(0)  # GeV
)