import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("MtdSuperClusOnly", Phase2C17I13M9)

# Load standard services and geometry
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")  # Needed for TrackingParticles
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Input file (must contain SimClusters + TrackingParticles)
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring([f"root://eosuser.cern.ch///eos/user/n/npalmeri/ntuples/MTD/PhotonReco/crab_MTDPhotonReco/CRAB_UserFiles/SingleGammaFlatPt0p1To10_Run4D110_aging1000_noPU_MTDPhotonReco/250424_080407/0000/step2_{i}.root" for i in range(1, 11)]),
    fileNames = cms.untracked.vstring([f"file:/eos/cms/store/relval/CMSSW_15_1_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_noPU-v1/2580000/03ec5b66-690c-415a-9602-362b351d2a08.root"]), #photon gun

)

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Load your producer
from SimFastTiming.MtdSimSuperClusterProducers.mtdSimSuperClusterProducer_cfi import mtdSimSuperClusterProducer
process.mtdSimSuperClusterProducer = mtdSimSuperClusterProducer.clone()
# process.superClusterSequence = cms.Sequence(process.mtdSimSuperClusterProducer)

# Output module (optional - to store results)
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("mtdSimSuperClusters.root"),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_mtdSimSuperClusterProducer_*_*"
    )
)

# Load MTD truth map associators (needed for the producer)
from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi import mtdSimLayerClusterToTPAssociatorByTrackId
from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi import mtdSimLayerClusterToTPAssociation

# clone and add label 
process.mtdSimLayerClusterToTPAssociatorByTrackId = mtdSimLayerClusterToTPAssociatorByTrackId.clone()
process.mtdSimLayerClusterToTPAssociation = mtdSimLayerClusterToTPAssociation.clone()

process.superClusterSequence = cms.Sequence(process.mtdSimLayerClusterToTPAssociatorByTrackId + process.mtdSimLayerClusterToTPAssociation + process.mtdSimSuperClusterProducer)

# process.superClusterSequence = cms.Sequence(process.mtdSimSuperClusterProducer)

# Execution path
process.p = cms.Path(process.superClusterSequence)
process.out_step = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p, process.out_step)
