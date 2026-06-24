import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9
process = cms.Process('mtdHarvesting',Phase2C22I13M9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtendedRun4D121Reco_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

# Input source
process.source = cms.Source("DQMRootSource",
    # fileNames = cms.untracked.vstring('file:step3_inDQM.root')
    # fileNames = cms.untracked.vstring('file:/eos/home-n/npalmeri/MTD/MTD_mergedCluster_validation/CMSSW_20_0_0_pre1/src/step3_inDQM.root')
    # # EDITED, pions
    # fileNames = cms.untracked.vstring('file:/eos/home-n/npalmeri/MTD/MTD_mergedCluster_validation/samples/CMSSW_20_0_0_pre1_mergedClusters/singlepi_0PU/step3_inDQM_pions.root')
    # EDITED, ttbar
    fileNames = cms.untracked.vstring('file:/eos/home-n/npalmeri/MTD/MTD_mergedCluster_validation/samples/CMSSW_20_0_0_pre1_mergedClusters/ttbar_200PU/step3_inDQM.root')
    # # CLEAN, pions
    # fileNames = cms.untracked.vstring('file:/eos/home-n/npalmeri/MTD/MTD_mergedCluster_validation/samples/CMSSW_20_0_0_pre1_vanilla/singlepi_0PU/step3_inDQM_pions.root')
    # # CLEAN, ttbar
    # fileNames = cms.untracked.vstring('file:/eos/user/t/tipaulet/MTDTask/mergedCluster/noemicleaned/latest/vanilla/CMSSW_20_0_0_pre1/sample/ttbarPU/step3_inDQM.root')
)

# Path and EndPath definitions

process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# --- PostProcessing

process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdEleIsoPostProcessor_cfi")
process.load("Validation.MtdValidation.Primary4DVertexPostProcessor_cfi")

process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.MtdEleIsoPostProcessor + process.Primary4DVertexPostProcessor)
process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.MtdTracksPostProcessor)

process.p = cms.Path( process.harvesting )

process.schedule = cms.Schedule( process.edmtome_step , process.p , process.dqmsave_step )
