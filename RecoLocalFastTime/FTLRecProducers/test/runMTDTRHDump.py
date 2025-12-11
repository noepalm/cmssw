import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("MTDTrackDump", Phase2C17I13M9)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/eos/home-n/npalmeri/MTD/MTD_supercluster/association_maps_merge/src/RecoLocalFastTime/FTLClusterizer/test/mtdMergedClusters_mergetest_numEvent1000.root')  
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.mtdTrackingRecHitDump = cms.EDAnalyzer("MTDTrackingRecHitDump",
    src = cms.InputTag("mtdTrackingRecHits")  
)

process.p = cms.Path(process.mtdTrackingRecHitDump)