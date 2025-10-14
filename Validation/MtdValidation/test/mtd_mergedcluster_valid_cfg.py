import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("MergedClusterValidationExample", Phase2C17I13M9)

# Load standard configurations
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.MessageLogger.cerr.threshold = 'INFO'

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:/eos/user/p/pakrap/MTD/CMSSW_SC_TEST/src/RecoLocalFastTime/FTLClusterizer/mergedcluster_reco_singlepi.root'
        'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/CMSSW_15_1_0_pre2/src/RecoLocalFastTime/FTLClusterizer/test/mtdMergedClusters_numEvent1000.root'
    )
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Output file for histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("mtd_recomergedcluster_validation.root")
)

# Validation analyzer
process.mtdMergedClusterValidation = cms.EDAnalyzer("MergedClusterValidationExample",
    mergedClusters = cms.InputTag("mtdMergedClusters", "FTLBarrel"),
    clusters = cms.InputTag("mtdClusters", "FTLBarrel")
)

# Path
process.p = cms.Path( 
    process.mtdMergedClusterValidation
)

print("Running MTD MergedCluster Validation...")