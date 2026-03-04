import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

options.register('inputFile',
                #  'file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_mergedclusterDev/src/SimFastTiming/MtdSimMergedClusterProducers/test/mtdSimMergedClusters_history_1k.root',
                #  'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/CMSSW_15_1_0_pre2/src/RecoLocalFastTime/FTLClusterizer/test/mtdMergedClusters_numEvent1000.root',
                #  'file:/eos/user/p/pakrap/MTD/CMSSW_SC_TEST/src/RecoLocalFastTime/FTLClusterizer/mergedcluster_reco_singlepi.root'
                 'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/CMSSW_15_1_0_pre2/src/RecoLocalFastTime/FTLClusterizer/test/mtdMergedClusters_forDev261125_numEvent1000.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input ROOT file")

# Parse command line arguments
options.parseArguments()

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("MergedClusterValidation", Phase2C17I13M9)


# Load standard configurations
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.MessageLogger.cerr.threshold = 'INFO'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFile)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file:" + options.outputFile)
)

# Validation analyzer
process.mtdMergedClusterValidation = cms.EDAnalyzer("MergedClusterValidation_dev",
    mergedClusters = cms.InputTag("mtdMergedClusters", "FTLBarrel"),
    clusters = cms.InputTag("mtdClusters", "FTLBarrel"),
    simMergedClusters = cms.InputTag("mtdSimMergedClusterProducer"),
    simLayerClusters = cms.InputTag("mix", "MergedMtdTruthLC")
)

# Path
process.p = cms.Path( 
    process.mtdMergedClusterValidation
)

print("Running MTD MergedCluster Validation...")
print(f"Input file: {options.inputFile}")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")