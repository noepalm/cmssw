import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

# Define custom options
options.register('inputFile',
                #  'file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_mergedclusterDev/src/SimFastTiming/MtdSimMergedClusterProducers/test/mtdSimMergedClusters_history_1k.root',
                 'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/CMSSW_15_1_0_pre2/src/RecoLocalFastTime/FTLClusterizer/test/mtdMergedClusters_numEvent1000.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input ROOT file")

# outputFile, maxEvents registered by default

# Parse command line arguments
options.parseArguments()

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
    fileNames = cms.untracked.vstring(options.inputFile)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

# Output file for histograms
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file:" + options.outputFile)
)

# process.mtdMergedClusters = cms.EDProducer("MTDMergedClusterProducer",
#     srcBarrel = cms.InputTag("mtdClusters", "FTLBarrel"),
#     BarrelMergedClusterName = cms.string("FTLBarrel"),
#     timeThreshold = cms.double(10.0),
#     energyThreshold = cms.double(1.0)
# )

# Validation analyzer
process.mtdSimMergedClusterValidation = cms.EDAnalyzer("SimMergedClusterValidationExample",
    simMergedClusters = cms.InputTag("mtdSimMergedClusterProducer"), #tag name? Not sure where it's specified
    simLayerClusters = cms.InputTag("mixData:MergedMtdTruthLC"),
    # mergedClusters = cms.InputTag("mtdMergedClusters", "FTLBarrel"),
    # clusters = cms.InputTag("mtdClusters", "FTLBarrel")
)

# Path
process.p = cms.Path(
    # process.mtdMergedClusters * 
    process.mtdSimMergedClusterValidation
)

print("Running MTD MergedCluster Validation...")
print(f"Input file: {options.inputFile}")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")