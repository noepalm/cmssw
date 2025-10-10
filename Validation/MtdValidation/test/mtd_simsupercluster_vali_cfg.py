import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

# Define custom options
options.register('inputFile',
                 'file:/eos/home-n/npalmeri/MTD/MTD_photonReco/CMSSW_15_1_0_pre2_superclusterDev/src/SimFastTiming/MtdSimSuperClusterProducers/test/mtdSimSuperClusters_history_1k.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input ROOT file")

# outputFile, maxEvents registered by default

# Parse command line arguments
options.parseArguments()

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("SimSuperClusterValidationExample", Phase2C17I13M9)

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

# process.mtdSuperClusters = cms.EDProducer("MTDSuperClusterProducer",
#     srcBarrel = cms.InputTag("mtdClusters", "FTLBarrel"),
#     BarrelSuperClusterName = cms.string("FTLBarrel"),
#     timeThreshold = cms.double(10.0),
#     energyThreshold = cms.double(1.0)
# )

# Validation analyzer
process.mtdSimSuperClusterValidation = cms.EDAnalyzer("SimSuperClusterValidationExample",
    simSuperClusters = cms.InputTag("mtdSimSuperClusterProducer"), #tag name? Not sure where it's specified
    simLayerClusters = cms.InputTag("mixData:MergedMtdTruthLC"),
    # superClusters = cms.InputTag("mtdSuperClusters", "FTLBarrel"),
    # clusters = cms.InputTag("mtdClusters", "FTLBarrel")
)

# Path
process.p = cms.Path(
    # process.mtdSuperClusters * 
    process.mtdSimSuperClusterValidation
)

print("Running MTD SuperCluster Validation...")
print(f"Input file: {options.inputFile}")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")