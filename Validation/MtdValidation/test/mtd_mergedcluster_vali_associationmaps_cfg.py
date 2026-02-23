import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

options.register('inputFile',
                #  'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/CMSSW_15_1_0_pre2/src/SimFastTiming/MtdAssociatorProducers/test/mtdMergedAssociationMaps_forDev031225.root',
                 'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/association_maps_validation/src/SimFastTiming/MtdAssociatorProducers/test/output.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input ROOT file")

# Parse command line arguments
options.parseArguments()

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("vali", Phase2C17I13M9)

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
process.mtdMergedClusterValidation = cms.EDAnalyzer("MergedClusterValidation_withAssociationMaps",
    mergedClusters = cms.InputTag("mtdMergedClusters", "FTLBarrel"),
    clusters = cms.InputTag("mtdClusters", "FTLBarrel"),
    simMergedClusters = cms.InputTag("mtdSimMergedClusterProducer"),
    simLayerClusters = cms.InputTag("mix", "MergedMtdTruthLC"),
    trackingParticles = cms.InputTag("mix", "MergedTrackTruth"),
    mergedRecoToSimMap = cms.InputTag("mtdRecoMergedClusterToSimMergedClusterAssociation", ""),
    mergedSimToRecoMap = cms.InputTag("mtdRecoMergedClusterToSimMergedClusterAssociation", ""),
    mergedSimToTPMap = cms.InputTag("mtdSimMergedClusterToTPAssociation", ""),
    mergedTPToSimMap = cms.InputTag("mtdSimMergedClusterToTPAssociation", ""),
)

# Path
process.p = cms.Path( 
    process.mtdMergedClusterValidation
)

print("Running MTD MergedCluster Validation...")
print(f"Input file: {options.inputFile}")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")