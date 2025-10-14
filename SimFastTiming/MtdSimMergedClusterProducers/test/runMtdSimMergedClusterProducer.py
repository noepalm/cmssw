import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

# Define custom options
options.register('useTopologicalClustering',
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Enable topological clustering in MergedCluster producer")

# Parse command line arguments
options.parseArguments()

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("MtdSuperClusOnly", Phase2C17I13M9)

# Load standard services and geometry
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")  # Needed for TrackingParticles
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

# Configure logging levels for your producer
process.MessageLogger.debugModules = ["*"]
process.MessageLogger.cerr.MtdSimMergedClusterProducer = cms.untracked.PSet(
    limit = cms.untracked.int32(-1),  # No limit on messages
    # Choose your debug level:
    # INFO: Shows LogInfo and above (basic info)
    # DEBUG: Shows LogDebug and above (more detailed)  
    # TRACE: Shows LogTrace and above (most detailed)
    reportEvery = cms.untracked.int32(1)
)
# Set overall threshold:
process.MessageLogger.cerr.threshold = cms.untracked.string('INFO')    # Basic info
# process.MessageLogger.cerr.threshold = cms.untracked.string('DEBUG')   # More detailed
# process.MessageLogger.cerr.threshold = cms.untracked.string('TRACE')   # Most detailed

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag  
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Input file (must contain SimClusters + TrackingParticles)
process.source = cms.Source("PoolSource",
    # fileNames = cms.untracked.vstring(options.inputFiles)
    fileNames = cms.untracked.vstring([f"file:/eos/home-n/npalmeri/ntuples/MTD/PhotonReco/crab_MTDPhotonReco/CRAB_UserFiles/SingleGammaFlatPt0p1To10_Run4D110_aging1000_noPU_MTDPhotonReco/250424_080407/0000/step2_{i}.root" for i in range(1, 11)]),
    skipEvents = cms.untracked.uint32(175)  # Skip event 0, start from event 1
)

# Number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Load your producer
from SimFastTiming.MtdSimMergedClusterProducers.mtdSimMergedClusterProducer_cfi import mtdSimMergedClusterProducer
process.mtdSimMergedClusterProducer = mtdSimMergedClusterProducer.clone(
    useTopologicalClustering = cms.bool(options.useTopologicalClustering)
)

# Output module (optional - to store results)
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_mtdSimMergedClusterProducer_*_*",
        "keep *_mix_*_*",  # Keep mix to have TrackingParticles, MtdSimLayerClusters
        "keep *_genParticles_*_*",  # keep GenParticles
        "keep *_mtdSimLayerClusterToTPAssociation_*_*",
        "keep *_mtdSimLayerClusterToTPAssociatorByTrackId_*_*",
    )
)

# Load MTD truth map associators (needed for the producer)
from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi import mtdSimLayerClusterToTPAssociatorByTrackId
from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi import mtdSimLayerClusterToTPAssociation

# clone and add label 
process.mtdSimLayerClusterToTPAssociatorByTrackId = mtdSimLayerClusterToTPAssociatorByTrackId.clone()
process.mtdSimLayerClusterToTPAssociation = mtdSimLayerClusterToTPAssociation.clone()

process.mergedClusterSequence = cms.Sequence(process.mtdSimLayerClusterToTPAssociatorByTrackId + process.mtdSimLayerClusterToTPAssociation + process.mtdSimMergedClusterProducer)

# process.mergedClusterSequence = cms.Sequence(process.mtdSimMergedClusterProducer)

# Execution path
process.p = cms.Path(process.mergedClusterSequence)
process.out_step = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p, process.out_step)

print("Running MTD MergedCluster Producer...")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")
print(f"Use topological clustering: {options.useTopologicalClustering}")
