"""
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

options.register('inputFile',
                 "file:/eos/user/p/pakrap/MTD/CMSSW_15_0_0_pre2/src/Validation/MtdValidation/29706.0_SinglePiFlatPt0p7To10+Run4D110/1000evt/step3.root",
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

# Load the merged cluster producer (reco)
process.load('RecoLocalFastTime.FTLClusterizer.mtdMergedClusterProducer_cfi')

# Load sim layer cluster to TP associator algorithm
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi')

# Load sim layer cluster to TP association producer
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi')

# Load sim merged cluster producer
process.load('SimFastTiming.MtdSimMergedClusterProducers.mtdSimMergedClusterProducer_cfi')

# Load the reco-to-sim merged cluster associator algorithm (ByHits implementation)
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits_cfi')

# Load the reco-to-sim merged cluster association map producer
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociation_cfi')

# Validation analyzer
process.mtdMergedClusterValidation = cms.EDAnalyzer("MergedClusterValidation_withAssociationMaps",
    mergedClusters = cms.InputTag("mtdMergedClusters", "FTLBarrel"),
    clusters = cms.InputTag("mtdClusters", "FTLBarrel"),
    simMergedClusters = cms.InputTag("mtdSimMergedClusterProducer"),
    simLayerClusters = cms.InputTag("mix", "MergedMtdTruthLC"),
    mergedRecoToSimMap = cms.InputTag("mtdRecoMergedClusterToSimMergedClusterAssociation", ""),
    mergedSimToRecoMap = cms.InputTag("mtdRecoMergedClusterToSimMergedClusterAssociation", ""),
)

# Path - order matters!
process.p = cms.Path( 
    process.mtdMergedClusters *                                           # 1. Create reco merged clusters
    process.mtdSimLayerClusterToTPAssociatorByTrackId *                  # 2. Sim layer cluster to TP associator
    process.mtdSimLayerClusterToTPAssociation *                          # 3. Sim layer cluster to TP association
    process.mtdSimMergedClusterProducer *                                # 4. Create sim merged clusters
    process.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits *     # 5. Reco-sim merged cluster associator
    process.mtdRecoMergedClusterToSimMergedClusterAssociation *          # 6. Create association maps
    process.mtdMergedClusterValidation                                   # 7. Run validation
)

print("Running MTD MergedCluster Validation...")
print(f"Input file: {options.inputFile}")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}") 

"""

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# Set up command line argument parsing
options = VarParsing('analysis')

options.register('inputFile',
                 "file:/eos/user/p/pakrap/MTD/CMSSW_15_0_0_pre2/src/Validation/MtdValidation/29706.0_SinglePiFlatPt0p7To10+Run4D110/1000evt/step3.root",
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input ROOT file")

options.parseArguments()

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("MergedClusterValidation", Phase2C17I13M9)

# Load standard configurations
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# Load CPE and full reconstruction sequence
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

# Configure mtdTrackingRecHits to use merged clusters
process.mtdTrackingRecHits.barrelClusters = cms.InputTag("mtdMergedClusters", "FTLBarrel")
process.mtdTrackingRecHits.endcapClusters = cms.InputTag("mtdMergedClusters", "FTLEndcap")

process.MessageLogger.cerr.threshold = 'INFO'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFile)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("file:" + options.outputFile)
)

# Load the merged cluster producer (reco)
process.load('RecoLocalFastTime.FTLClusterizer.mtdMergedClusterProducer_cfi')

# Load sim layer cluster to TP associator algorithm
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi')

# Load sim layer cluster to TP association producer
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi')

# Load sim merged cluster producer
process.load('SimFastTiming.MtdSimMergedClusterProducers.mtdSimMergedClusterProducer_cfi')

# Load the reco-to-sim merged cluster associator algorithm (ByHits implementation)
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits_cfi')

# Load the reco-to-sim merged cluster association map producer
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociation_cfi')

# Load DQM services
process.load('DQMServices.Core.DQMStore_cfi')

# Load MTD Tracks validation and configure for barrel-only
process.load("Validation.MtdValidation.mtdTracksValid_cfi")
process.mtdTracksValid.trackMinimumEtlEta = cms.double(10.0)  # Effectively disable ETL validation
process.mtdTracksValid.trackMaximumEtlEta = cms.double(10.0)  # Effectively disable ETL validation

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

# Clear mix digitizers
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# DQM output 
process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Path - order matters!
process.p = cms.Path( 
    process.mix *
    process.mtdMergedClusters *                                           # 1. Create reco merged clusters
    process.mtdTrackingRecHits *                                          # 2. Create tracking rec hits from merged clusters
    process.trackExtenderWithMTD *                                        # 3. Extend tracks with MTD
    process.mtdSimLayerClusterToTPAssociatorByTrackId *                   # 4. Sim layer cluster to TP associator
    process.mtdSimLayerClusterToTPAssociation *                           # 5. Sim layer cluster to TP association
    process.mtdSimMergedClusterProducer *                                 # 6. Create sim merged clusters
    process.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits *      # 7. Reco-sim merged cluster associator
    process.mtdRecoMergedClusterToSimMergedClusterAssociation *           # 8. Create association maps
    process.mtdTracksValid *                                              # 9. MTD Tracks validation (barrel only)
    process.mtdMergedClusterValidation                                    # 10. Merged cluster validation
)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(process.p, process.endjob_step, process.DQMoutput_step)

print("Running MTD MergedCluster Validation with Barrel Tracking Validation...")
print(f"Input file: {options.inputFile}")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")