import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("MtdMergedAssociationMapTest", Phase2C17I13M9)

# Set up command line argument parsing
options = VarParsing('analysis')

# Define custom options
options.register('useSimTopologicalClustering',
                 True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Enable topological clustering in SIM MergedCluster producer")

options.register('inputFile',
                 'file:/eos/home-n/npalmeri/MTD/MTD_supercluster/association_maps_validation/src/RecoLocalFastTime/FTLClusterizer/test/output.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Input ROOT file")

# Parse command line arguments
options.parseArguments()

# essential things
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

# process.MessageLogger.cerr.threshold = 'INFO'
# process.MessageLogger.cerr.threshold = 'DEBUG'
# process.MessageLogger.debugModules = ["*"]

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFile),
)


# Looking for type: reco::MtdRecoMergedClusterToSimMergedClusterAssociator
# Looking for module label: MtdRecoMergedClusterToSimMergedClusterAssociatorByHits

# first, create associator instance through producer
# MtdRecoMergedClusterToSimMergedClusterAssociatorByHitsProducer
from SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits_cfi import mtdRecoMergedClusterToSimMergedClusterAssociatorByHits
process.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits = mtdRecoMergedClusterToSimMergedClusterAssociatorByHits.clone()

# Load reco merged cluster => sim merged cluster producer
from SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociation_cfi import mtdRecoMergedClusterToSimMergedClusterAssociation

process.mtdRecoMergedClusterToSimMergedClusterAssociation = mtdRecoMergedClusterToSimMergedClusterAssociation.clone()


# And now we do the same for SimMergedCluster <-> TP producer
from SimFastTiming.MtdAssociatorProducers.mtdSimMergedClusterToTPAssociatorByTrackId_cfi import mtdSimMergedClusterToTPAssociatorByTrackId
process.mtdSimMergedClusterToTPAssociatorByTrackId = mtdSimMergedClusterToTPAssociatorByTrackId.clone()

from SimFastTiming.MtdAssociatorProducers.mtdSimMergedClusterToTPAssociation_cfi import mtdSimMergedClusterToTPAssociation
process.mtdSimMergedClusterToTPAssociation = mtdSimMergedClusterToTPAssociation.clone()

# # Load MTD truth map associators (needed for the producer)
# from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi import mtdSimLayerClusterToTPAssociatorByTrackId
# from SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi import mtdSimLayerClusterToTPAssociation

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    outputCommands = cms.untracked.vstring(
        #'keep *'
        'drop *',
        #'keep *_mtdRecHits_*_*',            # keep rec hits
        #'keep *_mtdClusters_*_*',           # keep original clusters
        #'keep *_mtdMergedClusters_*_*',      # keep MergedClusters
        #'keep *_genParticles_*_*',          # keep gen info if needed?
        #'keep *_simHits_*_*',               # keep sim hits if needed for validation?

        "keep *_genParticles_*_*",  # keep GenParticles
        "keep *_mtdSimLayerClusterToTPAssociation_*_*",
        "keep *_mix_*_*",  # Keep mix to have TrackingParticles, MtdSimLayerClusters
        # 'keep *_mix_FTLBarrel_*',
        # 'keep *_mix_FTLEndcap_*',
        'keep *_mtdRecoClusterToSimLayerClusterAssociation_*_*',
        'keep *_mtdSimLayerClusterToRecoClusterAssociation_*_*',
        'keep *_mtdSimLayerClusterToTPAssociation_*_*',
        "keep *_mtdSimLayerClusterToTPAssociatorByTrackId_*_*",
        'keep *_mtdRecoMergedClusterToSimMergedClusterAssociation_*_*',
        'keep *_mtdSimMergedClusterToRecoMergedClusterAssociation_*_*',
        'keep *_mtdSimMergedClusterToTPAssociation_*_*',
        "keep *_mtdSimMergedClusterToTPAssociatorByTrackId_*_*",
        'keep *_mtdRecHits_FTLBarrel_*',
        'keep *_mtdRecHits_FTLEndcap_*',
        'keep *_mtdUncalibratedRecHits_FTLBarrel_*',
        'keep *_mtdUncalibratedRecHits_FTLEndcap_*',
        'keep *_trackExtenderWithMTD_generalTrackmtdpos_*',
        'keep *_trackExtenderWithMTD_generalTracksigmatmtd_*',
        'keep *_trackExtenderWithMTD_generalTracktmtd_*',
        'keep *_mtdTrackQualityMVA_mtdQualMVA_*',
        'keep *_mtdClusters_FTLBarrel_*',
        'keep *_mtdClusters_FTLEndcap_*',
        'keep *_mtdTrackingRecHits_*_*',
        'keep *_mtdMergedClusters_*_*',  
        "keep *_mtdSimMergedClusterProducer_*_*",
    ),
    #SelectEvents = cms.untracked.PSet(
    #    SelectEvents = cms.vstring('p')
    #)
)

process.mergedClusterSequence = cms.Sequence(
    process.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits + 
    process.mtdRecoMergedClusterToSimMergedClusterAssociation +
    process.mtdSimMergedClusterToTPAssociatorByTrackId +
    process.mtdSimMergedClusterToTPAssociation
)

process.p = cms.Path(process.mergedClusterSequence)
process.out_step = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p, process.out_step)

print("Testing BTL MTDMergedClusterProducer with adjacent cluster algorithm...")
print(f"Output file: {options.outputFile}")
print(f"Max events: {options.maxEvents}")