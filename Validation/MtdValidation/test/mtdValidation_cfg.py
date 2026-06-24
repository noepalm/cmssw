import FWCore.ParameterSet.Config as cms


from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9
process = cms.Process('mtdValidation',Phase2C22I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi') # No pileup

process.load("Configuration.Geometry.GeometryExtendedRun4D121Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# For playback pileup mode
# process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
# process.load('Configuration.StandardSequences.Services_cff')
# Other statements
# process.mix.input.nbPileupEvents.averageNumber = cms.double(200.000000)
# process.mix.bunchspace = cms.int32(25)
# process.mix.minBunch = cms.int32(-3)
# process.mix.maxBunch = cms.int32(3)
# process.mix.input.fileNames = cms.untracked.vstring([]) # MinBias, from step3 confif file

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#Setup FWK for multithreaded
process.options.numberOfThreads = 96
process.options.numberOfStreams = 0
process.options.numberOfConcurrentLuminosityBlocks = 0
process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10),
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:/eos/home-n/npalmeri/MTD/MTD_mergedCluster_validation/samples/CMSSW_20_0_0_pre1_mergedClusters/ttbar_200PU/step3.root'
        'file:/eos/home-n/npalmeri/MTD/MTD_mergedCluster_validation/samples/CMSSW_20_0_0_pre1_mergedClusters/singlepi_0PU/step3_pions.root'
        # 'root://cms-xrd-global.cern.ch///store/relval/CMSSW_20_0_0_pre1/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_150X_mcRun4_realistic_v1_STD_D121_RegeneratedGS_PU-v1/2590000/795569f1-d53a-4649-8998-88766b5fdcea.root',
        # 'root://cms-xrd-global.cern.ch///store/relval/CMSSW_20_0_0_pre1/RelValSinglePiFlatPt0p7To10/GEN-SIM-DIGI-RAW/150X_mcRun4_realistic_v1_STD_RegeneratedGS_D121_noPU-v1/2590000/877b6071-9538-4fcd-a966-5850d2ccc754.root',
    )
)

# For playback pileup mode
# process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
# process.mix.playback = True

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

# --- BTL Validation
process.load("Validation.MtdValidation.btlSimHitsValid_cfi")
process.load("Validation.MtdValidation.btlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.btlLocalRecoValid_cfi")
btlValidation = cms.Sequence(process.btlSimHitsValid + process.btlDigiHitsValid + process.btlLocalRecoValid)

# --- ETL Validation
process.load("Validation.MtdValidation.etlSimHitsValid_cfi")
process.load("Validation.MtdValidation.etlDigiHitsValid_cfi")
process.load("Validation.MtdValidation.etlLocalRecoValid_cfi")
etlValidation = cms.Sequence(process.etlSimHitsValid + process.etlDigiHitsValid + process.etlLocalRecoValid)

# --- Reconstruction
process.load("RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi")

# --- Associators for Merged Clusters
process.load("SimFastTiming.MtdSimMergedClusterProducers.mtdSimMergedClusterProducer_cfi")
process.load("SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociation_cfi")
process.load("SimFastTiming.MtdAssociatorProducers.mtdSimMergedClusterToTPAssociation_cfi")
process.load('SimFastTiming.MtdAssociatorProducers.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits_cfi')
process.load('SimFastTiming.MtdAssociatorProducers.mtdSimMergedClusterToTPAssociatorByTrackId_cfi')

process.load("SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociation_cfi")
process.load("SimFastTiming.MtdAssociatorProducers.mtdSimLayerClusterToTPAssociatorByTrackId_cfi")
process.load("SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociation_cfi")
process.load("SimFastTiming.MtdAssociatorProducers.mtdRecoClusterToSimLayerClusterAssociatorByHits_cfi")

# --- Global Validation
process.load("Validation.MtdValidation.mtdTracksValid_cfi")
process.load("Validation.MtdValidation.vertices4DValid_cff")

# process.btlSimHitsValid.optionalPlots = True
# process.btlDigiHitsValid.optionalPlots = True
# process.etlDigiHitsValid.optionalPlots = True
# process.btlLocalRecoValid.optionalPlots = True
# process.etlLocalRecoValid.optionalPlots = True
# process.mtdTracksValid.optionalPlots = True
# process.vertices4DValid.optionalPlots = True

process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string("tree.root")
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

process.validation = cms.Sequence(btlValidation + etlValidation + process.vertices4DValid + process.mtdTracksValid + process.mtdMergedClusterValidation)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)


process.p = cms.Path( process.mix 
                      + process.mtdSimLayerClusterToTPAssociatorByTrackId 
                      + process.mtdSimLayerClusterToTPAssociation 
                      + process.mtdRecoClusterToSimLayerClusterAssociatorByHits
                      + process.mtdRecoClusterToSimLayerClusterAssociation
                      + process.mtdTrackingRecHits 
                      + process.mtdSimMergedClusterProducer 
                      + process.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits 
                      + process.mtdSimMergedClusterToTPAssociatorByTrackId 
                      + process.mtdSimMergedClusterToTPAssociation 
                      + process.mtdRecoMergedClusterToSimMergedClusterAssociation 
                      + process.validation )
# process.p = cms.Path( process.mix 
#                     + process.mtdClusters
#                     + process.mtdMergedClusters
#                     + process.mtdTrackingRecHits 
#                     + process.mtdSimMergedClusterProducer 
#                     + process.mtdRecoMergedClusterToSimMergedClusterAssociatorByHits 
#                     + process.mtdSimMergedClusterToTPAssociatorByTrackId 
#                     + process.mtdSimMergedClusterToTPAssociation 
#                     + process.mtdRecoMergedClusterToSimMergedClusterAssociation)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )
