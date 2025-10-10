import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process("TEST", Phase2C17I13M9)

# essential things
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtendedRun4D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

process.MessageLogger.cerr.threshold = 'INFO'
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "file:/eos/user/p/pakrap/MTD/CMSSW_15_0_0_pre2/src/Validation/MtdValidation/29706.0_SinglePiFlatPt0p7To10+Run4D110/1000evt/step3.root"
    )
)

process.mtdSuperClusters = cms.EDProducer("MTDSuperClusterProducer",
    btlBarrel = cms.InputTag("mtdClusters", "FTLBarrel"),
    btlSuperClusterInstance = cms.string("FTLBarrel"),
    timeThreshold = cms.double(10.0),
    energyThreshold = cms.double(1.0)
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('supercluster_reco_singlepi.root'),
    outputCommands = cms.untracked.vstring(
        #'keep *'
        'drop *',
        #'keep *_mtdRecHits_*_*',            # keep rec hits
        #'keep *_mtdClusters_*_*',           # keep original clusters
        #'keep *_mtdSuperClusters_*_*',      # keep SuperClusters
        #'keep *_genParticles_*_*',          # keep gen info if needed?
        #'keep *_simHits_*_*',               # keep sim hits if needed for validation?

        'keep *_mix_FTLBarrel_*',
        'keep *_mix_FTLEndcap_*',
        'keep *_mtdRecoClusterToSimLayerClusterAssociation_*_*',
        'keep *_mtdSimLayerClusterToRecoClusterAssociation_*_*',
        'keep *_mtdSimLayerClusterToTPAssociation_*_*',
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
        'keep *_mtdSuperClusters_*_*',  
    ),
    #SelectEvents = cms.untracked.PSet(
    #    SelectEvents = cms.vstring('p')
    #)
)

process.mtd_reco = cms.Path(process.mtdSuperClusters)
process.outpath = cms.EndPath(process.output)
    
print("Testing BTL MTDSuperClusterProducer with adjacent cluster algorithm...")
