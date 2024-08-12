is_PU = True

import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdValidation',Phase2C17I13M9)

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

# process.load("Configuration.Geometry.GeometryExtended2026D98Reco_cff")
process.load("Configuration.Geometry.GeometryExtended2026D110Reco_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
# process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '') # T21 previously
process.load('RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi')
process.load("Configuration.StandardSequences.Reconstruction_cff")


# max_events = 1 if is_PU else 1
max_events = 300 if is_PU else 2000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(max_events) )

# # Setup FWK for multithreaded
# process.options.numberOfThreads = 4
# process.options.numberOfStreams = 0
# process.options.numberOfConcurrentLuminosityBlocks = 0
# process.options.eventSetup.numberOfConcurrentIOVs = 1

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(10),
)

if is_PU:
    folder_name = "file:/eos/cms/store/relval/CMSSW_14_1_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/PU_140X_mcRun4_realistic_v3_STD_2026D110_PU-v1/2580000/"
else:
    folder_name = "file:/eos/cms/store/relval/CMSSW_14_1_0_pre4/RelValTTbar_14TeV/GEN-SIM-RECO/140X_mcRun4_realistic_v4_STD_2026D110_noPU-v1/2590000/"

# read file list from relval_200PU_files.txt
input_file = "relval_200PU_files.txt" if is_PU else "relval_NoPU_files.txt"
with open(input_file) as f:
    file_list = f.readlines()
file_list = [folder_name + f.strip() for f in file_list]

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        *file_list
    )
)

process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)

outname = "MTD_ntuples_200PU.root" if is_PU else "MTD_ntuples_NoPU.root"
process.TFileService = cms.Service("TFileService",
                                    fileName = cms.string(outname)
                                   )

# --- Global Validation
process.load("Validation.MtdValidation.vertices4DValid_cfi")

process.validation = cms.Sequence(process.vertices4DValid)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

process.p = cms.Path( process.mix + process.mtdTrackingRecHits + process.validation )
process.endjob_step = cms.EndPath(process.endOfProcess)
process.DQMoutput_step = cms.EndPath( process.DQMoutput )

process.schedule = cms.Schedule( process.p , process.endjob_step , process.DQMoutput_step )
