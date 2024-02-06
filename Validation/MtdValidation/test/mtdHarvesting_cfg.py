import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('mtdHarvesting',Phase2C17I13M9)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.EDMtoMEAtRunEnd_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.load("Configuration.Geometry.GeometryExtended2026D98Reco_cff")

# [DEBUG] MY FLAG
is_minbias = False
# if passed argument == "minbias", set is_minbias = True
import sys
if len(sys.argv) > 1:
    if sys.argv[1] == "-minbias":
        is_minbias = True
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(-1),
)

infile = 'file:step3_inDQM_TTbar.root'
if is_minbias:
    infile = 'file:step3_inDQM_MinBias.root'

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring('file:' + infile)
)

# Path and EndPath definitions

process.edmtome_step = cms.Path(process.EDMtoME)
process.dqmsave_step = cms.Path(process.DQMSaver)

# --- PostProcessing

process.load("Validation.MtdValidation.btlSimHitsPostProcessor_cfi")
process.load("Validation.MtdValidation.btlLocalRecoPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdTracksPostProcessor_cfi")
process.load("Validation.MtdValidation.MtdEleIsoPostProcessor_cfi")
process.load("Validation.MtdValidation.Primary4DVertexPostProcessor_cfi")

process.harvesting = cms.Sequence(process.btlSimHitsPostProcessor + process.btlLocalRecoPostProcessor + process.MtdTracksPostProcessor + process.MtdEleIsoPostProcessor + process.Primary4DVertexPostProcessor)

process.p = cms.Path( process.harvesting )

process.schedule = cms.Schedule( process.edmtome_step , process.p , process.dqmsave_step )
