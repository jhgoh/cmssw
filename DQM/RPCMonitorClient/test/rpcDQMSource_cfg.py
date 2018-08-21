# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --runUnscheduled --conditions 101X_dataRun2_Prompt_v11 --geometry DB:Extended -s RAW2DIGI,DQM --datatier DQMIO -n -1 --era Run2_2017 --eventcontent DQM --filein file:/xrootd/store/user/jhgoh/data/Run2018C_SingleMuon_RAW-RECO_ZMu-PromptReco-v3/FE016AAA-F88F-E811-8F6E-FA163EFA55B0.root --fileout step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('DQM',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/xrootd/store/user/jhgoh/data/Run2018C_SingleMuon_RAW-RECO_ZMu-PromptReco-v3/FE016AAA-F88F-E811-8F6E-FA163EFA55B0.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step3.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v11', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
#process.dqmoffline_step = cms.EndPath(process.DQMOffline)
process.dqmoffline_step = cms.EndPath(process.rpcTier0Source)
#process.dqmofflineOnPAT_step = cms.EndPath(process.DQMOffline)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.dqmoffline_step,process.dqmofflineOnPAT_step,process.DQMoutput_step)
process.schedule = cms.Schedule(process.raw2digi_step,process.dqmoffline_step,process.DQMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
