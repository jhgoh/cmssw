# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step4 --filetype DQM --mc --conditions auto:run2_mc --era Run2_2016 -s HARVESTING:@standardValidation+@standardDQM+@ExtraHLT+@miniAODValidation+@miniAODDQM+@L1TMuon -n 100 --filein file:step3_inDQM.root --fileout file:step4.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HARVESTING',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DQMSaverAtRunEnd_cff')
process.load('Configuration.StandardSequences.Harvesting_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring(
        'file:DQM/step3_000.root',
        'file:DQM/step3_001.root',
        'file:DQM/step3_002.root',
        'file:DQM/step3_003.root',
        'file:DQM/step3_004.root',
        'file:DQM/step3_005.root',
        'file:DQM/step3_006.root',
        'file:DQM/step3_007.root',
        'file:DQM/step3_008.root',
        'file:DQM/step3_009.root',
        'file:DQM/step3_010.root',
        'file:DQM/step3_011.root',
        'file:DQM/step3_012.root',
        'file:DQM/step3_013.root',
        'file:DQM/step3_014.root',
        'file:DQM/step3_015.root',
        'file:DQM/step3_016.root',
        'file:DQM/step3_017.root',
        'file:DQM/step3_018.root',
        'file:DQM/step3_019.root',
    ),
)

process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound'),
    fileMode = cms.untracked.string('FULLMERGE')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step4 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition
process.dqmSaver.saveAtJobEnd = True

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# Path and EndPath definitions
process.dqmHarvesting = cms.Path(process.rpcTier0Client)
process.dqmsave_step = cms.Path(process.DQMSaver)

# Schedule definition
process.schedule = cms.Schedule(process.dqmHarvesting,process.dqmsave_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
