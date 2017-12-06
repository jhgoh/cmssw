
import FWCore.ParameterSet.Config as cms

process = cms.Process('DIGI')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic50ns13TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load('RecoLocalMuon.RPCRecHit.rpcRecHits_cfi')
from RecoLocalMuon.RPCRecHit.rpcRecHits_cfi import *

process.load('L1Trigger.L1TMuonCPPF.simCppfDigis_cfi')
from L1Trigger.L1TMuonCPPF.simCppfDigis_cfi import *

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )


# Input source
readFiles = cms.untracked.vstring()
process.source = cms.Source(
	'PoolSource',
	fileNames = readFiles
	)

AWB_dir = '/afs/cern.ch/work/a/abrinke1/public/EMTF/MC/SingleMu_Pt1To1000_FlatRandomOneOverPt/'
readFiles.extend([
		'file:'+AWB_dir+'26CA310A-4164-E611-BE48-001E67248566.root'
])

process.options = cms.untracked.PSet(
)

process.treeOut = cms.OutputModule("PoolOutputModule",
				   fileName = cms.untracked.string("CPPF_MC_Tree_SingleMu_10k.root"),
				   outputCommands = cms.untracked.vstring(
		'drop *',
		'keep *_simMuonRPCDigis_*_*', 
		'keep *_simCppfDigis_*_*',
		)
)

# Additional output definition

# Other statements

# process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.rpcRecHits.rpcDigiLabel = 'simMuonRPCDigis'
process.simCppfDigis.recHitLabel = 'rpcRecHits'

# Path and EndPath definitions
process.rpcrechits_step = cms.Path(process.rpcRecHits)
process.simCppfDigis_step = cms.Path(process.simCppfDigis)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.treeOut_step = cms.EndPath(process.treeOut)


# Schedule definition
process.schedule = cms.Schedule(
    process.rpcrechits_step,
    process.simCppfDigis_step,
    process.endjob_step,
    process.treeOut_step )

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Customisation from command line
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
