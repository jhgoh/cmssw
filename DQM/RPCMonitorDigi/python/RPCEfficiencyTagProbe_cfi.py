import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
rpcEfficiencyTagProbe = DQMEDAnalyzer('RPCEfficiencyTagProbe',
    vertex = cms.InputTag("offlinePrimaryVertices"),
    muons = cms.InputTag("muons"),
    minMuonPt = cms.double(30),
    maxMuonAbsEta = cms.double(2.4),
    maxMuonRelIso = cms.double(0.25), # Loose isolation ~98% eff. (tight=0.15)
    minTrackPt = cms.double(10),
    maxTrackAbsEta = cms.double(2.1),
    doCheckSign = cms.bool(True),
    minDR = cms.double(0.1),
    minMass = cms.double(70),
    maxMass = cms.double(110),
    triggerObjects = cms.InputTag("hltTriggerSummaryAOD"),
    triggerResults = cms.InputTag("TriggerResults::HLT"),
    triggerPaths = cms.vstring("HLT_IsoMu27", "HLT_IsoMu30", "HLT_IsoMu24", "HLT_Mu50", "HLT_Mu55"),
    triggerModules = cms.vstring(""), ## Make it to be a pair with the trigger path if given
    rpcRecHits = cms.InputTag("rpcRecHits"),
    dtSegments = cms.InputTag("dt4DSegments"),
    cscSegments = cms.InputTag("cscSegments"),
)

rpcefficiency = cms.Sequence(rpcEfficiencyTagProbe)


