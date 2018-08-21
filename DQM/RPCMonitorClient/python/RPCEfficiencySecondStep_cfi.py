import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

rpcEfficiencySecond = DQMEDHarvester("RPCEfficiencySecond",
    SaveFile = cms.untracked.bool(False),
    NameFile = cms.untracked.string('/tmp/carrillo/RPCEfficiency.root'),
    debug = cms.untracked.bool(False),
)

rpcEfficiencyClientTagProbe = DQMEDHarvester("RPCEfficiencyClientTagProbe",
    minNExtrapolation = cms.untracked.uint32(100),
)

#rpcefficiencysecond = cms.Sequence(rpcEfficiencySecond)
rpcefficiencysecond = cms.Sequence(rpcEfficiencyClientTagProbe)


