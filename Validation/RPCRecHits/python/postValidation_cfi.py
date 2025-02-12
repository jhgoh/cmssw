import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

def efficSet(nameIn, titleIn, numeratorIn, denominatorIn, typeIn="eff"):
    pset = cms.PSet(name=cms.untracked.string(nameIn),
                    title=cms.untracked.string(titleIn),
                    numerator=cms.untracked.string(numeratorIn),
                    denominator=cms.untracked.string(denominatorIn),
                    type=cms.untracked.string(typeIn))
    return pset

rpcRecHitSimRecoClient = DQMEDHarvester("RPCRecHitValidClient",
    subDir = cms.string("RPC/RPCRecHitV"),
)

rpcRecHitPostValidation = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring("RPC/RPCRecHitV",),
    #subDirs = cms.untracked.vstring("RPC/RPCRecHitV/SimVsReco",
    #                                "RPC/RPCRecHitV/SimVsDTExt",
    #                                "RPC/RPCRecHitV/SimVsCSCExt"),
    efficiency = cms.vstring(),
    resolution = cms.vstring(),
    efficiencyProfileSets = cms.untracked.VPSet(
        efficSet("Efficiency/Effic_wheel", "Barrel SimHit to RecHit matching efficiency;Wheel",
                 "Occupancy/MatchBarrelOccupancy_wheel", "Occupancy/RefHitBarrelOccupancy_wheel"),
        efficSet("Efficiency/Effic_station", "Barrel SimHit to RecHit matching efficiency;Station",
                 "Occupancy/MatchBarrelOccupancy_station", "Occupancy/RefHitBarrelOccupancy_station"),
        efficSet("Efficiency/Effic_disk", "Endcap SimHit to RecHit matching efficiency;Disk",
                 "Occupancy/MatchEndcapOccupancy_disk", "Occupancy/RefHitEndcapOccupancy_disk"),
    ),
    resolutionSets = cms.untracked.VPSet(
    ),
    outputFileName = cms.untracked.string("")
)

rpcPointVsRecHitPostValidation = DQMEDHarvester("DQMGenericClient",
    subDirs = cms.untracked.vstring("RPC/RPCRecHitV/DTVsReco",
                                    "RPC/RPCRecHitV/CSCVsReco"),
#                                    "RPC/RPCRecHitV/TrackVsReco"),
    efficiency = cms.vstring(),
    resolution = cms.vstring(),
    efficiencyProfileSets = cms.untracked.VPSet(
        efficSet("Efficiency/Effic_wheel", "Barrel RPCPoint to RecHit matching efficiency;Wheel",
                 "Occupancy/MatchBarrelOccupancy_wheel", "Occupancy/RefHitBarrelOccupancy_wheel"),
        efficSet("Efficiency/Effic_station", "Barrel RPCPoint to RecHit matching efficiency;Station",
                 "Occupancy/MatchBarrelOccupancy_station", "Occupancy/RefHitBarrelOccupancy_station"),
        efficSet("Efficiency/Effic_disk", "Endcap RPCPoint to RecHit matching efficiency;Disk",
                 "Occupancy/MatchEndcapOccupancy_disk", "Occupancy/RefHitEndcapOccupancy_disk"),
    ),
    resolutionSets = cms.untracked.VPSet(
    ),
    outputFileName = cms.untracked.string("")
)

rpcRecHitPostValidation_step = cms.Sequence(rpcRecHitPostValidation+rpcRecHitSimRecoClient)
rpcPointVsRecHitPostValidation_step = cms.Sequence(rpcPointVsRecHitPostValidation)
