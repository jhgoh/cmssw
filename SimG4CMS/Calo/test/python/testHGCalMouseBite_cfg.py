import FWCore.ParameterSet.Config as cms
import os, sys, imp, re
import FWCore.ParameterSet.VarParsing as VarParsing

####################################################################
### SETUP OPTIONS
options = VarParsing.VarParsing('standard')
options.register('type',
                 "V18",
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,
                  "type of operations: V16, V17, V17Shift")

### get and parse the command line arguments
options.parseArguments()
print(options)

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9

process = cms.Process("HGCalMouseBiteTest",Phase2C17I13M9)

####################################################################
# Use the options
if (options.type == "V18"):
    geomFile = "Geometry.HGCalCommonData.testHGCalV18Reco_cff"
elif (options.type == "V16"):
    geomFile = "Configuration.Geometry.GeometryExtended2026D88_cff"
else:
    geomFile = "Configuration.Geometry.GeometryExtended2026D92_cff"

print("Geometry file: ", geomFile)

process.load("SimGeneral.HepPDTESSource.pdt_cfi")
process.load(geomFile)
process.load('FWCore.MessageService.MessageLogger_cfi')

if hasattr(process,'MessageLogger'):
    process.MessageLogger.HGCalGeom=dict()
    process.MessageLogger.HGCSim=dict()

process.load("IOMC.RandomEngine.IOMC_cff")
process.RandomNumberGeneratorService.generator.initialSeed = 456789

process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("FlatRandomEGunProducer",
                                   PGunParameters = cms.PSet(
                                       PartID = cms.vint32(14),
                                       MinEta = cms.double(-3.5),
                                       MaxEta = cms.double(3.5),
                                       MinPhi = cms.double(-3.14159265359),
                                       MaxPhi = cms.double(3.14159265359),
                                       MinE   = cms.double(9.99),
                                       MaxE   = cms.double(10.01)
                                   ),
                                   AddAntiParticle = cms.bool(False),
                                   Verbosity       = cms.untracked.int32(0),
                                   firstRun        = cms.untracked.uint32(1)
                               )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)


process.load("Geometry.HGCalCommonData.hgcalMouseBiteTester_cfi")

 
process.p1 = cms.Path(process.generator*process.hgcalMouseBiteTester)
