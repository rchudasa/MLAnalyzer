# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: Configuration/GenProduction/python/Guns/DiPi0ToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8_cfi.py --conditions 102X_upgrade2018_realistic_v15 -n 10 --mc --era Run2_2018 --eventcontent AODSIM --runUnscheduled --step GEN,SIM,DIGI,L1,DIGI2RAW,HLT:@fake,RAW2DIGI,L1Reco,RECO --geometry DB:Extended --datatier GEN-SIM-DIGI-RAW-RECO --beamspot Realistic25ns13TeVEarly2018Collision --fileout file:step_GEN2AODSIM.root --python_filename DiPi0ToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8_pythia8_noPU_v2_cfg.py --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('HLT',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic25ns13TeVEarly2018Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_Fake_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.genHToTauTauFilter = cms.EDFilter("GenHToTauTauFilter",
#    src       = cms.InputTag("genParticles"), #GenParticles collection as input
#    nHiggs    = cms.double(1),    #Number of pdgID=25 candidates 
#    tauPtCut  = cms.double(10.0), #at least a GenTau with this minimum pT 
#    tauEtaCut = cms.double(2.4),  #GenTau eta
#    taudRCut  = cms.double(0.4)   #GenTauTau cut
#)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/Guns/DiPi0ToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8_cfi.py nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    outputCommands = process.AODSIMEventContent.outputCommands+cms.untracked.vstring('keep *_simSiPixelDigis_*_*','keep *_simSiStripDigis_*_*','keep *_siPixelClusters_*_*','keep *_siStripClusters_*_*','keep *_g4SimHits_*_*','keep *_generalTracks_*_*','keep *_siStripMatchedRecHits_*_*','keep *_siPixelRecHits_*_*'),
    #outputCommands = process.AODSIMEventContent.outputCommands
    fileName = cms.untracked.string('file:aToTauTau_m3p6To15_pT30To200_ctau0To3_eta0To1p4_pythia8_noPU.root')
)

# Additional output definition

# Other statements
process.XMLFromDBSource.label = cms.string("Extended")
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v15_L1v1', '')

process.generator = cms.EDFilter("Pythia8PtGunV2",
    PGunParameters = cms.PSet(
        AddAntiParticle = cms.bool(True),
        MaxCTau = cms.double(3.0),
        MaxEta = cms.double(1.4),
        MaxMass = cms.double(15.0),
        MaxPhi = cms.double(3.14159265359),
        MaxPt = cms.double(200.0),
        MinCTau = cms.double(0.0),
        MinEta = cms.double(-1.4),
        MinMass = cms.double(3.6),
        MinPhi = cms.double(-3.14159265359),
        MinPt = cms.double(30.0),
        ParticleID = cms.vint32(25)
    ),
    PythiaParameters = cms.PSet(
        parameterSets = cms.vstring('processParameters'),
        processParameters = cms.vstring(
            '25:onMode = off', 
            '25:onIfAny = 15'
        )
    ),
    Verbosity = cms.untracked.int32(0),
    firstRun = cms.untracked.uint32(1),
    maxEventsToPrint = cms.untracked.int32(1),
    psethack = cms.string('Higgs pTgun'),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    pythiaPylistVerbosity = cms.untracked.int32(1)
)


# Path and EndPath definitions
#process.generation_step = cms.Path(process.pgen+process.genHToTauTauFilter)
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.AODSIMoutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

process.source.numberEventsInLuminosityBlock=cms.untracked.uint32(500)
#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
