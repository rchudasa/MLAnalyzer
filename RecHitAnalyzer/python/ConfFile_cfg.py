
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.register('skipEvents', 
    default=0, 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.int,
    info = "skipEvents")
options.register('processMode', 
    default='JetLevel', 
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process task: JetLevel or EventLevel")
options.register('processTask',
    default='tau_classification',
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.string,
    info = "process task: tau_classification")
options.register('processIsDebug',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "IsDebug: True or False")
options.register('processIsData',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "IsData: True or False")
options.register('processIsSignal',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "IsSignal: True or False")
options.register('processIsW',
    default=False,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "Is W plus jet: True or False")
options.register('processIsTrain',
    default=True,
    mult=VarParsing.VarParsing.multiplicity.singleton,
    mytype=VarParsing.VarParsing.varType.bool,
    info = "Is train sample: True or False")
options.parseArguments()

process = cms.Process("FEVTAnalyzer")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('106X_upgrade2018_realistic_v11_L1v1')
process.es_prefer_GlobalTag = cms.ESPrefer('PoolDBESSource','GlobalTag')

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents) 
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      options.inputFiles
      )
    , skipEvents = cms.untracked.uint32(options.skipEvents)
    )
print " >> Loaded",len(options.inputFiles),"input files from list."

process.load("MLAnalyzer.RecHitAnalyzer.RHAnalyzer_cfi")
process.fevt.isDebug = cms.bool(options.processIsDebug)
print " >> Debug mode:",(process.fevt.isDebug)
process.fevt.isData = cms.bool(options.processIsData)
print " >> Is data:",(process.fevt.isData)
process.fevt.mode = cms.string(options.processMode)
print " >> Processing as:",(process.fevt.mode)
process.fevt.task = cms.string(options.processTask)
print " >> Task:",(process.fevt.task)
process.fevt.isSignal = cms.bool(options.processIsSignal)
print " >> Is Signal:",(process.fevt.isSignal)
process.fevt.isW = cms.bool(options.processIsW)
print " >> Is W+jet:",(process.fevt.isW)
process.fevt.isTrain = cms.bool(options.processIsTrain)
print " >> Is Train:",(process.fevt.isTrain)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
    )

if (process.fevt.isData): hltpath = cms.vstring('HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*','HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*','HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*')
else : hltpath = cms.vstring('*')

process.hltFilter = cms.EDFilter("HLTHighLevel",
                                          eventSetupPathsKey = cms.string(''),
                                          TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                          HLTPaths = hltpath,
                                          andOr = cms.bool(True),
                                          throw = cms.bool(False)
                                          )

#process.SimpleMemoryCheck = cms.Service( "SimpleMemoryCheck", ignoreTotal = cms.untracked.int32(1) )
process.p = cms.Path(
  process.hltFilter*
#  process.patDefaultSequence*
  process.fevt
)
