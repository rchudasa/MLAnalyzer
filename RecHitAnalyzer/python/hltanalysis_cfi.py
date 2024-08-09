import FWCore.ParameterSet.Config as cms

from MLAnalyzer.RecHitAnalyzer.dummybranches_cff import *

hltanalysis = cms.EDAnalyzer('HLTBitAnalyzer',
    HLTProcessName = cms.string('HLT'),
    hltresults = cms.InputTag('TriggerResults::HLT'),
    l1results = cms.InputTag('gtStage2Digis'),
    UseTFileService = cms.untracked.bool(True),
    RunParameters = cms.PSet(
        isData = cms.untracked.bool(True)),
    #dummyBranches = dummy_branches_for_PbPb_2018_HLT,
    #l1dummyBranches = dummy_branches_for_PbPb_2018_L1,

    mctruth = cms.InputTag(''),
    genEventInfo = cms.InputTag(''),
    OfflinePrimaryVertices0 = cms.InputTag(''),
    simhits = cms.InputTag(''),
)
