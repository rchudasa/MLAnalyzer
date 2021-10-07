

import FWCore.ParameterSet.Config as cms 

from RecoMET.METProducers.METSignificanceParams_cfi import METSignificanceParams

fevt = cms.EDAnalyzer('RecHitAnalyzer'
    , transTrackBuilder              = cms.ESInputTag("TransientTrackBuilder")
    , reducedEBRecHitCollection      = cms.InputTag("reducedEcalRecHitsEB")
    , reducedEERecHitCollection      = cms.InputTag("reducedEcalRecHitsEE")
    , reducedHBHERecHitCollection    = cms.InputTag("reducedHcalRecHits:hbhereco")
    , JetCollection                  = cms.InputTag("slimmedJets")
    , genParticleCollection          = cms.InputTag("prunedGenParticles")
    , PhotonCollection               = cms.InputTag("slimmedPhotons")
    , PFJetCollection                = cms.InputTag("slimmedJets")
    , trackCollection                = cms.InputTag("generalTracks")
    , vertexCollection               = cms.InputTag("offlineSlimmedPrimaryVertices")
    , siPixelRecHitCollection        = cms.InputTag("siPixelRecHits")
    , siStripMatchedRecHitCollection = cms.InputTag("siStripMatchedRecHits", "matchedRecHit")
    , siStripRphiRecHits             = cms.InputTag("siStripMatchedRecHits", "rphiRecHit")
    , siStripStereoRecHits           = cms.InputTag("siStripMatchedRecHits", "stereoRecHit")
    #, metCollection                  = cms.InputTag("slimmedMETs")
    #, metPuppiCollection             = cms.InputTag("slimmedMETsPuppi")
    #, eleCollection                  = cms.InputTag("slimmedElectrons")
    , tauCollection                  = cms.InputTag("slimmedTaus")
    #, tauBoostedCollection           = cms.InputTag("slimmedTausBoosted")
    , triggerResultsTag              = cms.InputTag("TriggerResults", "", "HLT")
    , isDebug                        = cms.bool(False)
    , isData                         = cms.bool(False)
    , isSignal                       = cms.bool(True)
    , isW                            = cms.bool(False)
    , mode                           = cms.string("JetLevel")
    , task                           = cms.string("tau_classification")
    , rhoLabel                       = cms.InputTag("fixedGridRhoAll")

    # Jet level cfg
    , nJets     = cms.int32(-1)
    , minJetPt  = cms.double(20.)
    , maxJetEta = cms.double(2.4)
    , z0PVCut   = cms.double(0.1)

    # MET parameter
    , parameters = METSignificanceParams
 
    #granularity multiplier wrt ECAL maps for tracker and tracking RH images
    , granularityMultiPhi = cms.int32(5)
    , granularityMultiEta = cms.int32(5)
    )
