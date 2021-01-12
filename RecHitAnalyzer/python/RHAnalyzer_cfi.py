

import FWCore.ParameterSet.Config as cms 

fevt = cms.EDAnalyzer('RecHitAnalyzer'
    #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
    #, EBRecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEB')
    , reducedEBRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
    #, EERecHitCollection = cms.InputTag('ecalRecHit:EcalRecHitsEE')
    , reducedEERecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
    #, EBDigiCollection = cms.InputTag('simEcalDigis:ebDigis')
    #, selectedEBDigiCollection = cms.InputTag('selectDigi:selectedEcalEBDigiCollection')
    , reducedHBHERecHitCollection = cms.InputTag('reducedHcalRecHits:hbhereco')
    , genParticleCollection = cms.InputTag('genParticles')
    , gedPhotonCollection = cms.InputTag('gedPhotons')
    , ak4PFJetCollection = cms.InputTag('ak4PFJets')
    , genJetCollection = cms.InputTag('ak4GenJets')
    , trackRecHitCollection = cms.InputTag('generalTracks')
    , trackCollection = cms.InputTag("generalTracks")
   #, vertexCollection = cms.InputTag("offlinePrimaryVerticesWithBS")
    , vertexCollection = cms.InputTag("offlinePrimaryVertices")
    , siPixelRecHitCollection = cms.InputTag("siPixelRecHits")
    , siStripRecHitCollection = cms.VInputTag("siStripMatchedRecHits")
    , pfCollection        = cms.InputTag("particleFlow")
    , recoJetsForBTagging = cms.InputTag("ak4PFJetsCHS")
    , jetTagCollection    = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags")
    , tauCollection       = cms.InputTag("hpsPFTauProducer")
    , tauDiscriminator    = cms.InputTag("hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw")
    , ipTagInfoCollection = cms.InputTag("pfImpactParameterTagInfos")
    , mode                = cms.string("JetLevel")

    # Jet level cfg
    , nJets = cms.int32(-1)
    , minJetPt = cms.double(20.)
    , maxJetEta = cms.double(2.4)
    , z0PVCut  = cms.double(0.1)
    )
