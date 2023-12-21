// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
//
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//

#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
{
  //debug      = iConfig.getParameter<bool>("isDebug");  
  mode_      = iConfig.getParameter<std::string>("mode");
  task_      = iConfig.getParameter<std::string>("task");
  isMC_      = iConfig.getParameter<bool>("isMC");
  isSignal_  = iConfig.getParameter<bool>("isSignal");
  isW_       = iConfig.getParameter<bool>("isW");
  isBoostedTop_   = iConfig.getParameter<bool>("isBoostedTop");

  //EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  //EBDigiCollectionT_      = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection"));
  //EBDigiCollectionT_      = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  //EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHitCollection"));
  HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));
  TRKRecHitCollectionT_   = consumes<TrackingRecHitCollection>(iConfig.getParameter<edm::InputTag>("trackRecHitCollection"));

  if(isMC_){
  	genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  }  

  photonCollectionT_      = consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("gedPhotonCollection"));
  jetCollectionT_         = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("ak4PFJetCollection"));
  pfjetsToken_            = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("srcPfJets"));
  genJetCollectionT_      = consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetCollection"));
  trackCollectionT_       = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  transientTrackBuilderT_ = iConfig.getParameter<edm::ESInputTag>("transTrackBuilder");
  pfCollectionT_          = consumes<PFCollection>(iConfig.getParameter<edm::InputTag>("pfCollection"));

  pfCandidatesToken_      = consumes<edm::View<reco::Candidate> >(iConfig.getParameter<edm::InputTag>("srcPFCandidates"));

  vertexCollectionT_      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  secVertexCollectionT_   = consumes<reco::VertexCompositePtrCandidateCollection>(iConfig.getParameter<edm::InputTag>("secVertexCollection"));

  recoJetsT_              = consumes<edm::View<reco::Jet> >(iConfig.getParameter<edm::InputTag>("recoJetsForBTagging"));
  jetTagCollectionT_      = consumes<reco::JetTagCollection>(iConfig.getParameter<edm::InputTag>("jetTagCollection"));
  ipTagInfoCollectionT_   = consumes<std::vector<reco::CandIPTagInfo> > (iConfig.getParameter<edm::InputTag>("ipTagInfoCollection"));

  siPixelRecHitCollectionT_ = consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("siPixelRecHitCollection"));

  //siStripRecHitCollectionT_ = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<std::vector<edm::InputTag> >("siStripRecHitCollection"));
  siStripMatchedRecHitCollectionT_ = consumes<SiStripMatchedRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripMatchedRecHitCollection"));
  siStripRPhiRecHitCollectionT_    = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripRphiRecHits"));
  siStripUnmatchedRPhiRecHitCollectionT_ = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripUnmatchedRphiRecHits"));
  siStripStereoRecHitCollectionT_  = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripStereoRecHits"));
  siStripUnmatchedStereoRecHitCollectionT_  = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripUnmatchedStereoRecHits"));
 
  metCollectionT_           = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("metCollection"));

  tauCollectionT_           = consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("tauCollection"));
  tauDecayMode_             = consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tauDecayMode"));
  tauMVAIsolation_          = consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tauMVAIsolationRaw"));
  tauMuonRejection_         = consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tauMuonRejectionLoose"));
  tauElectronRejectionMVA6_ = consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tauElectronRejectionMVA6VLoose"));

  eleCollectionT_           = consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("eleCollection"));

  processName_              = iConfig.getUntrackedParameter<std::string>("processName","HLT");
  //triggerResultsToken_      = consumes<edm::TriggerResults> (iConfig.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")));
  triggerResultsToken_      = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerResultsTag"));
  
  jetSFType_                = iConfig.getParameter<std::string>("srcJetSF");
  jetResPtType_             = iConfig.getParameter<std::string>("srcJetResPt");
  jetResPhiType_            = iConfig.getParameter<std::string>("srcJetResPhi");
  rhoLabel_                 = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"));

  std::vector<edm::InputTag> srcLeptonsTags        = iConfig.getParameter< std::vector<edm::InputTag> >("srcLeptons");
  for(std::vector<edm::InputTag>::const_iterator it=srcLeptonsTags.begin();it!=srcLeptonsTags.end();it++) {
    lepTokens_.push_back( consumes<edm::View<reco::Candidate> >( *it ) );
  }

  metSigAlgo_               = new metsig::METSignificance(iConfig);

  //johnda add configuration
  minJetPt_  = iConfig.getParameter<double>("minJetPt");
  maxJetEta_ = iConfig.getParameter<double>("maxJetEta");
  z0PVCut_   = iConfig.getParameter<double>("z0PVCut");

  std::cout << " >> Mode set to " << mode_ << std::endl;
  if ( mode_ == "JetLevel" ) {
    doJets_ = true;
    nJets_ = iConfig.getParameter<int>("nJets");
    std::cout << "\t>> nJets set to " << nJets_ << std::endl;
  } else if ( mode_ == "EventLevel" ) {
    doJets_ = false;
  } else {
    std::cout << " >> Assuming EventLevel Config. " << std::endl;
    doJets_ = false;
  }



  // Initialize file writer
  // NOTE: initializing dynamic-memory histograms outside of TFileService
  // will cause memory leaks
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_sel = fs->make<TH1F>("h_sel", "isSelected;isSelected;Events", 2, 0., 2.);

  ///////////adjustable granularity stuff

  granularityMultiPhi[0]  = iConfig.getParameter<int>("granularityMultiPhi");
  granularityMultiEta[0]  = iConfig.getParameter<int>("granularityMultiEta");

  granularityMultiPhi[1] = 3;
  granularityMultiEta[1] = 3;

  for (unsigned int proj=0; proj<Nadjproj; proj++)
  {

    int totalMultiEta = granularityMultiEta[proj] * granularityMultiECAL;

    for (int i=0; i<eta_nbins_HBHE; i++)
    {
      double step=(eta_bins_HBHE[i+1]-eta_bins_HBHE[i])/totalMultiEta;
      for (int j=0; j<totalMultiEta; j++)
      {
        adjEtaBins[proj].push_back(eta_bins_HBHE[i]+step*j);
      }
    }
    adjEtaBins[proj].push_back(eta_bins_HBHE[eta_nbins_HBHE]);

    totalEtaBins[proj] = totalMultiEta*(eta_nbins_HBHE);
    totalPhiBins[proj] = granularityMultiPhi[proj] * granularityMultiECAL*HBHE_IPHI_NUM;

  }

  //////////// TTree //////////

  // These will be use to create the actual images
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  if ( doJets_ ) {
    branchesEvtSel_jet( RHTree, fs );
  } else {
    branchesEvtSel( RHTree, fs );
  }
  branchesEB           ( RHTree, fs );
  branchesEE           ( RHTree, fs );
  branchesHBHE         ( RHTree, fs );
  branchesECALatHCAL   ( RHTree, fs );
  branchesECALstitched ( RHTree, fs );
  branchesHCALatEBEE   ( RHTree, fs );
  branchesTracksAtEBEE(RHTree, fs);
  branchesTracksAtECALstitched( RHTree, fs);
  branchesPFCandsAtEBEE(RHTree, fs);
  branchesPFCandsAtECALstitched( RHTree, fs);
  //branchesTRKlayersAtEBEE(RHTree, fs);
  //branchesTRKlayersAtECAL(RHTree, fs);
  //branchesTRKvolumeAtEBEE(RHTree, fs);
  //branchesTRKvolumeAtECAL(RHTree, fs);
  branchesJetInfoAtECALstitched( RHTree, fs);
  branchesScalarInfo( RHTree, fs);  
  branchesTRKlayersAtECALstitched(RHTree, fs);

  // For FC inputs
  //RHTree->Branch("FC_inputs",      &vFC_inputs_);

} // constructor
//
RecHitAnalyzer::~RecHitAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete metSigAlgo_;  //FIXME
}
//
// member functions
//
// ------------ method called for each event  ------------
void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  nTotal++;
  using namespace edm;
 // ----- Apply event selection cuts ----- //

  bool passedSelection = false;
  if ( doJets_ ) {
    passedSelection = runEvtSel_jet( iEvent, iSetup );
  } else {
    passedSelection = runEvtSel( iEvent, iSetup );
  }

  if ( !passedSelection ) {
    if ( debug ) std::cout << "!!!!!!!!!!! DID NOT PASS EVENT/JET SELECTION !!!!!!!!!!!" << std::endl;
    h_sel->Fill( 0. );;
    return;
  }

  fillEB( iEvent, iSetup );
  fillEE( iEvent, iSetup );
  fillHBHE( iEvent, iSetup );
  fillECALatHCAL( iEvent, iSetup );
  fillECALstitched( iEvent, iSetup );
  fillHCALatEBEE( iEvent, iSetup );
  fillTracksAtEBEE( iEvent, iSetup );
  for (unsigned int i=0;i<Nproj;i++)
  {
    fillTracksAtECALstitched( iEvent, iSetup, i );
  }
  fillPFCandsAtEBEE( iEvent, iSetup );
  fillPFCandsAtECALstitched( iEvent, iSetup );
  //fillTRKlayersAtEBEE( iEvent, iSetup );
  //fillTRKlayersAtECAL( iEvent, iSetup );
  //fillTRKvolumeAtEBEE( iEvent, iSetup );
  //fillTRKvolumeAtECAL( iEvent, iSetup );
  fillJetInfoAtECALstitched( iEvent, iSetup );
  fillScalarInfo( iEvent, iSetup );
  for (unsigned int i=0;i<Nhitproj;i++)
  {
    fillTRKlayersAtECALstitched( iEvent, iSetup, i );
  }


  ////////////// 4-Momenta //////////
  //fillFC( iEvent, iSetup );

  // Fill RHTree
  RHTree->Fill();
  h_sel->Fill( 1. );
  nPassed++;

} // analyze()


// ------------ method called once each job just before starting event loop  ------------
void 
RecHitAnalyzer::beginJob()
{
  nTotal = 0;
  nPassed = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecHitAnalyzer::endJob() 
{
  std::cout << " selected: " << nPassed << "/" << nTotal << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

const reco::PFCandidate*
RecHitAnalyzer::getPFCand(edm::Handle<PFCollection> pfCands, float eta, float phi, float& minDr, bool debug ) {

  minDr = 10;
  const reco::PFCandidate* minDRCand = nullptr;
  
  for ( PFCollection::const_iterator iPFC = pfCands->begin();
        iPFC != pfCands->end(); ++iPFC ) {

    const reco::Track* thisTrk = iPFC->bestTrack();
    if ( !thisTrk ) continue;

    float thisdR = reco::deltaR( eta, phi, thisTrk->eta(), thisTrk->phi() );
    if (debug) std::cout << "\tthisdR: " << thisdR << " " << thisTrk->pt() << " " << iPFC->particleId() << std::endl;

    const reco::PFCandidate& thisPFCand = (*iPFC);
      
    if ( (thisdR < 0.01) && (thisdR <minDr) ) {
      minDr    = thisdR; 
      minDRCand = &thisPFCand;
    }
  }

  return minDRCand;  
}

const reco::Track*
RecHitAnalyzer::getTrackCand(edm::Handle<reco::TrackCollection> trackCands, float eta, float phi, float& minDr, bool debug ) {

  minDr = 10;
  const reco::Track* minDRCand = nullptr;
  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  for ( reco::TrackCollection::const_iterator iTk = trackCands->begin();
        iTk != trackCands->end(); ++iTk ) {
    if ( !(iTk->quality(tkQt_)) ) continue;  

    float thisdR = reco::deltaR( eta, phi, iTk->eta(),iTk->phi() );
    if (debug) std::cout << "\tthisdR: " << thisdR << " " << iTk->pt() << std::endl;

    const reco::Track& thisTrackCand = (*iTk);
      
    if ( (thisdR < 0.01) && (thisdR <minDr) ) {
      minDr    = thisdR; 
      minDRCand = &thisTrackCand;
    }
  }

  return minDRCand;  
}




int RecHitAnalyzer::getTruthLabel(const reco::PFJetRef& recJet, edm::Handle<reco::GenParticleCollection> genParticles, float dRMatch , bool debug ){
  if ( debug ) {
    std::cout << " Mathcing reco jetPt:" << recJet->pt() << " jetEta:" << recJet->eta() << " jetPhi:" << recJet->phi() << std::endl;
  }

  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
       iGen != genParticles->end();
       ++iGen) {

    // Do not want to match to the final particles in the shower
    if ( iGen->status() > 99 ) continue;
    
    // Only want to match to partons/leptons/bosons
    if ( iGen->pdgId() > 25 ) continue;

    float dR = reco::deltaR( recJet->eta(),recJet->phi(), iGen->eta(),iGen->phi() );

    if ( debug ) std::cout << " \t >> dR " << dR << " id:" << iGen->pdgId() << " status:" << iGen->status() << " nDaught:" << iGen->numberOfDaughters() << " pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;

    if ( dR > dRMatch ) continue; 
    if ( debug ) std::cout << " Matched pdgID " << iGen->pdgId() << std::endl;

    return iGen->pdgId();

  } // gen particles 

  return -99;
}


float RecHitAnalyzer::getBTaggingValue(const reco::PFJetRef& recJet, edm::Handle<edm::View<reco::Jet> >& recoJetCollection, edm::Handle<reco::JetTagCollection>& btagCollection, float dRMatch, bool debug ){

  // loop over jets
  for( edm::View<reco::Jet>::const_iterator jetToMatch = recoJetCollection->begin(); jetToMatch != recoJetCollection->end(); ++jetToMatch )
    {
      reco::Jet thisJet = *jetToMatch;
      float dR = reco::deltaR( recJet->eta(),recJet->phi(), thisJet.eta(),thisJet.phi() );
      if(dR > 0.1) continue;

      size_t idx = (jetToMatch - recoJetCollection->begin());
      edm::RefToBase<reco::Jet> jetRef = recoJetCollection->refAt(idx);

      if(debug) std::cout << "btag discriminator value = " << (*btagCollection)[jetRef] << std::endl;
      return (*btagCollection)[jetRef];
  
    }

  if(debug){
    std::cout << "ERROR  No btag match: " << std::endl;
    
    // loop over jets
    for( edm::View<reco::Jet>::const_iterator jetToMatch = recoJetCollection->begin(); jetToMatch != recoJetCollection->end(); ++jetToMatch )
      {
	const reco::Jet thisJet = *jetToMatch;
	std::cout << "\t Match attempt pt: " <<  thisJet.pt() << " vs " <<  recJet->pt()
		  << " eta: " << thisJet.eta() << " vs " << recJet->eta()
		  << "phi: "<< thisJet.phi() << " vs " << recJet->phi()
		  << std::endl;
	float dR = reco::deltaR( recJet->eta(),recJet->phi(), thisJet.eta(),thisJet.phi() );
	std::cout << "dR " << dR << std::endl;
      }
  }    

  return -99;
}


Measurement1D RecHitAnalyzer::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
  VertexDistanceXY dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

Measurement1D RecHitAnalyzer::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
  VertexDistance3D dist;
  reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
  reco::Vertex svtx(svcand.vertex(), csv);
  return dist.distance(svtx, pv);
}

float RecHitAnalyzer::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
  reco::Candidate::Vector p = sv.momentum();
  reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
  return p.Unit().Dot(d.Unit());
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
