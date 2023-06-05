#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
TH1D *h_tau_gen_pT;
TH1D *h_tau_gen_prongs;
TH1D *h_tau_jet_pT;
TH1D *h_tau_jet_E;
TH1D *h_tau_jet_eta;
TH1D *h_tau_jet_m0;
TH1D *h_tau_jet_nJet;
TH1D *h_tau_jet_isTau;
TH1D *h_tau_jet_dR;
TH1D *h_tau_goodvertices;
vector<float> v_jetIsTau;
vector<float> v_jetdR;
vector<float> v_goodvertices;
vector<float> v_taupT;
vector<float> v_tauDaughters;
vector<float> v_tauJetPdgId;

vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_pt_;
vector<float> v_tau_gen_pt_;
vector<float> v_tau_gen_prongs_;
vector<float> v_tau_jetPdgIds_;
vector<float> v_tau_jetIsTau_;
vector<float> v_tau_jetdR_;
vector<float> v_tau_goodvertices_;

vector<float> v_tau_subJetE_[nJets];
vector<float> v_tau_subJetPx_[nJets];
vector<float> v_tau_subJetPy_[nJets];
vector<float> v_tau_subJetPz_[nJets];


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                   , 100,  0., 500.);
  h_tau_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                           , 100,  0., 500.);
  h_tau_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_tau_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_tau_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  h_tau_jet_isTau      = fs->make<TH1D>("h_jet_isTau"      , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_tau_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{jet,#tau};dR_{jet,#tau};Jets"           ,  50,  0.,   1.);
  h_tau_goodvertices   = fs->make<TH1D>("h_goodvertices"   , "good vertices;good vertices;Jets"           ,  15,  0.,  75.);
  h_tau_gen_pT         = fs->make<TH1D>("h_gen_pT"         , "p_{T};p_{T}; Gen part"                      ,  30,  0., 300.);
  h_tau_gen_prongs     = fs->make<TH1D>("h_gen_prongs"     , "prongs; prongs; Gen part"                   ,  10,  0.,  10.);

  tree->Branch("jet_M",         &v_tau_jet_m0_);
  tree->Branch("jet_Pt",        &v_tau_jet_pt_);
  tree->Branch("jet_PdgIds",    &v_tau_jetPdgIds_);
  tree->Branch("jet_IsTau",     &v_tau_jetIsTau_);
  tree->Branch("gen_pt",        &v_tau_gen_pt_);
  tree->Branch("gen_Prongs",    &v_tau_gen_prongs_);
  tree->Branch("jet_dR",        &v_tau_jetdR_);
  tree->Branch("goodvertices", &v_tau_goodvertices_);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_tau_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_tau_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_tau_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_tau_subJetPz_[iJ]);
  }

} // branchesEvtSel_jet_dijet_tau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_tau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);

  vJetIdxs.clear();
  v_tau_jetPdgIds_.clear();
  v_jetIsTau.clear();
  v_jetdR.clear();
  v_goodvertices.clear();
  v_taupT.clear();
  v_tauDaughters.clear();
  v_tauJetPdgId.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_tau_jetFakePhoIdxs;
  */

  unsigned int nMatchedJets = 0;
  unsigned int goodVertices = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
 
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;

  if ( debug ) std::cout << "\t" << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;

  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

    unsigned int PdgId        = 0;
    float jetdR               = -99.;
    float taupT               = -99.;
    int tauDaughters          = -1;
    int tauPi0                = -1;
    bool JetIsTau             = false;

    reco::PFJetRef iJet( jets, iJ );
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    //if (debug ) std::cout << "\t\t"<< " Jet [" << iJ << "] ->  Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0;
   
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.4 ) continue;

      //if (debug ) std::cout << "\t\t"<< " Jet [" << iJ << "] ->  Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
      if ( iGen->pt() > 20 && (std::abs(iGen->pdgId()) == 11 || std::abs(iGen->pdgId()) == 13) ) break; //only clean jets (lepton veto) 
      if ( std::abs(iGen->pdgId()) == 12 || std::abs(iGen->pdgId()) == 14 || std::abs(iGen->pdgId()) == 16 ) continue;

      if (isSignal_ && !(std::abs(iGen->pdgId()) == 15 && iGen->status() == 2) ) continue;      // for drell yan and HiggsToTauTau
      if ( !isSignal_ && !isW_ && !( iGen->status() == 23 ) ) continue;                         //for QCD background
      if ( !isSignal_ &&  isW_ && !( iGen->status() == 71 ) ) continue;                         //only for W + jet background

      //if ( !(std::abs(iGen->pdgId()) == 6) ) continue;

      if ( debug ) std::cout << "\t\t\t" << " GEN particle " << iGenParticle << ", status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << ", nMoms: " <<iGen->numberOfMothers() << ", mother ID: " << iGen->mother()->pdgId() << ", pt: "<< iGen->pt() << ", eta: " <<iGen->eta() << ", phi: " <<iGen->phi() << ", dR: " << dR << std::endl;
      if ( debug && iGen->numberOfDaughters() ==1 ) std::cout << "\t\t\t" << " Daughter 1 ID " << iGen->daughter(0)->pdgId() << std::endl;
      if ( debug && iGen->numberOfDaughters() ==2 ) std::cout << "\t\t\t" << " Daughter 1 ID " << iGen->daughter(0)->pdgId() << " Daughter 2 ID " << iGen->daughter(1)->pdgId() << std::endl;

      
     if ( std::abs(iGen->pdgId()) == 15 ) {
        JetIsTau = true;
        PdgId = std::abs(iGen->pdgId());
        jetdR = dR;
        taupT = iGen->pt();
        tauDaughters = 0;
        tauPi0 = 0; 

        for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
          if ( debug ) std::cout << "\t\t\t\t" <<" Tau daughter [" << iDaughter << "] : "<<  std::abs(iGen->daughter(iDaughter)->pdgId()); 
          if ( debug ) std::cout << " charge : "<< iGen->daughter(iDaughter)->charge() << "  | pt : "<< iGen->daughter(iDaughter)->pt();
          if ( debug ) std::cout << " eta:" << iGen->daughter(iDaughter)->eta() << " |Energy:" << iGen->daughter(iDaughter)->energy() << std::endl;
          if ( abs(iGen->daughter(iDaughter)->pdgId()) == 111 ) tauPi0++;
          if ( iGen->daughter(iDaughter)->charge() == 0 ) continue;          
          tauDaughters++;
        }
   
        if ( debug ) std::cout << "\t\t\t"<<" Tau prongs = " << tauDaughters << " + Tau pi0 = " << tauPi0 << std::endl;

 	if (!isSignal_){  
        passedGenSel = false; // for backg
        break; 		      // for backg
	}

      } // if pDG ID ==15
	else if ( taupT < iGen->pt() ){
        PdgId = std::abs(iGen->pdgId());
        jetdR = dR;
        taupT = iGen->pt();
      }

      passedGenSel = true;
      ++iGenParticle;

    } // primary gen particles
    if (passedGenSel) { 
      ++nMatchedJets;
      vJetIdxs.push_back( iJ );
      v_tauJetPdgId.push_back( PdgId );
      v_taupT.push_back( taupT );
      v_tauDaughters.push_back( tauDaughters );
      v_jetdR.push_back( jetdR );
      v_goodvertices.push_back( goodVertices );
      v_jetIsTau.push_back( JetIsTau );

    }

  } // reco jets
  if ( debug ) std::cout << "\t" <<" Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << "\t" <<" Event contains a tau candidate" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_jet_nJet->Fill( vJetIdxs.size() );

  v_tau_jet_pt_.clear();
  v_tau_gen_pt_.clear();
  v_tau_gen_prongs_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jetIsTau_.clear();
  v_tau_jetdR_.clear();
  v_tau_goodvertices_.clear();
  v_tau_jetPdgIds_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_jet_eta->Fill( iJet->eta() );
    h_tau_jet_E->Fill( iJet->energy() );
    h_tau_jet_m0->Fill( iJet->mass() );
    h_tau_jet_isTau->Fill( v_jetIsTau[iJ] );
    h_tau_jet_dR->Fill( v_jetdR[iJ] );
    h_tau_goodvertices->Fill( v_goodvertices[iJ] );
    h_tau_gen_pT->Fill( v_taupT[iJ] );
    h_tau_gen_prongs->Fill( v_tauDaughters[iJ] );

    // Fill branches 
    v_tau_jet_pt_.push_back( iJet->pt() );
    v_tau_jet_m0_.push_back( iJet->mass() );
    v_tau_jetIsTau_.push_back( v_jetIsTau[iJ] );
    v_tau_jetdR_.push_back( v_jetdR[iJ] );
    v_tau_goodvertices_.push_back( v_goodvertices[iJ] );
    v_tau_gen_pt_.push_back( v_taupT[iJ] );
    v_tau_gen_prongs_.push_back( v_tauDaughters[iJ] );
    v_tau_jetPdgIds_.push_back(v_tauJetPdgId[iJ]);

    // Gen jet constituents
    v_tau_subJetE_[iJ].clear();
    v_tau_subJetPx_[iJ].clear();
    v_tau_subJetPy_[iJ].clear();
    v_tau_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      //if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_dijet_tau()
