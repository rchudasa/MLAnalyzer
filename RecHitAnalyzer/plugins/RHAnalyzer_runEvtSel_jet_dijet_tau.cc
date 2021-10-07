#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
TH1D *h_tau_gen_pT;
TH1D *h_tau_gen_prongs;
TH1D *h_tau_gen_pi0;
TH1D *h_tau_jet_pT;
TH1D *h_tau_jet_E;
TH1D *h_tau_jet_eta;
TH1D *h_tau_jet_m0;
TH1D *h_tau_jet_nJet;
TH1D *h_tau_jet_isTau;
TH1D *h_tau_jet_dR;
TH1D *h_tau_goodvertices;
TH1D *h_tau_loosedeeptau;
TH1D *h_tau_mediumdeeptau;
TH1D *h_tau_tightdeeptau;
vector<float> v_jetIsTau;
vector<float> v_jetdR;
vector<float> v_jetPdgIds;
vector<float> v_goodvertices;
vector<float> v_loosedeeptau;
vector<float> v_mediumdeeptau;
vector<float> v_tightdeeptau;
vector<float> v_taupT;
vector<float> v_tauDaughters;
vector<float> v_taupi0;

vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_pt_;
vector<float> v_tau_gen_pt_;
vector<float> v_tau_gen_prongs_;
vector<float> v_tau_gen_pi0_;
vector<float> v_tau_jetPdgIds_;
vector<float> v_tau_jetIsTau_;
vector<float> v_tau_jetdR_;
vector<float> v_tau_goodvertices_;
vector<float> v_tau_loosedeeptau_;
vector<float> v_tau_mediumdeeptau_;
vector<float> v_tau_tightdeeptau_;

//vector<float> v_tau_subJetE_[nJets];
//vector<float> v_tau_subJetPx_[nJets];
//vector<float> v_tau_subJetPy_[nJets];
//vector<float> v_tau_subJetPz_[nJets];


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_jet_E          = fs->make<TH1D>("h_jet_E"            , "E;E;Jets"                                   , 100,  0., 500.);
  h_tau_jet_pT         = fs->make<TH1D>("h_jet_pT"           , "p_{T};p_{T};Jets"                           , 100,  0., 500.);
  h_tau_jet_eta        = fs->make<TH1D>("h_jet_eta"          , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_tau_jet_nJet       = fs->make<TH1D>("h_jet_nJet"         , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_tau_jet_m0         = fs->make<TH1D>("h_jet_m0"           , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  h_tau_jet_isTau      = fs->make<TH1D>("h_jet_isTau"        , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_tau_jet_dR         = fs->make<TH1D>("h_jet_dR"           , "dR_{jet,#tau};dR_{jet,#tau};Jets"           ,  50,  0.,   1.);
  h_tau_goodvertices   = fs->make<TH1D>("h_goodvertices"     , "good vertices;good vertices;Jets"           ,  15,  0.,  75.);
  h_tau_gen_pT         = fs->make<TH1D>("h_gen_pT"           , "p_{T};p_{T}; Gen part"                      ,  30,  0., 300.);
  h_tau_gen_prongs     = fs->make<TH1D>("h_gen_prongs"       , "prongs; prongs; Gen part"                   ,  10,  0.,  10.);
  h_tau_gen_pi0        = fs->make<TH1D>("h_gen_pi0"          , "pi0; pi0; Gen part"                         ,  10,  0.,  10.);
  h_tau_loosedeeptau   = fs->make<TH1D>("h_tau_LooseDeepTau" , "LooseDeepTau; LooseDeepTau; Gen part"       ,   2,  0., 1.01);
  h_tau_mediumdeeptau  = fs->make<TH1D>("h_tau_MediumDeepTau", "MediumDeepTau; MediumDeepTau; Gen part"     ,   2,  0., 1.01);
  h_tau_tightdeeptau   = fs->make<TH1D>("h_tau_TightDeepTau" , "TightDeepTau; TightDeepTau; Gen part"       ,   2,  0., 1.01);

  tree->Branch("jet_M",         &v_tau_jet_m0_);
  tree->Branch("jet_Pt",        &v_tau_jet_pt_);
  tree->Branch("jet_PdgIds",    &v_tau_jetPdgIds_);
  tree->Branch("jet_IsTau",     &v_tau_jetIsTau_);
  tree->Branch("gen_pt",        &v_tau_gen_pt_);
  tree->Branch("gen_Prongs",    &v_tau_gen_prongs_);
  tree->Branch("gen_pi0",       &v_tau_gen_pi0_);
  tree->Branch("jet_dR",        &v_tau_jetdR_);
  tree->Branch("goodvertices",  &v_tau_goodvertices_);
  tree->Branch("LooseDeepTau",  &v_tau_loosedeeptau_);
  tree->Branch("MediumDeepTau", &v_tau_mediumdeeptau_);
  tree->Branch("TightDeepTau",  &v_tau_tightdeeptau_);

  //char hname[50];
  //for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
  //  sprintf(hname, "subJet%d_E", iJ);
  //  tree->Branch(hname,            &v_tau_subJetE_[iJ]);
  //  sprintf(hname, "subJet%d_Px", iJ);
  //  tree->Branch(hname,            &v_tau_subJetPx_[iJ]);
  //  sprintf(hname, "subJet%d_Py", iJ);
  //  tree->Branch(hname,            &v_tau_subJetPy_[iJ]);
  //  sprintf(hname, "subJet%d_Pz", iJ);
  //  tree->Branch(hname,            &v_tau_subJetPz_[iJ]);
  //}

} // branchesEvtSel_jet_dijet_tau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_tau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);


  vJetIdxs.clear();
  v_jetPdgIds.clear();
  v_jetIsTau.clear();
  v_jetdR.clear();
  v_goodvertices.clear();
  v_loosedeeptau.clear();
  v_mediumdeeptau.clear();
  v_tightdeeptau.clear();
  v_taupT.clear();
  v_tauDaughters.clear();
  v_taupi0.clear();

  unsigned int nMatchedJets = 0;
  unsigned int goodVertices = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
 
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << " good vertices in the event (PU) = " << goodVertices << std::endl;

  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;

  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {

    unsigned int PdgId        = 0;
    float jetdR               = -99.;
    float taupT               = -99.;
    int tauDaughters          = -1;
    int taupi0                = -1;
    bool JetIsTau             = false;
    float LooseDeepTau        = 0;
    float MediumDeepTau       = 0;
    float TightDeepTau        = 0;

    pat::Jet iJet = (*jets)[iJ];
    if ( std::abs(iJet.pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet.eta()) > maxJetEta_ ) continue;
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] ->  Pt: " << iJet.pt() << ", Eta: " << iJet.eta() << ", Phi: " << iJet.phi() << std::endl;
    if (!isData_) {
      bool passedGenSel = false;
      unsigned int iGenParticle = 0;
      for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
        float dR = reco::deltaR( iJet.eta(),iJet.phi(), iGen->eta(),iGen->phi() );
        if ( dR > 0.4 ) continue;

        if ( iGen->pt() > 20 && (std::abs(iGen->pdgId()) == 11 || std::abs(iGen->pdgId()) == 13) ) break; //only clean jets (lepton veto) 
        if ( std::abs(iGen->pdgId()) == 12 || std::abs(iGen->pdgId()) == 14 || std::abs(iGen->pdgId()) == 16 ) continue;

        if (  isSignal_ && !( std::abs(iGen->pdgId()) == 15 && iGen->status() == 2 ) ) continue;  //only for tau signal
        if ( !isSignal_ && !isW_ && !( iGen->status() == 23 ) ) continue;                         //for QCD background
        if ( !isSignal_ &&  isW_ && !( iGen->status() == 71 ) ) continue;                         //only for W + jet background

        if ( debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR: " << dR << std::endl;

        bool isHadronic = true;
        if ( std::abs(iGen->pdgId()) == 15 ) {
          for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
           if ( abs(iGen->daughter(iDaughter)->pdgId()) == 11 || abs(iGen->daughter(iDaughter)->pdgId()) == 13 ) isHadronic = false;
          }
          if (isSignal_ && !isHadronic) continue;
          JetIsTau = true;

          PdgId = std::abs(iGen->pdgId());
          jetdR = dR;
          taupT = iGen->pt();
          tauDaughters = 0;
          taupi0 = 0;
          for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
            if ( debug ) std::cout << "    Tau daughter [" << iDaughter << "] : "<<  std::abs(iGen->daughter(iDaughter)->pdgId()) << std::endl;
            if ( abs(iGen->daughter(iDaughter)->pdgId()) == 111 ) taupi0++;
            if ( iGen->daughter(iDaughter)->charge() == 0 ) continue;
            tauDaughters++;
          }
          if ( debug ) std::cout << "    Tau prongs = " << tauDaughters << " + Tau pi0 = " << taupi0 << std::endl;

          float recoTaudR = 0.4;
          for (const pat::Tau &tau : *taus) {
            //if(tau.pt() < 18.0) continue;
            //if(fabs(tau.eta()) > 2.3) continue;
            float TaudR = reco::deltaR( tau.eta(),tau.phi(), iGen->eta(),iGen->phi() );
            if ( TaudR > recoTaudR ) continue;
            recoTaudR = TaudR;
            LooseDeepTau  = tau.tauID("byLooseDeepTau2017v2p1VSjet");
            MediumDeepTau = tau.tauID("byMediumDeepTau2017v2p1VSjet");
            TightDeepTau  = tau.tauID("byTightDeepTau2017v2p1VSjet");
          }
          if ( debug ) std::cout << "   DeepTau loose  = " << LooseDeepTau << std::endl;
          if ( debug ) std::cout << "   DeepTau medium = " << MediumDeepTau << std::endl;
          if ( debug ) std::cout << "   DeepTau tight  = " << TightDeepTau << std::endl;

          if (!isSignal_){
            passedGenSel = false;  //only for background
            break;                 //only for background
          }

        } else if ( taupT < iGen->pt() ){
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
        v_jetPdgIds.push_back( PdgId );
        v_taupT.push_back( taupT );
        v_tauDaughters.push_back( tauDaughters );
        v_taupi0.push_back( taupi0 );
        v_jetdR.push_back( jetdR );
        v_goodvertices.push_back( goodVertices );
        v_loosedeeptau.push_back( LooseDeepTau );
        v_mediumdeeptau.push_back( MediumDeepTau );
        v_tightdeeptau.push_back( TightDeepTau );
        v_jetIsTau.push_back( JetIsTau );

      }
    } // End MC selection

  } // reco jets
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << " >> Event contains a tau candidate" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_jet_nJet->Fill( vJetIdxs.size() );

  v_tau_jet_pt_.clear();
  v_tau_gen_pt_.clear();
  v_tau_gen_prongs_.clear();
  v_tau_gen_pi0_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jetIsTau_.clear();
  v_tau_jetPdgIds_.clear();
  v_tau_jetdR_.clear();
  v_tau_goodvertices_.clear();
  v_tau_loosedeeptau_.clear();
  v_tau_mediumdeeptau_.clear();
  v_tau_tightdeeptau_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    pat::Jet iJet = (*jets)[iJ];

    // Fill histograms 
    h_tau_jet_pT->Fill( std::abs(iJet.pt()) );
    h_tau_jet_eta->Fill( iJet.eta() );
    h_tau_jet_E->Fill( iJet.energy() );
    h_tau_jet_m0->Fill( iJet.mass() );
    h_tau_jet_isTau->Fill( v_jetIsTau[iJ] );
    h_tau_jet_dR->Fill( v_jetdR[iJ] );
    h_tau_goodvertices->Fill( v_goodvertices[iJ] );
    h_tau_loosedeeptau->Fill( v_loosedeeptau[iJ] );
    h_tau_mediumdeeptau->Fill( v_mediumdeeptau[iJ] );
    h_tau_tightdeeptau->Fill( v_tightdeeptau[iJ] );
    h_tau_gen_pT->Fill( v_taupT[iJ] );
    h_tau_gen_prongs->Fill( v_tauDaughters[iJ] );
    h_tau_gen_pi0->Fill( v_taupi0[iJ] );

    // Fill branches 
    v_tau_jet_pt_.push_back( iJet.pt() );
    v_tau_jet_m0_.push_back( iJet.mass() );
    v_tau_jetIsTau_.push_back( v_jetIsTau[iJ] );
    v_tau_jetPdgIds_.push_back( v_jetPdgIds[iJ] );
    v_tau_jetdR_.push_back( v_jetdR[iJ] );
    v_tau_goodvertices_.push_back( v_goodvertices[iJ] );
    v_tau_loosedeeptau_.push_back( v_loosedeeptau[iJ] );
    v_tau_mediumdeeptau_.push_back( v_mediumdeeptau[iJ] );
    v_tau_tightdeeptau_.push_back( v_tightdeeptau[iJ] );
    v_tau_gen_pt_.push_back( v_taupT[iJ] );
    v_tau_gen_prongs_.push_back( v_tauDaughters[iJ] );
    v_tau_gen_pi0_.push_back( v_taupi0[iJ] );

    // Gen jet constituents
    //v_tau_subJetE_[iJ].clear();
    //v_tau_subJetPx_[iJ].clear();
    //v_tau_subJetPy_[iJ].clear();
    //v_tau_subJetPz_[iJ].clear();
    //unsigned int nConstituents = iJet.numberOfDaughters();
    //if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    //for ( unsigned int j = 0; j < nConstituents; j++ ) {
    //  if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << iJet.daughter(j).energy() << " px:" << iJet.daughter(j).px() << " py:" << iJet.daughter(j).py() << " pz:" << iJet.daughter(j).pz() << std::endl;
    //  v_tau_subJetE_[iJ].push_back( iJet.daughter(j).energy() );
    //  v_tau_subJetPx_[iJ].push_back( iJet.daughter(j).px() );
    //  v_tau_subJetPy_[iJ].push_back( iJet.daughter(j).py() );
    //  v_tau_subJetPz_[iJ].push_back( iJet.daughter(j).pz() );
    //}
  }

} // fillEvtSel_jet_dijet_tau()
