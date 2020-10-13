#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 2; //TODO: use cfg level nJets_
TH1D *h_tau_jet_pT;
TH1D *h_tau_jet_E;
TH1D *h_tau_jet_eta;
TH1D *h_tau_jet_m0;
TH1D *h_tau_jet_ma;
TH1D *h_tau_jet_nJet;
TH1D *h_tau_jet_isDiTau;
TH1D *h_tau_jet_dR;
TH1D *h_tau_jet_TaudR;
TH1D *h_tau_jet_Tau1dR;
TH1D *h_tau_jet_Tau2dR;
TH1D *h_tau_jet_NrecoTaus;
TH1D *h_tau_jet_recoTau1dR;
TH1D *h_tau_jet_recoTau2dR;
TH1D *h_tau_jet_n1dR;
TH1D *h_tau_jet_n2dR;
vector<float> v_jetIsDiTau;
vector<float> v_jetdR;
vector<float> v_ma;
vector<float> v_jetTaudR;
vector<float> v_jetTau1dR;
vector<float> v_jetTau2dR;
vector<float> v_jetNrecoTaus;
vector<float> v_jetrecoTau1dR;
vector<float> v_jetrecoTau2dR;
vector<float> v_jetn1dR;
vector<float> v_jetn2dR;

vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_ma_;
vector<float> v_tau_jet_pt_;
vector<float> v_tau_jetPdgIds_;
vector<float> v_tau_jetIsDiTau_;
vector<float> v_tau_jetdR_;
vector<float> v_tau_jetTaudR_;
vector<float> v_tau_jetTau1dR_;
vector<float> v_tau_jetTau2dR_;
vector<float> v_tau_jetNrecoTaus_;
vector<float> v_tau_jetrecoTau1dR_;
vector<float> v_tau_jetrecoTau2dR_;
vector<float> v_tau_jetn1dR_;
vector<float> v_tau_jetn2dR_;

vector<float> v_tau_subJetE_[nJets];
vector<float> v_tau_subJetPx_[nJets];
vector<float> v_tau_subJetPy_[nJets];
vector<float> v_tau_subJetPz_[nJets];

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                     , 100,  0., 500.);
  h_tau_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                             , 100,  0., 800.);
  h_tau_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                       , 100, -5.,   5.);
  h_tau_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                     ,  10,  0.,  10.);
  h_tau_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                 , 100,  0., 100.);
  h_tau_jet_ma         = fs->make<TH1D>("h_jet_ma"         , "m_{a};m_{a};Jets"                     , 100,  0., 100.);
  h_tau_jet_isDiTau    = fs->make<TH1D>("h_jet_isDiTau"    , "nIsDiTau;nIsDiTau;Jets"               ,  10,  0.,  10.);
  h_tau_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{a,j};dR_{a,j};Jets"               ,  50,  0.,  0.5);
  h_tau_jet_TaudR      = fs->make<TH1D>("h_jet_TaudR"      , "dR_{#tau,#tau};dR_{#tau,#tau};Jets"   ,  50,  0.,  0.5);
  h_tau_jet_Tau1dR     = fs->make<TH1D>("h_jet_Tau1dR"     , "dR_{#tau_{1},j};dR_{#tau_{1},j};Jets" ,  50,  0.,  0.5);
  h_tau_jet_Tau2dR     = fs->make<TH1D>("h_jet_Tau2dR"     , "dR_{#tau_{2},j};dR_{#tau_{2},j};Jets" ,  50,  0.,  0.5);
  h_tau_jet_NrecoTaus  = fs->make<TH1D>("h_jet_NrecoTaus"  , "N#tau^{RECO};N#tau^{RECO};Jets"       ,   5, 0.,   5.);
  h_tau_jet_recoTau1dR = fs->make<TH1D>("h_jet_recoTau1dR" , "dR_{#tau_{1}^{RECO},j};dR_{#tau_{1}^{RECO},j};Jets" ,  50,  0.,  0.5);
  h_tau_jet_recoTau2dR = fs->make<TH1D>("h_jet_recoTau2dR" , "dR_{#tau_{2}^{RECO},j};dR_{#tau_{2}^{RECO},j};Jets" ,  25,  0.,  0.5);
  h_tau_jet_n1dR       = fs->make<TH1D>("h_jet_n1dR"       , "dR_{#eta_{1},j};dR_{#eta_{1},j};Jets" ,  25,  0.,  0.5);
  h_tau_jet_n2dR       = fs->make<TH1D>("h_jet_n2dR"       , "dR_{#eta_{2},j};dR_{#eta_{2},j};Jets" ,  25,  0.,  0.5);

  tree->Branch("jetM",       &v_tau_jet_m0_);
  tree->Branch("jetPt",      &v_tau_jet_pt_);
  tree->Branch("jetPdgIds",  &v_tau_jetPdgIds_);
  tree->Branch("jetdR",      &v_tau_jetdR_);
  tree->Branch("jetIsDiTau", &v_tau_jetIsDiTau_);
  tree->Branch("a_m",        &v_tau_jet_ma_);
  tree->Branch("TaudR",      &v_tau_jetTaudR_);
  tree->Branch("Tau1dR",     &v_tau_jetTau1dR_);
  tree->Branch("Tau2dR",     &v_tau_jetTau2dR_);
  tree->Branch("NrecoTaus",  &v_tau_jetNrecoTaus_);
  tree->Branch("recoTau1dR", &v_tau_jetrecoTau1dR_);
  tree->Branch("recoTau2dR", &v_tau_jetrecoTau2dR_);
  tree->Branch("n1dR",       &v_tau_jetn1dR_);
  tree->Branch("n2dR",       &v_tau_jetn2dR_);

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
  float dR;

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  vJetIdxs.clear();
  v_tau_jetPdgIds_.clear();
  v_jetIsDiTau.clear();
  v_jetdR.clear();
  v_ma.clear();
  v_jetTaudR.clear();
  v_jetTau1dR.clear();
  v_jetTau2dR.clear();
  v_jetNrecoTaus.clear();
  v_jetrecoTau1dR.clear();
  v_jetrecoTau2dR.clear();
  v_jetn1dR.clear();
  v_jetn2dR.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_tau_jetFakePhoIdxs;
  */

  unsigned int nMatchedJets = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    bool passedGenSel = false;
    bool MatchedPseudoScalar = false;
    unsigned int iGenParticle = 0; 
    unsigned int NTau1Daughters = 0;
    unsigned int NTau2Daughters = 0;
    unsigned int nMatchedTaus = 0;
    float a_mass = -99.;
    float tausdR = -99.;
    float tau1dR = -99.;
    float tau2dR = -99.;
    float recotau1dR = -99.;
    float recotau2dR = -99.;
    float n1dR   = -99.;
    float n2dR   = -99.;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      //if ( iJ == 0 && debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;
      if ( iGen->status() != 22 ) continue;
      if ( abs(iGen->pdgId()) != 25 ) continue;
      if ( iGen->numberOfDaughters() != 2 ) continue;
      if ( abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;
      passedGenSel = true;
      //if ( iGen->mother()->pdgId() != 25 ) continue;
      //if ( iGen->status() != 2 && iGen->numberOfMothers() == 1 ) continue; // pythia8: 23 = outgoing particle from hard scatter
      ++iGenParticle;
      if ( debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.4 ) continue;
      ++nMatchedJets;
      if ( abs(iGen->pdgId()) == 25 ) a_mass = iGen->mass();
      if (abs(iGen->daughter(0)->pdgId()) == 15 && abs(iGen->daughter(1)->pdgId()) == 15) MatchedPseudoScalar = true;
      if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] matched particle [" << iGenParticle << "] -> pdgId: " << std::abs(iGen->pdgId()) << " | dR: " << dR << "| Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
      tausdR = reco::deltaR( iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
      tau1dR = reco::deltaR( iJet->eta(),iJet->phi()                          , iGen->daughter(0)->eta(),iGen->daughter(0)->phi() );
      tau2dR = reco::deltaR( iJet->eta(),iJet->phi()                          , iGen->daughter(1)->eta(),iGen->daughter(1)->phi() );
      if ( debug ) std::cout << "   >>>>>> Taus dR = " << tausdR << " , tau1 dR = " << tau1dR << " , tau2 dR = " << tau2dR << std::endl;
      NTau1Daughters = iGen->daughter(0)->numberOfDaughters();
      NTau2Daughters = iGen->daughter(1)->numberOfDaughters();
      if ( debug ) std::cout << "    >>>>>> # Tau 1 daughters = " << NTau1Daughters << ",  # Tau 2 daughters = "<< NTau2Daughters << std::endl;
      for (unsigned int iDaughter = 0; iDaughter != NTau1Daughters; ++iDaughter){
        if ( debug ) std::cout << "   >>>>>> Tau 1 Decay = " << iGen->daughter(0)->daughter(iDaughter)->pdgId() << std::endl;
        if ( abs(iGen->daughter(0)->daughter(iDaughter)->pdgId()) == 16 ) n1dR   = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(0)->daughter(iDaughter)->eta(), iGen->daughter(0)->daughter(iDaughter)->phi() );
      }
      for (unsigned int iDaughter = 0; iDaughter != NTau2Daughters; ++iDaughter){
        if ( debug ) std::cout << "   >>>>>> Tau 2 Decay = " << iGen->daughter(1)->daughter(iDaughter)->pdgId() << std::endl;
        if ( abs(iGen->daughter(1)->daughter(iDaughter)->pdgId()) == 16 ) n2dR   = reco::deltaR( iJet->eta(),iJet->phi(), iGen->daughter(1)->daughter(iDaughter)->eta(), iGen->daughter(1)->daughter(iDaughter)->phi() );
      }
      if (debug ) std::cout << "   >>>>>> n1 dR = " << n1dR << " , n2 dR = " << n2dR << std::endl;
      vJetIdxs.push_back( iJ );
      v_tau_jetPdgIds_.push_back( std::abs(iGen->pdgId()) );
    } // primary gen particles

    //Lookin at RecoTaus
    if (passedGenSel){
      if (taus->size() == 0) {
        if (debug ) std::cout << "   !!!!!!!!!!  NO RECO TAUS IN THIS EVENT  !!!!!!!!!!"<< std::endl;
      }
      for ( unsigned iT(0); iT != taus->size(); ++iT ) {
        reco::PFTauRef iTau( taus, iT );
        float recotaudR = reco::deltaR( iJet->eta(),iJet->phi(), iTau->eta(),iTau->phi() );
        if ( recotaudR < 0.5 && nMatchedTaus == 0 ) {
          if ( debug ) std::cout << "Reco Tau [" << iT << "]  matched jet [" << iJ << "] : dR = " << recotaudR << std::endl;
          recotau1dR = recotaudR;
          ++nMatchedTaus;
        } else if ( recotaudR < 0.5 && nMatchedTaus == 1 ) {
          if ( debug ) std::cout << "Reco Tau [" << iT << "]  matched jet [" << iJ << "] : dR = " << recotaudR << std::endl;
          if (recotaudR < recotau1dR) {
            recotau2dR = recotau1dR;
            recotau1dR = recotaudR;
          } else recotau2dR = recotaudR;
          ++nMatchedTaus;
        } else if ( debug && recotaudR < 0.5 && nMatchedTaus > 1 ) {
          std::cout << "   !!!!!!!!!!  FOUND MORE THAN 2 TAUS INSIDE JET CONE OF 0.5 !!!!!!!!!!"<< std::endl;
          if (recotaudR < recotau2dR && recotaudR < recotau1dR) { 
            if (recotau1dR < recotau2dR) recotau2dR = recotau1dR;
            recotau1dR = recotaudR;
          } else if (recotaudR < recotau2dR && recotaudR > recotau1dR) recotau2dR = recotaudR; 
          ++nMatchedTaus;
        } else if ( debug ) {
          std::cout << "   !!!!!!!!!!  NO MATCH FOR Reco Tau [" << iT << "]  with jet [" << iJ << "] : dR = " << recotaudR << std::endl;
        }
      }
    }
    v_jetdR.push_back( dR );
    v_jetTaudR.push_back( tausdR );
    v_ma.push_back( a_mass );
    v_jetTau1dR.push_back( tau1dR );
    v_jetTau2dR.push_back( tau2dR );
    v_jetn1dR.push_back( n1dR );
    v_jetn2dR.push_back( n2dR );
    v_jetNrecoTaus.push_back( nMatchedTaus );
    v_jetrecoTau1dR.push_back( recotau1dR );
    v_jetrecoTau2dR.push_back( recotau2dR );

    if ( MatchedPseudoScalar ){
      v_jetIsDiTau.push_back(true);
    }
    else {
      v_jetIsDiTau.push_back(false);
    }

  } // reco jets
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl; 

  // Check jet multiplicity
  if ( nMatchedJets != 2 ) return false;
  if ( debug ) std::cout << " >> has_jet_dijet_tau: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_jet_nJet->Fill( vJetIdxs.size() );

  v_tau_jet_pt_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jet_ma_.clear();
  v_tau_jetIsDiTau_.clear();
  v_tau_jetdR_.clear();
  v_tau_jetTaudR_.clear();
  v_tau_jetTau1dR_.clear();
  v_tau_jetTau2dR_.clear();
  v_tau_jetNrecoTaus_.clear();
  v_tau_jetrecoTau1dR_.clear();
  v_tau_jetrecoTau2dR_.clear();
  v_tau_jetn1dR_.clear();
  v_tau_jetn2dR_.clear();
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_jet_eta->Fill( iJet->eta() );
    h_tau_jet_E->Fill( iJet->energy() );
    h_tau_jet_m0->Fill( iJet->mass() );
    h_tau_jet_isDiTau->Fill( v_jetIsDiTau[iJ] );
    h_tau_jet_dR->Fill( v_jetdR[iJ] );
    h_tau_jet_ma->Fill( v_ma[iJ] );
    h_tau_jet_TaudR->Fill( v_jetTaudR[iJ] );
    h_tau_jet_Tau1dR->Fill( v_jetTau1dR[iJ] );
    h_tau_jet_Tau2dR->Fill( v_jetTau2dR[iJ] );
    h_tau_jet_NrecoTaus->Fill( v_jetNrecoTaus[iJ] );
    h_tau_jet_recoTau1dR->Fill( v_jetrecoTau1dR[iJ] );
    h_tau_jet_recoTau2dR->Fill( v_jetrecoTau2dR[iJ] );
    h_tau_jet_n1dR->Fill( v_jetn1dR[iJ] );
    h_tau_jet_n2dR->Fill( v_jetn2dR[iJ] );

    // Fill branches 
    v_tau_jet_pt_.push_back( iJet->pt() );
    v_tau_jet_m0_.push_back( iJet->mass() );
    v_tau_jet_ma_.push_back( v_ma[iJ] );
    v_tau_jetIsDiTau_.push_back( v_jetIsDiTau[iJ] );
    v_tau_jetdR_.push_back( v_jetdR[iJ] );
    v_tau_jetTaudR_.push_back( v_jetTaudR[iJ] );
    v_tau_jetTau1dR_.push_back( v_jetTau1dR[iJ] );
    v_tau_jetTau2dR_.push_back( v_jetTau2dR[iJ] );
    v_tau_jetNrecoTaus_.push_back( v_jetNrecoTaus[iJ] );
    v_tau_jetrecoTau1dR_.push_back( v_jetrecoTau1dR[iJ] );
    v_tau_jetrecoTau2dR_.push_back( v_jetrecoTau2dR[iJ] );
    v_tau_jetn1dR_.push_back( v_jetn1dR[iJ] );
    v_tau_jetn2dR_.push_back( v_jetn2dR[iJ] );

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
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_dijet_tau()
