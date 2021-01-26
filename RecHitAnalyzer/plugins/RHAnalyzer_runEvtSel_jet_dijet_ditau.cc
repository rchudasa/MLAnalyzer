#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets      = 50; //TODO: use cfg level nJets_
TH1D *h_tau_att_jet_pT;
TH1D *h_tau_att_jet_E;
TH1D *h_tau_att_jet_eta;
TH1D *h_tau_att_jet_m0;
TH1D *h_tau_att_jet_ma;
TH1D *h_tau_att_jet_pta;
TH2D *h_tau_att_jet_a_m_pt;
TH1D *h_tau_att_jet_nJet;
TH1D *h_tau_att_jet_isDiTau;
TH1D *h_tau_att_jet_dR;
TH1D *h_tau_att_jet_TaudR;
TH1D *h_tau_att_jet_Tau1dR;
TH1D *h_tau_att_jet_Tau2dR;
TH1D *h_tau_att_jet_NrecoTaus;
TH1D *h_tau_att_jet_NGenTaus;
TH1D *h_tau_att_jet_recoTau1dR;
TH1D *h_tau_att_jet_recoTau2dR;
TH1D *h_tau_att_jet_n1dR;
TH1D *h_tau_att_jet_n2dR;
vector<float> v_att_tau_Idxs;
vector<float> v_att_tau_combs;
vector<float> v_att_jetIsDiTau;
vector<float> v_att_ma;
vector<float> v_att_pta;
vector<float> v_att_jetTaudR;
vector<float> v_att_jetTau1dR;
vector<float> v_att_jetTau2dR;
vector<float> v_att_jetNGenTaus;
vector<float> v_att_jetNrecoTaus;
vector<float> v_att_jetrecoTau1dR;
vector<float> v_att_jetrecoTau2dR;
vector<float> v_att_jetn1dR;
vector<float> v_att_jetn2dR;

vector<float> v_att_tau_jet_m0_;
vector<float> v_att_tau_jet_ma_;
vector<float> v_att_tau_jet_pta_;
vector<float> v_att_tau_jet_pt_;
vector<float> v_att_tau_jetIsDiTau_;
vector<float> v_att_tau_jetTaudR_;
vector<float> v_att_tau_jetTau1dR_;
vector<float> v_att_tau_jetTau2dR_;
vector<float> v_att_tau_jetNGenTaus_;
vector<float> v_att_tau_jetNrecoTaus_;
vector<float> v_att_tau_jetrecoTau1dR_;
vector<float> v_att_tau_jetrecoTau2dR_;
vector<float> v_att_tau_jetn1dR_;
vector<float> v_att_tau_jetn2dR_;

vector<float> v_att_tau_subJetE_[nJets];
vector<float> v_att_tau_subJetPx_[nJets];
vector<float> v_att_tau_subJetPy_[nJets];
vector<float> v_att_tau_subJetPz_[nJets];

TLorentzVector DiTau(Float_t tau1_pt, Float_t tau1_eta, Float_t tau1_phi, Float_t tau1_mass, Float_t tau2_pt, Float_t tau2_eta, Float_t tau2_phi, Float_t tau2_mass){
   TLorentzVector t1, t2, DiTau_Candidate;
    t1.SetPtEtaPhiM(tau1_pt, tau1_eta, tau1_phi, tau1_mass);
    t2.SetPtEtaPhiM(tau2_pt, tau2_eta, tau2_phi, tau2_mass);
    DiTau_Candidate = t1 + t2;
    return DiTau_Candidate;
}

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_ditau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_att_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                   , 100,  0., 800.);
  h_tau_att_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                           , 100,  0., 800.);
  h_tau_att_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_tau_att_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_tau_att_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  h_tau_att_jet_isDiTau    = fs->make<TH1D>("h_jet_isDiTau"    , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_tau_att_jet_TaudR      = fs->make<TH1D>("h_jet_TaudR"      , "dR_{#tau,#tau};dR_{#tau,#tau};Jets"         ,  50,  0.,   1.);

  tree->Branch("jetM",       &v_att_tau_jet_m0_);
  tree->Branch("jetIsDiTau", &v_att_tau_jetIsDiTau_);
  tree->Branch("jetpT",      &v_att_tau_jet_pt_);
  tree->Branch("TaudR",      &v_att_tau_jetTaudR_);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_att_tau_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_att_tau_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_att_tau_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_att_tau_subJetPz_[iJ]);
  }

} // branchesEvtSel_jet_dijet_ditau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_ditau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  edm::Handle<reco::PFTauDiscriminator> MVAIsolationLoose;
  iEvent.getByToken(tauMVAIsolationLoose_, MVAIsolationLoose);
  edm::Handle<reco::PFTauDiscriminator> MuonRejectionLoose;
  iEvent.getByToken(tauMuonRejectionLoose_, MuonRejectionLoose);
  edm::Handle<reco::PFTauDiscriminator> ElectronRejectionMVA6Loose;
  iEvent.getByToken(tauElectronRejectionMVA6Loose_, ElectronRejectionMVA6Loose);

  vJetIdxs.clear();
  v_att_tau_Idxs.clear();
  v_att_tau_combs.clear();
  v_att_jetIsDiTau.clear();
  v_att_jetTaudR.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_att_tau_jetFakePhoIdxs;
  */

  bool IsSignal             = true;
  //bool IsSignal             = false;


  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  if ( debug ) std::cout << " MVAIsolationLoose size is " << MVAIsolationLoose->size() << std::endl;

  //Lookin at RecoTaus
  if ( taus->size() < 2 ) {
    if (debug ) std::cout << "   !!!!!!!!!!  ONLY " << taus->size() << " TAUS IN THIS EVENT  !!!!!!!!!!"<< std::endl;
    return false;
  }
  if ( debug ) std::cout << " TAUS IN THE EVENT = " << taus->size() << " | Selection requires tau discrimator = NAME" << std::endl;
  for ( unsigned iT1(0); iT1 != taus->size(); ++iT1 ) {
    reco::PFTauRef iTau1( taus, iT1 );
    if (!((*MuonRejectionLoose)[iTau1])) continue;
    if (!((*ElectronRejectionMVA6Loose)[iTau1])) continue;
    if (!((*MVAIsolationLoose)[iTau1])) continue;
    if (debug ) std::cout << " TAU PASSED SELECTION "<< std::endl;
    unsigned int tau_combinations = 0;
    for ( unsigned iT2(0); iT2 != taus->size(); ++iT2 ) {
      if ( iT2 == iT1 ) continue;
      reco::PFTauRef iTau2( taus, iT2 );
      float recotaudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iTau2->eta(),iTau2->phi() );
      if ( recotaudR < 0.5 ) continue;
      float ditau_mass = DiTau(iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass(), iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()).M();
      if ( debug ) std::cout << " Tau pair " << iT1 << " + " << iT2 << " | ditau mass : " << ditau_mass << " GeV" << std::endl;
      //if ( IsSignal ) {
      //  if ( ditau_mass < 80 || ditau_mass > 160 ) {
      //    if ( debug ) std::cout << " TAU PAIR NOT IN MASS RANGE = " << v_att_tau_Idxs.size() << std::endl;
      //    continue;
      //  }
      //} else {
      //  if ( ditau_mass > 80 && ditau_mass < 160 ) continue;
      //}
      ++tau_combinations;
    }
    if ( tau_combinations == 0 ) continue;
    v_att_tau_Idxs.push_back( iT1 );
    v_att_tau_combs.push_back( tau_combinations );
  }
  
  if ( v_att_tau_Idxs.size() == 0 ) return false;
  if ( debug ) std::cout << " SELECTED TAUS = " << v_att_tau_Idxs.size() << std::endl;

  unsigned int nMatchedJets = 0;
  bool IsDitau = false;
  // Loop over jets
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << " ] -> Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;

    //loop over selected taus
    float nMatchedTaus = 0;
    float TaupT     = -99;
    float TaudR     = -99;
    //int TauIdx      = -99;
    for ( unsigned iT(0); iT != v_att_tau_Idxs.size(); ++iT ) {
      reco::PFTauRef iTau( taus, v_att_tau_Idxs[iT] );
      float recotaudR = reco::deltaR( iJet->eta(),iJet->phi(), iTau->eta(),iTau->phi() );
      if ( recotaudR > 0.5 ) continue;
      if ( iTau->pt() > TaupT ) {
        TaupT = iTau->pt();
        //TauIdx   = iT;
        TaudR   = recotaudR;
      }
      ++nMatchedTaus;
    } // end tau loop
    //filling jet variables
    if ( nMatchedTaus == 0 ) continue;
    if ( nMatchedTaus > 1 ) IsDitau = true;
    if ( debug ) std::cout << " JET MATCHED A TAU CANDIDATE " << nMatchedJets << std::endl;
    ++nMatchedJets;
    vJetIdxs.push_back( iJ );
    v_att_jetTaudR.push_back( TaudR );
    v_att_jetIsDiTau.push_back( IsDitau );
  } // end jet loop
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << " >> has_jet_dijet_ditau: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet_ditau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_ditau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_att_jet_nJet->Fill( vJetIdxs.size() );

  v_att_tau_jet_pt_.clear();
  v_att_tau_jet_m0_.clear();
  v_att_tau_jetIsDiTau_.clear();
  v_att_tau_jetTaudR_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_att_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_att_jet_eta->Fill( iJet->eta() );
    h_tau_att_jet_E->Fill( iJet->energy() );
    h_tau_att_jet_m0->Fill( iJet->mass() );
    h_tau_att_jet_isDiTau->Fill( v_att_jetIsDiTau[iJ] );
    h_tau_att_jet_TaudR->Fill( v_att_jetTaudR[iJ] );

    // Fill branches 
    v_att_tau_jet_pt_.push_back( iJet->pt() );
    v_att_tau_jet_m0_.push_back( iJet->mass() );
    v_att_tau_jetIsDiTau_.push_back( v_att_jetIsDiTau[iJ] );
    v_att_tau_jetTaudR_.push_back( v_att_jetTaudR[iJ] );

    // Gen jet constituents
    v_att_tau_subJetE_[iJ].clear();
    v_att_tau_subJetPx_[iJ].clear();
    v_att_tau_subJetPy_[iJ].clear();
    v_att_tau_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_att_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_att_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_att_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_att_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_dijet_ditau()
