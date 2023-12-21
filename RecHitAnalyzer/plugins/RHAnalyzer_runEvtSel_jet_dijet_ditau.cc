#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;
using namespace trigger;

const unsigned nJets      = 50; //TODO: use cfg level nJets_
TH1D *h_tau_att_genHiggs_M;
TH1D *h_tau_att_genTau1_pT;
TH1D *h_tau_att_genTau2_pT;
TH1D *h_tau_att_jet_pT;
TH1D *h_tau_att_jet_E;
TH1D *h_tau_att_jet_eta;
TH1D *h_tau_att_jet_m0;
TH1D *h_tau_att_jet_nJet;
TH1D *h_tau_att_jet_isSignal;
TH1D *h_tau_att_jet_dR;
TH1D *h_tau_att_dRtautau;
TH1D *h_tau_att_Mvis;
TH1D *h_tau_att_Mtautau;
TH1D *h_tau_att_Mtautau_svFit;
TH1D *h_tau_att_Mth_svFit;
TH1D *h_tau_att_pTh_svFit;
TH1D *h_tau_att_mth;
TH1D *h_tau_att_mcoll;
TH1D *h_tau_att_pTvis;
TH1D *h_tau_att_pTh;
TH1D *h_tau_att_tau_pT;
TH1D *h_tau_att_tau_mva;
TH1D *h_tau_att_tau_dm;
TH1D *h_tau_att_pfMET;
TH1D *h_tau_att_dphillmet;
TH1D *h_tau_att_dphill;

float v_att_genHiggs_M;
float v_att_genTau1_pT;
float v_att_genTau2_pT;
float v_att_genTau3_pT;
float v_att_genTau4_pT;
float v_att_dRtautau;
float v_att_pfMET;
float v_att_dphillmet;
float v_att_dphill;
float v_att_Mvis;
float v_att_pTvis;
float v_att_mth;
float v_att_mcoll;
float v_att_Mtautau;
float v_att_Mtautau_svFit;
float v_att_Mth_svFit;
float v_att_pTh_svFit;
float v_att_pTh;
vector<float> v_att_gentau_Idxs;
vector<float> v_att_tau_Idxs;
vector<float> v_att_tau_combs;
vector<float> v_att_tau_pT;
vector<float> v_att_tau_mva;
vector<float> v_att_tau_dm;
vector<float> v_att_jetIsSignal;
vector<float> v_att_jetdR;

vector<float> v_att_genHiggs_M_;
vector<float> v_att_genTau1_pT_;
vector<float> v_att_genTau2_pT_;
vector<float> v_att_Mvis_;
vector<float> v_att_pTvis_;
vector<float> v_att_Mtautau_;
vector<float> v_att_Mtautau_svFit_;
vector<float> v_att_Mth_svFit_;
vector<float> v_att_pTh_svFit_;
vector<float> v_att_pTh_;
vector<float> v_att_tau_pT_;
vector<float> v_att_tau_mva_;
vector<float> v_att_tau_dm_;
vector<float> v_att_mth_;
vector<float> v_att_mcoll_;
vector<float> v_att_dRtautau_;
vector<float> v_att_pfMET_;
vector<float> v_att_dphillmet_;
vector<float> v_att_dphill_;
vector<float> v_att_tau_jet_m0_;
vector<float> v_att_tau_jet_pt_;
vector<float> v_att_tau_jetIsSignal_;
vector<float> v_att_tau_jetdR_;

vector<float> v_att_tau_subJetE_[nJets];
vector<float> v_att_tau_subJetPx_[nJets];
vector<float> v_att_tau_subJetPy_[nJets];
vector<float> v_att_tau_subJetPz_[nJets];

TLorentzVector SetTau(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
  TLorentzVector Tau_Candidate;
  Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
  return Tau_Candidate;
}

TLorentzVector SetMET(Float_t met, Float_t metphi){
  TLorentzVector Met;
  Met.SetPtEtaPhiM(met, 0, metphi, 0.);
  return Met;
}

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_ditau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_att_genHiggs_M     = fs->make<TH1D>("h_genHiggs_M"   , "m^{gen H};m^{gen H};Events"                     ,  10,  120, 130);
  h_tau_att_genTau1_pT     = fs->make<TH1D>("h_genTau1_pT"   , "p_{T}^{gen #tau_{1}};p_{T}^{gen #tau_{1}};Events" ,  50,  0., 150.);
  h_tau_att_genTau2_pT     = fs->make<TH1D>("h_genTau2_pT"   , "p_{T}^{gen #tau_{2}};p_{T}^{gen #tau_{2}};Events" ,  50,  0., 150.);
  h_tau_att_jet_E          = fs->make<TH1D>("h_jet_E"        , "E;E;Jets"                                       ,  50,  0., 500.);
  h_tau_att_jet_pT         = fs->make<TH1D>("h_jet_pT"       , "p_{T};p_{T};Jets"                               ,  30,  0., 300.);
  h_tau_att_jet_eta        = fs->make<TH1D>("h_jet_eta"      , "#eta;#eta;Jets"                                 ,  30, -3.,   3.);
  h_tau_att_jet_m0         = fs->make<TH1D>("h_jet_m0"       , "m_{jet};m_{jet};Jets"                           ,  50,  0.,  50.);
  h_tau_att_jet_isSignal   = fs->make<TH1D>("h_jet_isSignal" , "IsSignal;IsSignal;Jets"                         ,   2,  0.,   2.);
  h_tau_att_jet_dR         = fs->make<TH1D>("h_jet_dR"       , "dR_{jet,#tau};dR_{jet,#tau};Jets"               ,  40,  0.,  0.4);
  h_tau_att_jet_nJet       = fs->make<TH1D>("h_jet_nJet"     , "nJet;nJet;Events"                               ,  10,  0.,  10.);
  h_tau_att_dRtautau       = fs->make<TH1D>("h_dRtautau"     , "dR_{#tau,#tau};dR_{#tau,#tau};Events"           ,  50,  0.,   5.);
  h_tau_att_Mvis           = fs->make<TH1D>("h_Mvis"         , "m^{vis};m^{vis};Events"                         ,  40,  0., 200.);
  h_tau_att_pTvis          = fs->make<TH1D>("h_pTvis"        , "p_{T}^{vis};p_{T}^{vis};Events"                 ,  40,  0., 200.);
  h_tau_att_Mtautau        = fs->make<TH1D>("h_Mtautau"      , "m^{#tau#tau};m^{#tau#tau};Events"               ,  40,  0., 400.);
  h_tau_att_Mtautau_svFit  = fs->make<TH1D>("h_Mtautau_svFit", "m^{#tau#tau}(SVFit);m^{#tau#tau}(SVFit);Events" ,  40,  0., 400.);
  h_tau_att_Mth_svFit      = fs->make<TH1D>("h_Mth_svFit"    , "m_{T}^{H}(SVFit);m_{T}^{H}(SVFit);Events"       ,  40,  0., 400.);
  h_tau_att_pTh_svFit      = fs->make<TH1D>("h_pTh_svFit"    , "p_{T}^{H}(SVFit);p_{T}^{H}(SVFit);Events"       ,  40,  0., 400.);
  h_tau_att_pTh            = fs->make<TH1D>("h_pTh"          , "p_{T}^{#tau#tau};p_{T}^{#tau#tau};Events"       ,  40,  0., 200.);
  h_tau_att_tau_pT         = fs->make<TH1D>("h_tau_pT"       , "p_{T}^{#tau};p_{T}^{#tau};Jets"                 ,  40,  0., 200.);
  h_tau_att_tau_mva        = fs->make<TH1D>("h_tau_mva"      , "#tau MVA ID;#tau MVA ID;Jets"                   ,  21, -1.,   1.);
  h_tau_att_tau_dm         = fs->make<TH1D>("h_tau_dm"       , "#tau Decay Mode;#tau Decay Mode;Jets"           ,  21, -10,  10.);
  h_tau_att_mth            = fs->make<TH1D>("h_mth"          , "m_{T}^{H};m_{T}^{H};Events"                     ,  40,  0., 400.);
  h_tau_att_mcoll          = fs->make<TH1D>("h_mcoll"        , "m_{coll}^{H};m_{coll}^{H};Events"               ,  40,  0., 400.);
  h_tau_att_pfMET          = fs->make<TH1D>("h_pfMET"        , "p_{T}^{miss};p_{T}^{miss};Events"               ,  40,  0., 200.);
  h_tau_att_dphill         = fs->make<TH1D>("h_dphill"       , "#Delta#phi_{l,l};#Delta#phi_{l,l};Events"       ,  40,  0., 3.1416);
  h_tau_att_dphillmet      = fs->make<TH1D>("h_dphillmet"    , "#Delta#phi_{#tau#tau,pfMET};#Delta #phi_{#tau#tau,pfMET};Events" ,40,0.,3.1416);

  tree->Branch("GenHiggsM",        &v_att_genHiggs_M_);
  tree->Branch("GenTau1pT",        &v_att_genTau1_pT_);
  tree->Branch("GenTau2pT",        &v_att_genTau2_pT_);
  tree->Branch("jetM",             &v_att_tau_jet_m0_);
  tree->Branch("jetIsSignal",      &v_att_tau_jetIsSignal_);
  tree->Branch("jetpT",            &v_att_tau_jet_pt_);
  tree->Branch("dR",               &v_att_tau_jetdR_);
  tree->Branch("TaudR",            &v_att_dRtautau_);
  tree->Branch("TauMvis",          &v_att_Mvis_);
  tree->Branch("TauMtautau",       &v_att_Mtautau_);
  tree->Branch("TauMtautau_svFit", &v_att_Mtautau_svFit_);
  tree->Branch("TauMth_svFit",     &v_att_Mth_svFit_);
  tree->Branch("TaupTh_svFit",     &v_att_pTh_svFit_);
  tree->Branch("TaupTvis",         &v_att_pTvis_);
  tree->Branch("TaupTh",           &v_att_pTh_);
  tree->Branch("TaupT",            &v_att_tau_pT_);
  tree->Branch("TauMVA",           &v_att_tau_mva_);
  tree->Branch("TauDM",            &v_att_tau_dm_);
  tree->Branch("Taumth",           &v_att_mth_);
  tree->Branch("Taumcoll",         &v_att_mcoll_);
  tree->Branch("pfMET",            &v_att_pfMET_);
  tree->Branch("dphillmet",        &v_att_dphillmet_);
  tree->Branch("dphill",           &v_att_dphill_);

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
  if(isMC_){
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );
  }
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<edm::View<reco::Jet> > Jets;
  iEvent.getByToken( pfjetsToken_, Jets );

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  edm::Handle<reco::PFMETCollection> pfmet;
  iEvent.getByToken(metCollectionT_, pfmet);

  edm::Handle<reco::PFTauDiscriminator> DecayMode;
  iEvent.getByToken(tauDecayMode_, DecayMode);
  edm::Handle<reco::PFTauDiscriminator> MVAIsolation;
  iEvent.getByToken(tauMVAIsolation_, MVAIsolation);
  edm::Handle<reco::PFTauDiscriminator> MuonRejection;
  iEvent.getByToken(tauMuonRejection_, MuonRejection);
  edm::Handle<reco::PFTauDiscriminator> ElectronRejectionMVA6;
  iEvent.getByToken(tauElectronRejectionMVA6_, ElectronRejectionMVA6);

  JME::JetResolution resPtObj            = JME::JetResolution::get(iSetup, jetResPtType_);
  JME::JetResolution resPhiObj           = JME::JetResolution::get(iSetup, jetResPhiType_);
  JME::JetResolutionScaleFactor resSFObj = JME::JetResolutionScaleFactor::get(iSetup, jetSFType_);

  edm::Handle<PFCollection> pfCandsH_;
  iEvent.getByToken( pfCollectionT_, pfCandsH_ );

  edm::Handle<reco::CandidateView> pfCandidates;
  iEvent.getByToken( pfCandidatesToken_, pfCandidates );

  std::vector< edm::Handle<reco::CandidateView> > leptons;
  for ( std::vector<edm::EDGetTokenT<edm::View<reco::Candidate> > >::const_iterator srcLeptons_i = lepTokens_.begin(); srcLeptons_i != lepTokens_.end(); ++srcLeptons_i ) {
    edm::Handle<reco::CandidateView> leptons_i;
    iEvent.getByToken(*srcLeptons_i, leptons_i);
    leptons.push_back( leptons_i );
  }

  edm::Handle<double> rho;
  iEvent.getByToken(rhoLabel_, rho);
  
  vJetIdxs.clear();
  v_att_genHiggs_M    = -1;
  v_att_genTau1_pT    = -1;
  v_att_genTau2_pT    = -1;
  v_att_genTau3_pT    = -1;
  v_att_genTau4_pT    = -1;
  v_att_dRtautau      = -1;
  v_att_pfMET         = -1;
  v_att_dphillmet     = -1;
  v_att_dphill        = -1;
  v_att_Mvis          = -1;
  v_att_pTvis         = -1;
  v_att_Mtautau       = -1;
  v_att_Mtautau_svFit = -1;
  v_att_Mth_svFit     = -1;
  v_att_pTh_svFit     = -1;
  v_att_pTh           = -1;
  v_att_mth           = -1;
  v_att_mcoll         = -1;
  
  v_att_gentau_Idxs.clear();
  v_att_tau_Idxs.clear();
  v_att_tau_combs.clear();
  v_att_jetIsSignal.clear();
  v_att_jetdR.clear();
  v_att_tau_pT.clear();
  v_att_tau_mva.clear();
  v_att_tau_dm.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_att_tau_jetFakePhoIdxs;
  */

  bool IsSignal             = true;
  //bool IsMC                 = false;
  //if (iEvent.isRealData()) IsMC = false;
  //TAU SELECTION
  float tau_sel_mvaID       = -2;
  float tau_sel_pT          = 20;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;

  //Lookin at RecoTaus
  if ( taus->size() < 2 ) {
    if (debug ) std::cout << "   !!!!!!!!!!  ONLY " << taus->size() << " TAUS IN THIS EVENT  !!!!!!!!!!"<< std::endl;
    return false;
  }
  if ( debug ) std::cout << " TAUS IN THE EVENT = " << taus->size() << std::endl;
  unsigned best_tau_1 = 99;
  unsigned best_tau_2 = 99;
  float best_tau_pt_1 = -1;
  float best_tau_pt_2 = -1;
  float best_tau_mva_1 = -1;
  float best_tau_mva_2 = -1;
  float genHiggs_mass = -1;
  float genTau1_pT = -1;
  float genTau2_pT = -1;
  float genTau3_pT = -1;
  float genTau4_pT = -1;
  if ( isMC_ ) {
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      if ( IsSignal ) {
        if ( abs(iGen->pdgId()) != 35 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 25 || iGen->daughter(1)->pdgId() != 25 ) continue;
        if ( abs(iGen->daughter(0)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(0)->daughter(1)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->daughter(1)->pdgId()) != 15 ) continue;
        TLorentzVector GenTau1  = SetTau(iGen->daughter(0)->daughter(0)->pt(), iGen->daughter(0)->daughter(0)->eta(), iGen->daughter(0)->daughter(0)->phi(), iGen->daughter(0)->daughter(0)->mass());
        TLorentzVector GenTau2  = SetTau(iGen->daughter(0)->daughter(1)->pt(), iGen->daughter(0)->daughter(1)->eta(), iGen->daughter(0)->daughter(1)->phi(), iGen->daughter(0)->daughter(1)->mass());
        TLorentzVector GenTau3  = SetTau(iGen->daughter(1)->daughter(0)->pt(), iGen->daughter(1)->daughter(0)->eta(), iGen->daughter(1)->daughter(0)->phi(), iGen->daughter(1)->daughter(0)->mass());
        TLorentzVector GenTau4  = SetTau(iGen->daughter(1)->daughter(1)->pt(), iGen->daughter(1)->daughter(1)->eta(), iGen->daughter(1)->daughter(1)->phi(), iGen->daughter(1)->daughter(1)->mass());
        TLorentzVector GenHiggs = GenTau1 + GenTau2 + GenTau3 + GenTau4;
        genHiggs_mass = GenHiggs.M();
        if (GenTau1.Pt() > GenTau2.Pt() ) {
          genTau1_pT = GenTau1.Pt();
          genTau2_pT = GenTau2.Pt();
        }else {
          genTau1_pT = GenTau2.Pt();
          genTau2_pT = GenTau1.Pt();
        }
        if (GenTau3.Pt() > GenTau4.Pt() ) {
          genTau3_pT = GenTau3.Pt();
          genTau4_pT = GenTau4.Pt();
        } else {
          genTau3_pT = GenTau4.Pt();
          genTau4_pT = GenTau3.Pt();
        }
      } else {
        if ( abs(iGen->pdgId()) != 25 || iGen->numberOfDaughters() != 2 || iGen->daughter(0)->pdgId() != 15 || iGen->daughter(1)->pdgId() != 15 ) continue;
        TLorentzVector GenTau1  = SetTau(iGen->daughter(0)->pt(), iGen->daughter(0)->eta(), iGen->daughter(0)->phi(), iGen->daughter(0)->mass());
        TLorentzVector GenTau2  = SetTau(iGen->daughter(1)->pt(), iGen->daughter(1)->eta(), iGen->daughter(1)->phi(), iGen->daughter(1)->mass());
        TLorentzVector GenHiggs = GenTau1 + GenTau2;
        genHiggs_mass = GenHiggs.M();
        if (GenTau1.Pt() > GenTau2.Pt() ) {
          genTau1_pT = GenTau1.Pt();
          genTau2_pT = GenTau1.Pt();
        } else {
          genTau1_pT = GenTau2.Pt();
          genTau2_pT = GenTau1.Pt();
        }
      }
    }
  }
  for ( unsigned iT1(0); iT1 != taus->size(); ++iT1 ) {
    reco::PFTauRef iTau1( taus, iT1 );
    if ( iTau1->pt() < tau_sel_pT ) continue;
    //if (!((*MuonRejection)[iTau1])) continue;
    //if (!((*ElectronRejectionMVA6)[iTau1])) continue;
    if ( (*MVAIsolation)[iTau1] < tau_sel_mvaID ) continue;  
    if ( debug ) std::cout << " TAU PASSED SELECTION "<< std::endl;
    if ( isMC_ ) {
      bool skip_tau = true;
      edm::Handle<reco::GenParticleCollection> genParticles;
      iEvent.getByToken( genParticleCollectionT_, genParticles );
      for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
        if ( abs(iGen->pdgId()) != 15 ) continue;
        float gentaudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iGen->eta(), iGen->phi() );
        if ( iGen->mother()->pdgId() != 25 || iGen->mother()->numberOfDaughters() != 2 || abs(iGen->mother()->daughter(0)->pdgId()) != 15 || abs(iGen->mother()->daughter(1)->pdgId()) != 15 ) continue;
        float gendRtautau = reco::deltaR( iGen->mother()->daughter(0)->eta(), iGen->mother()->daughter(0)->phi(), iGen->mother()->daughter(1)->eta(), iGen->mother()->daughter(1)->phi() );
        if ( gentaudR < 0.4 && gendRtautau < 0.4 ){
          skip_tau = false;
          break;
        }
      }
      if ( skip_tau ) {
        if ( debug ) std::cout << " RECO TAU DO NOT MATCH A GEN TAU" << std::endl;
        continue;
      }
    } 

    unsigned int tau_combinations = 0;
    for ( unsigned iT2(0); iT2 != taus->size(); ++iT2 ) {
      if ( iT2 == iT1 ) continue;
      reco::PFTauRef iTau2( taus, iT2 );
      if ( iTau2->pt() < tau_sel_pT ) continue;
      //if (!((*MuonRejection)[iTau2])) continue;
      //if (!((*ElectronRejectionMVA6)[iTau2])) continue;
      if ( (*MVAIsolation)[iTau2] < tau_sel_mvaID ) continue;  
      if ( isMC_ ) {
        bool skip_tau = true;
        edm::Handle<reco::GenParticleCollection> genParticles;
        iEvent.getByToken( genParticleCollectionT_, genParticles );
        for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
          //if ( abs(iGen->pdgId()) != 15 ) continue;
          if ( abs(iGen->pdgId()) != 15 || iGen->numberOfMothers() < 1 || iGen->mother()->pdgId() != 25) continue;
          float gentaudR = reco::deltaR( iTau2->eta(),iTau2->phi(), iGen->eta(), iGen->phi() );
          if ( iGen->mother()->pdgId() != 25 || iGen->mother()->numberOfDaughters() != 2 || abs(iGen->mother()->daughter(0)->pdgId()) != 15 || abs(iGen->mother()->daughter(1)->pdgId()) != 15 ) continue;
          float gendRtautau = reco::deltaR( iGen->mother()->daughter(0)->eta(), iGen->mother()->daughter(0)->phi(), iGen->mother()->daughter(1)->eta(), iGen->mother()->daughter(1)->phi() );
          if ( gentaudR < 0.4 && gendRtautau < 0.4 ){
            skip_tau = false;
            break;
          }
        }
        if ( skip_tau ) {
          if ( debug ) std::cout << " RECO TAU DO NOT MATCH A GEN TAU" << std::endl;
          continue;
        }
      }
      // tau dR selection
      float recotaudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iTau2->eta(),iTau2->phi() );
      if ( recotaudR < 0.5 ) continue;
      if ( debug ) std::cout << " RECO TAU dR = " << recotaudR << std::endl;
      ++tau_combinations;
    
     /* 
     // PAIR SELECTION BASED ON pT
     if ( ( iT2 != best_tau_1 && iT1 != best_tau_2 ) && ( best_tau_pt_1 < iTau1->pt() || best_tau_pt_2 < iTau2->pt() ) ) {
        if ( iTau1->pt() > iTau2->pt() ) {
          best_tau_2    = iT2;
          best_tau_pt_2 = iTau2->pt();
          best_tau_1    = iT1;
          best_tau_pt_1 = iTau2->pt();
        } else {
          best_tau_2    = iT1;
          best_tau_pt_2 = iTau1->pt();
          best_tau_1    = iT2;
          best_tau_pt_1 = iTau2->pt();
        }
      }*/
      // PAIR SELECTION BASED ON MVA ID
      if ( (iT2 != best_tau_1 && iT1 != best_tau_2) && ( best_tau_mva_1 < (*MVAIsolation)[iTau1] || best_tau_mva_2 < (*MVAIsolation)[iTau2] ) ) {
        if ( (*MVAIsolation)[iTau1] > (*MVAIsolation)[iTau2] ) {
          best_tau_2     = iT2;
          best_tau_mva_2 = (*MVAIsolation)[iTau2];
          best_tau_1     = iT1;
          best_tau_mva_1 = (*MVAIsolation)[iTau1];
        } else {
          best_tau_2     = iT1;
          best_tau_mva_2 = (*MVAIsolation)[iTau1];
          best_tau_1     = iT2;
          best_tau_mva_1 = (*MVAIsolation)[iTau2];
        }
      }
    }
    if ( tau_combinations == 0 ) continue;
    v_att_tau_Idxs.push_back( iT1 );
    v_att_tau_combs.push_back( tau_combinations );
  } //end iT1
  if ( v_att_tau_Idxs.size() < 2 ) return false;
  if ( debug ) std::cout << " SELECTED TAUS = " << v_att_tau_Idxs.size() << std::endl;

  // Pair selection
  reco::PFTauRef iTau1( taus, best_tau_1 );
  reco::PFTauRef iTau2( taus, best_tau_2 );
  float taudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iTau2->eta(),iTau2->phi() );
  TLorentzVector Tau1  = SetTau(iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()); 
  TLorentzVector Tau2  = SetTau(iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass());
  TLorentzVector DiTau = Tau1+Tau2;
  float dphill  = abs(Tau1.DeltaPhi(Tau2));  
  float diMvis  = DiTau.M();
  float dipTvis = DiTau.Pt();
  if ( debug ) std::cout << " Tau pair " << best_tau_1 << " ("<< iTau1->pt() << " GeV) + " << best_tau_2 << " (" << iTau2->pt() << " GeV) | ditau Mvis : " << diMvis << " GeV | dR = " << taudR << std::endl;

  
  float pfMET    = (pfmet->front()).et();
  float pfMETphi = (pfmet->front()).phi();
  TLorentzVector MET = SetMET(pfMET,pfMETphi);
  if ( debug ) std::cout << " PF MET = " << pfMET << std::endl;
  // define MET for SVFit
  double measuredMETx = (pfmet->front()).px();
  double measuredMETy = (pfmet->front()).py();
  if ( debug ) std::cout << " PF MET X = " << measuredMETx << " | PF MET Y = " << measuredMETy << std::endl;
  float dphillmet = abs(DiTau.DeltaPhi(MET));  
  float ditau_mth = sqrt( 2. * dipTvis * pfMET * ( 1. - cos (dphillmet) ));

  TLorentzVector Higgs = Tau1+Tau2+MET;
  float ditau_M   = Higgs.M();
  float ditau_pT  = Higgs.Pt();

  float the1    = 2.*atan(exp(-iTau1->eta()));
  float the2    = 2.*atan(exp(-iTau2->eta()));
  float pTmiss1 = ( (sin(iTau1->phi())*measuredMETx) - (cos(iTau1->phi())*measuredMETy) ) / (sin(the2)*( (sin(iTau1->phi())*cos(iTau2->phi())) - (cos(iTau1->phi())*sin(iTau2->phi())) ) );
  float pTmiss2 = ( (sin(iTau2->phi())*measuredMETx) - (cos(iTau2->phi())*measuredMETy) ) / (sin(the1)*( (sin(iTau2->phi())*cos(iTau1->phi())) - (cos(iTau2->phi())*sin(iTau1->phi())) ) );
  float mcoll = ditau_M / sqrt( (pTmiss1/(pTmiss1+iTau1->pt())) + (pTmiss2/(pTmiss2+iTau2->pt())) );

  // define MET covariance
  double sumPtUnclestered = 0;
  const reco::METCovMatrix cov = metSigAlgo_->getCovariance( *Jets, leptons, pfCandidates, *rho, resPtObj, resPhiObj, resSFObj, iEvent.isRealData(), sumPtUnclestered);
  TMatrixD covMET(2, 2);
  covMET[0][0] = cov[0][0]; 
  covMET[1][1] = cov[1][1]; 
  covMET[0][1] = cov[0][1]; 
  covMET[1][0] = cov[1][0]; 
  if ( debug ) std::cout << " MET COV xx = " << covMET[0][0] << std::endl;
  if ( debug ) std::cout << " MET COV yy = " << covMET[1][1] << std::endl;
  if ( debug ) std::cout << " MET COV xy = " << covMET[0][1] << std::endl;
  if ( debug ) std::cout << " MET COV yx = " << covMET[1][0] << std::endl;

  std::vector<MeasuredTauLepton> measuredTauLeptons;
  if ( iTau1->pt() > iTau2->pt() ) {
    measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()) );
    measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()) );
  }
  else {
    measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()) );
    measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()) );
  }

  FastMTT aFastMTTAlgo;
  aFastMTTAlgo.run(measuredTauLeptons,measuredMETx,measuredMETy,covMET);
  LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
  float svFitMass = ttP4.M(); // return value is in units of GeV
  float svFitMt = ttP4.Mt();
  float svFitPt = ttP4.Pt();
  //float svFitEta = ttP4.Eta();
  //float svFitPhi = ttP4.Phi();
  std::cout << " Found mass = " << svFitMass << std::endl;

  /*
  double higgsMass_1stRun = -1;
  double higgsMth_1stRun = -1;
  double higgsMass_2ndRun = -1;
  double higgsMth_2ndRun = -1;
 
  int verbosity_svFit = 0;
  //if ( debug ) verbosity_svFit = 1;
  ClassicSVfit svFitAlgo(verbosity_svFit);
  #ifdef USE_SVFITTF
  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  //svFitAlgo.setHadTauTF(hadTauTF);
  //svFitAlgo.enableHadTauTF();
  #endif
  svFitAlgo.addLogM_fixed(true, 5.);
  //svFitAlgo.setHistogramAdapter(new classic_svFit::DiTauSystemHistogramAdapter());
  svFitAlgo.setDiTauMassConstraint(-1.0);
  //svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  //svFitAlgo.addLogM_dynamic(true, "(m/1000.)*15.");
  //svFitAlgo.setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
  //svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_1stRun = svFitAlgo.isValidSolution();
  if ( isValidSolution_1stRun ) {
    higgsMass_1stRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMassErr();
    higgsMth_1stRun  = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
    if (debug) std::cout << "found valid solution: mass = " << higgsMass_1stRun << " , transverse mass = " << higgsMth_1stRun <<std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
   
  // re-run with mass constraint
  double massContraint = 125.06;
  std::cout << "Testing integration with ditau mass constraint set to " << massContraint << std::endl;
  svFitAlgo.setLikelihoodFileName("testClassicSVfit_withMassContraint.root");
  svFitAlgo.setDiTauMassConstraint(massContraint);
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_2ndRun = svFitAlgo.isValidSolution();
  if ( isValidSolution_2ndRun ) {
    higgsMass_2ndRun = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getMassErr();
    higgsMth_2ndRun  = static_cast<HistogramAdapterDiTau*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
    if (debug) std::cout << "found 2nd valid solution: mass = " << higgsMass_2ndRun << " , transverse mass = " << higgsMth_2ndRun <<std::endl;
  } else {
    higgsMass_2ndRun = ditau_M;
    higgsMth_2ndRun  = ditau_mth;
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  */
  
  //Signal Region selection
  if ( svFitMass > 50 && svFitMass < 130 && ditau_mth < 50 ) IsSignal = true;
  //if (!IsSignal) return false;
  
  // Study of trigger bit
  //hltConfig_.init(iEvent,iSetup)
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByToken(triggerResultsToken_,hltresults);
  if (!hltresults.isValid()) {
    std::cout << "!!! Error in getting TriggerResults product from Event !!!" << std::endl;
    return false;
  }
  int ntrigs = hltresults->size();
  edm::TriggerNames const& triggerNames = iEvent.triggerNames(*hltresults);
  for (int itrig = 0; itrig != ntrigs; ++itrig) {
    TString trigName = triggerNames.triggerName(itrig);
    bool accept = hltresults->accept(itrig);
    if (!(accept)) continue;
    std::cout << " Accept " << trigName << std::endl;
  }

  float tau_pT  = -99;
  float tau_mva = -99;
  float tau_dm  = -99;
  unsigned int nMatchedJets = 0;
  // Loop over jets
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    //if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] -> Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;

    float jetdR  = -99;
    float taudR1 = reco::deltaR( iJet->eta(),iJet->phi(), iTau1->eta(),iTau1->phi() );
    float taudR2 = reco::deltaR( iJet->eta(),iJet->phi(), iTau2->eta(),iTau2->phi() );
    if ( taudR1 < taudR2 ) {
      jetdR   = taudR1;
      tau_pT  = iTau1->pt();
      tau_mva = ((*MVAIsolation)[iTau1]);
      tau_dm  = ((*DecayMode)[iTau1]);
    } else {
      jetdR = taudR2; 
      tau_pT = iTau2->pt();
      tau_mva = ((*MVAIsolation)[iTau2]);
      tau_dm  = ((*DecayMode)[iTau2]);
    }
    if ( jetdR > 0.4 ) continue;

    //filling jet variables
    if ( debug ) std::cout << " JET [" << iJ << "] MATCHED A TAU CANDIDATE: dR = " << jetdR << std::endl;
    ++nMatchedJets;
    vJetIdxs.push_back( iJ );
    v_att_jetdR.push_back( jetdR );
    v_att_jetIsSignal.push_back( IsSignal );
    v_att_tau_pT.push_back( tau_pT );
    v_att_tau_mva.push_back( tau_mva );
    v_att_tau_dm.push_back( tau_dm );
  } // end jet loop
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  v_att_genHiggs_M    = genHiggs_mass;
  v_att_genTau1_pT    = genTau1_pT;
  v_att_genTau2_pT    = genTau2_pT;
  v_att_genTau3_pT    = genTau3_pT;
  v_att_genTau4_pT    = genTau4_pT;
  v_att_dRtautau      = taudR;
  v_att_pfMET         = pfMET;
  v_att_dphillmet     = dphillmet;
  v_att_dphill        = dphill;
  v_att_Mvis          = diMvis;
  v_att_pTvis         = dipTvis;
  v_att_Mtautau       = ditau_M;
  v_att_Mtautau_svFit = svFitMass;
  v_att_pTh           = ditau_pT;
  v_att_pTh_svFit     = svFitPt;
  v_att_mth           = ditau_mth;
  v_att_mcoll         = mcoll;
  v_att_Mth_svFit     = svFitMt;

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
  v_att_tau_jetIsSignal_.clear();
  v_att_tau_jetdR_.clear();
  v_att_dRtautau_.clear();
  v_att_Mvis_.clear();
  v_att_pTvis_.clear();
  v_att_Mtautau_.clear();
  v_att_Mtautau_svFit_.clear();
  v_att_Mth_svFit_.clear();
  v_att_pTh_svFit_.clear();
  v_att_pTh_.clear();
  v_att_tau_pT_.clear();
  v_att_tau_mva_.clear();
  v_att_tau_dm_.clear();
  v_att_mth_.clear();
  v_att_mcoll_.clear();
  v_att_pfMET_.clear();
  v_att_dphillmet_.clear();
  v_att_dphill_.clear();
  v_att_genHiggs_M_.clear();
  v_att_genTau1_pT_.clear();
  v_att_genTau2_pT_.clear();
  
  h_tau_att_genHiggs_M->Fill( v_att_genHiggs_M );
  h_tau_att_genTau1_pT->Fill( v_att_genTau1_pT );
  h_tau_att_genTau2_pT->Fill( v_att_genTau2_pT );
  h_tau_att_genTau1_pT->Fill( v_att_genTau3_pT );
  h_tau_att_genTau2_pT->Fill( v_att_genTau4_pT );
  h_tau_att_dRtautau->Fill( v_att_dRtautau );
  h_tau_att_Mvis->Fill( v_att_Mvis );
  h_tau_att_pTvis->Fill( v_att_pTvis );
  h_tau_att_Mtautau->Fill( v_att_Mtautau );
  h_tau_att_Mtautau_svFit->Fill( v_att_Mtautau_svFit );
  h_tau_att_Mth_svFit->Fill( v_att_Mth_svFit );
  h_tau_att_pTh_svFit->Fill( v_att_pTh_svFit );
  h_tau_att_pTh->Fill( v_att_pTh );
  h_tau_att_mth->Fill( v_att_mth );
  h_tau_att_mcoll->Fill( v_att_mcoll );
  h_tau_att_pfMET->Fill( v_att_pfMET );
  h_tau_att_dphillmet->Fill( v_att_dphillmet );
  h_tau_att_dphill->Fill( v_att_dphill );

  v_att_genHiggs_M_.push_back( v_att_genHiggs_M );
  v_att_genTau1_pT_.push_back( v_att_genTau1_pT );
  v_att_genTau2_pT_.push_back( v_att_genTau2_pT );
  v_att_genTau1_pT_.push_back( v_att_genTau3_pT );
  v_att_genTau2_pT_.push_back( v_att_genTau4_pT );
  v_att_dRtautau_.push_back( v_att_dRtautau );
  v_att_Mvis_.push_back( v_att_Mvis );
  v_att_pTvis_.push_back( v_att_pTvis );
  v_att_Mtautau_.push_back( v_att_Mtautau );
  v_att_Mtautau_svFit_.push_back( v_att_Mtautau_svFit );
  v_att_Mth_svFit_.push_back( v_att_Mth_svFit );
  v_att_pTh_svFit_.push_back( v_att_pTh_svFit );
  v_att_pTh_.push_back( v_att_pTh );
  v_att_mth_.push_back( v_att_mth );
  v_att_mcoll_.push_back( v_att_mcoll );
  v_att_pfMET_.push_back( v_att_pfMET );
  v_att_dphillmet_.push_back( v_att_dphillmet );
  v_att_dphill_.push_back( v_att_dphill );


  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_att_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_att_jet_eta->Fill( iJet->eta() );
    h_tau_att_jet_E->Fill( iJet->energy() );
    h_tau_att_jet_m0->Fill( iJet->mass() );
    h_tau_att_jet_isSignal->Fill( v_att_jetIsSignal[iJ] );
    h_tau_att_jet_dR->Fill( v_att_jetdR[iJ] );
    h_tau_att_tau_pT->Fill( v_att_tau_pT[iJ] );
    h_tau_att_tau_mva->Fill( v_att_tau_mva[iJ] );
    h_tau_att_tau_dm->Fill( v_att_tau_dm[iJ] );

    // Fill branches 
    v_att_tau_jet_pt_.push_back( iJet->pt() );
    v_att_tau_jet_m0_.push_back( iJet->mass() );
    v_att_tau_jetIsSignal_.push_back( v_att_jetIsSignal[iJ] );
    v_att_tau_jetdR_.push_back( v_att_jetdR[iJ] );
    v_att_tau_pT_.push_back( v_att_tau_pT[iJ] );
    v_att_tau_mva_.push_back( v_att_tau_mva[iJ] );
    v_att_tau_dm_.push_back( v_att_tau_dm[iJ] );

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
