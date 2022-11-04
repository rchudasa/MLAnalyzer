#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 4; //TODO: use cfg level nJets_
TH1D *h_ele_zee_mr_jet_pT;
TH1D *h_ele_zee_mr_jet_E;
TH1D *h_ele_zee_mr_jet_eta;
TH1D *h_ele_zee_mr_jet_m0;
TH1D *h_ele_zee_mr_jet_ma;
TH1D *h_ele_zee_mr_jet_pta;
TH2D *h_ele_zee_mr_jet_a_m_pt;
TH1D *h_ele_zee_mr_jet_nJet;
TH1D *h_ele_zee_mr_jet_isDiEle;
TH1D *h_ele_zee_mr_jet_dR;
TH1D *h_ele_zee_mr_jet_EledR;
TH1D *h_ele_zee_mr_jet_Ele1dR;
TH1D *h_ele_zee_mr_jet_Ele2dR;
TH1D *h_ele_zee_mr_jet_Ele1pT;
TH1D *h_ele_zee_mr_jet_Ele2pT;
TH1D *h_ele_zee_mr_jet_NrecoEle;
TH1D *h_ele_zee_mr_jet_NGenEle;
TH1D *h_ele_zee_mr_jet_recoEle1dR;
TH1D *h_ele_zee_mr_jet_recoEle2dR;
vector<float> v_zee_mr_jetIsDiEle;
vector<float> v_zee_mr_aJetdR;
vector<float> v_zee_mr_ma;
vector<float> v_zee_mr_pta;
vector<float> v_zee_mr_jetEledR;
vector<float> v_zee_mr_jetEle1dR;
vector<float> v_zee_mr_jetEle2dR;
vector<float> v_zee_mr_jetEle1pT;
vector<float> v_zee_mr_jetEle2pT;
vector<float> v_zee_mr_jetNGenEle;
vector<float> v_zee_mr_jetNrecoEle;
vector<float> v_zee_mr_jetrecoEle1dR;
vector<float> v_zee_mr_jetrecoEle2dR;

vector<float> v_zee_mr_ele_jet_m0_;
vector<float> v_zee_mr_ele_jet_ma_;
vector<float> v_zee_mr_ele_jet_pta_;
vector<float> v_zee_mr_ele_jet_pt_;
vector<float> v_zee_mr_ele_jetPdgIds_;
vector<float> v_zee_mr_ele_jetIsDiEle_;
vector<float> v_zee_mr_ele_aJetdR_;
vector<float> v_zee_mr_ele_jetEledR_;
vector<float> v_zee_mr_ele_jetEle1dR_;
vector<float> v_zee_mr_ele_jetEle2dR_;
vector<float> v_zee_mr_ele_jetEle1pT_;
vector<float> v_zee_mr_ele_jetEle2pT_;
vector<float> v_zee_mr_ele_jetNGenEle_;
vector<float> v_zee_mr_ele_jetNrecoEle_;
vector<float> v_zee_mr_ele_jetrecoEle1dR_;
vector<float> v_zee_mr_ele_jetrecoEle2dR_;

vector<float> v_zee_mr_ele_subJetE_[nJets];
vector<float> v_zee_mr_ele_subJetPx_[nJets];
vector<float> v_zee_mr_ele_subJetPy_[nJets];
vector<float> v_zee_mr_ele_subJetPz_[nJets];



// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_DY_ele_classification ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_ele_zee_mr_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                      , 100,  0., 800.);
  h_ele_zee_mr_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                              , 100,  0., 800.);
  h_ele_zee_mr_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                                , 100, -5.,   5.);
  h_ele_zee_mr_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                              ,  10,  0.,  10.);
  h_ele_zee_mr_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                          , 100,  0., 100.);
  h_ele_zee_mr_jet_a_m_pt     = fs->make<TH2D>("h_a_m_pT"         , "m^{a} vs p_{T}^{a};m^{a} vs p_{T}^{a};Jets"    ,  59, 0.1,  6., 26, 20., 150.);
  h_ele_zee_mr_jet_ma         = fs->make<TH1D>("h_jet_ma"         , "m^{a};m^{a};Jets"                              ,  59, 0.1,  6.);
  h_ele_zee_mr_jet_pta        = fs->make<TH1D>("h_jet_pta"        , "p_{T}^{a};p_{T}^{a};Jets"                      ,  26, 20., 150.);
  h_ele_zee_mr_jet_isDiEle    = fs->make<TH1D>("h_jet_isDiEle"    , "nIsDiEle;nIsDiEle;Jets"                        ,  10,  0.,  10.);
  h_ele_zee_mr_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{a,j};dR_{a,j};Jets"                        ,  50,  0.,  0.5);
  h_ele_zee_mr_jet_EledR      = fs->make<TH1D>("h_jet_EledR"      , "dR_{e,e};dR_{e,e};Jets"                        ,  50,  0.,   1.);
  h_ele_zee_mr_jet_Ele1dR     = fs->make<TH1D>("h_jet_Ele1dR"     , "dR_{e_{1},j};dR_{e_{1},j};Jets"                ,  50,  0.,  0.5);
  h_ele_zee_mr_jet_Ele2dR     = fs->make<TH1D>("h_jet_Ele2dR"     , "dR_{e_{2},j};dR_{e_{2},j};Jets"                ,  50,  0.,  0.5);
  h_ele_zee_mr_jet_Ele1pT     = fs->make<TH1D>("h_jet_Ele1pT"     , "p_{T}^{e_{1}};p_{T}^{e_{1}};Jets"              ,  50,  0.,  100);
  h_ele_zee_mr_jet_Ele2pT     = fs->make<TH1D>("h_jet_Ele2pT"     , "p_{T}^{e_{2}};p_{T}^{e_{2}};Jets"              ,  50,  0.,  100);
  h_ele_zee_mr_jet_NGenEle  = fs->make<TH1D>("h_jet_NGenEle"      , "# e^{RECO};# e^{RECO};Jets"                    ,   5,  0.,   5.);
  h_ele_zee_mr_jet_NrecoEle  = fs->make<TH1D>("h_jet_NrecoEle"    , "# e^{RECO};# e^{RECO};Jets"                    ,   5,  0.,   5.);
  h_ele_zee_mr_jet_recoEle1dR = fs->make<TH1D>("h_jet_recoEle1dR" , "dR_{e_{1}^{RECO},j};dR_{e_{1}^{RECO},j};Jets"  ,  50,  0.,  0.5);
  h_ele_zee_mr_jet_recoEle2dR = fs->make<TH1D>("h_jet_recoEle2dR" , "dR_{e_{2}^{RECO},j};dR_{e_{2}^{RECO},j};Jets"  ,  25,  0.,  0.5);

  tree->Branch("jetM",       &v_zee_mr_ele_jet_m0_);
  tree->Branch("jetPt",      &v_zee_mr_ele_jet_pt_);
  tree->Branch("jetPdgIds",  &v_zee_mr_ele_jetPdgIds_);
  tree->Branch("aJetdR",     &v_zee_mr_ele_aJetdR_);
  tree->Branch("jetIsDiEle", &v_zee_mr_ele_jetIsDiEle_);
  tree->Branch("a_m",        &v_zee_mr_ele_jet_ma_);
  tree->Branch("a_pt",       &v_zee_mr_ele_jet_pta_);
  tree->Branch("jetpT",      &v_zee_mr_ele_jet_pt_);
  tree->Branch("EledR",      &v_zee_mr_ele_jetEledR_);
  tree->Branch("Ele1dR",     &v_zee_mr_ele_jetEle1dR_);
  tree->Branch("Ele2dR",     &v_zee_mr_ele_jetEle2dR_);
  tree->Branch("Ele1pT",     &v_zee_mr_ele_jetEle1pT_);
  tree->Branch("Ele2pT",     &v_zee_mr_ele_jetEle2pT_);
  tree->Branch("NGenEle",   &v_zee_mr_ele_jetNGenEle_);
  tree->Branch("NrecoEle",  &v_zee_mr_ele_jetNrecoEle_);
  tree->Branch("recoEle1dR", &v_zee_mr_ele_jetrecoEle1dR_);
  tree->Branch("recoEle2dR", &v_zee_mr_ele_jetrecoEle2dR_);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_zee_mr_ele_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_zee_mr_ele_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_zee_mr_ele_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_zee_mr_ele_subJetPz_[iJ]);
  }

} // branchesEvtSel_jet_DY_ele_classification()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_DY_ele_classification( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::GsfElectronCollection> ele;    //TODO
  iEvent.getByToken(eleCollectionT_, ele);         //TODO

  vJetIdxs.clear();
  v_zee_mr_ele_jetPdgIds_.clear();
  v_zee_mr_jetIsDiEle.clear();
  v_zee_mr_aJetdR.clear();
  v_zee_mr_ma.clear();
  v_zee_mr_pta.clear();
  v_zee_mr_jetEledR.clear();
  v_zee_mr_jetEle1dR.clear();
  v_zee_mr_jetEle2dR.clear();
  v_zee_mr_jetEle1pT.clear();
  v_zee_mr_jetEle2pT.clear();
  v_zee_mr_jetNGenEle.clear();
  v_zee_mr_jetNrecoEle.clear();
  v_zee_mr_jetrecoEle1dR.clear();
  v_zee_mr_jetrecoEle2dR.clear();

  unsigned int nMatchedJets = 0;
  unsigned int nMatchedRecoEle = 0;
  unsigned int aPdgId           = 0;
  bool MatchedPseudoScalar = false;
  float a_mass = -99.;
  float a_pt   = -99.;
  float dRa    = -99.;
  float eledR =  99.;
  float ele1dR = -99.;
  float ele2dR = -99.;
  float ele1pT = -99.;
  float ele2pT = -99.;
  float recoele1dR = -99.;
  float recoele2dR = -99.;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  
    // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] => Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    unsigned int nMatchedGenParticles = 0;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0; 

  
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    if ( abs(iGen->pdgId()) != 11 ) continue;
    if ( iGen->numberOfMothers() != 1 ) continue;
    if ( iGen->mother()->pdgId() != 23)continue;
    ++iGenParticle;
    float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
    if ( debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR = "<< dR << std::endl ;
    if ( dR > 0.4 ) continue;
         
    if ( debug ) std::cout << "       MOTHER => status: " << iGen->mother()->status() << ", id: " << iGen->mother()->pdgId() << ", nDaught: " << iGen->mother()->numberOfDaughters() << " | pt: "<< iGen->mother()->pt() << " eta: " <<iGen->mother()->eta() << " phi: " <<iGen->mother()->phi() << " mass: " <<iGen->mother()->mass() << std::endl;
    
    if (nMatchedGenParticles == 0) {
      if ( debug ) std::cout << "       nMatchedGenParticles: " << nMatchedGenParticles << std::endl;
      aPdgId = std::abs(iGen->mother()->pdgId());
      if ( abs(iGen->mother()->pdgId()) == 23 ) {
	MatchedPseudoScalar = true;
      } else {
	MatchedPseudoScalar = false;
      }
      if (!MatchedPseudoScalar) {
	if ( debug ) std::cout << " MOTHER is not the pseudoscalar a"  << std::endl;
	continue;
      }
      a_mass = iGen->mother()->mass();
      a_pt   = iGen->mother()->pt();
      dRa = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->eta(),iGen->mother()->phi() );
      if (abs(iGen->mother()->daughter(0)->pdgId()) == 11 && abs(iGen->mother()->daughter(1)->pdgId()) == 11) { //TODO 
	eledR = reco::deltaR( iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
	if ( iGen->mother()->daughter(0)->pt() > iGen->mother()->daughter(1)->pt() ) {
	  ele1pT = iGen->mother()->daughter(0)->pt();
	  ele2pT = iGen->mother()->daughter(1)->pt();
	  ele1dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi() );
	  ele2dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
	} else {
	  ele1pT = iGen->mother()->daughter(1)->pt();
	  ele2pT = iGen->mother()->daughter(0)->pt();
	  ele1dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
	  ele2dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi() );
	} // end else pt2 > pt1
      }
    } //matched gen particle = 0
    if ( eledR > 0.4 ) {
        if ( debug ) std::cout << " !! Electrons are not merged: gen dR_ee = " << eledR << " !! " << std::endl;
        //continue;
      } 
     if ( !MatchedPseudoScalar ) continue;
     ++nMatchedGenParticles; 
  } // primary gen particles
  
   if ( nMatchedGenParticles > 0 ) passedGenSel = true;
    if (passedGenSel) { 
      ++nMatchedJets;
       
  
     std::cout<< "filling the branches" << std::endl;
     std::cout<< "mass:" << a_mass << std::endl;
 
     
      //Lookin at RecoEle
      if (debug ) std::cout << "  Looking at RECO electrons: "<< std::endl;
      if (ele->size() == 0) {
        if (debug ) std::cout << "   !!!!!!!!!!  NO RECO ELE IN THIS EVENT  !!!!!!!!!!"<< std::endl;
      }
      for ( unsigned iT(0); iT != ele->size(); ++iT ) {
        //edm::Ref<reco::GsfElectronCollection> iEle( ele, iT );
        reco::GsfElectronRef iEle( ele, iT );
        float recoeledR = reco::deltaR( iJet->eta(),iJet->phi(), iEle->eta(),iEle->phi() );
        if ( recoeledR < 0.4 && nMatchedRecoEle == 0 ) {
          if ( debug ) std::cout << "   Reco Ele [" << iT << "]  matched jet [" << iJ << "] : dR = " << recoeledR << std::endl;
          recoele1dR = recoeledR;
          ++nMatchedRecoEle;
        } else if ( recoeledR < 0.4 && nMatchedRecoEle == 1 ) {
          if ( debug ) std::cout << "   Reco Ele [" << iT << "]  matched jet [" << iJ << "] : dR = " << recoeledR << std::endl;
          if (recoeledR < recoele1dR) {
            recoele2dR = recoele1dR;
            recoele1dR = recoeledR;
          } else recoele2dR = recoeledR;
          ++nMatchedRecoEle;
        } else if ( debug && recoeledR < 0.4 && nMatchedRecoEle > 1 ) {
          std::cout << "   !!!!!!!!!!  FOUND MORE THAN 2 ELE INSIDE JET CONE OF 0.4 !!!!!!!!!!"<< std::endl;
          if (recoeledR < recoele2dR && recoeledR < recoele1dR) { 
            if (recoele1dR < recoele2dR) recoele2dR = recoele1dR;
            recoele1dR = recoeledR;
          } else if (recoeledR < recoele2dR && recoeledR > recoele1dR) recoele2dR = recoeledR; 
          ++nMatchedRecoEle;
        } else if ( debug ) {
          std::cout << "   !!!!!!!!!!  NO MATCH FOR Reco Ele [" << iT << "]  with jet [" << iJ << "] : dR = " << recoeledR << std::endl;
        }
      }

      vJetIdxs.push_back( iJ );
      v_zee_mr_ele_jetPdgIds_.push_back( aPdgId );
      v_zee_mr_aJetdR.push_back( dRa );
      v_zee_mr_jetEledR.push_back( eledR );
      v_zee_mr_ma.push_back( a_mass );
      v_zee_mr_pta.push_back( a_pt );
      v_zee_mr_jetEle1dR.push_back( ele1dR );
      v_zee_mr_jetEle2dR.push_back( ele2dR );
      v_zee_mr_jetEle1pT.push_back( ele1pT );
      v_zee_mr_jetEle2pT.push_back( ele2pT );
      v_zee_mr_jetNGenEle.push_back( nMatchedGenParticles );
      v_zee_mr_jetNrecoEle.push_back( nMatchedRecoEle );
      v_zee_mr_jetrecoEle1dR.push_back( recoele1dR );
      v_zee_mr_jetrecoEle2dR.push_back( recoele2dR );
      v_zee_mr_jetIsDiEle.push_back( MatchedPseudoScalar );

    }

  } // reco jets
  if ( debug ) std::cout << " Matched GEN particle-jet pairs " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << " >> has_jet_dijet_ele_classification: passed" << std::endl;
  return true;

} // runEvtSel_jet_DY_ele_classification()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_DY_ele_classification ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_ele_zee_mr_jet_nJet->Fill( vJetIdxs.size() );

  v_zee_mr_ele_jet_pt_.clear();
  v_zee_mr_ele_jet_m0_.clear();
  v_zee_mr_ele_jet_ma_.clear();
  v_zee_mr_ele_jet_pta_.clear();
  v_zee_mr_ele_jetIsDiEle_.clear();
  v_zee_mr_ele_aJetdR_.clear();
  v_zee_mr_ele_jetEledR_.clear();
  v_zee_mr_ele_jetEle1dR_.clear();
  v_zee_mr_ele_jetEle2dR_.clear();
  v_zee_mr_ele_jetEle1pT_.clear();
  v_zee_mr_ele_jetEle2pT_.clear();
  v_zee_mr_ele_jetNGenEle_.clear();
  v_zee_mr_ele_jetNrecoEle_.clear();
  v_zee_mr_ele_jetrecoEle1dR_.clear();
  v_zee_mr_ele_jetrecoEle2dR_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_ele_zee_mr_jet_pT->Fill( std::abs(iJet->pt()) );
    h_ele_zee_mr_jet_eta->Fill( iJet->eta() );
    h_ele_zee_mr_jet_E->Fill( iJet->energy() );
    h_ele_zee_mr_jet_m0->Fill( iJet->mass() );
    h_ele_zee_mr_jet_isDiEle->Fill( v_zee_mr_jetIsDiEle[iJ] );
    h_ele_zee_mr_jet_dR->Fill( v_zee_mr_aJetdR[iJ] );
    h_ele_zee_mr_jet_ma->Fill( v_zee_mr_ma[iJ] );
    h_ele_zee_mr_jet_pta->Fill( v_zee_mr_pta[iJ] );
    h_ele_zee_mr_jet_a_m_pt->Fill( v_zee_mr_ma[iJ], v_zee_mr_pta[iJ] );
    h_ele_zee_mr_jet_EledR->Fill( v_zee_mr_jetEledR[iJ] );
    h_ele_zee_mr_jet_Ele1dR->Fill( v_zee_mr_jetEle1dR[iJ] );
    h_ele_zee_mr_jet_Ele2dR->Fill( v_zee_mr_jetEle2dR[iJ] );
    h_ele_zee_mr_jet_Ele1pT->Fill( v_zee_mr_jetEle1pT[iJ] );
    h_ele_zee_mr_jet_Ele2pT->Fill( v_zee_mr_jetEle2pT[iJ] );
    h_ele_zee_mr_jet_NGenEle->Fill( v_zee_mr_jetNGenEle[iJ] );
    h_ele_zee_mr_jet_NrecoEle->Fill( v_zee_mr_jetNrecoEle[iJ] );
    h_ele_zee_mr_jet_recoEle1dR->Fill( v_zee_mr_jetrecoEle1dR[iJ] );
    h_ele_zee_mr_jet_recoEle2dR->Fill( v_zee_mr_jetrecoEle2dR[iJ] );

    // Fill branches 
    v_zee_mr_ele_jet_pt_.push_back( iJet->pt() );
    v_zee_mr_ele_jet_m0_.push_back( iJet->mass() );
    v_zee_mr_ele_jet_ma_.push_back( v_zee_mr_ma[iJ] );
    v_zee_mr_ele_jet_pta_.push_back( v_zee_mr_pta[iJ] );
    v_zee_mr_ele_jetIsDiEle_.push_back( v_zee_mr_jetIsDiEle[iJ] );
    v_zee_mr_ele_aJetdR_.push_back( v_zee_mr_aJetdR[iJ] );
    v_zee_mr_ele_jetEledR_.push_back( v_zee_mr_jetEledR[iJ] );
    v_zee_mr_ele_jetEle1dR_.push_back( v_zee_mr_jetEle1dR[iJ] );
    v_zee_mr_ele_jetEle2dR_.push_back( v_zee_mr_jetEle2dR[iJ] );
    v_zee_mr_ele_jetEle1pT_.push_back( v_zee_mr_jetEle1pT[iJ] );
    v_zee_mr_ele_jetEle2pT_.push_back( v_zee_mr_jetEle2pT[iJ] );
    v_zee_mr_ele_jetNGenEle_.push_back( v_zee_mr_jetNGenEle[iJ] );
    v_zee_mr_ele_jetNrecoEle_.push_back( v_zee_mr_jetNrecoEle[iJ] );
    v_zee_mr_ele_jetrecoEle1dR_.push_back( v_zee_mr_jetrecoEle1dR[iJ] );
    v_zee_mr_ele_jetrecoEle2dR_.push_back( v_zee_mr_jetrecoEle2dR[iJ] );

    // Gen jet constituents
    v_zee_mr_ele_subJetE_[iJ].clear();
    v_zee_mr_ele_subJetPx_[iJ].clear();
    v_zee_mr_ele_subJetPy_[iJ].clear();
    v_zee_mr_ele_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_zee_mr_ele_subJetE_[iJ].push_back( subJet->energy() );
      v_zee_mr_ele_subJetPx_[iJ].push_back( subJet->px() );
      v_zee_mr_ele_subJetPy_[iJ].push_back( subJet->py() );
      v_zee_mr_ele_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_DY_ele_classification()
