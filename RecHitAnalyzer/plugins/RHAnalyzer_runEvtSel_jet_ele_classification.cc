#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 4; //TODO: use cfg level nJets_
TH1D *h_ele_gen_pT;
TH1D *h_ele_jet_pT;
TH1D *h_ele_jet_E;
TH1D *h_ele_jet_eta;
TH1D *h_ele_jet_m0;
TH1D *h_ele_jet_nJet;
TH1D *h_ele_jet_isEle;
TH1D *h_ele_jet_dR;
TH1D *h_ele_goodvertices;

vector<float> v_jetIsEle;
vector<float> v_eledR;
vector<float> v_elegoodvertices;
vector<float> v_elepT;
vector<float> v_eleJetPt;
vector<float> v_eleJetEta;
vector<float> v_eleJetE;
vector<float> v_eleJetM0;
vector<float> v_eleJetPdgId;
vector<float> v_jetNrecoEle;
vector<float> v_jetrecoEle1dR;
vector<float> v_jetrecoEle2dR;

vector<float> v_ele_jet_m0_;
vector<float> v_ele_jet_pt_;
vector<float> v_ele_jet_eta_;
vector<float> v_ele_jet_e_;
vector<float> v_ele_gen_pt_;
vector<float> v_ele_jetPdgIds_;
vector<float> v_ele_jetIsEle_;
vector<float> v_ele_jetdR_;
vector<float> v_ele_goodvertices_;
vector<float> v_ele_jetNrecoEle_;
vector<float> v_ele_jetrecoEle1dR_;
vector<float> v_ele_jetrecoEle2dR_;

/*vector<float> v_ele_subJetE_[nJets];
vector<float> v_ele_subJetPx_[nJets];
vector<float> v_ele_subJetPy_[nJets];
vector<float> v_ele_subJetPz_[nJets];
*/

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_ele_classification ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_ele_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                      , 100,  0., 800.);
  h_ele_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                              , 100,  0., 800.);
  h_ele_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                                , 100, -5.,   5.);
  h_ele_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                              ,  10,  0.,  10.);
  h_ele_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                          , 100,  0., 100.);
  h_ele_jet_isEle      = fs->make<TH1D>("h_jet_isEle"      , "nIsEle;nIsEle;Jets"                            ,  10,  0.,  10.);
  h_ele_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{a,j};dR_{a,j};Jets"                        ,  50,  0.,  0.5);
  h_ele_goodvertices   = fs->make<TH1D>("h_goodvertices"   , "good vertices;good vertices;Jets"           ,  15,  0.,  75.);
  h_ele_gen_pT         = fs->make<TH1D>("h_gen_pT"         , "p_{T};p_{T}; Gen part"                      ,  30,  0., 300.);
  

  tree->Branch("jet_M",         &v_ele_jet_m0_);
  tree->Branch("jet_Pt",        &v_ele_jet_pt_);
  tree->Branch("jet_Eta",       &v_ele_jet_eta_);
  tree->Branch("jet_Energy",    &v_ele_jet_e_);
  tree->Branch("jet_PdgIds",    &v_ele_jetPdgIds_);
  tree->Branch("jet_IsEle",     &v_ele_jetIsEle_);
  tree->Branch("gen_pt",        &v_ele_gen_pt_);
  tree->Branch("jet_dR",        &v_ele_jetdR_);
  tree->Branch("goodvertices",  &v_ele_goodvertices_);
  tree->Branch("NrecoEle",  &v_ele_jetNrecoEle_);
  tree->Branch("recoEle1dR", &v_ele_jetrecoEle1dR_);
  tree->Branch("recoEle2dR", &v_ele_jetrecoEle2dR_);

 /* char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_ele_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_ele_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_ele_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_ele_subJetPz_[iJ]);
  }*/

} // branchesEvtSel_jet_ele_classification()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_ele_classification( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::GsfElectronCollection> ele;    //TODO
  iEvent.getByToken(eleCollectionT_, ele);         //TODO
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);
  
  vJetIdxs.clear();
  v_ele_jetPdgIds_.clear();
  v_jetIsEle.clear();
  v_eledR.clear();
  v_elegoodvertices.clear();
  v_elepT.clear();
  v_eleJetPt.clear();
  v_eleJetEta.clear();
  v_eleJetE.clear();
  v_eleJetM0.clear();

  v_eleJetPdgId.clear();
  v_jetNrecoEle.clear();
  v_jetrecoEle1dR.clear();
  v_jetrecoEle2dR.clear();

  unsigned int nMatchedJets = 0;
  unsigned int nMatchedRecoEle = 0;
  unsigned int goodVertices = 0;
  unsigned int PdgId        = 0;
  float jetdR               = -99.;
  float elepT               = -99.;
  bool JetIsEle             = false;
  float recoele1dR = -99.;
  float recoele2dR = -99.;

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
  

    
    reco::PFJetRef iJet( jets, iJ );
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    if (debug ) std::cout << "\t\t" << "Jet [" << iJ << "] => Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0; 
  
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
    if ( dR > 0.4 ) continue; 
  
    //if ( std::abs(iGen->pdgId()) == 12 || std::abs(iGen->pdgId()) == 14 || std::abs(iGen->pdgId()) == 16 ) continue; 
    if (isSignal_ && !(std::abs(iGen->pdgId()) == 11 && iGen->status() == 23 && iGen->numberOfMothers() == 1 ) ) continue;      // for drell yan to e+e-
    //if (isSignal_ && !(std::abs(iGen->pdgId()) == 11) ) continue;      // for drell yan to e+e-
    //if (!isSignal_ && !( iGen->status() == 23 ) ) continue;                         //for QCD background
 
    ++iGenParticle;

    if ( debug ) std::cout << "\t\t\t" <<"GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR = "<< dR << std::endl;         
    if ( debug ) std::cout << "\t\t\t" <<"MOTHER => status: " << iGen->mother()->status() << ", id: " << iGen->mother()->pdgId() << ", nDaught: " << iGen->mother()->numberOfDaughters() << " | pt: "<< iGen->mother()->pt() << " eta: " <<iGen->mother()->eta() << " phi: " <<iGen->mother()->phi() << " mass: " <<iGen->mother()->mass() << std::endl;
    

     if ( std::abs(iGen->pdgId()) == 11 ) {
        JetIsEle = true;
        PdgId = std::abs(iGen->pdgId());
        jetdR = dR;
        elepT = iGen->pt();
   
 	if (!isSignal_){  
        passedGenSel = false; // for backg
        break; 		      // for backg
	}

      } // if pDG ID ==11
	else if ( elepT < iGen->pt() ){
        PdgId = std::abs(iGen->pdgId());
        jetdR = dR;
        elepT = iGen->pt();
      }

      passedGenSel = true;
      ++iGenParticle;
  } // primary gen particles
  
    if (passedGenSel) { 
      ++nMatchedJets;
      
      
      //Lookin at RecoEle
      if (debug ) std::cout << "\t\t\t\t" << "  Looking at RECO electrons: "<< std::endl;
      if (ele->size() == 0) {
        if (debug ) std::cout << "\t\t\t\t" << "   !!!!!!!!!!  NO RECO ELE IN THIS EVENT  !!!!!!!!!!"<< std::endl;
      }
      for ( unsigned iT(0); iT != ele->size(); ++iT ) {
        //edm::Ref<reco::GsfElectronCollection> iEle( ele, iT );
        reco::GsfElectronRef iEle( ele, iT );
        float recoeledR = reco::deltaR( iJet->eta(),iJet->phi(), iEle->eta(),iEle->phi() );
        if ( recoeledR < 0.4 && nMatchedRecoEle == 0 ) {
          if ( debug ) std::cout << "\t\t\t\t" << "   Reco Ele [" << iT << "]  matched jet [" << iJ << "] : dR = " << recoeledR << " electron pT: " << iEle->pt() << "  eta:" << iEle->eta() << "  phi:"<< iEle->phi() <<  std::endl;
          recoele1dR = recoeledR;
          ++nMatchedRecoEle;
        } else if ( recoeledR < 0.4 && nMatchedRecoEle == 1 ) {
          if ( debug ) std::cout << "\t\t\t\t" << "   Reco Ele [" << iT << "]  matched jet [" << iJ << "] : dR = " << recoeledR << " electron pT: " << iEle->pt() << "  eta:" << iEle->eta() << "  phi:"<< iEle->phi() <<  std::endl;
          if (recoeledR < recoele1dR) {
            recoele2dR = recoele1dR;
            recoele1dR = recoeledR;
          } else recoele2dR = recoeledR;
          ++nMatchedRecoEle;
        } else if ( debug && recoeledR < 0.4 && nMatchedRecoEle > 1 ) {
          std::cout << "\t\t\t\t" << "   !!!!!!!!!!  FOUND MORE THAN 2 ELE INSIDE JET CONE OF 0.4 !!!!!!!!!!"<< std::endl;
          if (recoeledR < recoele2dR && recoeledR < recoele1dR) { 
            if (recoele1dR < recoele2dR) recoele2dR = recoele1dR;
            recoele1dR = recoeledR;
          } else if (recoeledR < recoele2dR && recoeledR > recoele1dR) recoele2dR = recoeledR; 
          ++nMatchedRecoEle;
        } else if ( debug ) {
          std::cout << "\t\t\t\t" << "   !!!!!!!!!!  NO MATCH FOR Reco Ele [" << iT << "]  with jet [" << iJ << "] : dR = " << recoeledR << std::endl;
        }
      }
      
      vJetIdxs.push_back( iJ );
      v_eleJetPdgId.push_back( PdgId );
      v_elepT.push_back( elepT );
      v_eleJetPt.push_back(iJet->pt());
      v_eleJetEta.push_back(iJet->eta());
      v_eleJetE.push_back(iJet->energy());
      v_eleJetM0.push_back(iJet->mass());
      
      std::cout << iJ << "First loop gen PDGID  :" << PdgId << " pt:" << elepT << std::endl;

      v_eledR.push_back( jetdR );
      v_elegoodvertices.push_back( goodVertices );
      v_jetIsEle.push_back( JetIsEle ); 
      v_jetNrecoEle.push_back( nMatchedRecoEle );
      v_jetrecoEle1dR.push_back( recoele1dR );
      v_jetrecoEle2dR.push_back( recoele2dR );

    } //passed gen level selection

  } // reco jets
  if ( debug ) std::cout << "\t" << " Matched GEN particle-jet pairs " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << "\t" << " >> has_jet_dijet_ele_classification: passed" << std::endl;
  return true;

} // runEvtSel_jet_ele_classification()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_ele_classification ( const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<int> passedJetIdxs, std::vector<int> failedJetIdxs ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  //h_ele_jet_nJet->Fill( vJetIdxs.size() );
  h_ele_jet_nJet->Fill( passedJetIdxs.size() );

  v_ele_jet_pt_.clear();
  v_ele_jet_eta_.clear();
  v_ele_jet_e_.clear();
  v_ele_jet_m0_.clear();
  v_ele_gen_pt_.clear();
  v_ele_jetIsEle_.clear();
  v_ele_jetdR_.clear();
  v_ele_goodvertices_.clear();
  v_ele_jetPdgIds_.clear();
  v_ele_jetNrecoEle_.clear();
  v_ele_jetrecoEle1dR_.clear();
  v_ele_jetrecoEle2dR_.clear();

  if(passedJetIdxs.size() == failedJetIdxs.size()){
   std::cout << "Passed Jet Index :" << passedJetIdxs[0] << "  Failed Jet Index :" << failedJetIdxs[0] << std::endl;
  }
  //std::cout << "************ filling up the jet branches after clearing vectors**** " << "passed jet indexes size :" << passedJetIdxs.size() << std::endl;
 
  if(passedJetIdxs.size()>0){
  if((passedJetIdxs.size())!=v_eleJetPdgId.size()) std::cout << " ++++++ there are failed jets due to HE cut ++++++++++++++++++"<< std::endl;
}
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    std::cout << iJ << "Second loop gen PDGID  :" << v_eleJetPdgId[iJ] << " pt:" << v_elepT[iJ] << std::endl;
    //std::cout << iJ << " ======================Jet eta is greater than 5 :" << v_eleJetEta[iJ] << "  " << iJet->eta() << std::endl<<std::endl;
    // Fill histograms 
    h_ele_jet_pT->Fill( iJet->pt());
    h_ele_jet_eta->Fill( iJet->eta());
    h_ele_jet_E->Fill( iJet->energy());
    h_ele_jet_m0->Fill( iJet->mass());
    h_ele_jet_isEle->Fill( v_jetIsEle[iJ] );
    h_ele_jet_dR->Fill( v_eledR[iJ] );
    h_ele_goodvertices->Fill( v_elegoodvertices[iJ] );
    h_ele_gen_pT->Fill( v_elepT[iJ] );

    // Fill branches 
    v_ele_jet_pt_.push_back( iJet->pt() );
    v_ele_jet_eta_.push_back( iJet->eta() );
    v_ele_jet_e_.push_back( iJet->energy() );
    v_ele_jet_m0_.push_back( iJet->mass() );
    v_ele_jetIsEle_.push_back( v_jetIsEle[iJ] );
    v_ele_jetdR_.push_back( v_eledR[iJ] );
    v_ele_goodvertices_.push_back( v_elegoodvertices[iJ] );
    v_ele_gen_pt_.push_back( v_elepT[iJ] );
    v_ele_jetPdgIds_.push_back(v_eleJetPdgId[iJ]);
    v_ele_jetNrecoEle_.push_back( v_jetNrecoEle[iJ] );
    v_ele_jetrecoEle1dR_.push_back( v_jetrecoEle1dR[iJ] );
    v_ele_jetrecoEle2dR_.push_back( v_jetrecoEle2dR[iJ] );

  }
  
  //for ( int passedJetIdx : passedJetIdxs ) {
  /*for ( auto & passedJetIdx : passedJetIdxs ) {
    std::cout << " in ele classification code passed jet index is :" << passedJetIdx << std::endl;
    reco::PFJetRef iJet( jets, passedJetIdxs[passedJetIdx] );

    std::cout << passedJetIdx << " ======================Jet eta is greater than 5 :" << v_eleJetEta[0] << "  " << iJet->eta() << std::endl;
    // Fill histograms 
    /*h_ele_jet_pT->Fill( v_eleJetPt[passedJetIdx] );
    h_ele_jet_eta->Fill( v_eleJetEta[passedJetIdx] );
    h_ele_jet_E->Fill( v_eleJetE[passedJetIdx] );
    h_ele_jet_m0->Fill( v_eleJetM0[passedJetIdx] );
    h_ele_jet_isEle->Fill( v_jetIsEle[passedJetIdx] );
    h_ele_jet_dR->Fill( v_eledR[passedJetIdx] );
    h_ele_goodvertices->Fill( v_elegoodvertices[passedJetIdx] );
    h_ele_gen_pT->Fill( v_elepT[passedJetIdx] );

    // Fill branches 
    std::cout << passedJetIdx << " ======================Jet eta is greater than 5 :" << v_eleJetEta[passedJetIdx] << "  " << iJet->eta() << std::endl;
    v_ele_jet_pt_.push_back( v_eleJetPt[passedJetIdx] );
    v_ele_jet_eta_.push_back( v_eleJetEta[passedJetIdx] );
    v_ele_jet_e_.push_back( v_eleJetE[passedJetIdx] );
    v_ele_jet_m0_.push_back( v_eleJetM0[passedJetIdx] );
    v_ele_jetIsEle_.push_back( v_jetIsEle[passedJetIdx] );
    v_ele_jetdR_.push_back( v_eledR[passedJetIdx] );
    v_ele_goodvertices_.push_back( v_elegoodvertices[passedJetIdx] );
    v_ele_gen_pt_.push_back( v_elepT[passedJetIdx] );
    v_ele_jetPdgIds_.push_back(v_eleJetPdgId[passedJetIdx]);
    v_ele_jetNrecoEle_.push_back( v_jetNrecoEle[passedJetIdx] );
    v_ele_jetrecoEle1dR_.push_back( v_jetrecoEle1dR[passedJetIdx] );
    v_ele_jetrecoEle2dR_.push_back( v_jetrecoEle2dR[passedJetIdx] );

    // Gen jet constituents
    v_ele_subJetE_[passedJetIdx].clear();
    v_ele_subJetPx_[passedJetIdx].clear();
    v_ele_subJetPy_[passedJetIdx].clear();
    v_ele_subJetPz_[passedJetIdx].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    //if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
     // if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_ele_subJetE_[passedJetIdx].push_back( subJet->energy() );
      v_ele_subJetPx_[passedJetIdx].push_back( subJet->px() );
      v_ele_subJetPy_[passedJetIdx].push_back( subJet->py() );
      v_ele_subJetPz_[passedJetIdx].push_back( subJet->pz() );
    }
  }*/

} // fillEvtSel_jet_ele_classification()
