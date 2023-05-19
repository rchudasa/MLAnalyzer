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
vector<float> v_eleJetPdgId;

vector<float> v_ele_jet_m0_;
vector<float> v_ele_jet_pt_;
vector<float> v_ele_gen_pt_;
vector<float> v_ele_jetPdgIds_;
vector<float> v_ele_jetIsEle_;
vector<float> v_ele_jetdR_;
vector<float> v_ele_goodvertices_;

vector<float> v_ele_subJetE_[nJets];
vector<float> v_ele_subJetPx_[nJets];
vector<float> v_ele_subJetPy_[nJets];
vector<float> v_ele_subJetPz_[nJets];


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
  tree->Branch("jet_PdgIds",    &v_ele_jetPdgIds_);
  tree->Branch("jet_IsEle",     &v_ele_jetIsEle_);
  tree->Branch("gen_pt",        &v_ele_gen_pt_);
  tree->Branch("jet_dR",        &v_ele_jetdR_);
  tree->Branch("goodvertices",  &v_ele_goodvertices_);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_ele_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_ele_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_ele_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_ele_subJetPz_[iJ]);
  }

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
  v_eleJetPdgId.clear();

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


  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  
    // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
  
    unsigned int PdgId        = 0;
    float jetdR               = -99.;
    float elepT               = -99.;
    bool JetIsEle             = false;
    
    reco::PFJetRef iJet( jets, iJ );
    if (debug ) std::cout <<"Jet [" << iJ << "] => Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0; 

  
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
    float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );//QCD
    if ( dR > 0.4 ) continue; //QCD
    
    if ( abs(iGen->pdgId()) != 11 ) continue; //DY
    //if ( iGen->numberOfMothers() != 1 ) continue; //DY
    //if ( iGen->mother()->pdgId() != 23)continue; //DY
    //if(abs(iGen->pdgId()) == 21 || abs(iGen->pdgId())<9)continue;//QCD
    //if ( !( iGen->status() == 23 ) ) continue; //QCD 
    ++iGenParticle;
    //float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );//DY
    if ( debug ) std::cout << "\t" <<"GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR = "<< dR << std::endl;         
    if ( debug ) std::cout << "\t" <<"MOTHER => status: " << iGen->mother()->status() << ", id: " << iGen->mother()->pdgId() << ", nDaught: " << iGen->mother()->numberOfDaughters() << " | pt: "<< iGen->mother()->pt() << " eta: " <<iGen->mother()->eta() << " phi: " <<iGen->mother()->phi() << " mass: " <<iGen->mother()->mass() << std::endl;
    

     if ( std::abs(iGen->pdgId()) == 11 ) {
        JetIsEle = true;
        PdgId = std::abs(iGen->pdgId());
        jetdR = dR;
        elepT = iGen->pt();
   
 	if (!isSignal_){  
        passedGenSel = false; // for backg
        break; 		      // for backg
	}

      } // if pDG ID ==15
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
      vJetIdxs.push_back( iJ );
      v_eleJetPdgId.push_back( PdgId );
      v_elepT.push_back( elepT );
      v_eledR.push_back( jetdR );
      v_elegoodvertices.push_back( goodVertices );
      v_jetIsEle.push_back( JetIsEle ); 

    } //passed gen level selection

  } // reco jets
  if ( debug ) std::cout << " Matched GEN particle-jet pairs " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << " >> has_jet_dijet_ele_classification: passed" << std::endl;
  return true;

} // runEvtSel_jet_ele_classification()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_ele_classification ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_ele_jet_nJet->Fill( vJetIdxs.size() );

  v_ele_jet_pt_.clear();
  v_ele_gen_pt_.clear();
  v_ele_jet_m0_.clear();
  v_ele_jetIsEle_.clear();
  v_ele_jetdR_.clear();
  v_ele_goodvertices_.clear();
  v_ele_jetPdgIds_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_ele_jet_pT->Fill( std::abs(iJet->pt()) );
    h_ele_jet_eta->Fill( iJet->eta() );
    h_ele_jet_E->Fill( iJet->energy() );
    h_ele_jet_m0->Fill( iJet->mass() );
    h_ele_jet_isEle->Fill( v_jetIsEle[iJ] );
    h_ele_jet_dR->Fill( v_eledR[iJ] );
    h_ele_goodvertices->Fill( v_elegoodvertices[iJ] );
    h_ele_gen_pT->Fill( v_elepT[iJ] );

    // Fill branches 
    v_ele_jet_pt_.push_back( iJet->pt() );
    v_ele_jet_m0_.push_back( iJet->mass() );
    v_ele_jetIsEle_.push_back( v_jetIsEle[iJ] );
    v_ele_jetdR_.push_back( v_eledR[iJ] );
    v_ele_goodvertices_.push_back( v_elegoodvertices[iJ] );
    v_ele_gen_pt_.push_back( v_elepT[iJ] );
    v_ele_jetPdgIds_.push_back(v_eleJetPdgId[iJ]);

    // Gen jet constituents
    v_ele_subJetE_[iJ].clear();
    v_ele_subJetPx_[iJ].clear();
    v_ele_subJetPy_[iJ].clear();
    v_ele_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_ele_subJetE_[iJ].push_back( subJet->energy() );
      v_ele_subJetPx_[iJ].push_back( subJet->px() );
      v_ele_subJetPy_[iJ].push_back( subJet->py() );
      v_ele_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_ele_classification()
