#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_

TH1D *h_top_goodvertices;

TH1D *h_top_jet_nJet;
TH1D *h_top_jet_isTau;
TH1D *h_top_jet_pT;
TH1D *h_top_jet_eta;
TH1D *h_top_jet_E;
TH1D *h_top_jet_m0;
TH1D *h_top_jet_dR;

vector<float> v_top_goodvertices_;
vector<float> v_top_jetIsTau_;
vector<float> v_top_jet_pt_;
vector<float> v_top_jet_eta_;
vector<float> v_top_jet_energy_;
vector<float> v_top_jet_m0_;
vector<float> v_top_jetdR_;

std::vector<std::vector<int> > seljet_genpart_collid;
std::vector<std::vector<int> > seljet_genpart_pdgid;
std::vector<std::vector<int> > seljet_genpart_charge;

std::vector<std::vector<float> > seljet_genpart_px;
std::vector<std::vector<float> > seljet_genpart_py;
std::vector<std::vector<float> > seljet_genpart_pz;
std::vector<std::vector<float> > seljet_genpart_energy;

std::vector<std::vector<int> > seljet_genpart_status;

std::vector<std::vector<int> > seljet_genpart_motherpdgid;
std::vector<std::vector<int> > seljet_genpart_dau1pdgid;
std::vector<std::vector<int> > seljet_genpart_dau2pdgid;


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_top ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_top_goodvertices   = fs->make<TH1D>("h_goodvertices"   , "good vertices;good vertices;Jets"           ,  15,  0.,  75.);
  h_top_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_top_jet_isTau      = fs->make<TH1D>("h_jet_isTau"      , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_top_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                           , 100,  0., 500.);
  h_top_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_top_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                   , 100,  0., 500.); 
  h_top_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  h_top_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{jet,#top};dR_{jet,#top};Jets"           ,  50,  0.,   1.);

  tree->Branch("goodvertices",  &v_top_goodvertices_);
  tree->Branch("jet_IsTau",     &v_top_jetIsTau_);
  tree->Branch("jet_Pt",        &v_top_jet_pt_);
  tree->Branch("jet_Eta",       &v_top_jet_eta_);
  tree->Branch("jet_Energy",    &v_top_jet_energy_);
  tree->Branch("jet_M",         &v_top_jet_m0_);
  tree->Branch("jet_dR",        &v_top_jetdR_);
 
  tree->Branch("seljet_genpart_collid", &seljet_genpart_collid);
  tree->Branch("seljet_genpart_pdgid", &seljet_genpart_pdgid);
  tree->Branch("seljet_genpart_charge", &seljet_genpart_charge);

  tree->Branch("seljet_genpart_px", &seljet_genpart_px);
  tree->Branch("seljet_genpart_py", &seljet_genpart_py);
  tree->Branch("seljet_genpart_pz", &seljet_genpart_pz);
  tree->Branch("seljet_genpart_energy", &seljet_genpart_energy);

  tree->Branch("seljet_genpart_status", &seljet_genpart_status);

  tree->Branch("seljet_genpart_motherpdgid", &seljet_genpart_motherpdgid);
  tree->Branch("seljet_genpart_dau1pdgid", &seljet_genpart_dau1pdgid);
  tree->Branch("seljet_genpart_dau2pdgid", &seljet_genpart_dau2pdgid);
} // branchesEvtSel_jet_dijet_top()


// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_top( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);

  vJetIdxs.clear();

  int nJet = 0;
  
  std::vector<TLorentzVector> had_tops,bdau,wdau;
  if (isttbar_) { //is a ttbar sample
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
  
    int id = iGen->pdgId();
    if(abs(id) != 6 || iGen->numberOfDaughters()!=2) continue;
    int iw=-1;
    int ib=-1;
    if (abs(iGen->daughter(0)->pdgId())==24 && abs(iGen->daughter(1)->pdgId())==5)
    {
      iw=0;ib=1;
    }
    else
    {
      if(abs(iGen->daughter(1)->pdgId())==24 && abs(iGen->daughter(0)->pdgId())==5)
      {
        iw=1;ib=0;
      }
      else continue;
    }
    const reco::Candidate *d = iGen->daughter(iw);
    const reco::Candidate *b = iGen->daughter(ib);
    while(d->numberOfDaughters() == 1) d = d->daughter(0);
    if(!(abs(d->daughter(0)->pdgId()) < 10 && abs(d->daughter(1)->pdgId()) < 10)) continue;
    TLorentzVector the_top,the_w,the_b;
    the_top.SetPtEtaPhiE(iGen->pt(),iGen->eta(),iGen->phi(),iGen->energy());
    the_w.SetPtEtaPhiE(d->pt(),d->eta(),d->phi(),d->energy());
    the_b.SetPtEtaPhiE(b->pt(),b->eta(),b->phi(),b->energy());
    had_tops.push_back(the_top);
    wdau.push_back(the_w);
    bdau.push_back(the_b);
  } //gen particle loop


  // Loop over jets
  for ( unsigned ihad=0;ihad<had_tops.size();ihad++)
  {
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ )
    {
      reco::PFJetRef iJet( jets, iJ );
      TLorentzVector vjet;
      vjet.SetPtEtaPhiE(iJet->pt(),iJet->eta(),iJet->phi(),iJet->energy());

      if ( std::abs(iJet->pt()) < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_) continue;
      if (had_tops[ihad].DeltaR(vjet)>0.8) continue;
      if (wdau[ihad].DeltaR(vjet)>0.8) continue;
      if (bdau[ihad].DeltaR(vjet)>0.8) continue;

      if ( debug ) std::cout << " >> jet[" << iJ << "]Pt:" << iJet->pt() << " jetE:" << iJet->energy() << " jetM:" << iJet->mass() << std::endl;

      vJetIdxs.push_back(iJ);

      nJet++;
      break;// This should allow two hardonic tops
    } // jets
    if ( (nJets_ > 0) && (nJet >= nJets_) ) break;
  } // hadronic tops
  } // isTTbar
 
  if ( debug ) {
    for(int thisJetIdx : vJetIdxs)
      std::cout << " >> vJetIdxs:" << thisJetIdx << std::endl;
  }

  if ( (nJets_ > 0) && (nJet != nJets_) ){
    if ( debug ) std::cout << " Fail jet multiplicity:  " << nJet << " < " << nJets_ << std::endl;
    return false;
  }

  if ( vJetIdxs.size() == 0){
    if ( debug ) std::cout << " No passing jets...  " << std::endl;
    return false;
  }

  if ( debug ) std::cout << " >> has_jet_dijet: passed" << std::endl;
  return true;
   
  unsigned int nMatchedJets = 0;

 

} // runEvtSel_jet_dijet_top()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_top ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);
  
  unsigned int goodVertices = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
 
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  h_top_jet_nJet->Fill( vJetIdxs.size() );

  v_top_goodvertices_.clear();
  v_top_jetIsTau_.clear();
  v_top_jet_pt_.clear();
  v_top_jet_eta_.clear();
  v_top_jet_energy_.clear();
  v_top_jet_m0_.clear();
  v_top_jetdR_.clear();

  
  seljet_genpart_collid.clear();
  seljet_genpart_pdgid.clear();
  seljet_genpart_charge.clear();

  seljet_genpart_px.clear();
  seljet_genpart_py.clear();
  seljet_genpart_pz.clear();
  seljet_genpart_energy.clear();

  seljet_genpart_status.clear();

  seljet_genpart_motherpdgid.clear();
  seljet_genpart_dau1pdgid.clear();
  seljet_genpart_dau2pdgid.clear();
  
  for (int iJ : vJetIdxs){

    reco::PFJetRef thisJet( jets, iJ );
    if ( debug ) std::cout << " >> Jet[" << iJ << "] Pt:" << thisJet->pt() << std::endl;

    // Fill histograms 
    h_top_jet_isTau->Fill( 0 );
    h_top_jet_pT->Fill( std::abs(thisJet->pt()) );
    h_top_jet_eta->Fill( thisJet->eta() );
    h_top_jet_E->Fill( thisJet->energy() );
    h_top_jet_m0->Fill( thisJet->mass() );
   
    v_top_jetIsTau_.push_back( 0 ); 
    v_top_jet_pt_.push_back( thisJet->pt() );
    v_top_jet_eta_.push_back( thisJet->eta() );
    v_top_jet_energy_.push_back( thisJet->energy() );
    v_top_jet_m0_.push_back( thisJet->mass() );
    
   
    TLorentzVector TLVJet(thisJet->px(),thisJet->py(),thisJet->pz(),thisJet->energy());
    double cosTheta = TLVJet.CosTheta();
      if (cosTheta*cosTheta >=0)
        TLVJet.SetPx(0.0001);

    std::vector<int> genpart_collid;
    std::vector<int> genpart_pdgid;
    std::vector<int> genpart_charge;

    std::vector<float> genpart_px;
    std::vector<float> genpart_py;
    std::vector<float> genpart_pz;
    std::vector<float> genpart_energy;

    std::vector<int> genpart_status;

    std::vector<int> genpart_motherpdgid;
    std::vector<int> genpart_dau1pdgid;
    std::vector<int> genpart_dau2pdgid;

    std::vector<reco::GenParticle>::const_iterator genpartIterator      = (genParticles.product())->begin();
    std::vector<reco::GenParticle>::const_iterator genpartIteratorEnd   = (genParticles.product())->end();
    for ( ; genpartIterator != genpartIteratorEnd; genpartIterator++)
    {

      TLorentzVector TLVgenpart(genpartIterator->px(),genpartIterator->py(),genpartIterator->pz(),genpartIterator->energy());
      cosTheta = TLVgenpart.CosTheta();
      if (cosTheta*cosTheta >=0)
        TLVgenpart.SetPx(0.0001);
      
      
      if (TLVJet.DeltaR(TLVgenpart)<0.8)
      { 
      	float dR = TLVJet.DeltaR(TLVgenpart);
      	v_top_jetdR_.push_back(dR);
        genpart_collid.push_back(genpartIterator->collisionId());
        genpart_pdgid.push_back(genpartIterator->pdgId());
        genpart_charge.push_back(genpartIterator->charge());

        genpart_px.push_back(genpartIterator->px());
        genpart_py.push_back(genpartIterator->py());
        genpart_pz.push_back(genpartIterator->pz());
        genpart_energy.push_back(genpartIterator->energy());

        genpart_status.push_back(genpartIterator->status());

        if (genpartIterator->numberOfMothers()>0)
        {
          genpart_motherpdgid.push_back(genpartIterator->mother(0)->pdgId());
        }
        else
        {
          genpart_motherpdgid.push_back(-9999);
        }

        switch (genpartIterator->numberOfDaughters())
        {
          case 0:
            genpart_dau1pdgid.push_back(-9999);
            genpart_dau2pdgid.push_back(-9999);
          break;

          case 1:
            genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
            genpart_dau2pdgid.push_back(-9999);
          break;

          default:
            genpart_dau1pdgid.push_back(genpartIterator->daughter(0)->pdgId());
            genpart_dau2pdgid.push_back(genpartIterator->daughter(1)->pdgId());
          break;
        }//switch
      }//DeltaR condition
    }//genpart loop
    seljet_genpart_collid.push_back(genpart_collid);
    seljet_genpart_pdgid.push_back(genpart_pdgid);
    seljet_genpart_charge.push_back(genpart_charge);

    seljet_genpart_px.push_back(genpart_px);
    seljet_genpart_py.push_back(genpart_py);
    seljet_genpart_pz.push_back(genpart_pz);
    seljet_genpart_energy.push_back(genpart_energy);

    seljet_genpart_status.push_back(genpart_status);

    seljet_genpart_motherpdgid.push_back(genpart_motherpdgid);
    seljet_genpart_dau1pdgid.push_back(genpart_dau1pdgid);
    seljet_genpart_dau2pdgid.push_back(genpart_dau2pdgid);    
  }

} // fillEvtSel_jet_dthisJet_top()
