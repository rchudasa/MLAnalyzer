#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
vector<float> v_tau_goodvertices_;
vector<float> v_tau_jetIsTau_;

vector<float> v_tau_jet_pt_;
vector<float> v_tau_jet_eta_;
vector<float> v_tau_jet_e_;
vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_Px_;
vector<float> v_tau_jet_Py_;
vector<float> v_tau_jet_Pz_;
vector<float> v_tau_genRecoJetdR_;
vector<float> v_tau_gen_pt_;
vector<float> v_tau_gen_eta_;
vector<float> v_tau_gen_pdgId_;
vector<float> v_tau_gen_prongs_;

//vector<float> v_tau_jetPFCandE_[nJets];
//vector<float> v_tau_jetPFCandPx_[nJets];
//vector<float> v_tau_jetPFCandPy_[nJets];
//vector<float> v_tau_jetPFCandPz_[nJets];

vector<vector<float>> v_tau_jetPFCandE_;
vector<vector<float>> v_tau_jetPFCandPx_;
vector<vector<float>> v_tau_jetPFCandPy_;
vector<vector<float>> v_tau_jetPFCandPz_;
vector<vector<int>> v_tau_jetPFCandType_;


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("goodvertices", &v_tau_goodvertices_);
  tree->Branch("jet_IsTau",     &v_tau_jetIsTau_);
  tree->Branch("jet_Pt",        &v_tau_jet_pt_);
  tree->Branch("jet_Eta",       &v_tau_jet_eta_);
  tree->Branch("jet_Energy",    &v_tau_jet_e_);
  tree->Branch("jet_M",         &v_tau_jet_m0_);
  tree->Branch("jet_Px",        &v_tau_jet_Px_);
  tree->Branch("jet_Py",        &v_tau_jet_Py_);
  tree->Branch("jet_Pz",        &v_tau_jet_Pz_);
  tree->Branch("jet_GendR",     &v_tau_genRecoJetdR_);
  tree->Branch("gen_pt",        &v_tau_gen_pt_);
  tree->Branch("gen_eta",       &v_tau_gen_eta_);
  tree->Branch("gen_pdgId",     &v_tau_gen_pdgId_);
  tree->Branch("gen_Prongs",    &v_tau_gen_prongs_);  
  
  tree->Branch("jetPFCandE",        &v_tau_jetPFCandE_);
  tree->Branch("jetPFCandPx",       &v_tau_jetPFCandPx_);
  tree->Branch("jetPFCandPy",       &v_tau_jetPFCandPy_);
  tree->Branch("jetPFCandPz",       &v_tau_jetPFCandPz_);
  tree->Branch("jetPFCandType",     &v_tau_jetPFCandType_);

 /* char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "jetPFCand%d_E", iJ);
    tree->Branch(hname,            &v_tau_jetPFCandE_[iJ]);
    sprintf(hname, "jetPFCand%d_Px", iJ);
    tree->Branch(hname,            &v_tau_jetPFCandPx_[iJ]);
    sprintf(hname, "jetPFCand%d_Py", iJ);
    tree->Branch(hname,            &v_tau_jetPFCandPy_[iJ]);
    sprintf(hname, "jetPFCand%d_Pz", iJ);
    tree->Branch(hname,            &v_tau_jetPFCandPz_[iJ]);
  }*/


} // branchesEvtSel_jet_dijet_tau()

// Define struct to handle mapping for gen pho<->matched reco photons<->matched presel photons
struct jet_tau_map {
  unsigned int idx;
  std::vector<unsigned int> matchedRecoJetIdxs;
  std::vector<unsigned int> matchedRecoTauIdxs;
};
std::vector<jet_tau_map> vTauJets;

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

  unsigned int nMatchedJets = 0;
  
  //identify gen hadronic taus
  float dR;
  std::vector<unsigned int> vGenTauIdxs;
  
  bool passedGenSel = false;
  unsigned int iGenParticle = 0;
  
  for (unsigned int iG=0; iG < genParticles->size(); iG++) {
    
    reco::GenParticleRef iGen( genParticles, iG );
    
    if ( iGen->pt() > 20 && (std::abs(iGen->pdgId()) == 11 || std::abs(iGen->pdgId()) == 13) ) break; //only clean jets (lepton veto) 
    if ( std::abs(iGen->pdgId()) == 12 || std::abs(iGen->pdgId()) == 14 || std::abs(iGen->pdgId()) == 16 ) continue;
    
    if (isSignal_ && !(std::abs(iGen->pdgId()) == 15 && iGen->status() == 2) ) continue;      // for drell yan and HiggsToTauTau
    if ( !isSignal_ && !isW_ && !( iGen->status() == 23 ) ) continue;                         //for QCD background
    if ( !isSignal_ &&  isW_ && !( iGen->status() == 71 ) ) continue;                         //only for W + jet background
    
    vGenTauIdxs.push_back(iG);
    
  } //genparticles
  
  if ( debug ) std::cout << " >> vGenTauIdxs.size: " << vGenTauIdxs.size() << std::endl;
  if ( vGenTauIdxs.empty() ) return false;
  
  
  ////////// Build gen Tau-jet mapping //////////
  
  // Create mapping between gen tau<->matched jets
  // For each gen tau, find "reco" jets matched to it,
  
  float minDR = 100.;
  float minDR_fpt = -10.;
  int minDR_idx = -1;
  vTauJets.clear();
  
  // Loop over valid gen pho idxs
  for ( auto& iG : vGenTauIdxs ) {
    
    reco::GenParticleRef iGenTau( genParticles, iG );
    if ( debug ) std::cout << " >> genTau[" << iG << "]" << " pt:" << iGenTau->pt() << " eta:" << iGenTau->eta() << std::endl;
    
    std::vector<unsigned int> vMatchedRecoJetIdxs;
    
    
    // Do dR match to closest reco jets
    minDR = 100.;
    minDR_fpt = -10;
    minDR_idx = -1;
    
    for ( unsigned int iJ = 0; iJ < jets->size(); iJ++ ) {
      
      reco::PFJetRef iJet( jets, iJ );
      
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      
      
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGenTau->eta(),iGenTau->phi() );
      if ( dR > minDR ) continue;
      
      minDR = dR;
      minDR_idx = iJ;
      minDR_fpt = iJet->pt()/iGenTau->pt();
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << " pt:" << iJet->pt() << " eta:" << iJet->eta() << std::endl;
      
    } // reco jets
    
    // Require minimum dR to declare match
    // Protects against matching to PU
    // minDR only needs to be generous enough so that one of the gen taus match to a reco jets for analysis
    if ( minDR > 0.4 ) continue;
    
    // Declare reco jet matching to gen tau: only store unique reco idxs
    if ( std::find(vMatchedRecoJetIdxs.begin(), vMatchedRecoJetIdxs.end(), minDR_idx) != vMatchedRecoJetIdxs.end() ) continue;
    vMatchedRecoJetIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << " f_pt(reco/gen):" << minDR_fpt << std::endl;
 
    // store the matched jets ID to a vector 
    for ( auto& iJ : vMatchedRecoJetIdxs ) {   
      vJetIdxs.push_back( iJ );

      if ( debug ) std::cout << " ----> matched jet [" << iJ << "] " <<  std::endl;
    }

    // Store this mappin 
    jet_tau_map iTau_obj = { iG, vMatchedRecoJetIdxs };
    vTauJets.push_back( iTau_obj );
                    
                          
    if ( debug ) std::cout << "\t iG" << iG << "  Jet size: " << vMatchedRecoJetIdxs.size() << std::endl;   
    
    if ( debug ) std::cout << "\t\t\t" << " GEN particle " << iGenParticle << ", status: " << iGenTau->status() << ", id: " << iGenTau->pdgId() << ", nDaught: " << iGenTau->numberOfDaughters() << ", nMoms: " <<iGenTau->numberOfMothers() << ", mother ID: " << iGenTau->mother()->pdgId() << ", pt: "<< iGenTau->pt() << ", eta: " <<iGenTau->eta() << ", phi: " <<iGenTau->phi() << ", dR: " << dR << std::endl;
    if ( debug && iGenTau->numberOfDaughters() ==1 ) std::cout << "\t\t\t" << " Daughter 1 ID " << iGenTau->daughter(0)->pdgId() << std::endl;
    if ( debug && iGenTau->numberOfDaughters() ==2 ) std::cout << "\t\t\t" << " Daughter 1 ID " << iGenTau->daughter(0)->pdgId() << " Daughter 2 ID " << iGenTau->daughter(1)->pdgId() << std::endl;
    
    
  }//gentau
  

  if(vJetIdxs.empty()) return false;

  if ( debug ) std::cout << "\t" <<" Event contains a tau candidate" << std::endl;
  return true;
  
} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);
  
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);

  v_tau_goodvertices_.clear();
  v_tau_jetIsTau_.clear();
  
  v_tau_jet_pt_.clear();
  v_tau_jet_eta_.clear();
  v_tau_jet_e_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jet_Px_.clear();
  v_tau_jet_Py_.clear();
  v_tau_jet_Pz_.clear();
  v_tau_genRecoJetdR_.clear();
  v_tau_gen_pt_.clear();
  v_tau_gen_eta_.clear();
  v_tau_gen_pdgId_.clear();
  v_tau_gen_prongs_.clear();
  
  v_tau_jetPFCandE_.clear();
  v_tau_jetPFCandPx_.clear();
  v_tau_jetPFCandPy_.clear();
  v_tau_jetPFCandPz_.clear();
  v_tau_jetPFCandType_.clear();


  unsigned int goodVertices = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
 
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;

  if ( debug ) std::cout << "\t" << " JETS IN THE EVENT = " << vTauJets.size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  
  v_tau_goodvertices_.push_back(goodVertices);
 
  int iJ = 0; 
  for ( auto const& ii: vTauJets ) {
    
    // Skip electrons which fails HE edge cut
    if(ii.matchedRecoJetIdxs.empty())continue;
    if ( std::find(vJetIdxs.begin(), vJetIdxs.end(), ii.matchedRecoJetIdxs[0]) == vJetIdxs.end()) continue;
    
    reco::GenParticleRef iGen( genParticles, ii.idx );
    reco::PFJetRef iJet( jets, ii.matchedRecoJetIdxs[0] );
    
    iJ  = ii.matchedRecoJetIdxs[0];
    
    if ( debug )  std::cout << " --------------------------------- Filling branches --------------------------------- " << std::endl;
    if ( debug )  std::cout << " Gen pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;       
    if ( debug )  std::cout << " Jet pt: "<< iJet->pt() << " eta: " <<iJet->eta() << " phi: " <<iJet->phi() << std::endl;       
    
    v_tau_jetIsTau_.push_back(1);
    
    
    // Fill branches 
    v_tau_jet_pt_.push_back( iJet->pt() );
    v_tau_jet_eta_.push_back( iJet->eta() );
    v_tau_jet_e_.push_back( iJet->energy() ); 
    v_tau_jet_m0_.push_back( iJet->mass() );
    v_tau_jet_Px_.push_back( iJet->px() );
    v_tau_jet_Py_.push_back( iJet->py() );
    v_tau_jet_Pz_.push_back( iJet->pz() );    
    
    TLorentzVector TLVJet(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
    TLorentzVector TLVGen(iGen->px(),iGen->py(),iGen->pz(),iGen->energy());
    
    v_tau_genRecoJetdR_.push_back(TLVJet.DeltaR(TLVGen) );
    
    v_tau_gen_pt_.push_back(iGen->pt());
    v_tau_gen_eta_.push_back(iGen->eta());
    v_tau_gen_pdgId_.push_back(iGen->pdgId());
    
    int tauDaughters          = 0;
    int tauPi0 = 0;
    
    for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
      if ( debug ) std::cout << "\t\t\t\t" <<" Tau daughter [" << iDaughter << "] : "<<  std::abs(iGen->daughter(iDaughter)->pdgId()); 
      if ( debug ) std::cout << " charge : "<< iGen->daughter(iDaughter)->charge() << "  | pt : "<< iGen->daughter(iDaughter)->pt();
      if ( debug ) std::cout << " eta:" << iGen->daughter(iDaughter)->eta() << " |Energy:" << iGen->daughter(iDaughter)->energy() << std::endl;
      if ( abs(iGen->daughter(iDaughter)->pdgId()) == 111 ) tauPi0++;
      if ( iGen->daughter(iDaughter)->charge() == 0 ) continue;          
      tauDaughters++;
    }
    
    if ( debug ) std::cout << "\t\t\t"<<" Tau prongs = " << tauDaughters << " + Tau pi0 = " << tauPi0 << std::endl;      
    v_tau_gen_prongs_.push_back(tauDaughters);
    
    // jet constituents
    vector<float> pfcand_px;
    vector<float> pfcand_py;
    vector<float> pfcand_pz;
    vector<float> pfcand_energy;
    vector<int> pfcand_type;
 
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr jetPFCand = iJet->getPFConstituent( j );
      //if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;
      std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;

      pfcand_px.push_back(jetPFCand->px());
      pfcand_py.push_back(jetPFCand->py());
      pfcand_pz.push_back(jetPFCand->pz());
      pfcand_energy.push_back(jetPFCand->energy());
      pfcand_type.push_back((int)jetPFCand->particleId());
    }//jet constituents loop
    iJ++;
   
      v_tau_jetPFCandE_.push_back( pfcand_energy );
      v_tau_jetPFCandPx_.push_back( pfcand_px );
      v_tau_jetPFCandPy_.push_back( pfcand_py );
      v_tau_jetPFCandPz_.push_back( pfcand_pz );
      v_tau_jetPFCandType_.push_back( pfcand_type );

  }//taujets loop
  
} // fillEvtSel_jet_dijet_tau()
