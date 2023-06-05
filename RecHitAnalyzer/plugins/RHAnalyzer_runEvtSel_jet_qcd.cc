#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 4; //TODO: use cfg level nJets_

vector<float> v_goodvertices_;
vector<float> v_qcdJetIsEle_;
vector<float> v_jet_pt_;
vector<float> v_jet_eta_;
vector<float> v_jet_e_;
vector<float> v_jet_m0_;
vector<float> v_genRecoJetdR_;
vector<float> v_gen_pt_;
vector<float> v_gen_eta_;
vector<float> v_gen_pdgId_;


// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_qcd ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("goodvertices",  &v_goodvertices_);
  tree->Branch("jet_IsEle",     &v_qcdJetIsEle_);
  tree->Branch("jet_Pt",        &v_jet_pt_);
  tree->Branch("jet_Eta",       &v_jet_eta_);
  tree->Branch("jet_Energy",    &v_jet_e_);
  tree->Branch("jet_M",         &v_jet_m0_);
  tree->Branch("jet_GendR",     &v_genRecoJetdR_);
  tree->Branch("gen_pt",        &v_gen_pt_);
  tree->Branch("gen_eta",       &v_gen_eta_);
  tree->Branch("gen_pdgId",     &v_gen_pdgId_);

} // branchesEvtSel_jet_qcd()

// Define struct to handle mapping for gen pho<->matched reco photons<->matched presel photons
struct jet_gen_map {
  unsigned int idx;
  std::vector<unsigned int> matchedRecoJetIdxs;
};
std::vector<jet_gen_map> vGens;

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_qcd( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  //if(iEvent.id().event()>657412){

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  vJetIdxs.clear();
    
  // Identify a parton jet
  float dR;
  std::vector<unsigned int> vGenIdxs;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
    
    reco::GenParticleRef iGen( genParticles, iG );
    if (!isSignal_ && !( iGen->status() == 23 ) ) continue;     
    
    vGenIdxs.push_back( iG );
    
  } // genParticles
  
  if ( debug ) std::cout << " >> vGenIdxs.size: " << vGenIdxs.size() << std::endl;
  if ( vGenIdxs.empty() ) return false;
 
  
  
  ////////// Build gen-jet mapping //////////
  
  // Create mapping between gen <->matched jets
  // For each gen partons, find "reco" jets matched to it,
  
  float minDR = 100.;
  float minDR_fpt = -10.;
  int minDR_idx = -1;
  vGens.clear();
  
  // Loop over valid gen pho idxs
  for ( auto& iG : vGenIdxs ) {
    
    reco::GenParticleRef iGenEle( genParticles, iG );
    if ( debug ) std::cout << " >> genParton[" << iG << "]" << " pt:" << iGenEle->pt() << " eta:" << iGenEle->eta() << std::endl;
    
    std::vector<unsigned int> vMatchedRecoJetIdxs;
    
    // Do dR match to closest reco photon
    minDR = 100.;
    minDR_fpt = -10;
    minDR_idx = -1;
    for ( unsigned int iJ = 0; iJ < jets->size(); iJ++ ) {
      
      reco::PFJetRef iJet( jets, iJ );
      
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      
      
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGenEle->eta(),iGenEle->phi() );
      if ( dR > minDR ) continue;
      
      minDR = dR;
      minDR_idx = iJ;
      minDR_fpt = iJet->pt()/iGenEle->pt();      
    } // reco jets
    
    
    // Require minimum dR to declare match
    // Protects against matching to PU
    // minDR only needs to be generous enough so that one of the gen electrons match to a reco jets for analysis
    if ( minDR > 0.4 ) continue;
    
    // Declare reco jet matching to gen electron: only store unique reco idxs
    if ( std::find(vMatchedRecoJetIdxs.begin(), vMatchedRecoJetIdxs.end(), minDR_idx) != vMatchedRecoJetIdxs.end() ) continue;
    vMatchedRecoJetIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << " f_pt(reco/gen):" << minDR_fpt << std::endl;


     
    // Check if matched reco jets also matches to electron:
    for ( auto& iJ : vMatchedRecoJetIdxs ) {
      
      reco::PFJetRef iJet( jets, iJ );      
      vJetIdxs.push_back( iJ );
      
      if ( debug ) std::cout << " ----> matched jet [" << iJ << "] " << " jet pt: " << iJet->pt() << " eta: " << iJet->eta() << std::endl; 
    } // matched reco jets
    
   // Store this mapping
   if(vMatchedRecoJetIdxs.empty()) continue;
   jet_gen_map iEle_obj = { iG, vMatchedRecoJetIdxs };
   vGens.push_back( iEle_obj );
   
  } //gen electrons
  //}

  
  if ( vJetIdxs.size() == 0){
    if ( debug ) std::cout << " No passing jets...  " << std::endl;
    return false;
  }

  if ( debug ) std::cout << "\t" << " >> has_jet_dijet_qcd: passed" << std::endl;
  return true;

} // runEvtSel_jet_qcd()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_qcd ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );
  
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);
  
  v_goodvertices_.clear();
  v_qcdJetIsEle_.clear();
  v_jet_pt_.clear();
  v_jet_eta_.clear();
  v_jet_e_.clear();
  v_jet_m0_.clear();
  v_gen_pt_.clear();
  v_gen_eta_.clear();
  v_gen_pdgId_.clear();
  v_genRecoJetdR_.clear();

  unsigned int goodVertices = 0;
  
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;
      
  if ( debug )  std::cout << " --------------------------------- vGen size" << vGens.size() << " --------------------------------- " << std::endl;

  v_goodvertices_.push_back(goodVertices);
  
  for ( auto const& ii: vGens ) {
    
    // Skip electrons which fails HE edge cut
    if(ii.matchedRecoJetIdxs.empty())continue;
    if ( std::find(vJetIdxs.begin(), vJetIdxs.end(), ii.matchedRecoJetIdxs[0]) == vJetIdxs.end()) continue;
    
    reco::GenParticleRef iGen( genParticles, ii.idx );
    reco::PFJetRef iJet( jets, ii.matchedRecoJetIdxs[0] );
    
    if ( debug )  std::cout << " --------------------------------- Filling branches --------------------------------- " << std::endl;
    if ( debug )  std::cout << " Gen pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;       
    if ( debug )  std::cout << " Jet pt: "<< iJet->pt() << " eta: " <<iJet->eta() << " phi: " <<iJet->phi() << std::endl;       
    
    v_qcdJetIsEle_.push_back(0);
    v_jet_pt_.push_back( iJet->pt() );
    v_jet_eta_.push_back( iJet->eta() );
    v_jet_e_.push_back( iJet->energy() );
    v_jet_m0_.push_back( iJet->mass() );
    v_gen_pt_.push_back(iGen->pt());
    v_gen_eta_.push_back(iGen->eta());
    v_gen_pdgId_.push_back(iGen->pdgId());
 
    TLorentzVector TLVJet(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
    TLorentzVector TLVGen(iGen->px(),iGen->py(),iGen->pz(),iGen->energy());
       
    v_genRecoJetdR_.push_back(TLVJet.DeltaR(TLVGen) );
  }//gen loop
  


} // fillEvtSel_jet_qcd()
