#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 4; //TODO: use cfg level nJets_

vector<float> v_ele_goodvertices_;
vector<float> v_jetIsEle_;


vector<float> v_ele_jet_pt_;
vector<float> v_ele_jet_eta_;
vector<float> v_ele_jet_e_;
vector<float> v_ele_jet_m0_;
vector<float> v_ele_genRecoJetdR_;
vector<float> v_ele_gen_pt_;
vector<float> v_ele_gen_eta_;
vector<float> v_ele_gen_pdgId_;
vector<float> v_ele_NrecoEle_;
vector<float> v_ele_genRecoEledR_;
vector<float> v_ele_jetEledR_;



// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_ele_classification ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("goodvertices",  &v_ele_goodvertices_);
  tree->Branch("jet_IsEle",     &v_jetIsEle_);
  tree->Branch("jet_Pt",        &v_ele_jet_pt_);
  tree->Branch("jet_Eta",       &v_ele_jet_eta_);
  tree->Branch("jet_Energy",    &v_ele_jet_e_);
  tree->Branch("jet_M",         &v_ele_jet_m0_);
  tree->Branch("jet_GendR",     &v_ele_genRecoJetdR_);
  tree->Branch("jet_EledR",     &v_ele_jetEledR_);
  tree->Branch("gen_pt",        &v_ele_gen_pt_);
  tree->Branch("gen_eta",       &v_ele_gen_eta_);
  tree->Branch("gen_pdgId",     &v_ele_gen_pdgId_);
  tree->Branch("NrecoEle",      &v_ele_NrecoEle_);
  tree->Branch("gen_EledR",     &v_ele_genRecoEledR_);





} // branchesEvtSel_jet_ele_classification()

// Define struct to handle mapping for gen pho<->matched reco photons<->matched presel photons
struct jet_ele_map {
  unsigned int idx;
  std::vector<unsigned int> matchedRecoJetIdxs;
  std::vector<unsigned int> matchedRecoEleIdxs;
};
std::vector<jet_ele_map> vEles;

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_ele_classification( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  
  //if(iEvent.id().event()>657412){

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::GsfElectronCollection> ele;    //TODO
  iEvent.getByToken(eleCollectionT_, ele);         //TODO
  
  vJetIdxs.clear();
  
  unsigned int nMatchedRecoEle = 0;
  float recoele1dR = -99.;
  float recoele2dR = -99.;
    
  // Identify particle gun gen electrons from the event
  float dR;
  std::vector<unsigned int> vGenEleIdxs;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {
    
    reco::GenParticleRef iGen( genParticles, iG );
    
    if (isSignal_ && !(std::abs(iGen->pdgId()) == 11 && iGen->status() == 23 && iGen->numberOfMothers() == 1 ) ) continue; // for drell yan to e+e-
    if (!isSignal_ && !( iGen->status() == 23 ) ) continue;     
    
    vGenEleIdxs.push_back( iG );
    
  } // genParticles
  
  if ( debug ) std::cout << " >> vGenEleIdxs.size: " << vGenEleIdxs.size() << std::endl;
  if ( vGenEleIdxs.empty() ) return false;
 
  
  
  ////////// Build gen electron-jet mapping //////////
  
  // Create mapping between gen pho<->matched jets<->matched electrons
  // For each gen electron, find "reco" jets matched to it,
  // then check if that reco jet belongs to an electron or not
  float minDR = 100.;
  float minDR_fpt = -10.;
  int minDR_idx = -1;
  vEles.clear();
  
  // Loop over valid gen pho idxs
  for ( auto& iG : vGenEleIdxs ) {
    
    reco::GenParticleRef iGenEle( genParticles, iG );
    if ( debug ) std::cout << " >> genPho[" << iG << "]" << " pt:" << iGenEle->pt() << " eta:" << iGenEle->eta() << std::endl;
    
    std::vector<unsigned int> vMatchedRecoJetIdxs;
    std::vector<unsigned int> vMatchedRecoJetEleIdxs;
    std::vector<unsigned int> vMatchedRecoEleIdxs;
    
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
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << " pt:" << iJet->pt() << " eta:" << iJet->eta() << std::endl;
      
    } // reco jets
    
    
    // Require minimum dR to declare match
    // Protects against matching to PU
    // minDR only needs to be generous enough so that one of the gen electrons match to a reco jets for analysis
    if ( minDR > 0.4 ) continue;
    
    // Declare reco jet matching to gen electron: only store unique reco idxs
    if ( std::find(vMatchedRecoJetIdxs.begin(), vMatchedRecoJetIdxs.end(), minDR_idx) != vMatchedRecoJetIdxs.end() ) continue;
    vMatchedRecoJetIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << " f_pt(reco/gen):" << minDR_fpt << std::endl;

     int minDREle_idx=-1;
     int minDREle1_idx=-1;
     int minDREle2_idx=-1;
     
    // Check if matched reco jets also matches to electron:
    for ( auto& iJ : vMatchedRecoJetIdxs ) {
      
      reco::PFJetRef iJet( jets, iJ );
      
      //Lookin at RecoEle
      if (debug ) std::cout << "\t\t\t\t" << "  Looking at RECO electrons: "<< std::endl;
      if (debug && ele->size() == 0 ) std::cout << "\t\t\t\t" << "   !!!!!!!!!!  NO RECO ELE IN THIS EVENT  !!!!!!!!!!"<< std::endl;
      if (ele->size() == 0) continue;
  
      minDREle_idx  = -1;    
      minDREle1_idx = -1;    
      minDREle2_idx = -1;    
      int eleIdx = -1;    
      for ( unsigned iT(0); iT != ele->size(); ++iT ) {
	
        reco::GsfElectronRef iEle( ele, iT );
        float recoeledR = reco::deltaR( iJet->eta(),iJet->phi(), iEle->eta(),iEle->phi() );
        minDREle_idx = iT;
        if ( recoeledR < 0.4 && nMatchedRecoEle == 0 ) {
          if ( debug ) std::cout << "\t\t\t\t" << "   Reco Ele [" << iT << "]  matched jet [" << iJ << "] : dR = " << recoeledR << " electron pT: " << iEle->pt() << "  eta:" << iEle->eta() << "  phi:"<< iEle->phi() <<  std::endl;
          recoele1dR = recoeledR;
          minDREle1_idx = iT;
          ++nMatchedRecoEle;
        } 
        
        else if ( recoeledR < 0.4 && nMatchedRecoEle == 1 ) {
          if ( debug ) std::cout << "\t\t\t\t" << "  second Reco Ele [" << iT << "]  matched jet [" << iJ << "] : dR = " << recoeledR << " electron pT: " << iEle->pt() << "  eta:" << iEle->eta() << "  phi:"<< iEle->phi() <<  std::endl;
          if (recoeledR < recoele1dR) {
            recoele2dR = recoele1dR;
            minDREle2_idx = minDREle1_idx;
            recoele1dR = recoeledR;
            minDREle1_idx = minDREle_idx;
          } else { recoele2dR = recoeledR; minDREle2_idx = minDREle_idx;}
          ++nMatchedRecoEle;
        } 
        
        else if ( recoeledR < 0.4 && nMatchedRecoEle > 1 ) {
          std::cout << "\t\t\t\t" << "   !!!!!!!!!!  FOUND MORE THAN 2 ELE INSIDE JET CONE OF 0.4 !!!!!!!!!!"<< std::endl;
          if (recoeledR < recoele2dR && recoeledR < recoele1dR) { 
            if (recoele1dR < recoele2dR) recoele2dR = recoele1dR;
            recoele1dR = recoeledR;
          } else if (recoeledR < recoele2dR && recoeledR > recoele1dR) recoele2dR = recoeledR; 
          ++nMatchedRecoEle;
        } 
        
        else if ( debug ) {
          std::cout << "\t\t\t\t" << "   !!!!!!!!!!  NO MATCH FOR Reco Ele [" << iT << "]  with jet [" << iJ << "] : dR = " << recoeledR << std::endl;
        }
      } //electrons
      
      if(nMatchedRecoEle==0) continue;
      if ( debug )std::cout << " There are " << nMatchedRecoEle << " matched electrons" <<  std::endl;
      if ( debug )std::cout << " matched jet [" << iJ << "] " << "Electron 1 [" << minDREle1_idx << "]"<< "  Electron 2 [" << minDREle2_idx << "]"<<  std::endl;

      if(minDREle1_idx >-1 && minDREle2_idx > -1 )continue;
      if(minDREle1_idx ==-1 && minDREle2_idx==-1 )continue;
      if(minDREle1_idx ==-1 && minDREle2_idx>-1 )eleIdx = minDREle2_idx;
       if(minDREle2_idx ==-1 && minDREle1_idx>-1 )eleIdx = minDREle1_idx;
      
      vMatchedRecoEleIdxs.push_back(eleIdx);
      
      vJetIdxs.push_back( iJ );
      vMatchedRecoJetEleIdxs.push_back(iJ);
      
      if ( debug ) std::cout << " ----> matched jet [" << iJ << "] " << "Electron Idx [" << eleIdx << "]"<<  std::endl;
      if ( debug ) std::cout << " >> presel jet pt: " << iJet->pt() << " eta: " << iJet->eta() << std::endl;


      
    } // matched reco jets
    
   // Store this mappin 
   if(vMatchedRecoJetEleIdxs.empty() || vMatchedRecoEleIdxs.empty()) continue;
   jet_ele_map iEle_obj = { iG, vMatchedRecoJetEleIdxs, vMatchedRecoEleIdxs };
   vEles.push_back( iEle_obj );
   
      
  if ( debug ) std::cout << "\t iG" << iG << "  Jet-electron size: " << vMatchedRecoJetEleIdxs.size() << "  vMatchedRecoEleIdxs size:" <<  vMatchedRecoEleIdxs.size() << std::endl;   
  } //gen electrons
  //}


  if(vJetIdxs.empty()) return false;

  if ( debug ) std::cout << "\t" << " >> has_jet_dijet_ele_classification: passed" << std::endl;
  return true;

} // runEvtSel_jet_ele_classification()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_ele_classification ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );
  
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::GsfElectronCollection> ele;   
  iEvent.getByToken(eleCollectionT_, ele);      

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);
  
  v_ele_goodvertices_.clear();
  v_jetIsEle_.clear();
  v_ele_jet_pt_.clear();
  v_ele_jet_eta_.clear();
  v_ele_jet_e_.clear();
  v_ele_jet_m0_.clear();
  v_ele_gen_pt_.clear();
  v_ele_gen_eta_.clear();
  v_ele_gen_pdgId_.clear();
  v_ele_genRecoJetdR_.clear();
  v_ele_NrecoEle_.clear();
  v_ele_genRecoEledR_.clear();
  v_ele_jetEledR_.clear(); 

  unsigned int goodVertices = 0;
  
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;
      
  if ( debug )  std::cout << " --------------------------------- vEle size" << vEles.size() << " --------------------------------- " << std::endl;

  v_ele_goodvertices_.push_back(goodVertices);
  
  int nRecoEle = 0; 
  for ( auto const& ii: vEles ) {
    
    // Skip electrons which fails HE edge cut
    if(ii.matchedRecoJetIdxs.empty() || ii.matchedRecoEleIdxs.empty())continue;
    if ( std::find(vJetIdxs.begin(), vJetIdxs.end(), ii.matchedRecoJetIdxs[0]) == vJetIdxs.end()) continue;
    
    reco::GenParticleRef iGen( genParticles, ii.idx );
    reco::PFJetRef iJet( jets, ii.matchedRecoJetIdxs[0] );
    reco::GsfElectronRef iEle( ele, ii.matchedRecoEleIdxs[0] );
    
    if ( debug )  std::cout << " --------------------------------- Filling branches --------------------------------- " << std::endl;
    if ( debug )  std::cout << " Gen pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;       
    if ( debug )  std::cout << " Jet pt: "<< iJet->pt() << " eta: " <<iJet->eta() << " phi: " <<iJet->phi() << std::endl;       
    if ( debug )  std::cout << " Ele pt: "<< iEle->pt() << " eta: " <<iEle->eta() << " phi: " <<iEle->phi() << std::endl;       
    
    v_jetIsEle_.push_back(1);
    v_ele_jet_pt_.push_back( iJet->pt() );
    v_ele_jet_eta_.push_back( iJet->eta() );
    v_ele_jet_e_.push_back( iJet->energy() );
    v_ele_jet_m0_.push_back( iJet->mass() );
    v_ele_gen_pt_.push_back(iGen->pt());
    v_ele_gen_eta_.push_back(iGen->eta());
    v_ele_gen_pdgId_.push_back(iGen->pdgId());
 
    TLorentzVector TLVJet(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
    TLorentzVector TLVGen(iGen->px(),iGen->py(),iGen->pz(),iGen->energy());
    TLorentzVector TLVEle(iEle->px(),iEle->py(),iEle->pz(),iEle->energy());   
       
    v_ele_genRecoJetdR_.push_back(TLVJet.DeltaR(TLVGen) );
    v_ele_genRecoEledR_.push_back(TLVEle.DeltaR(TLVGen) );
    v_ele_jetEledR_.push_back(TLVEle.DeltaR(TLVJet) );

    nRecoEle++;
  }//gen ele loop
  
  v_ele_NrecoEle_.push_back(nRecoEle);

} // fillEvtSel_jet_ele_classification()
