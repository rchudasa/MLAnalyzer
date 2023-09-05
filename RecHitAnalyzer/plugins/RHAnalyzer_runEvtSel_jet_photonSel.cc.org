#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

vector<float> v_PhoGoodvertices_;
vector<float> v_PhoIsEle_;
vector<float> v_pho_pt_;
vector<float> v_pho_eta_;
vector<float> v_pho_e_;
vector<float> v_pho_m0_;
vector<float> v_genRecoPhodR_;
vector<float> v_genPho_pt_;
vector<float> v_genPho_eta_;
vector<float> v_genPho_pdgId_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_photonSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("goodvertices",  &v_PhoGoodvertices_);
  tree->Branch("jet_IsEle",     &v_PhoIsEle_);
  tree->Branch("jet_Pt",        &v_pho_pt_);
  tree->Branch("jet_Eta",       &v_pho_eta_);
  tree->Branch("jet_Energy",    &v_pho_e_);
  tree->Branch("jet_M",         &v_pho_m0_);
  tree->Branch("jet_GendR",     &v_genRecoPhodR_);
  tree->Branch("gen_pt",        &v_genPho_pt_);
  tree->Branch("gen_eta",       &v_genPho_eta_);
  tree->Branch("gen_pdgId",     &v_genPho_pdgId_);

}

// Define struct to handle mapping for gen pho<->matched reco photons<->matched presel photons
struct pho_map {
  unsigned int idx;
  std::vector<unsigned int> matchedRecoPhoIdxs;
  std::vector<unsigned int> matchedRecoJetPhoIdxs;
};
std::vector<pho_map> vGamma;

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_photonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::Photon> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  ////////// Gen-level validation //////////

  // Identify particle gen photons from the gamma+jet events
  float dR;
  std::vector<unsigned int> vGenPhoIdxs;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );

    if ( iGen->pdgId() != 22 ) continue; //photon
    if ( iGen->status() != 23 ) continue;

    vGenPhoIdxs.push_back( iG );

  } // genParticles
  if ( debug ) std::cout << " >> vGenPhoIdxs.size: " << vGenPhoIdxs.size() << std::endl;
  if ( vGenPhoIdxs.empty() ) return false;

  ////////// Build gen pho-reco photon mapping /////////
  // Create mapping between gen pho<->matched reco photons<->matched jet
  // For each gen pho, find "reco" photons matched to it,
  // then check if that reco photon passes photon preselection criteria
  
  int phoIdx = -1;
  float minDR = 100.;
  float minDR_fpt = -10.;
  int minDR_idx = -1;
  
  float minJetDR = 100.;
  int minJetDR_idx = -1;
  
  vGamma.clear();
  // Loop over valid gen pho idxs
  for ( auto& iG : vGenPhoIdxs ) {

    reco::GenParticleRef iGenPho( genParticles, iG );
    if ( debug ) std::cout << " >> genPho[" << iG << "]" << " pt:" << iGenPho->pt() << " eta:" << iGenPho->eta() << std::endl;

    std::vector<unsigned int> vMatchedRecoPhoIdxs;
    std::vector<unsigned int> vMatchedJetPhoIdxs;

    // Do dR match to closest reco photon
    minDR = 100.;
    minDR_fpt = -10;
    minDR_idx = -1;
    phoIdx = -1;
    //for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {
    for (auto iP = photons->begin(); iP != photons->end(); ++iP) {

	phoIdx++;

      // Definition of a "reco" photon--highly subject to interpretation
      if ( iP->pt() < minJetPt_ ) continue; 
      if ( std::abs(iP->eta()) > maxJetEta_ ) continue; 

      dR = reco::deltaR( iP->eta(),iP->phi(), iGenPho->eta(),iGenPho->phi() );
      if ( dR > minDR ) continue;

      minDR = dR;
      minDR_idx = phoIdx;
      minDR_fpt = iP->pt()/iGenPho->pt();
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << " pt:" << iP->pt() << " eta:" << iP->eta() << std::endl;

    } // reco photons


    // Require minimum dR to declare match
    // Protects against matching to PU, although not a major issue since these will likely fail preselection
    // minDR only needs to be generous enough so that one of the gen photons match to a reco photon for analysis
    if ( minDR > 0.04 ) continue;

    // Declare reco photon matching to gen pho: only store unique reco idxs
    if ( std::find(vMatchedRecoPhoIdxs.begin(), vMatchedRecoPhoIdxs.end(), minDR_idx) != vMatchedRecoPhoIdxs.end() ) continue;
    vMatchedRecoPhoIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << " f_pt(reco/gen):" << minDR_fpt << std::endl;

   /* minJetDR = 100.;
    minJetDR_idx = -1;

    for ( auto& iP : vMatchedRecoPhoIdxs ) {
 
      reco::PhotonRef iRecoPho( photons, iP );
      
      //Lookin at RecoJet
      
      for ( unsigned int iJ = 0; iJ < jets->size(); iJ++ ) {
        reco::PFJetRef iJet( jets, iJ );    
        
         if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
         if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
         
         dR = reco::deltaR( iJet->eta(),iJet->phi(), iP->eta(),iP->phi() ); 
         if ( dR > minJetDR ) continue;
           minJetDR = dR;
           minJetDR_idx = iJ;
           if ( debug ) std::cout << "   >> minJetDR_idx:" << minJetDR_idx << " " << minDR << " pt:" << iJet->pt() << " eta:" << iJet->eta() << std::endl;             
      } //jets
    
      if(minJetDR>0.4) continue;
      
      if ( debug ) std::cout << "   >> minJetDR_idx:" << minJetDR_idx << " " << minDR << " pt:" << iJet->pt() << " eta:" << iJet->eta() << std::endl; 
      
      vJetIdxs.push_back( iP );
      vMatchedJetPhoIdxs.push_back( minJetDR_idx );
      
    } // matched photon loop


    // Store this mapping
    if(vMatchedRecoPhoIdxs.empty())continue;
    pho_map iPho_obj = { iG, vMatchedRecoPhoIdxs, vMatchedJetPhoIdxs };
    vGamma.push_back( iPho_obj );*/
    
    if ( debug ) std::cout << " >> Matched gen-photon size: " << vMatchedRecoPhoIdxs.size() << std::endl;

  } // gen phos



  if(vJetIdxs.empty()) return false;

  // Photon gets passed to cropping routine.
  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  return true;

} // runPhotonSel()

// Fill branches ___________________________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_photonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  edm::Handle<reco::Photon> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);

  v_PhoGoodvertices_.clear();
  v_PhoIsEle_.clear();
  v_pho_pt_.clear();
  v_pho_eta_.clear();
  v_pho_e_.clear();
  v_pho_m0_.clear();
  v_genPho_pt_.clear();
  v_genPho_eta_.clear();
  v_genPho_pdgId_.clear();
  v_genRecoPhodR_.clear();

  unsigned int goodVertices = 0;
  
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;
      
  if ( debug )  std::cout << " --------------------------------- vGen size" << vGamma.size() << " --------------------------------- " << std::endl;
   
   v_PhoGoodvertices_.push_back(goodVertices); 
  
  ////////// Store gen-level pho kinematics //////////

  std::vector<unsigned int> matchedRecoPhoIdxs;
  std::vector<unsigned int> matchedRecoJetPhoIdxs;
  
  
  /*for ( auto const& iPho : vGamma ) {

    // Skip phos which are not valid for regression
    if ( iPho.matchedRecoPhoIdxs.empty() ) continue;

    //if matched preselcted photon does not pass zero suppression ( they are:regressed photon IDx) then continue
    if ( std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), iPho.matchedRecoPhoIdxs[0]) == vRegressPhoIdxs_.end() ) continue;

    reco::GenParticleRef iGen( genParticles, iPho.idx );
    reco::PhotonRef iRecoPho( photons, iPho.matchedRecoPhoIdxs[0] );
    
    v_PhoIsEle_.push_back(0);
    v_pho_pt_.push_back( iRecoPho->pt() );
    v_pho_eta_.push_back( iRecoPho->eta() );
    v_pho_e_.push_back( iRecoPho->energy() );
    v_pho_m0_.push_back( iRecoPho->mass() );
    v_genPho_pt_.push_back(iGen->pt());
    v_genPho_eta_.push_back(iGen->eta());
    v_genPho_pdgId_.push_back(iGen->pdgId());
 
    TLorentzVector TLVPho(iRecoPho->px(),iRecoPho->py(),iRecoPho->pz(),iRecoPho->energy());
    TLorentzVector TLVGen(iGen->px(),iGen->py(),iGen->pz(),iGen->energy());
       
    v_genRecoPhodR_.push_back(TLVPho.DeltaR(TLVGen) );
    
  } // gen phos*/

} // fillPhotonSel()
