#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 4; //TODO: use cfg level nJets_

vector<float> v_goodvertices_;
vector<float> v_qcdJetIsEle_;
vector<float> v_qcdJetIsTau_;

vector<float> v_jet_pt_;
vector<float> v_jet_eta_;
vector<float> v_jet_e_;
vector<float> v_jet_m0_;
vector<float> v_jet_Px_;
vector<float> v_jet_Py_;
vector<float> v_jet_Pz_;

vector<float> v_genRecoJetdR_;
vector<float> v_gen_pt_;
vector<float> v_gen_eta_;
vector<float> v_gen_pdgId_;

vector<vector<float>> v_jetPFCandE_;
vector<vector<float>> v_jetPFCandPx_;
vector<vector<float>> v_jetPFCandPy_;
vector<vector<float>> v_jetPFCandPz_;
vector<vector<int>> v_jetPFCandType_;

vector<vector<float>> v_jetSV_PtRel_;
vector<vector<float>> v_jetSV_ERel_;
vector<vector<float>> v_jetSV_PhiRel_;
vector<vector<float>> v_jetSV_EtaRel_;
vector<vector<float>> v_jetSV_DeltaR_;
vector<vector<float>> v_jetSV_Pt_;
vector<vector<float>> v_jetSV_Eta_;
vector<vector<float>> v_jetSV_Phi_;
vector<vector<float>> v_jetSV_Mass_;

vector<vector<float>> v_jetSV_ntracks_;
vector<vector<float>> v_jetSV_chi2_;
vector<vector<float>> v_jetSV_ndf_;
vector<vector<float>> v_jetSV_normchi2_;

vector<vector<float>> v_jetSV_dxy_;
vector<vector<float>> v_jetSV_dxyerr_;
vector<vector<float>> v_jetSV_dxysig_;

vector<vector<float>> v_jetSV_d3d_;
vector<vector<float>> v_jetSV_d3derr_;
vector<vector<float>> v_jetSV_d3dsig_;
vector<vector<float>> v_jetSV_costhetasvpv_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_qcd ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("goodvertices",  &v_goodvertices_);
  tree->Branch("jet_IsEle",     &v_qcdJetIsEle_);
  tree->Branch("jet_IsTau",     &v_qcdJetIsTau_);
  tree->Branch("jet_Pt",        &v_jet_pt_);
  tree->Branch("jet_Eta",       &v_jet_eta_);
  tree->Branch("jet_Energy",    &v_jet_e_);
  tree->Branch("jet_M",         &v_jet_m0_);
  tree->Branch("jet_Px",        &v_jet_Px_);
  tree->Branch("jet_Py",        &v_jet_Py_);
  tree->Branch("jet_Pz",        &v_jet_Pz_); 
  
  tree->Branch("jet_GendR",     &v_genRecoJetdR_);
  tree->Branch("gen_pt",        &v_gen_pt_);
  tree->Branch("gen_eta",       &v_gen_eta_);
  tree->Branch("gen_pdgId",     &v_gen_pdgId_);
  
  tree->Branch("jetPFCandE",        &v_jetPFCandE_);
  tree->Branch("jetPFCandPx",       &v_jetPFCandPx_);
  tree->Branch("jetPFCandPy",       &v_jetPFCandPy_);
  tree->Branch("jetPFCandPz",       &v_jetPFCandPz_);
  tree->Branch("jetPFCandType",     &v_jetPFCandType_);

  tree->Branch("jetSV_PtRel",     &v_jetSV_PtRel_);
  tree->Branch("jetSV_ERel",      &v_jetSV_ERel_);
  tree->Branch("jetSV_PhiRel",    &v_jetSV_PhiRel_);
  tree->Branch("jetSV_EtaRel",  &v_jetSV_EtaRel_);
  tree->Branch("jetSV_DeltaR",  &v_jetSV_DeltaR_);
  tree->Branch("jetSV_Pt",      &v_jetSV_Pt_);
  tree->Branch("jetSV_Eta",     &v_jetSV_Eta_);
  tree->Branch("jetSV_Phi",     &v_jetSV_Phi_);
  tree->Branch("jetSV_Mass",    &v_jetSV_Mass_);

  tree->Branch("jetSV_ntracks",   &v_jetSV_ntracks_);
  tree->Branch("jetSV_chi2",      &v_jetSV_chi2_);
  tree->Branch("jetSV_ndf",       &v_jetSV_ndf_);
  tree->Branch("jetSV_normchi2",  &v_jetSV_normchi2_);
 
  tree->Branch("jetSV_dxy",       &v_jetSV_dxy_);
  tree->Branch("jetSV_dxyerr",    &v_jetSV_dxyerr_);
  tree->Branch("jetSV_dxysig",    &v_jetSV_dxysig_);

  tree->Branch("jetSV_d3d",       &v_jetSV_d3d_);
  tree->Branch("jetSV_d3derr",    &v_jetSV_d3derr_);
  tree->Branch("jetSV_d3dsig",    &v_jetSV_d3dsig_);
  tree->Branch("jetSV_costhetasvpv",  &v_jetSV_costhetasvpv_);

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
    if (abs(iGen->pdgId()) ==11 || abs(iGen->pdgId()) == 13  )continue; 
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
    
    reco::GenParticleRef iGen( genParticles, iG );
    if ( debug ) std::cout << " >> genParton[" << iG << "]" << " PDG ID: " << iGen->pdgId() << "  pt:" << iGen->pt() << " eta:" << iGen->eta() << std::endl;
    
    std::vector<unsigned int> vMatchedRecoJetIdxs;
    
    // Do dR match to closest reco photon
    minDR = 100.;
    minDR_fpt = -10;
    minDR_idx = -1;
    for ( unsigned int iJ = 0; iJ < jets->size(); iJ++ ) {
      
      reco::PFJetRef iJet( jets, iJ );
      
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      
      
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR > minDR ) continue;
      
      minDR = dR;
      minDR_idx = iJ;
      minDR_fpt = iJet->pt()/iGen->pt();      
    } // reco jets
    
    
    // Require minimum dR to declare match
    // Protects against matching to PU
    // minDR only needs to be generous enough so that one of the gen partons match to a reco jets for analysis
    if ( minDR > 0.4 ) continue;
    
    // Declare reco jet matching to gen parton: only store unique reco idxs
    if ( std::find(vMatchedRecoJetIdxs.begin(), vMatchedRecoJetIdxs.end(), minDR_idx) != vMatchedRecoJetIdxs.end() ) continue;
    vMatchedRecoJetIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << " f_pt(reco/gen):" << minDR_fpt << std::endl;


     
    // Check if matched reco jets also matches to parton:
    for ( auto& iJ : vMatchedRecoJetIdxs ) {
      
      reco::PFJetRef iJet( jets, iJ );      
      vJetIdxs.push_back( iJ );
      
      if ( debug ) std::cout << " ----> matched jet [" << iJ << "] " << " jet pt: " << iJet->pt() << " eta: " << iJet->eta() << std::endl; 
    } // matched reco jets
    
   // Store this mapping
   if(vMatchedRecoJetIdxs.empty()) continue;
   jet_gen_map iJet_obj = { iG, vMatchedRecoJetIdxs };
   vGens.push_back( iJet_obj );
   
  } //gen partons
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
  
  edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
  iEvent.getByToken(secVertexCollectionT_, secVertices);
  
  v_goodvertices_.clear();
  v_qcdJetIsEle_.clear();
  v_qcdJetIsTau_.clear();
  v_jet_pt_.clear();
  v_jet_eta_.clear();
  v_jet_e_.clear();
  v_jet_m0_.clear();
  v_jet_Px_.clear();
  v_jet_Py_.clear();
  v_jet_Pz_.clear();
  v_gen_pt_.clear();
  v_gen_eta_.clear();
  v_gen_pdgId_.clear();
  v_genRecoJetdR_.clear();
  
  v_jetPFCandE_.clear();
  v_jetPFCandPx_.clear();
  v_jetPFCandPy_.clear();
  v_jetPFCandPz_.clear();
  v_jetPFCandType_.clear();

  v_jetSV_PtRel_.clear();
  v_jetSV_ERel_.clear();
  v_jetSV_PhiRel_.clear();
  v_jetSV_EtaRel_.clear();
  v_jetSV_DeltaR_.clear();
  v_jetSV_Pt_.clear();
  v_jetSV_Eta_.clear();
  v_jetSV_Phi_.clear();
  v_jetSV_Mass_.clear();
  
  v_jetSV_ntracks_.clear();
  v_jetSV_chi2_.clear();
  v_jetSV_ndf_.clear();
  v_jetSV_normchi2_.clear();
  
  v_jetSV_dxy_.clear();
  v_jetSV_dxyerr_.clear();
  v_jetSV_dxysig_.clear();
  
  v_jetSV_d3d_.clear();
  v_jetSV_d3derr_.clear();
  v_jetSV_d3dsig_.clear();
  v_jetSV_costhetasvpv_.clear();

  unsigned int goodVertices = 0;
  
  if (vertices.isValid())
    if (vertices->size() > 0)
      for (auto v : *vertices)
        if (v.ndof() >= 4 && !v.isFake())
          ++goodVertices;
  if ( debug ) std::cout << "\t" << " good vertices in the event (PU) = " << goodVertices << std::endl;
      
  if ( debug )  std::cout << " --------------------------------- vGen size" << vGens.size() << " --------------------------------- " << std::endl;

  v_goodvertices_.push_back(goodVertices);
  
 
  int iJ = 0; 

  for ( auto const& ii: vGens ) {
    
    // Skip jets which fails HE edge cut
    if(ii.matchedRecoJetIdxs.empty())continue;
    if ( std::find(vJetIdxs.begin(), vJetIdxs.end(), ii.matchedRecoJetIdxs[0]) == vJetIdxs.end()) continue;
    
    reco::GenParticleRef iGen( genParticles, ii.idx );
    reco::PFJetRef iJet( jets, ii.matchedRecoJetIdxs[0] );
    
    iJ  = ii.matchedRecoJetIdxs[0];
       
    if ( debug )  std::cout << " --------------------------------- Filling branches --------------------------------- " << std::endl;
    if ( debug )  std::cout << " Gen pdgId: "<< iGen->pdgId() << " pt " << iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;       
    if ( debug )  std::cout << " Jet pt: "<< iJet->pt() << " eta: " <<iJet->eta() << " phi: " <<iJet->phi() << std::endl;       
    
    v_qcdJetIsEle_.push_back(0);
    v_qcdJetIsTau_.push_back(0);
    v_jet_pt_.push_back( iJet->pt() );
    v_jet_eta_.push_back( iJet->eta() );
    v_jet_e_.push_back( iJet->energy() );
    v_jet_m0_.push_back( iJet->mass() );
    v_jet_Px_.push_back( iJet->px() );
    v_jet_Py_.push_back( iJet->py() );
    v_jet_Pz_.push_back( iJet->pz() );  
    v_gen_pt_.push_back(iGen->pt());
    v_gen_eta_.push_back(iGen->eta());
    v_gen_pdgId_.push_back(iGen->pdgId());
 
    TLorentzVector TLVJet(iJet->px(),iJet->py(),iJet->pz(),iJet->energy());
    TLorentzVector TLVGen(iGen->px(),iGen->py(),iGen->pz(),iGen->energy());
       
    v_genRecoJetdR_.push_back(TLVJet.DeltaR(TLVGen) );
    
    
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
      //std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;

      pfcand_px.push_back(jetPFCand->px());
      pfcand_py.push_back(jetPFCand->py());
      pfcand_pz.push_back(jetPFCand->pz());
      pfcand_energy.push_back(jetPFCand->energy());
      pfcand_type.push_back((int)jetPFCand->particleId());
    }//jet constituents loop
    iJ++;
   
    v_jetPFCandE_.push_back( pfcand_energy );
    v_jetPFCandPx_.push_back( pfcand_px );
    v_jetPFCandPy_.push_back( pfcand_py );
    v_jetPFCandPz_.push_back( pfcand_pz );
    v_jetPFCandType_.push_back( pfcand_type );

    std::vector<const reco::VertexCompositePtrCandidate*> jetSVs;
    for (const auto &sv : *secVertices){
      if (reco::deltaR(sv, *iJet) < 0.4) {
	jetSVs.push_back(&sv);
      }
    }

    // sort by dxy significance
    const auto &pv = vertices->at(0);
    std::sort(jetSVs.begin(), jetSVs.end(), [&](const reco::VertexCompositePtrCandidate *sv1, const reco::VertexCompositePtrCandidate *sv2){
	return vertexDxy(*sv1, pv).significance() > vertexDxy(*sv2, pv).significance();
      });


    vector<float> sv_PtRel;
    vector<float> sv_ERel;
    vector<float> sv_PhiRel;
    vector<float> sv_EtaRel;
    vector<float> sv_DeltaR;
    vector<float> sv_Pt;
    vector<float> sv_Eta;
    vector<float> sv_Phi;
    vector<float> sv_Mass;
  
    vector<float> sv_ntracks;
    vector<float> sv_chi2;
    vector<float> sv_ndf;
    vector<float> sv_normchi2;
  
    vector<float> sv_dxy;
    vector<float> sv_dxyerr;
    vector<float> sv_dxysig;
  
    vector<float> sv_d3d;
    vector<float> sv_d3derr;
    vector<float> sv_d3dsig;
    vector<float> sv_costhetasvpv;
  
  
    float etasign = iJet->eta()>0 ? 1 : -1;
  
    for (const auto *sv : jetSVs){
    
      std::cout << "=================================================Secondary vertex pt:" << sv->pt() << " eta:" << sv->eta() << "  phi:" << sv->phi() << std::endl;
    
      sv_PtRel.push_back(sv->pt()/iJet->pt());
      sv_ERel.push_back(sv->energy()/iJet->energy());
      sv_PhiRel.push_back(reco::deltaPhi(*sv, *iJet));
      sv_EtaRel.push_back(etasign * (sv->eta() - iJet->eta()));
      sv_DeltaR.push_back(reco::deltaR(*sv, *iJet));
      sv_Pt.push_back(sv->pt());
      sv_Eta.push_back(sv->eta());
      sv_Phi.push_back(sv->phi());
      sv_Mass.push_back(sv->mass());
    
      //sv properties
      sv_ntracks.push_back(sv->numberOfDaughters());	
      sv_chi2.push_back(sv->vertexChi2());	
      sv_ndf.push_back(sv->vertexNdof());
      sv_normchi2.push_back(catchInfs(sv->vertexNormalizedChi2()));		
    
    
      const auto &dxy = vertexDxy(*sv, pv);
      sv_dxy.push_back(dxy.value());
      sv_dxyerr.push_back(dxy.error());
      sv_dxysig.push_back(dxy.significance());	
    
    
      const auto &d3d = vertexD3d(*sv, pv);
      sv_d3d.push_back(d3d.value());
      sv_d3derr.push_back(d3d.error());
      sv_d3dsig.push_back(d3d.significance());
      sv_costhetasvpv.push_back(vertexDdotP(*sv, pv));
    
    }
  
    v_jetSV_PtRel_.push_back(sv_PtRel);
    v_jetSV_ERel_.push_back(sv_ERel);
    v_jetSV_PhiRel_.push_back(sv_PhiRel);
    v_jetSV_EtaRel_.push_back(sv_EtaRel);
    v_jetSV_DeltaR_.push_back(sv_DeltaR);
    v_jetSV_Pt_.push_back(sv_Pt);
    v_jetSV_Eta_.push_back(sv_Eta);
    v_jetSV_Phi_.push_back(sv_Phi);
    v_jetSV_Mass_.push_back(sv_Mass);
  
    v_jetSV_ntracks_.push_back(sv_ntracks);
    v_jetSV_chi2_.push_back(sv_chi2);
    v_jetSV_ndf_.push_back(sv_ndf);
    v_jetSV_normchi2_.push_back(sv_normchi2);
  
    v_jetSV_dxy_.push_back(sv_dxy);
    v_jetSV_dxyerr_.push_back(sv_dxyerr);
    v_jetSV_dxysig_.push_back(sv_dxysig);
  
    v_jetSV_d3d_.push_back(sv_d3d);
    v_jetSV_d3derr_.push_back(sv_d3derr);
    v_jetSV_d3dsig_.push_back(sv_d3dsig);
    v_jetSV_costhetasvpv_.push_back(sv_costhetasvpv);
  }//gen loop
  


} // fillEvtSel_jet_qcd()
