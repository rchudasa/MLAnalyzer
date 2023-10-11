#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

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

vector<vector<float>> v_tau_jetSV_PtRel_;
vector<vector<float>> v_tau_jetSV_ERel_;
vector<vector<float>> v_tau_jetSV_PhiRel_;
vector<vector<float>> v_tau_jetSV_EtaRel_;
vector<vector<float>> v_tau_jetSV_DeltaR_;
vector<vector<float>> v_tau_jetSV_Pt_;
vector<vector<float>> v_tau_jetSV_Eta_;
vector<vector<float>> v_tau_jetSV_Phi_;
vector<vector<float>> v_tau_jetSV_Mass_;

vector<vector<float>> v_tau_jetSV_ntracks_;
vector<vector<float>> v_tau_jetSV_chi2_;
vector<vector<float>> v_tau_jetSV_ndf_;
vector<vector<float>> v_tau_jetSV_normchi2_;

vector<vector<float>> v_tau_jetSV_dxy_;
vector<vector<float>> v_tau_jetSV_dxyerr_;
vector<vector<float>> v_tau_jetSV_dxysig_;

vector<vector<float>> v_tau_jetSV_d3d_;
vector<vector<float>> v_tau_jetSV_d3derr_;
vector<vector<float>> v_tau_jetSV_d3dsig_;
vector<vector<float>> v_tau_jetSV_costhetasvpv_;

vector<float> v_tau_deeptau_;
vector<float> v_tau_loosedeeptau_;
vector<float> v_tau_mediumdeeptau_;
vector<float> v_tau_tightdeeptau_;



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

  tree->Branch("jetSV_PtRel",     &v_tau_jetSV_PtRel_);
  tree->Branch("jetSV_ERel",      &v_tau_jetSV_ERel_);
  tree->Branch("jetSV_PhiRel",    &v_tau_jetSV_PhiRel_);
  tree->Branch("jetSV_EtaRel",  &v_tau_jetSV_EtaRel_);
  tree->Branch("jetSV_DeltaR",  &v_tau_jetSV_DeltaR_);
  tree->Branch("jetSV_Pt",      &v_tau_jetSV_Pt_);
  tree->Branch("jetSV_Eta",     &v_tau_jetSV_Eta_);
  tree->Branch("jetSV_Phi",     &v_tau_jetSV_Phi_);
  tree->Branch("jetSV_Mass",    &v_tau_jetSV_Mass_);

  tree->Branch("jetSV_ntracks",   &v_tau_jetSV_ntracks_);
  tree->Branch("jetSV_chi2",      &v_tau_jetSV_chi2_);
  tree->Branch("jetSV_ndf",       &v_tau_jetSV_ndf_);
  tree->Branch("jetSV_normchi2",  &v_tau_jetSV_normchi2_);
 
  tree->Branch("jetSV_dxy",       &v_tau_jetSV_dxy_);
  tree->Branch("jetSV_dxyerr",    &v_tau_jetSV_dxyerr_);
  tree->Branch("jetSV_dxysig",    &v_tau_jetSV_dxysig_);

  tree->Branch("jetSV_d3d",       &v_tau_jetSV_d3d_);
  tree->Branch("jetSV_d3derr",    &v_tau_jetSV_d3derr_);
  tree->Branch("jetSV_d3dsig",    &v_tau_jetSV_d3dsig_);
  tree->Branch("jetSV_costhetasvpv",  &v_tau_jetSV_costhetasvpv_);

  tree->Branch("DeepTau",       &v_tau_deeptau_);
  tree->Branch("LooseDeepTau",  &v_tau_loosedeeptau_);
  tree->Branch("MediumDeepTau", &v_tau_mediumdeeptau_);
  tree->Branch("TightDeepTau",  &v_tau_tightdeeptau_);

} // branchesEvtSel_jet_dijet_tau()

float catchInfs(const float& in, float replace_value=0){
    if(std::isinf(in) || std::isnan(in))
      return replace_value;
    else if(in < -1e32 || in > 1e32)
      return replace_value;
    return in;
}

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

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<pat::TauCollection> taus;
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
   
    if ( debug ) std::cout << "\t\t\t" << " GEN particle " <<  ", status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << ", nMoms: " <<iGen->numberOfMothers() << ", mother ID: " << iGen->mother()->pdgId() << ", pt: "<< iGen->pt() << ", eta: " <<iGen->eta() << ", phi: " <<iGen->phi() << std::endl; 
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
      
      pat::Jet iJet = (*jets)[iJ];
      
      if ( std::abs(iJet.pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet.eta()) > maxJetEta_ ) continue;
      
      
      dR = reco::deltaR( iJet.eta(),iJet.phi(), iGenTau->eta(),iGenTau->phi() );
      if ( dR > minDR ) continue;
      
      minDR = dR;
      minDR_idx = iJ;
      minDR_fpt = iJet.pt()/iGenTau->pt();
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << " pt:" << iJet.pt() << " eta:" << iJet.eta() << std::endl;
      
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

  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vertexCollectionT_, vertices);

  edm::Handle<reco::VertexCompositePtrCandidateCollection> secVertices;
  iEvent.getByToken(secVertexCollectionT_, secVertices);

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

  v_tau_jetSV_PtRel_.clear();
  v_tau_jetSV_ERel_.clear();
  v_tau_jetSV_PhiRel_.clear();
  v_tau_jetSV_EtaRel_.clear();
  v_tau_jetSV_DeltaR_.clear();
  v_tau_jetSV_Pt_.clear();
  v_tau_jetSV_Eta_.clear();
  v_tau_jetSV_Phi_.clear();
  v_tau_jetSV_Mass_.clear();
  
  v_tau_jetSV_ntracks_.clear();
  v_tau_jetSV_chi2_.clear();
  v_tau_jetSV_ndf_.clear();
  v_tau_jetSV_normchi2_.clear();
  
  v_tau_jetSV_dxy_.clear();
  v_tau_jetSV_dxyerr_.clear();
  v_tau_jetSV_dxysig_.clear();
  
  v_tau_jetSV_d3d_.clear();
  v_tau_jetSV_d3derr_.clear();
  v_tau_jetSV_d3dsig_.clear();
  v_tau_jetSV_costhetasvpv_.clear();
 
  v_tau_deeptau_.clear();
  v_tau_loosedeeptau_.clear();
  v_tau_mediumdeeptau_.clear();
  v_tau_tightdeeptau_.clear();
 
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

    pat::Jet iJet = ( *jets)[ii.matchedRecoJetIdxs[0]];
    
    iJ  = ii.matchedRecoJetIdxs[0];
    
    if ( debug )  std::cout << " --------------------------------- Filling branches --------------------------------- " << std::endl;
    if ( debug )  std::cout << " Gen pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;       
    if ( debug )  std::cout << " Jet pt: "<< iJet.pt() << " eta: " <<iJet.eta() << " phi: " <<iJet.phi() << std::endl;       
     
    
    // Fill branches 
    v_tau_jet_pt_.push_back( iJet.pt() );
    v_tau_jet_eta_.push_back( iJet.eta() );
    v_tau_jet_e_.push_back( iJet.energy() ); 
    v_tau_jet_m0_.push_back( iJet.mass() );
    v_tau_jet_Px_.push_back( iJet.px() );
    v_tau_jet_Py_.push_back( iJet.py() );
    v_tau_jet_Pz_.push_back( iJet.pz() );    
    
    TLorentzVector TLVJet(iJet.px(),iJet.py(),iJet.pz(),iJet.energy());
    TLorentzVector TLVGen(iGen->px(),iGen->py(),iGen->pz(),iGen->energy());
    
    v_tau_genRecoJetdR_.push_back(TLVJet.DeltaR(TLVGen) );
    
    v_tau_gen_pt_.push_back(iGen->pt());
    v_tau_gen_eta_.push_back(iGen->eta());
    v_tau_gen_pdgId_.push_back(iGen->pdgId());
    
    int tauDaughters          = 0;
    int tauPi0 = 0;
    bool JetIsTau             = false;

    if ( std::abs(iGen->pdgId()) == 15 ) {
      JetIsTau = true;
      for (unsigned int iDaughter = 0; iDaughter != iGen->numberOfDaughters(); ++iDaughter ){
	if ( debug ) std::cout << "\t\t\t\t" <<" Tau daughter [" << iDaughter << "] : "<<  std::abs(iGen->daughter(iDaughter)->pdgId()); 
	if ( debug ) std::cout << " charge : "<< iGen->daughter(iDaughter)->charge() << "  | pt : "<< iGen->daughter(iDaughter)->pt();
	if ( debug ) std::cout << " eta:" << iGen->daughter(iDaughter)->eta() << " |Energy:" << iGen->daughter(iDaughter)->energy() << std::endl;
	if ( abs(iGen->daughter(iDaughter)->pdgId()) == 111 ) tauPi0++;
	if ( iGen->daughter(iDaughter)->charge() == 0 ) continue;          
	tauDaughters++;
      }
    }
    v_tau_jetIsTau_.push_back(JetIsTau);
    
    
    if ( debug ) std::cout << "\t\t\t"<<" Tau prongs = " << tauDaughters << " + Tau pi0 = " << tauPi0 << std::endl;      
    v_tau_gen_prongs_.push_back(tauDaughters);
    
    // jet constituents
    vector<float> pfcand_px;
    vector<float> pfcand_py;
    vector<float> pfcand_pz;
    vector<float> pfcand_energy;
    vector<int> pfcand_type;
    
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet.getPFConstituents();
    //unsigned int nConstituents = iJet.getPFConstituents().size();
    const auto& pfCands = iJet.getJetConstituents();
    
    //std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    //if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    //for ( unsigned int j = 0; j < nConstituents; j++ ) {
    //const reco::PFCandidatePtr jetPFCand = iJet.getPFConstituent( j );
    //if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;
    //std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << jetPFCand->energy() << " px:" << jetPFCand->px() << " py:" << jetPFCand->py() << " pz:" << jetPFCand->pz() << std::endl;
    
    for (const auto& cand : pfCands) {
      const auto* jetPFCand = dynamic_cast<const pat::PackedCandidate*>(&(*cand));
      
      pfcand_px.push_back(jetPFCand->px());
      pfcand_py.push_back(jetPFCand->py());
      pfcand_pz.push_back(jetPFCand->pz());
      pfcand_energy.push_back(jetPFCand->energy());
      pfcand_type.push_back((int)jetPFCand->pdgId());
    }//jet constituents loop
    iJ++;
   
    v_tau_jetPFCandE_.push_back( pfcand_energy );
    v_tau_jetPFCandPx_.push_back( pfcand_px );
    v_tau_jetPFCandPy_.push_back( pfcand_py );
    v_tau_jetPFCandPz_.push_back( pfcand_pz );
    v_tau_jetPFCandType_.push_back( pfcand_type );
    
    std::vector<const reco::VertexCompositePtrCandidate*> jetSVs;
    for (const auto &sv : *secVertices){
      if (reco::deltaR(sv, iJet) < 0.4) {
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
    
    
    float etasign = iJet.eta()>0 ? 1 : -1;
    
    for (const auto *sv : jetSVs){
      
      std::cout << "=================================================Secondary vertex pt:" << sv->pt() << " eta:" << sv->eta() << "  phi:" << sv->phi() << std::endl;
      
      sv_PtRel.push_back(sv->pt()/iJet.pt());
      sv_ERel.push_back(sv->energy()/iJet.energy());
      sv_PhiRel.push_back(reco::deltaPhi(*sv, iJet));
      sv_EtaRel.push_back(etasign * (sv->eta() - iJet.eta()));
      sv_DeltaR.push_back(reco::deltaR(*sv, iJet));
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
    
    v_tau_jetSV_PtRel_.push_back(sv_PtRel);
    v_tau_jetSV_ERel_.push_back(sv_ERel);
    v_tau_jetSV_PhiRel_.push_back(sv_PhiRel);
    v_tau_jetSV_EtaRel_.push_back(sv_EtaRel);
    v_tau_jetSV_DeltaR_.push_back(sv_DeltaR);
    v_tau_jetSV_Pt_.push_back(sv_Pt);
    v_tau_jetSV_Eta_.push_back(sv_Eta);
    v_tau_jetSV_Phi_.push_back(sv_Phi);
    v_tau_jetSV_Mass_.push_back(sv_Mass);
    
    v_tau_jetSV_ntracks_.push_back(sv_ntracks);
    v_tau_jetSV_chi2_.push_back(sv_chi2);
    v_tau_jetSV_ndf_.push_back(sv_ndf);
    v_tau_jetSV_normchi2_.push_back(sv_normchi2);
    
    v_tau_jetSV_dxy_.push_back(sv_dxy);
    v_tau_jetSV_dxyerr_.push_back(sv_dxyerr);
    v_tau_jetSV_dxysig_.push_back(sv_dxysig);
    
    v_tau_jetSV_d3d_.push_back(sv_d3d);
    v_tau_jetSV_d3derr_.push_back(sv_d3derr);
    v_tau_jetSV_d3dsig_.push_back(sv_d3dsig);
    v_tau_jetSV_costhetasvpv_.push_back(sv_costhetasvpv);
    
    float DeepTau             = -1;
    float LooseDeepTau        = -1;
    float MediumDeepTau       = -1;
    float TightDeepTau        = -1; 
   
      if ( isSignal_ ){
        float recoTaudR = 0.4;
        for (const pat::Tau &tau : *taus) {
          if(tau.pt() < 18.0) continue;
          float TaudR = reco::deltaR( tau.eta(),tau.phi(), iJet.eta(), iJet.phi() );
          if ( TaudR > recoTaudR ) continue;
          recoTaudR = TaudR;
          DeepTau       = tau.tauID("byDeepTau2017v2p1VSjetraw");
          LooseDeepTau  = tau.tauID("byLooseDeepTau2017v2p1VSjet");
          MediumDeepTau = tau.tauID("byMediumDeepTau2017v2p1VSjet");
          TightDeepTau  = tau.tauID("byTightDeepTau2017v2p1VSjet");
        }
        //if ( debug ) std::cout << "   DeepTau  = " << DeepTau << std::endl;
        std::cout << "   DeepTau  = " << DeepTau << std::endl;
      } //if signal and not training

        v_tau_deeptau_.push_back( DeepTau );
        v_tau_loosedeeptau_.push_back( LooseDeepTau );
        v_tau_mediumdeeptau_.push_back( MediumDeepTau );
        v_tau_tightdeeptau_.push_back( TightDeepTau );
 
  }//taujets loop
  
} // fillEvtSel_jet_dijet_tau()

