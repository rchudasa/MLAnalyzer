#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <algorithm>

using std::vector;
using namespace trigger;

const unsigned nJets      = 50; //TODO: use cfg level nJets_

int  nKinTau_;

TH1F * hNpassed_kin;
TH1F * hNpassed_MVA;
TH1F * hNpassed_mGG;
TH1F * hNpassed_nRecoPho;
TH1F * hNpassed_img;

//gen variables
vector<float> v_att_genHiggs_M_;
vector<float> v_att_genPS_M_;
vector<float> v_att_genTau_pT_;

//jet variables
int v_att_tau_njet_;
vector<float> v_att_tau_jet_m0_;
vector<float> v_att_tau_jet_pt_;
vector<float> v_att_tau_jetIsSignal_;
vector<float> v_att_tau_jetdR_;

//single tau variables
int v_att_tau_ntau_;

vector<float> v_att_tau_pT_;
vector<float> v_att_tau_mva_;
vector<float> v_att_tau_dm_;

//ditau variables
int v_att_tau_combs_;
vector<float> v_att_dRtautau_;
vector<float> v_att_Mvis_;
vector<float> v_att_pTvis_;
vector<float> v_att_Mtautau_;
vector<float> v_att_Mtautau_svFit_;
vector<float> v_att_Mth_svFit_;
vector<float> v_att_pTh_svFit_;
vector<float> v_att_pTh_;
vector<float> v_att_mth_;
vector<float> v_att_mcoll_;
vector<float> v_att_pfMET_;
vector<float> v_att_dphillmet_;
vector<float> v_att_dphill_;


TLorentzVector SetTau(Float_t tau_pt, Float_t tau_eta, Float_t tau_phi, Float_t tau_mass){
  TLorentzVector Tau_Candidate;
  Tau_Candidate.SetPtEtaPhiM(tau_pt, tau_eta, tau_phi, tau_mass);
  return Tau_Candidate;
}

TLorentzVector SetMET(Float_t met, Float_t metphi){
  TLorentzVector Met;
  Met.SetPtEtaPhiM(met, 0, metphi, 0.);
  return Met;
}

struct jet_tau_obj {
  unsigned int recoJetIdxs;
  unsigned int recoTauIdxs;
};

struct jet_tau_allCut_obj {
  unsigned int jetIdxs;
  unsigned int tauIdxs;
};
std::vector<jet_tau_allCut_obj> vJetTauCut; 
std::vector<jet_tau_allCut_obj> vJetTauFrameCropped; 

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_ditau ( TTree* tree, edm::Service<TFileService> &fs ) {

  tree->Branch("nKinTau",        &nKinTau_);
  //gen variables
  tree->Branch("GenHiggsM",        &v_att_genHiggs_M_);
  tree->Branch("GenPSM",           &v_att_genPS_M_);
  tree->Branch("GenTaupT",         &v_att_genTau_pT_);
  //tau variables
  tree->Branch("nJets",            &v_att_tau_njet_);
  tree->Branch("jetM",             &v_att_tau_jet_m0_);
  tree->Branch("jetIsSignal",      &v_att_tau_jetIsSignal_);
  tree->Branch("jetpT",            &v_att_tau_jet_pt_);
  tree->Branch("dR",               &v_att_tau_jetdR_);
  //single tau variables 
  tree->Branch("nTaus",            &v_att_tau_ntau_);
  tree->Branch("TaupT",            &v_att_tau_pT_);
  tree->Branch("TauMVA",           &v_att_tau_mva_);
  tree->Branch("TauDM",            &v_att_tau_dm_);
  //ditau variables
  tree->Branch("TauPairs",         &v_att_tau_combs_);
  tree->Branch("TaudR",            &v_att_dRtautau_);
  tree->Branch("TauMvis",          &v_att_Mvis_);
  tree->Branch("TaupTvis",         &v_att_pTvis_);
  tree->Branch("TauMtautau",       &v_att_Mtautau_);
  tree->Branch("TauMtautau_svFit", &v_att_Mtautau_svFit_);
  tree->Branch("TauMth_svFit",     &v_att_Mth_svFit_);
  tree->Branch("TaupTh_svFit",     &v_att_pTh_svFit_);  
  tree->Branch("TaupTh",           &v_att_pTh_);
  tree->Branch("Taumth",           &v_att_mth_);
  tree->Branch("Taumcoll",         &v_att_mcoll_);  
  tree->Branch("pfMET",            &v_att_pfMET_);
  tree->Branch("dphillmet",        &v_att_dphillmet_);
  tree->Branch("dphill",           &v_att_dphill_);
  

  hNpassed_kin      = fs->make<TH1F>("hNpassed_kin", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_MVA      = fs->make<TH1F>("hNpassed_MVA", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_mGG      = fs->make<TH1F>("hNpassed_mGG", "isPassed;isPassed;N", 2, 0., 2);
  hNpassed_nRecoPho = fs->make<TH1F>("hNpassed_nRecoPho", "isPassed;isPassed;N", 2, 0., 2);
  
  
} // branchesEvtSel_jet_dijet_ditau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_ditau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  edm::Handle<edm::View<reco::Jet> > Jets;
  iEvent.getByToken( pfjetsToken_, Jets );
  
  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);
  
  edm::Handle<reco::PFMETCollection> pfmet;
  iEvent.getByToken(metCollectionT_, pfmet);
  
  edm::Handle<reco::PFTauDiscriminator> DecayMode;
  iEvent.getByToken(tauDecayMode_, DecayMode);
  edm::Handle<reco::PFTauDiscriminator> MVAIsolation;
  iEvent.getByToken(tauMVAIsolation_, MVAIsolation);
  edm::Handle<reco::PFTauDiscriminator> MuonRejection;
  iEvent.getByToken(tauMuonRejection_, MuonRejection);
  edm::Handle<reco::PFTauDiscriminator> ElectronRejectionMVA6;
  iEvent.getByToken(tauElectronRejectionMVA6_, ElectronRejectionMVA6);
  
  JME::JetResolution resPtObj            = JME::JetResolution::get(iSetup, jetResPtType_);
  JME::JetResolution resPhiObj           = JME::JetResolution::get(iSetup, jetResPhiType_);
  JME::JetResolutionScaleFactor resSFObj = JME::JetResolutionScaleFactor::get(iSetup, jetSFType_);
  
  edm::Handle<PFCollection> pfCandsH_;
  iEvent.getByToken( pfCollectionT_, pfCandsH_ );

  edm::Handle<reco::CandidateView> pfCandidates;
  iEvent.getByToken( pfCandidatesToken_, pfCandidates );

  std::vector< edm::Handle<reco::CandidateView> > leptons;
  for ( std::vector<edm::EDGetTokenT<edm::View<reco::Candidate> > >::const_iterator srcLeptons_i = lepTokens_.begin(); srcLeptons_i != lepTokens_.end(); ++srcLeptons_i ) {
    edm::Handle<reco::CandidateView> leptons_i;
    iEvent.getByToken(*srcLeptons_i, leptons_i);
    leptons.push_back( leptons_i );
  }
  
  edm::Handle<double> rho;
  iEvent.getByToken(rhoLabel_, rho);
  
  vJetIdxs.clear();
  
  
  bool IsSignal             = true;
  //bool IsMC                 = false;
  //if (iEvent.isRealData()) IsMC = false;
  //TAU SELECTION
  float tau_sel_mvaID       = -2;
  float tau_sel_pT          = 20;
  
  //vDiTaus.clear();


  ////////// Apply selection //////////
  
  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  //if ( debug ) std::cout << " Tau collection size:" << taus->size() << std::endl;
  
  std::vector<unsigned int> vMatchedRecoJetIdxs;
  std::vector<unsigned int> vMatchedRecoTauIdxs;
     
  
  // Loop over jets
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
     
    if ( iJet->pt() < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
     
    //count number of "reco" taus
     
    //Lookin at RecoTaus passing pt, eta cut and ateast 2 taus
    for ( unsigned iT1(0); iT1 != taus->size(); ++iT1 ) {
      reco::PFTauRef iTau1( taus, iT1 );
       
      if ( iTau1->pt() < tau_sel_pT ) continue;
      if ( abs(iTau1->eta()) >= 2.4 ) continue;
       
      float jet_taudR = reco::deltaR( iJet->eta(),iJet->phi(), iTau1->eta(),iTau1->phi() );
       
      if( jet_taudR > 0.4) continue;
       
      std ::cout<< "Jet matched with tau jetID " << iJ << " tau ID:" << iT1 << " tau-Jet dR: " << jet_taudR << " jet pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi();
      std:: cout<< " tau pT: " << iTau1->pt() << " eta:" << iTau1->eta() << " phi "<< iTau1->phi() << std::endl;

      vMatchedRecoJetIdxs.push_back(iJ);
      vMatchedRecoTauIdxs.push_back(iT1);
    }
  }//jet loop
   
  if ( vMatchedRecoTauIdxs.size() < 2 ) {
    hNpassed_kin->Fill(0);
    return false;
  }
  nKinTau_ = vMatchedRecoTauIdxs.size();
   
  //std::cout << "number of taus after pt cut:" << vMatchedRecoTauIdxs.size() << std::endl;
  hNpassed_kin->Fill(1.);
   
   
  std::vector<jet_tau_obj> vJetTau;
   
  vJetTau.clear();
   
  // Ensure two taus passing basic MVA identification cuts such as electron & muon rejection   
  for ( unsigned int j = 0; j < vMatchedRecoJetIdxs.size(); j++ ) { 
    //std::cout << " jet ID:" << vMatchedRecoJetIdxs[j] << "   tau ID:" << vMatchedRecoTauIdxs[j]<< std::endl;
     
    reco::PFTauRef iTau1( taus, vMatchedRecoTauIdxs[j] );
    if (!((*MuonRejection)[iTau1])) continue; 
    if (!((*ElectronRejectionMVA6)[iTau1])) continue;
    if ( (*MVAIsolation)[iTau1] < tau_sel_mvaID ) continue;  
     
    jet_tau_obj Jet_tau_obj = { vMatchedRecoJetIdxs[j],vMatchedRecoTauIdxs[j] };
    vJetTau.push_back( Jet_tau_obj );
     
  } //Tau MVA
   
  //if ( debug ) std::cout << " Taus passing MVA cut:" << vJetTau.size() << std::endl;
  if ( vJetTau.size() < 2 ) {
    hNpassed_MVA->Fill(0);
    return false;
  }
  hNpassed_MVA->Fill(1);
   


  std::vector<std::pair<unsigned int, unsigned int>> vecOf_jetTauPair; 
     
  vecOf_jetTauPair.clear();
   
  // apply mass cuts on combination of two taus
  for ( unsigned int j = 0; j < vJetTau.size()-1; j++ ) {
    reco::PFTauRef iTau1( taus, vJetTau[j].recoTauIdxs );
    reco::PFJetRef iJet1( jets, vJetTau[j].recoJetIdxs );
     
    for ( unsigned int k = j+1; k < vJetTau.size(); k++ ) {
     
      if ( k <= j ) continue;
      reco::PFTauRef iTau2( taus, vJetTau[k].recoTauIdxs );
      reco::PFJetRef iJet2( jets, vJetTau[k].recoJetIdxs );
       
      float recotaudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iTau2->eta(),iTau2->phi() );
       
      if(recotaudR > 0.5){
	 
	std::cout << "================         In Combination loop        ====================" << std::endl;
	std::cout<< " j: "<< j << " Tau1 pt:" << iTau1->pt() << "  eta:" << iTau1->eta() << "  k:" << k << "  tau2 pt: " << iTau2->pt() << "  eta:" << iTau2->eta() << std::endl;
	 
	TLorentzVector Tau1  = SetTau(iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()); 
	TLorentzVector Tau2  = SetTau(iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass());
	TLorentzVector DiTau = Tau1+Tau2;
	float dphill  = abs(Tau1.DeltaPhi(Tau2));  
	float diMvis  = DiTau.M();
	float dipTvis = DiTau.Pt();
	//if ( debug ) std::cout << " Tau pair " << " ("<< iTau1->pt() << " GeV) + " << " (" << iTau2->pt() << " GeV) | ditau Mvis : " << diMvis << " GeV | dR = " << taudR << std::endl;
	 
	 
	float pfMET    = (pfmet->front()).et();
	float pfMETphi = (pfmet->front()).phi();
	TLorentzVector MET = SetMET(pfMET,pfMETphi);
	//if ( debug ) std::cout << " PF MET = " << pfMET << std::endl;
	 
	//define MET for SVFit
	double measuredMETx = (pfmet->front()).px();
	double measuredMETy = (pfmet->front()).py();
	//if ( debug ) std::cout << " PF MET X = " << measuredMETx << " | PF MET Y = " << measuredMETy << std::endl;
	float dphillmet = abs(DiTau.DeltaPhi(MET));  
	float ditau_mth = sqrt( 2. * dipTvis * pfMET * ( 1. - cos (dphillmet) ));
	 
	TLorentzVector Higgs = Tau1+Tau2+MET;
	float ditau_M   = Higgs.M();
	float ditau_pT  = Higgs.Pt();
	 
	float the1    = 2.*atan(exp(-iTau1->eta()));
	float the2    = 2.*atan(exp(-iTau2->eta()));
	float pTmiss1 = ( (sin(iTau1->phi())*measuredMETx) - (cos(iTau1->phi())*measuredMETy) ) / (sin(the2)*( (sin(iTau1->phi())*cos(iTau2->phi())) - (cos(iTau1->phi())*sin(iTau2->phi())) ) );
	float pTmiss2 = ( (sin(iTau2->phi())*measuredMETx) - (cos(iTau2->phi())*measuredMETy) ) / (sin(the1)*( (sin(iTau2->phi())*cos(iTau1->phi())) - (cos(iTau2->phi())*sin(iTau1->phi())) ) );
	float mcoll = ditau_M / sqrt( (pTmiss1/(pTmiss1+iTau1->pt())) + (pTmiss2/(pTmiss2+iTau2->pt())) );
	 
	// define MET covariance
	double sumPtUnclestered = 0;
	const reco::METCovMatrix cov = metSigAlgo_->getCovariance( *Jets, leptons, pfCandidates, *rho, resPtObj, resPhiObj, resSFObj, iEvent.isRealData(), sumPtUnclestered);
	TMatrixD covMET(2, 2);
	covMET[0][0] = cov[0][0]; 
	covMET[1][1] = cov[1][1]; 
	covMET[0][1] = cov[0][1]; 
	covMET[1][0] = cov[1][0]; 
	//if ( debug ) std::cout << " MET COV xx = " << covMET[0][0] << std::endl;
	//if ( debug ) std::cout << " MET COV yy = " << covMET[1][1] << std::endl;
	//if ( debug ) std::cout << " MET COV xy = " << covMET[0][1] << std::endl;
	//if ( debug ) std::cout << " MET COV yx = " << covMET[1][0] << std::endl;
	 
	std::vector<MeasuredTauLepton> measuredTauLeptons;
	if ( iTau1->pt() > iTau2->pt() ) {
	  measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()) );
	  measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()) );
	}
	else {
	  measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()) );
	  measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()) );
	}
	 
	FastMTT aFastMTTAlgo;
	aFastMTTAlgo.run(measuredTauLeptons,measuredMETx,measuredMETy,covMET);
	LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
	float svFitMass = ttP4.M(); // return value is in units of GeV
	float svFitMt = ttP4.Mt();
	float svFitPt = ttP4.Pt();
	//float svFitEta = ttP4.Eta();
	//float svFitPhi = ttP4.Phi();
	std::cout << " Found mass = " << svFitMass << std::endl;
	 
	//if (!( svFitMass > 50 && svFitMass < 130 && ditau_mth < 50 )) return false;
	//if ( svFitMass > 50 && svFitMass < 130 && ditau_mth < 50 ){
	if ( svFitMass > 50 && svFitMass < 150 ){ //HIG-19-010 approval presentation
	  //if (!IsSignal) return false;
	 
	  //make a pair of Jet and Tau index here for the tau pairs passing delta R cut
	  vecOf_jetTauPair.push_back(std::make_pair(vJetTau[j].recoJetIdxs, vJetTau[j].recoTauIdxs));
	  vecOf_jetTauPair.push_back(std::make_pair(vJetTau[k].recoJetIdxs, vJetTau[k].recoTauIdxs));
	} //mass cuts
	 
      } // if recotau dr > 0.5
    } //tau loop inside
  } //tau loop outside
  
   
  //if ( debug) std::cout << " Tau pt :" << iTau1->pt() << "  " << iTau1->eta() << "  tau2: " << iTau2->pt() << "  " << iTau2->eta() << std::endl;   
   
   

   
  // there will be duplicate entries of Jet-tau indexes here. therefore remove the duplicate ones
  // Sorting the vector of pairs
  std::sort(vecOf_jetTauPair.begin(), vecOf_jetTauPair.end());
   
  // Removing consecutive duplicates
  auto last = std::unique(vecOf_jetTauPair.begin(), vecOf_jetTauPair.end());
    
  // Erasing redundant elements
  vecOf_jetTauPair.erase(last, vecOf_jetTauPair.end());
  
  vJetTauCut.clear(); 

  unsigned int nMatchedJets = 0;
    
  // Accessing elements of the vector of pairs
  for (const auto& pair : vecOf_jetTauPair) {
    jet_tau_allCut_obj Jet_tau_allCut_obj = { pair.first, pair.second };
    vJetTauCut.push_back( Jet_tau_allCut_obj );
     
    //std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
    vJetIdxs.push_back( pair.first );
    nMatchedJets++;
 
  }
   
  
  if(vJetTauCut.empty() || vJetIdxs.empty()){
   hNpassed_mGG->Fill(0);
   return false;
   }
   
   hNpassed_mGG->Fill(1);

  
   
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;
  if ( debug ) std::cout << " >> has_jet_dijet_ditau: passed" << std::endl;
  return true;
     
} // runEvtSel_jet_dijet_ditau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_ditau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {


  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);
  
  edm::Handle<edm::View<reco::Jet> > Jets;
  iEvent.getByToken( pfjetsToken_, Jets );
  
  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);
  
  edm::Handle<reco::PFMETCollection> pfmet;
  iEvent.getByToken(metCollectionT_, pfmet);
  
  edm::Handle<reco::PFTauDiscriminator> DecayMode;
  iEvent.getByToken(tauDecayMode_, DecayMode);
  edm::Handle<reco::PFTauDiscriminator> MVAIsolation;
  iEvent.getByToken(tauMVAIsolation_, MVAIsolation);
  edm::Handle<reco::PFTauDiscriminator> MuonRejection;
  iEvent.getByToken(tauMuonRejection_, MuonRejection);
  edm::Handle<reco::PFTauDiscriminator> ElectronRejectionMVA6;
  iEvent.getByToken(tauElectronRejectionMVA6_, ElectronRejectionMVA6);
  
  JME::JetResolution resPtObj            = JME::JetResolution::get(iSetup, jetResPtType_);
  JME::JetResolution resPhiObj           = JME::JetResolution::get(iSetup, jetResPhiType_);
  JME::JetResolutionScaleFactor resSFObj = JME::JetResolutionScaleFactor::get(iSetup, jetSFType_);
  
  edm::Handle<PFCollection> pfCandsH_;
  iEvent.getByToken( pfCollectionT_, pfCandsH_ );

  edm::Handle<reco::CandidateView> pfCandidates;
  iEvent.getByToken( pfCandidatesToken_, pfCandidates );

  std::vector< edm::Handle<reco::CandidateView> > leptons;
  for ( std::vector<edm::EDGetTokenT<edm::View<reco::Candidate> > >::const_iterator srcLeptons_i = lepTokens_.begin(); srcLeptons_i != lepTokens_.end(); ++srcLeptons_i ) {
    edm::Handle<reco::CandidateView> leptons_i;
    iEvent.getByToken(*srcLeptons_i, leptons_i);
    leptons.push_back( leptons_i );
  }
  
  edm::Handle<double> rho;
  iEvent.getByToken(rhoLabel_, rho);

  //gen variables
  v_att_genHiggs_M_.clear();
  v_att_genPS_M_.clear();
  v_att_genTau_pT_.clear();

  //jet variables
  v_att_tau_jet_m0_.clear();    
  v_att_tau_jet_pt_.clear();
  v_att_tau_jetIsSignal_.clear();
  v_att_tau_jetdR_.clear();
  
  //single tau variables
  v_att_tau_pT_.clear();
  v_att_tau_mva_.clear();
  v_att_tau_dm_.clear();
  
  //ditau variables  
  v_att_dRtautau_.clear();
  v_att_Mvis_.clear();
  v_att_pTvis_.clear();
  v_att_Mtautau_.clear();
  v_att_Mtautau_svFit_.clear();
  v_att_Mth_svFit_.clear();
  v_att_pTh_svFit_.clear();
  v_att_pTh_.clear();
  v_att_mth_.clear();
  v_att_mcoll_.clear();
  v_att_pfMET_.clear();
  v_att_dphillmet_.clear();
  v_att_dphill_.clear();
  
  vJetTauFrameCropped.clear(); 
    
  int numJets=0; int numTaus=0; int tauPairs=0;
  
  v_att_tau_combs_ = 0;  v_att_tau_njet_ = 0; v_att_tau_ntau_ =0; 
    

  for ( unsigned int j = 0; j < vJetTauCut.size(); j++ ) {

    
    
    // Skip jets and taus which fails HE edge cut
    //if(vJetTauCut[j].tauIdxs.empty() || vJetTauCut[j].jetIdxs.empty())continue;
    if ( std::find(vJetIdxs.begin(), vJetIdxs.end(), vJetTauCut[j].jetIdxs) == vJetIdxs.end()) continue;
  
    reco::PFTauRef iTau1( taus, vJetTauCut[j].tauIdxs );
    reco::PFJetRef iJet1( jets, vJetTauCut[j].jetIdxs );
       
    //std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
    //std::cout << " after HBHE cut jet ID:" << vJetTauCut[j].jetIdxs << "  tau ID:" << vJetTauCut[j].tauIdxs << std::endl;
    
    //std ::cout<< " jet pt: " << iJet1->pt() << ", Eta: " << iJet1->eta() << ", Phi: " << iJet1->phi();
    //std:: cout<< " tau pT: " << iTau1->pt() << " eta:" << iTau1->eta() << " phi "<< iTau1->phi() << std::endl;
    
    // fill jet variables
    float jetTaudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iJet1->eta(),iJet1->phi() );  
    v_att_tau_jet_m0_.push_back(iJet1->mass());
    v_att_tau_jet_pt_.push_back( iJet1->pt() );
    v_att_tau_jetdR_.push_back( jetTaudR );
     
    // fill tau variables
    v_att_tau_pT_.push_back( iTau1->pt() );
    v_att_tau_mva_.push_back(((*MVAIsolation)[iTau1]) );
    v_att_tau_dm_.push_back( ((*DecayMode)[iTau1]) );
    
    jet_tau_allCut_obj Jet_tau_FrameCropped_obj = { vJetTauCut[j].jetIdxs, vJetTauCut[j].tauIdxs };
    vJetTauFrameCropped.push_back( Jet_tau_FrameCropped_obj );
    v_att_tau_njet_++;
    v_att_tau_ntau_++;
  }
    


  for ( unsigned int j = 0; j < vJetTauFrameCropped.size()-1; j++ ) {
    reco::PFTauRef iTau1( taus, vJetTauFrameCropped[j].tauIdxs );
    reco::PFJetRef iJet1( jets, vJetTauFrameCropped[j].jetIdxs );
        
         
    for ( unsigned int k = j+1; k < vJetTauFrameCropped.size(); k++ ) {
     
      if ( k <= j ) continue;
      reco::PFTauRef iTau2( taus, vJetTauFrameCropped[k].tauIdxs );
      reco::PFJetRef iJet2( jets, vJetTauFrameCropped[k].jetIdxs );
       
      float recotaudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iTau2->eta(),iTau2->phi() );
             
      std ::cout<< " jet pt: " << iJet1->pt() << ", Eta: " << iJet1->eta() << ", Phi: " << iJet1->phi();
      std:: cout<< " tau pT: " << iTau1->pt() << " eta:" << iTau1->eta() << " phi "<< iTau1->phi() << std::endl;

      std ::cout<< " jet 2 pt: " << iJet2->pt() << ", Eta: " << iJet2->eta() << ", Phi: " << iJet2->phi();
      std:: cout<< " tau 2 pT: " << iTau2->pt() << " eta:" << iTau2->eta() << " phi "<< iTau2->phi() << std::endl;
         	 
      float taudR = reco::deltaR( iTau1->eta(),iTau1->phi(), iTau2->eta(),iTau2->phi() );
      TLorentzVector Tau1  = SetTau(iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()); 
      TLorentzVector Tau2  = SetTau(iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass());
      TLorentzVector DiTau = Tau1+Tau2;
      float dphill  = abs(Tau1.DeltaPhi(Tau2));  
      float diMvis  = DiTau.M();
      float dipTvis = DiTau.Pt();
      //if ( debug ) std::cout << " Tau pair " << " ("<< iTau1->pt() << " GeV) + " << " (" << iTau2->pt() << " GeV) | ditau Mvis : " << diMvis << " GeV | dR = " << taudR << std::endl;
	 
	 
      float pfMET    = (pfmet->front()).et();
      float pfMETphi = (pfmet->front()).phi();
      TLorentzVector MET = SetMET(pfMET,pfMETphi);
      //if ( debug ) std::cout << " PF MET = " << pfMET << std::endl;
	 
      //define MET for SVFit
      double measuredMETx = (pfmet->front()).px();
      double measuredMETy = (pfmet->front()).py();
      //if ( debug ) std::cout << " PF MET X = " << measuredMETx << " | PF MET Y = " << measuredMETy << std::endl;
      float dphillmet = abs(DiTau.DeltaPhi(MET));  
      float ditau_mth = sqrt( 2. * dipTvis * pfMET * ( 1. - cos (dphillmet) ));
	 
      TLorentzVector Higgs = Tau1+Tau2+MET;
      float ditau_M   = Higgs.M();
      float ditau_pT  = Higgs.Pt();
	 
      float the1    = 2.*atan(exp(-iTau1->eta()));
      float the2    = 2.*atan(exp(-iTau2->eta()));
      float pTmiss1 = ( (sin(iTau1->phi())*measuredMETx) - (cos(iTau1->phi())*measuredMETy) ) / (sin(the2)*( (sin(iTau1->phi())*cos(iTau2->phi())) - (cos(iTau1->phi())*sin(iTau2->phi())) ) );
      float pTmiss2 = ( (sin(iTau2->phi())*measuredMETx) - (cos(iTau2->phi())*measuredMETy) ) / (sin(the1)*( (sin(iTau2->phi())*cos(iTau1->phi())) - (cos(iTau2->phi())*sin(iTau1->phi())) ) );
      float mcoll = ditau_M / sqrt( (pTmiss1/(pTmiss1+iTau1->pt())) + (pTmiss2/(pTmiss2+iTau2->pt())) );
	 
      // define MET covariance
      double sumPtUnclestered = 0;
      const reco::METCovMatrix cov = metSigAlgo_->getCovariance( *Jets, leptons, pfCandidates, *rho, resPtObj, resPhiObj, resSFObj, iEvent.isRealData(), sumPtUnclestered);
      TMatrixD covMET(2, 2);
      covMET[0][0] = cov[0][0]; 
      covMET[1][1] = cov[1][1]; 
      covMET[0][1] = cov[0][1]; 
      covMET[1][0] = cov[1][0]; 
      //if ( debug ) std::cout << " MET COV xx = " << covMET[0][0] << std::endl;
      //if ( debug ) std::cout << " MET COV yy = " << covMET[1][1] << std::endl;
      //if ( debug ) std::cout << " MET COV xy = " << covMET[0][1] << std::endl;
      //if ( debug ) std::cout << " MET COV yx = " << covMET[1][0] << std::endl;
	 
      std::vector<MeasuredTauLepton> measuredTauLeptons;
      if ( iTau1->pt() > iTau2->pt() ) {
	measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()) );
	measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()) );
      }
      else {
	measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau2->pt(), iTau2->eta(), iTau2->phi(), iTau2->mass()) );
	measuredTauLeptons.push_back( MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  iTau1->pt(), iTau1->eta(), iTau1->phi(), iTau1->mass()) );
      }
	 
      FastMTT aFastMTTAlgo;
      aFastMTTAlgo.run(measuredTauLeptons,measuredMETx,measuredMETy,covMET);
      LorentzVector ttP4 = aFastMTTAlgo.getBestP4();
      float svFitMass = ttP4.M(); // return value is in units of GeV
      float svFitMt = ttP4.Mt();
      float svFitPt = ttP4.Pt();
      //float svFitEta = ttP4.Eta();
      //float svFitPhi = ttP4.Phi();
      //std::cout << " Found mass = " << svFitMass << std::endl;
	 
	 
      v_att_dRtautau_.push_back(taudR);
      v_att_Mvis_.push_back(diMvis);
      v_att_pTvis_.push_back(dipTvis);
      v_att_Mtautau_.push_back(ditau_M);
      v_att_Mtautau_svFit_.push_back(ditau_M);
      v_att_Mth_svFit_.push_back(svFitMass);
      v_att_pTh_svFit_.push_back(svFitPt);
      v_att_pTh_.push_back(ditau_pT);
      v_att_mth_.push_back(ditau_mth);
      v_att_mcoll_.push_back(mcoll);
      v_att_pfMET_.push_back(pfMET);
      v_att_dphillmet_.push_back(dphillmet);
      v_att_dphill_.push_back(dphill);
  
      v_att_tau_combs_++;

    } //tau loop inside
  } //tau loop outside
  //v_att_tau_combs_ = tauPairs;
   
  
} // fillEvtSel_jet_dijet_ditau()
