#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 50; //TODO: use cfg level nJets_
TH1D *h_tau_mr_jet_pT;
TH1D *h_tau_mr_jet_E;
TH1D *h_tau_mr_jet_eta;
TH1D *h_tau_mr_jet_m0;
TH1D *h_tau_mr_jet_ma;
TH1D *h_tau_mr_jet_pta;
TH2D *h_tau_mr_jet_a_m_pt;
TH1D *h_tau_mr_jet_nJet;
TH1D *h_tau_mr_jet_isDiTau;
TH1D *h_tau_mr_jet_dR;
TH1D *h_tau_mr_jet_TaudR;
TH1D *h_tau_mr_jet_Tau1dR;
TH1D *h_tau_mr_jet_Tau2dR;
TH1D *h_tau_mr_jet_Tau1pT;
TH1D *h_tau_mr_jet_Tau2pT;
TH1D *h_tau_mr_jet_NrecoTaus;
TH1D *h_tau_mr_jet_NGenTaus;
TH1D *h_tau_mr_jet_recoTau1dR;
TH1D *h_tau_mr_jet_recoTau2dR;
TH1D *h_tau_mr_jet_n1dR;
TH1D *h_tau_mr_jet_n2dR;
vector<float> v_mr_jetIsDiTau;
vector<float> v_mr_jetadR;
vector<float> v_mr_ma;
vector<float> v_mr_pta;
vector<float> v_mr_jetTaudR;
vector<float> v_mr_jetTau1dR;
vector<float> v_mr_jetTau2dR;
vector<float> v_mr_jetTau1pT;
vector<float> v_mr_jetTau2pT;
vector<float> v_mr_jetNGenTaus;
vector<float> v_mr_jetNrecoTaus;
vector<float> v_mr_jetrecoTau1dR;
vector<float> v_mr_jetrecoTau2dR;
vector<float> v_mr_jetn1dR;
vector<float> v_mr_jetn2dR;

vector<float> v_mr_tau_jet_m0_;
vector<float> v_mr_tau_jet_ma_;
vector<float> v_mr_tau_jet_pta_;
vector<float> v_mr_tau_jet_pt_;
vector<float> v_mr_tau_jetPdgIds_;
vector<float> v_mr_tau_jetIsDiTau_;
vector<float> v_mr_tau_jetadR_;
vector<float> v_mr_tau_jetTaudR_;
vector<float> v_mr_tau_jetTau1dR_;
vector<float> v_mr_tau_jetTau2dR_;
vector<float> v_mr_tau_jetTau1pT_;
vector<float> v_mr_tau_jetTau2pT_;
vector<float> v_mr_tau_jetNGenTaus_;
vector<float> v_mr_tau_jetNrecoTaus_;
vector<float> v_mr_tau_jetrecoTau1dR_;
vector<float> v_mr_tau_jetrecoTau2dR_;
vector<float> v_mr_tau_jetn1dR_;
vector<float> v_mr_tau_jetn2dR_;

vector<float> v_mr_tau_subJetE_[nJets];
vector<float> v_mr_tau_subJetPx_[nJets];
vector<float> v_mr_tau_subJetPy_[nJets];
vector<float> v_mr_tau_subJetPz_[nJets];


int sum(vector <int> dist) {
    return std::accumulate(dist.begin(), dist.end(), 0);
}

float max_element(vector <float> dist) {
    float max = 0;
    int s = dist.size();
    for (int i = 0; i < s; i++) {
        float el = dist[i];
        if (max < el){max = el;}
    }
    return max;
}

vector <float> get_inverse_pdf(vector <int> dist) {
    vector <float> invpdf(dist.size());
    float sum_hist = sum(dist);
    int s = dist.size();
    for (int i = 0; i < s; i++) {
        if (dist[i] != 0 ) {invpdf[i] = sum_hist / dist[i];}
        else {invpdf[i] = 1;}
    }
    float max_invpdf = max_element(invpdf);
    for (int i = 0; i < s; i++) {
        invpdf[i] = invpdf[i] / max_invpdf;
    }
    return invpdf;
}

float lookup_mass_invpdf(float Mgen, vector <float> M_bins, vector <float> M_invpdf) {
    int ipt = 0;
    int s1 = M_bins.size();
    int s2 = M_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
        ipt = ib;
        if (ib + 1 >  s2 - 1) { break; }
        if (Mgen <= M_bins[ib]) { break; }
    }
    if (debug) std::cout << "   mass gen = " << Mgen << " | bin = " << ipt << " | mass bin = " << M_bins[ipt] << " | inv mass bin = " << M_invpdf[ipt] << std::endl;
    return M_invpdf[ipt];
}

float lookup_pt_invpdf(int pTgen, vector <int> pT_bins, vector <float> pT_invpdf) {
    int ipt = 0;
    int s1 = pT_bins.size();
    int s2 = pT_invpdf.size();
    for (int ib = 0; ib < s1; ib++) {
        ipt = ib;
        if (ib + 1 >  s2 - 1) { break; }
        if (pTgen <= pT_bins[ib]) { break; }
    }
    if (debug) std::cout << "   pT gen = " << pTgen << " | bin = " << ipt << " | pt bin = " << pT_bins[ipt] << " | inv pt bin = " << pT_invpdf[ipt] << std::endl;
    return pT_invpdf[ipt];
}

float lookup_invpdf(float Mgen, vector <float> M_bins, float pTgen, vector <int> pT_bins, vector <float> invpdf) {
    unsigned int ibin = 0;
    unsigned int ibinx = 0;
    unsigned int ibiny = 0;
    unsigned int m1  = M_bins.size();
    unsigned int pt1 = pT_bins.size();
    unsigned int inv = invpdf.size();
    bool found_mass = false;
    bool found_end  = false;
    if (debug) std::cout << "   Nbin = " << inv << " | Nbinx = " << m1 << " , Nbiny = " << pt1 << std::endl;
    for (unsigned int ibx = 0; ibx < m1; ibx++) {
        if (found_mass || found_end) { break; }
        for (unsigned int iby = 0; iby < pt1; iby++) {
            ibin = (ibx*pt1) + iby;
            ibinx = ibx;
            ibiny = iby;
            if ( ((ibx*pt1) + iby + 1) >  (inv - 1) ) { 
                found_end = true;
                break; 
            }
            if ( (Mgen  <= M_bins[ibx]) && (pTgen <= pT_bins[iby]) ) {
                found_mass = true;
                break; 
            }
        }
    }
    if (debug) std::cout << "   mass gen   = " << Mgen  << " | bin = " << ibinx << " | mass bin = " << M_bins[ibinx]  << std::endl;
    if (debug) std::cout << "   pt gen     = " << pTgen << " | bin = " << ibiny << " | pt bin   = " << pT_bins[ibiny] << std::endl;
    if (debug) std::cout << "   Global bin = " << ibin <<  " | inv mass bin = " << invpdf[ibin] << std::endl;
    return invpdf[ibin];
}

float get_rand_el(vector <int> dist) {
    int randomIndex = rand() % dist.size();
      return dist[randomIndex];
}

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau_massregression ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_mr_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                   , 100,  0., 800.);
  h_tau_mr_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                           , 100,  0., 800.);
  h_tau_mr_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_tau_mr_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_tau_mr_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  //h_tau_mr_jet_a_m_pt     = fs->make<TH2D>("h_a_m_pT"         , "m^{a} vs p_{T}^{a};m^{a} vs p_{T}^{a};Jets" ,  26, 3.6,  14, 18, 30., 120.);
  //h_tau_mr_jet_ma         = fs->make<TH1D>("h_jet_ma"         , "m^{a};m^{a};Jets"                           ,  26, 3.6,  14);
  //h_tau_mr_jet_pta        = fs->make<TH1D>("h_jet_pta"        , "p_{T}^{a};p_{T}^{a};Jets"                   ,  18, 30., 120.);
  h_tau_mr_jet_a_m_pt     = fs->make<TH2D>("h_a_m_pT"         , "m^{a} vs p_{T}^{a};m^{a} vs p_{T}^{a};Jets" ,  26, 3.6,  14, 24, 30., 150.);
  h_tau_mr_jet_ma         = fs->make<TH1D>("h_jet_ma"         , "m^{a};m^{a};Jets"                           ,  26, 3.6,  14);
  h_tau_mr_jet_pta        = fs->make<TH1D>("h_jet_pta"        , "p_{T}^{a};p_{T}^{a};Jets"                   ,  24, 30., 150.);
  h_tau_mr_jet_isDiTau    = fs->make<TH1D>("h_jet_isDiTau"    , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_tau_mr_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{a,j};dR_{a,j};Jets"                     ,  50,  0.,  0.5);
  h_tau_mr_jet_TaudR      = fs->make<TH1D>("h_jet_TaudR"      , "dR_{#tau,#tau};dR_{#tau,#tau};Jets"         ,  50,  0.,   1.);
  h_tau_mr_jet_Tau1dR     = fs->make<TH1D>("h_jet_Tau1dR"     , "dR_{#tau_{1},j};dR_{#tau_{1},j};Jets"       ,  50,  0.,  0.5);
  h_tau_mr_jet_Tau2dR     = fs->make<TH1D>("h_jet_Tau2dR"     , "dR_{#tau_{2},j};dR_{#tau_{2},j};Jets"       ,  50,  0.,  0.5);
  h_tau_mr_jet_Tau1pT     = fs->make<TH1D>("h_jet_Tau1pT"     , "p_{T}^{#tau_{1}};p_{T}^{#tau_{1}};Jets"     ,  50,  0.,  100);
  h_tau_mr_jet_Tau2pT     = fs->make<TH1D>("h_jet_Tau2pT"     , "p_{T}^{#tau_{2}};p_{T}^{#tau_{2}};Jets"     ,  50,  0.,  100);
  h_tau_mr_jet_NGenTaus  = fs->make<TH1D>("h_jet_NGenTaus"    , "N#tau^{RECO};N#tau^{RECO};Jets"             ,   5,  0.,   5.);
  h_tau_mr_jet_NrecoTaus  = fs->make<TH1D>("h_jet_NrecoTaus"  , "N#tau^{RECO};N#tau^{RECO};Jets"             ,   5,  0.,   5.);
  h_tau_mr_jet_recoTau1dR = fs->make<TH1D>("h_jet_recoTau1dR" , "dR_{#tau_{1}^{RECO},j};dR_{#tau_{1}^{RECO},j};Jets" ,  50,  0.,  0.5);
  h_tau_mr_jet_recoTau2dR = fs->make<TH1D>("h_jet_recoTau2dR" , "dR_{#tau_{2}^{RECO},j};dR_{#tau_{2}^{RECO},j};Jets" ,  25,  0.,  0.5);
  h_tau_mr_jet_n1dR       = fs->make<TH1D>("h_jet_n1dR"       , "dR_{#eta_{1},j};dR_{#eta_{1},j};Jets"       ,  25,  0.,  0.5);
  h_tau_mr_jet_n2dR       = fs->make<TH1D>("h_jet_n2dR"       , "dR_{#eta_{2},j};dR_{#eta_{2},j};Jets"       ,  25,  0.,  0.5);

  tree->Branch("jetM",       &v_mr_tau_jet_m0_);
  tree->Branch("jetPt",      &v_mr_tau_jet_pt_);
  tree->Branch("jetPdgIds",  &v_mr_tau_jetPdgIds_);
  tree->Branch("jetadR",     &v_mr_tau_jetadR_);
  tree->Branch("jetIsDiTau", &v_mr_tau_jetIsDiTau_);
  tree->Branch("a_m",        &v_mr_tau_jet_ma_);
  tree->Branch("a_pt",       &v_mr_tau_jet_pta_);
  tree->Branch("jetpT",      &v_mr_tau_jet_pt_);
  tree->Branch("TaudR",      &v_mr_tau_jetTaudR_);
  tree->Branch("Tau1dR",     &v_mr_tau_jetTau1dR_);
  tree->Branch("Tau2dR",     &v_mr_tau_jetTau2dR_);
  tree->Branch("Tau1pT",     &v_mr_tau_jetTau1pT_);
  tree->Branch("Tau2pT",     &v_mr_tau_jetTau2pT_);
  tree->Branch("NGenTaus",   &v_mr_tau_jetNGenTaus_);
  tree->Branch("NrecoTaus",  &v_mr_tau_jetNrecoTaus_);
  tree->Branch("recoTau1dR", &v_mr_tau_jetrecoTau1dR_);
  tree->Branch("recoTau2dR", &v_mr_tau_jetrecoTau2dR_);
  tree->Branch("n1dR",       &v_mr_tau_jetn1dR_);
  tree->Branch("n2dR",       &v_mr_tau_jetn2dR_);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_mr_tau_subJetPz_[iJ]);
  }

} // branchesEvtSel_jet_dijet_tau_massregression()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_tau_massregression( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  vJetIdxs.clear();
  v_mr_tau_jetPdgIds_.clear();
  v_mr_jetIsDiTau.clear();
  v_mr_jetadR.clear();
  v_mr_ma.clear();
  v_mr_pta.clear();
  v_mr_jetTaudR.clear();
  v_mr_jetTau1dR.clear();
  v_mr_jetTau2dR.clear();
  v_mr_jetTau1pT.clear();
  v_mr_jetTau2pT.clear();
  v_mr_jetNGenTaus.clear();
  v_mr_jetNrecoTaus.clear();
  v_mr_jetrecoTau1dR.clear();
  v_mr_jetrecoTau2dR.clear();
  v_mr_jetn1dR.clear();
  v_mr_jetn2dR.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_mr_tau_jetFakePhoIdxs;
  */

  unsigned int nMatchedJets = 0;
  unsigned int nMatchedRecoTaus = 0;
  unsigned int aPdgId           = 0;
  bool MatchedPseudoScalar = false;
  float a_mass = -99.;
  float a_pt   = -99.;
  float dRa    = -99.;
  float tausdR =  99.;
  float tau1dR = -99.;
  float tau2dR = -99.;
  float tau1pT = -99.;
  float tau2pT = -99.;
  float recotau1dR = -99.;
  float recotau2dR = -99.;
  float n1dR   = -99.;
  float n2dR   = -99.;

  bool unbiasing = false;
  vector <int> pT_bins   = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150};
  vector <float> m_bins = {4.0, 4.4, 4.800000000000001, 5.2, 5.6, 6.0, 6.4, 6.800000000000001, 7.2, 7.6, 8.0, 8.4, 8.8, 9.200000000000001, 9.6, 10.0, 10.4, 10.8, 11.200000000000001, 11.6, 12.0, 12.4, 12.8, 13.200000000000001, 13.6, 14.0};
  vector <int> occ = {
//  35,   40,   45,   50,   55,   60,   65,   70,   75,   80,   85,   90,   95,  100,  105,  110,  115,  120,  125,  130,  135,  140,  145,  150
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, // 3.6 -  4.0
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  4.4
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  4.8
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  5.2
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  5.6
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  6.0
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  6.4
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  6.8
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  7.2
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  7.6
   100,  114,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  8.0 #!!!!
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  8.4
   100,  100,  104,  100,  104,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  8.8
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  9.2
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    ->  9.6
   100,  100,  100,  152,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 10.0 #!!!!
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 10.4
   100,  100,  100,  100,  130,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 10.8 #!!!!
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 11.2
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 11.6
   100,  100,  100,  100,  100,  398,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 12.0 #!!!!
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 12.4
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 12.8 #!!!!
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 13.2
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100, //    -> 13.6
   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100  //    -> 14.0
};
  //vector <float> m_invpdf  = get_inverse_pdf(m_occ);
  //vector <float> pT_invpdf = get_inverse_pdf(pT_occ);
  vector <float> invpdf    = get_inverse_pdf(occ);

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] -> Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    unsigned int nMatchedGenParticles = 0;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0; 
    unsigned int NTau1Daughters = 0;
    unsigned int NTau2Daughters = 0;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      if ( abs(iGen->pdgId()) != 15 ) continue;
      if ( iGen->status() != 2 ) continue;
      if (iGen->numberOfMothers() < 1) continue;
      ++iGenParticle;
      float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.4 ) continue;

      if ( iGen->numberOfMothers() != 1 ) continue;
      aPdgId = std::abs(iGen->mother()->pdgId());
      if ( abs(iGen->mother()->pdgId()) == 25 && iGen->mother()->mass() < 15 && iGen->mother()->mass() > 3.5) {
        MatchedPseudoScalar = true;
      }
      else continue;

      if (nMatchedGenParticles == 0) {

        a_mass = iGen->mother()->mass();
        a_pt   = iGen->mother()->pt();
        dRa = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->eta(),iGen->mother()->phi() );
        if ( iGen->mother()->numberOfDaughters() == 2 ){ 
          if (abs(iGen->mother()->daughter(0)->pdgId()) == 15 && abs(iGen->mother()->daughter(1)->pdgId()) == 15){ 
            tausdR = reco::deltaR( iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
            if ( tausdR > 0.4 ) continue;

            if ( debug ) std::cout << "   TAUS MOTHER  -> status: " << iGen->mother()->status() << ", id: " << iGen->mother()->pdgId() << ", nDaught: " << iGen->mother()->numberOfDaughters() << " | pt: "<< iGen->mother()->pt() << " eta: " <<iGen->mother()->eta() << " phi: " <<iGen->mother()->phi() << " mass: " <<iGen->mother()->mass() << std::endl;

            if ( unbiasing ) {
              // 2D
              float rand_sampler  = rand() / float(RAND_MAX);
              float pT_gen        = iGen->mother()->pt();
              float m_gen         = iGen->mother()->mass();
              if ( debug ) std::cout << "   Testing with pT = " << pT_gen << " and m = " << m_gen << std::endl;
              float wgt           = lookup_invpdf(m_gen, m_bins, pT_gen, pT_bins, invpdf);
              if (debug) std::cout << "      wgt " << wgt  << " | rand_sampler " << rand_sampler << std::endl;
              if (rand_sampler > wgt) continue;
            }

            if ( iGen->mother()->daughter(0)->pt() > iGen->mother()->daughter(1)->pt() ) {
              tau1pT = iGen->mother()->daughter(0)->pt();
              tau2pT = iGen->mother()->daughter(1)->pt();
              tau1dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi() );
              tau2dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
              
              if ( debug ) std::cout << "   >>>>>> Taus dR = " << tausdR << " , tau1 dR = " << tau1dR << " , tau2 dR = " << tau2dR << std::endl;
              NTau1Daughters = iGen->mother()->daughter(0)->numberOfDaughters();
              NTau2Daughters = iGen->mother()->daughter(1)->numberOfDaughters();
              if ( debug ) std::cout << "    >>>>>> # Tau 1 daughters = " << NTau1Daughters << ",  # Tau 2 daughters = "<< NTau2Daughters << std::endl;
              for (unsigned int iDaughter = 0; iDaughter != NTau1Daughters; ++iDaughter){
                if ( debug ) std::cout << "   >>>>>> Tau 1 Decay = " << iGen->mother()->daughter(0)->daughter(iDaughter)->pdgId() << std::endl;
                if ( abs(iGen->mother()->daughter(0)->daughter(iDaughter)->pdgId()) == 16 ) n1dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(0)->daughter(iDaughter)->eta(), iGen->mother()->daughter(0)->daughter(iDaughter)->phi() );
              }
              for (unsigned int iDaughter = 0; iDaughter != NTau2Daughters; ++iDaughter){
                if ( debug ) std::cout << "   >>>>>> Tau 2 Decay = " << iGen->mother()->daughter(1)->daughter(iDaughter)->pdgId() << std::endl;
                if ( abs(iGen->mother()->daughter(1)->daughter(iDaughter)->pdgId()) == 16 ) n2dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(1)->daughter(iDaughter)->eta(), iGen->mother()->daughter(1)->daughter(iDaughter)->phi() );
              }
            } else {
              tau1pT = iGen->mother()->daughter(1)->pt();
              tau2pT = iGen->mother()->daughter(0)->pt();
              tau1dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
              tau2dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi() );

              if ( debug ) std::cout << "   >>>>>> Taus dR = " << tausdR << " , tau1 dR = " << tau1dR << " , tau2 dR = " << tau2dR << std::endl;
              NTau1Daughters = iGen->mother()->daughter(1)->numberOfDaughters();
              NTau2Daughters = iGen->mother()->daughter(0)->numberOfDaughters();
              if ( debug ) std::cout << "    >>>>>> # Tau 1 daughters = " << NTau1Daughters << ",  # Tau 2 daughters = "<< NTau2Daughters << std::endl;
              for (unsigned int iDaughter = 0; iDaughter != NTau1Daughters; ++iDaughter){
                if ( debug ) std::cout << "   >>>>>> Tau 1 Decay = " << iGen->mother()->daughter(1)->daughter(iDaughter)->pdgId() << std::endl;
                if ( abs(iGen->mother()->daughter(1)->daughter(iDaughter)->pdgId()) == 16 ) n1dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(1)->daughter(iDaughter)->eta(), iGen->mother()->daughter(1)->daughter(iDaughter)->phi() );
              }
              for (unsigned int iDaughter = 0; iDaughter != NTau2Daughters; ++iDaughter){
                if ( debug ) std::cout << "   >>>>>> Tau 2 Decay = " << iGen->mother()->daughter(0)->daughter(iDaughter)->pdgId() << std::endl;
                if ( abs(iGen->mother()->daughter(0)->daughter(iDaughter)->pdgId()) == 16 ) n2dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->daughter(0)->daughter(iDaughter)->eta(), iGen->mother()->daughter(0)->daughter(iDaughter)->phi() );
              }
            } // end else pt2 > pt1
          }
        }
        //}
        if (debug ) std::cout << "   >>>>>> n1 dR = " << n1dR << " , n2 dR = " << n2dR << std::endl;
      }

      if ( debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << " | dR = "<< dR << std::endl ;

      if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] matched particle [" << iGenParticle << "] -> pdgId: " << std::abs(iGen->pdgId()) << " | dR: " << dR << "| Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;

      ++nMatchedGenParticles;
    } // primary gen particles
    if ( nMatchedGenParticles > 0 ) passedGenSel = true;
    if (passedGenSel) { 
      ++nMatchedJets;

      //Lookin at RecoTaus
      if (taus->size() == 0) {
        if (debug ) std::cout << "   !!!!!!!!!!  NO RECO TAUS IN THIS EVENT  !!!!!!!!!!"<< std::endl;
      }
      for ( unsigned iT(0); iT != taus->size(); ++iT ) {
        reco::PFTauRef iTau( taus, iT );
        float recotaudR = reco::deltaR( iJet->eta(),iJet->phi(), iTau->eta(),iTau->phi() );
        if ( recotaudR < 0.4 && nMatchedRecoTaus == 0 ) {
          if ( debug ) std::cout << "Reco Tau [" << iT << "]  matched jet [" << iJ << "] : dR = " << recotaudR << std::endl;
          recotau1dR = recotaudR;
          ++nMatchedRecoTaus;
        } else if ( recotaudR < 0.4 && nMatchedRecoTaus == 1 ) {
          if ( debug ) std::cout << "Reco Tau [" << iT << "]  matched jet [" << iJ << "] : dR = " << recotaudR << std::endl;
          if (recotaudR < recotau1dR) {
            recotau2dR = recotau1dR;
            recotau1dR = recotaudR;
          } else recotau2dR = recotaudR;
          ++nMatchedRecoTaus;
        } else if ( debug && recotaudR < 0.4 && nMatchedRecoTaus > 1 ) {
          std::cout << "   !!!!!!!!!!  FOUND MORE THAN 2 TAUS INSIDE JET CONE OF 0.4 !!!!!!!!!!"<< std::endl;
          if (recotaudR < recotau2dR && recotaudR < recotau1dR) { 
            if (recotau1dR < recotau2dR) recotau2dR = recotau1dR;
            recotau1dR = recotaudR;
          } else if (recotaudR < recotau2dR && recotaudR > recotau1dR) recotau2dR = recotaudR; 
          ++nMatchedRecoTaus;
        } else if ( debug ) {
          std::cout << "   !!!!!!!!!!  NO MATCH FOR Reco Tau [" << iT << "]  with jet [" << iJ << "] : dR = " << recotaudR << std::endl;
        }
      }

      vJetIdxs.push_back( iJ );
      v_mr_tau_jetPdgIds_.push_back( aPdgId );
      v_mr_jetadR.push_back( dRa );
      v_mr_jetTaudR.push_back( tausdR );
      v_mr_ma.push_back( a_mass );
      v_mr_pta.push_back( a_pt );
      v_mr_jetTau1dR.push_back( tau1dR );
      v_mr_jetTau2dR.push_back( tau2dR );
      v_mr_jetTau1pT.push_back( tau1pT );
      v_mr_jetTau2pT.push_back( tau2pT );
      v_mr_jetn1dR.push_back( n1dR );
      v_mr_jetn2dR.push_back( n2dR );
      v_mr_jetNGenTaus.push_back( nMatchedGenParticles );
      v_mr_jetNrecoTaus.push_back( nMatchedRecoTaus );
      v_mr_jetrecoTau1dR.push_back( recotau1dR );
      v_mr_jetrecoTau2dR.push_back( recotau2dR );
      v_mr_jetIsDiTau.push_back( MatchedPseudoScalar );

    }

  } // reco jets
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 1 ) return false;

  if ( debug ) std::cout << " >> has_jet_dijet_tau_massregression: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau_massregression()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau_massregression ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_mr_jet_nJet->Fill( vJetIdxs.size() );

  v_mr_tau_jet_pt_.clear();
  v_mr_tau_jet_m0_.clear();
  v_mr_tau_jet_ma_.clear();
  v_mr_tau_jet_pta_.clear();
  v_mr_tau_jetIsDiTau_.clear();
  v_mr_tau_jetadR_.clear();
  v_mr_tau_jetTaudR_.clear();
  v_mr_tau_jetTau1dR_.clear();
  v_mr_tau_jetTau2dR_.clear();
  v_mr_tau_jetTau1pT_.clear();
  v_mr_tau_jetTau2pT_.clear();
  v_mr_tau_jetNGenTaus_.clear();
  v_mr_tau_jetNrecoTaus_.clear();
  v_mr_tau_jetrecoTau1dR_.clear();
  v_mr_tau_jetrecoTau2dR_.clear();
  v_mr_tau_jetn1dR_.clear();
  v_mr_tau_jetn2dR_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_mr_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_mr_jet_eta->Fill( iJet->eta() );
    h_tau_mr_jet_E->Fill( iJet->energy() );
    h_tau_mr_jet_m0->Fill( iJet->mass() );
    h_tau_mr_jet_isDiTau->Fill( v_mr_jetIsDiTau[iJ] );
    h_tau_mr_jet_dR->Fill( v_mr_jetadR[iJ] );
    h_tau_mr_jet_ma->Fill( v_mr_ma[iJ] );
    h_tau_mr_jet_pta->Fill( v_mr_pta[iJ] );
    h_tau_mr_jet_a_m_pt->Fill( v_mr_ma[iJ], v_mr_pta[iJ] );
    h_tau_mr_jet_TaudR->Fill( v_mr_jetTaudR[iJ] );
    h_tau_mr_jet_Tau1dR->Fill( v_mr_jetTau1dR[iJ] );
    h_tau_mr_jet_Tau2dR->Fill( v_mr_jetTau2dR[iJ] );
    h_tau_mr_jet_Tau1pT->Fill( v_mr_jetTau1pT[iJ] );
    h_tau_mr_jet_Tau2pT->Fill( v_mr_jetTau2pT[iJ] );
    h_tau_mr_jet_NGenTaus->Fill( v_mr_jetNGenTaus[iJ] );
    h_tau_mr_jet_NrecoTaus->Fill( v_mr_jetNrecoTaus[iJ] );
    h_tau_mr_jet_recoTau1dR->Fill( v_mr_jetrecoTau1dR[iJ] );
    h_tau_mr_jet_recoTau2dR->Fill( v_mr_jetrecoTau2dR[iJ] );
    h_tau_mr_jet_n1dR->Fill( v_mr_jetn1dR[iJ] );
    h_tau_mr_jet_n2dR->Fill( v_mr_jetn2dR[iJ] );

    // Fill branches 
    v_mr_tau_jet_pt_.push_back( iJet->pt() );
    v_mr_tau_jet_m0_.push_back( iJet->mass() );
    v_mr_tau_jet_ma_.push_back( v_mr_ma[iJ] );
    v_mr_tau_jet_pta_.push_back( v_mr_pta[iJ] );
    v_mr_tau_jetIsDiTau_.push_back( v_mr_jetIsDiTau[iJ] );
    v_mr_tau_jetadR_.push_back( v_mr_jetadR[iJ] );
    v_mr_tau_jetTaudR_.push_back( v_mr_jetTaudR[iJ] );
    v_mr_tau_jetTau1dR_.push_back( v_mr_jetTau1dR[iJ] );
    v_mr_tau_jetTau2dR_.push_back( v_mr_jetTau2dR[iJ] );
    v_mr_tau_jetTau1pT_.push_back( v_mr_jetTau1pT[iJ] );
    v_mr_tau_jetTau2pT_.push_back( v_mr_jetTau2pT[iJ] );
    v_mr_tau_jetNGenTaus_.push_back( v_mr_jetNGenTaus[iJ] );
    v_mr_tau_jetNrecoTaus_.push_back( v_mr_jetNrecoTaus[iJ] );
    v_mr_tau_jetrecoTau1dR_.push_back( v_mr_jetrecoTau1dR[iJ] );
    v_mr_tau_jetrecoTau2dR_.push_back( v_mr_jetrecoTau2dR[iJ] );
    v_mr_tau_jetn1dR_.push_back( v_mr_jetn1dR[iJ] );
    v_mr_tau_jetn2dR_.push_back( v_mr_jetn2dR[iJ] );

    // Gen jet constituents
    v_mr_tau_subJetE_[iJ].clear();
    v_mr_tau_subJetPx_[iJ].clear();
    v_mr_tau_subJetPy_[iJ].clear();
    v_mr_tau_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_mr_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_mr_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_mr_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_mr_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_dijet_tau_massregression()
