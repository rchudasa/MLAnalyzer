#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 4; //TODO: use cfg level nJets_
TH1D *h_tau_jet_pT;
TH1D *h_tau_jet_E;
TH1D *h_tau_jet_eta;
TH1D *h_tau_jet_m0;
TH1D *h_tau_jet_ma;
TH1D *h_tau_jet_pta;
TH2D *h_tau_jet_a_m_pt;
TH1D *h_tau_jet_nJet;
TH1D *h_tau_jet_isDiTau;
TH1D *h_tau_jet_dR;
TH1D *h_tau_jet_TaudR;
TH1D *h_tau_jet_Tau1dR;
TH1D *h_tau_jet_Tau2dR;
TH1D *h_tau_jet_NrecoTaus;
TH1D *h_tau_jet_NGenTaus;
TH1D *h_tau_jet_recoTau1dR;
TH1D *h_tau_jet_recoTau2dR;
TH1D *h_tau_jet_n1dR;
TH1D *h_tau_jet_n2dR;
vector<float> v_jetIsDiTau;
vector<float> v_jetadR;
vector<float> v_ma;
vector<float> v_pta;
vector<float> v_jetTaudR;
vector<float> v_jetTau1dR;
vector<float> v_jetTau2dR;
vector<float> v_jetNGenTaus;
vector<float> v_jetNrecoTaus;
vector<float> v_jetrecoTau1dR;
vector<float> v_jetrecoTau2dR;
vector<float> v_jetn1dR;
vector<float> v_jetn2dR;

vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_ma_;
vector<float> v_tau_jet_pta_;
vector<float> v_tau_jet_pt_;
vector<float> v_tau_jetPdgIds_;
vector<float> v_tau_jetIsDiTau_;
vector<float> v_tau_jetadR_;
vector<float> v_tau_jetTaudR_;
vector<float> v_tau_jetTau1dR_;
vector<float> v_tau_jetTau2dR_;
vector<float> v_tau_jetNGenTaus_;
vector<float> v_tau_jetNrecoTaus_;
vector<float> v_tau_jetrecoTau1dR_;
vector<float> v_tau_jetrecoTau2dR_;
vector<float> v_tau_jetn1dR_;
vector<float> v_tau_jetn2dR_;

vector<float> v_tau_subJetE_[nJets];
vector<float> v_tau_subJetPx_[nJets];
vector<float> v_tau_subJetPy_[nJets];
vector<float> v_tau_subJetPz_[nJets];


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
        invpdf[i] = invpdf[i] / max_invpdf;}
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
    if (debug) std::cout << "mass gen = " << Mgen << " | bin = " << ipt << " | mass bin = " << M_bins[ipt] << " | inv mass bin = " << M_invpdf[ipt] << std::endl;
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
    if (debug) std::cout << "pT gen = " << pTgen << " | bin = " << ipt << " | pt bin = " << pT_bins[ipt] << " | inv pt bin = " << pT_invpdf[ipt] << std::endl;
    return pT_invpdf[ipt];
}

float lookup_invpdf(float Mgen, vector <float> M_bins, int pTgen, vector <int> pT_bins, vector <float> invpdf) {
    int ibin = 0;
    int m1  = M_bins.size();
    int pt1 = pT_bins.size();
    int inv = invpdf.size();
    for (int ibx = 0; ibx < m1; ibx++) {
        for (int iby = 0; iby < pt1; iby++) {
        ibin = (ibx*pt1)+ iby + 1;
        if ((ibx*pt1)+ iby + 1 >  inv - 1) { break; }
        if (Mgen  <= M_bins[ibx]) { break; }
        if (pTgen <= pT_bins[iby]) { break; }
    }
    std::cout << "mass gen = " << Mgen  << " | bin = " << ibx << " | mass bin = " << M_bins[ibx]  << std::endl;
    std::cout << "pt gen   = " << pTgen << " | bin = " << iby << " | pt bin   = " << pT_bins[iby] << std::endl;
    std::cout << "Global bin = " << ibin << " | inv mass bin = " << invpdf[bin] << std::endl;
    return invpdf[ibin];
}

float get_rand_el(vector <int> dist) {
    int randomIndex = rand() % dist.size();
      return dist[randomIndex];
}

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_jet_E          = fs->make<TH1D>("h_jet_E"          , "E;E;Jets"                                   , 100,  0., 800.);
  h_tau_jet_pT         = fs->make<TH1D>("h_jet_pT"         , "p_{T};p_{T};Jets"                           , 100,  0., 800.);
  h_tau_jet_eta        = fs->make<TH1D>("h_jet_eta"        , "#eta;#eta;Jets"                             , 100, -5.,   5.);
  h_tau_jet_nJet       = fs->make<TH1D>("h_jet_nJet"       , "nJet;nJet;Events"                           ,  10,  0.,  10.);
  h_tau_jet_m0         = fs->make<TH1D>("h_jet_m0"         , "m_{jet};m_{jet};Jets"                       , 100,  0., 100.);
  h_tau_jet_a_m_pt     = fs->make<TH2D>("h_a_m_pT"         , "m^{a} vs p_{T}^{a};m^{a} vs p_{T}^{a};Jets" ,  24, 3.6,  15.,  30, 20., 200.);
  h_tau_jet_ma         = fs->make<TH1D>("h_jet_ma"         , "m^{a};m^{a};Jets"                           ,  24, 3.6,  15.);
  h_tau_jet_pta        = fs->make<TH1D>("h_jet_pta"        , "p_{T}^{a};p_{T}^{a};Jets"                   ,  30,  20., 200.);
  h_tau_jet_isDiTau    = fs->make<TH1D>("h_jet_isDiTau"    , "nIsDiTau;nIsDiTau;Jets"                     ,  10,  0.,  10.);
  h_tau_jet_dR         = fs->make<TH1D>("h_jet_dR"         , "dR_{a,j};dR_{a,j};Jets"                     ,  50,  0.,  0.5);
  h_tau_jet_TaudR      = fs->make<TH1D>("h_jet_TaudR"      , "dR_{#tau,#tau};dR_{#tau,#tau};Jets"         ,  50,  0.,   1.);
  h_tau_jet_Tau1dR     = fs->make<TH1D>("h_jet_Tau1dR"     , "dR_{#tau_{1},j};dR_{#tau_{1},j};Jets"       ,  50,  0.,  0.5);
  h_tau_jet_Tau2dR     = fs->make<TH1D>("h_jet_Tau2dR"     , "dR_{#tau_{2},j};dR_{#tau_{2},j};Jets"       ,  50,  0.,  0.5);
  h_tau_jet_NGenTaus  = fs->make<TH1D>("h_jet_NGenTaus"    , "N#tau^{RECO};N#tau^{RECO};Jets"             ,   5,  0.,   5.);
  h_tau_jet_NrecoTaus  = fs->make<TH1D>("h_jet_NrecoTaus"  , "N#tau^{RECO};N#tau^{RECO};Jets"             ,   5,  0.,   5.);
  h_tau_jet_recoTau1dR = fs->make<TH1D>("h_jet_recoTau1dR" , "dR_{#tau_{1}^{RECO},j};dR_{#tau_{1}^{RECO},j};Jets" ,  50,  0.,  0.5);
  h_tau_jet_recoTau2dR = fs->make<TH1D>("h_jet_recoTau2dR" , "dR_{#tau_{2}^{RECO},j};dR_{#tau_{2}^{RECO},j};Jets" ,  25,  0.,  0.5);
  h_tau_jet_n1dR       = fs->make<TH1D>("h_jet_n1dR"       , "dR_{#eta_{1},j};dR_{#eta_{1},j};Jets" ,  25,  0.,  0.5);
  h_tau_jet_n2dR       = fs->make<TH1D>("h_jet_n2dR"       , "dR_{#eta_{2},j};dR_{#eta_{2},j};Jets" ,  25,  0.,  0.5);

  tree->Branch("jetM",       &v_tau_jet_m0_);
  tree->Branch("jetPt",      &v_tau_jet_pt_);
  tree->Branch("jetPdgIds",  &v_tau_jetPdgIds_);
  tree->Branch("jetadR",     &v_tau_jetadR_);
  tree->Branch("jetIsDiTau", &v_tau_jetIsDiTau_);
  tree->Branch("a_m",        &v_tau_jet_ma_);
  tree->Branch("a_pt",       &v_tau_jet_pta_);
  tree->Branch("jetpT",      &v_tau_jet_pt_);
  tree->Branch("TaudR",      &v_tau_jetTaudR_);
  tree->Branch("Tau1dR",     &v_tau_jetTau1dR_);
  tree->Branch("Tau2dR",     &v_tau_jetTau2dR_);
  tree->Branch("NGenTaus",   &v_tau_jetNGenTaus_);
  tree->Branch("NrecoTaus",  &v_tau_jetNrecoTaus_);
  tree->Branch("recoTau1dR", &v_tau_jetrecoTau1dR_);
  tree->Branch("recoTau2dR", &v_tau_jetrecoTau2dR_);
  tree->Branch("n1dR",       &v_tau_jetn1dR_);
  tree->Branch("n2dR",       &v_tau_jetn2dR_);

  char hname[50];
  for ( unsigned iJ = 0; iJ != nJets; iJ++ ) {
    sprintf(hname, "subJet%d_E", iJ);
    tree->Branch(hname,            &v_tau_subJetE_[iJ]);
    sprintf(hname, "subJet%d_Px", iJ);
    tree->Branch(hname,            &v_tau_subJetPx_[iJ]);
    sprintf(hname, "subJet%d_Py", iJ);
    tree->Branch(hname,            &v_tau_subJetPy_[iJ]);
    sprintf(hname, "subJet%d_Pz", iJ);
    tree->Branch(hname,            &v_tau_subJetPz_[iJ]);
  }

} // branchesEvtSel_jet_dijet_tau()

// Run jet selection _____________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_tau( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken( genParticleCollectionT_, genParticles );

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauCollectionT_, taus);

  vJetIdxs.clear();
  v_tau_jetPdgIds_.clear();
  v_jetIsDiTau.clear();
  v_jetadR.clear();
  v_ma.clear();
  v_pta.clear();
  v_jetTaudR.clear();
  v_jetTau1dR.clear();
  v_jetTau2dR.clear();
  v_jetNGenTaus.clear();
  v_jetNrecoTaus.clear();
  v_jetrecoTau1dR.clear();
  v_jetrecoTau2dR.clear();
  v_jetn1dR.clear();
  v_jetn2dR.clear();

  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_tau_jetFakePhoIdxs;
  */

  unsigned int nMatchedJets = 0;
  unsigned int nMatchedGenParticles = 0;
  unsigned int nMatchedRecoTaus = 0;
  unsigned int aPdgId           = 0;
  bool MatchedPseudoScalar = false;
  float a_mass = -99.;
  float a_pt   = -99.;
  float dRa    = -99.;
  float tausdR = -99.;
  float tau1dR = -99.;
  float tau2dR = -99.;
  float recotau1dR = -99.;
  float recotau2dR = -99.;
  float n1dR   = -99.;
  float n2dR   = -99.;

  //vector <int> pT_bins  = {40, 60, 80, 100, 120, 140, 160, 180, 200};
  //vector <int> pT_occ   = {25081, 113920, 158016, 167288, 172264, 172800, 172829, 174200, 173920};
  //vector <float> m_bins = {4.74, 5.88, 7.02, 8.16, 9.3, 10.44, 11.58, 12.72, 13.86, 15};
  //vector <int> m_occ    = {147101, 147957, 142590, 137724, 134015, 128584, 125929, 124459, 121719, 120240};
  vector <int> pT_bins  = {26.0,  32.0,  38.0,  44.0,  50.0,  56.0,  62.0,  68.0,  74.0,  80.0,  86.0,  92.0,  98.0,  104.0,  110.0,  116.0,  122.0,  128.0,  134.0,  140.0,  146.0,  152.0,  158.0,  164.0,  170.0,  176.0,  182.0,  188.0,  194.0,  200.0};
  vector <float> m_bins = {4.075,  4.55,  5.025,  5.5,  5.975,  6.45,  6.925,  7.4,  7.875,  8.35,  8.825,  9.3,  9.775,  10.25,  10.725,  11.2,  11.675,  12.15,  12.625,  13.1,  13.575,  14.05,  14.525,  15};
  vector <int> occ      = {52.0,  380.0,  1020.0,  1257.0,  1665.0,  1836.0,  1981.0,  2130.0,  2144.0,  2112.0,  2103.0,  2163.0,  2183.0,  2160.0,  2153.0,  2213.0,  2139.0,  2076.0,  2251.0,  2308.0,  2147.0,  2184.0,  2241.0,  2039.0,  2148.0,  2271.0,  2189.0,  2148.0,  2249.0,  2089.0,  58.0,  478.0,  1116.0,  1434.0,  1878.0,  2093.0,  2005.0,  2277.0,  2370.0,  2332.0,  2490.0,  2433.0,  2505.0,  2354.0,  2428.0,  2488.0,  2555.0,  2512.0,  2544.0,  2498.0,  2545.0,  2571.0,  2371.0,  2592.0,  2462.0,  2425.0,  2458.0,  2522.0,  2512.0,  2505.0,  46.0,  443.0,  1109.0,  1468.0,  1921.0,  2060.0,  2246.0,  2354.0,  2316.0,  2265.0,  2417.0,  2383.0,  2463.0,  2459.0,  2541.0,  2499.0,  2454.0,  2403.0,  2450.0,  2602.0,  2285.0,  2503.0,  2564.0,  2435.0,  2554.0,  2554.0,  2409.0,  2393.0,  2629.0,  2383.0,  48.0,  458.0,  1059.0,  1583.0,  1724.0,  1944.0,  2114.0,  2356.0,  2288.0,  2296.0,  2387.0,  2319.0,  2309.0,  2518.0,  2425.0,  2395.0,  2376.0,  2390.0,  2619.0,  2599.0,  2349.0,  2421.0,  2474.0,  2552.0,  2502.0,  2469.0,  2382.0,  2462.0,  2559.0,  2425.0,  23.0,  324.0,  945.0,  1456.0,  1822.0,  2024.0,  2268.0,  2286.0,  2353.0,  2476.0,  2417.0,  2397.0,  2388.0,  2648.0,  2559.0,  2462.0,  2446.0,  2550.0,  2315.0,  2393.0,  2414.0,  2428.0,  2307.0,  2572.0,  2485.0,  2596.0,  2481.0,  2495.0,  2328.0,  2502.0,  33.0,  341.0,  933.0,  1443.0,  1753.0,  1837.0,  2139.0,  2239.0,  2193.0,  2381.0,  2348.0,  2362.0,  2452.0,  2368.0,  2590.0,  2514.0,  2422.0,  2486.0,  2517.0,  2407.0,  2426.0,  2517.0,  2481.0,  2512.0,  2543.0,  2534.0,  2514.0,  2539.0,  2371.0,  2459.0,  26.0,  290.0,  834.0,  1459.0,  1802.0,  1904.0,  2147.0,  2205.0,  2340.0,  2344.0,  2412.0,  2452.0,  2405.0,  2341.0,  2492.0,  2431.0,  2423.0,  2438.0,  2333.0,  2398.0,  2622.0,  2498.0,  2521.0,  2375.0,  2452.0,  2457.0,  2521.0,  2512.0,  2305.0,  2429.0,  22.0,  205.0,  853.0,  1416.0,  1781.0,  1948.0,  2092.0,  2245.0,  2187.0,  2399.0,  2419.0,  2357.0,  2460.0,  2443.0,  2477.0,  2516.0,  2519.0,  2432.0,  2492.0,  2425.0,  2466.0,  2557.0,  2453.0,  2471.0,  2485.0,  2499.0,  2476.0,  2464.0,  2460.0,  2460.0,  12.0,  176.0,  668.0,  1306.0,  1723.0,  1960.0,  2120.0,  2294.0,  2316.0,  2248.0,  2324.0,  2388.0,  2477.0,  2275.0,  2419.0,  2502.0,  2479.0,  2441.0,  2421.0,  2429.0,  2449.0,  2479.0,  2577.0,  2556.0,  2523.0,  2331.0,  2402.0,  2444.0,  2616.0,  2383.0,  10.0,  152.0,  562.0,  1180.0,  1667.0,  1946.0,  2132.0,  2101.0,  2193.0,  2332.0,  2451.0,  2449.0,  2330.0,  2409.0,  2481.0,  2499.0,  2480.0,  2349.0,  2439.0,  2421.0,  2437.0,  2536.0,  2413.0,  2550.0,  2418.0,  2480.0,  2379.0,  2528.0,  2449.0,  2594.0,  5.0,  123.0,  534.0,  1086.0,  1491.0,  1824.0,  2149.0,  2350.0,  2275.0,  2385.0,  2297.0,  2278.0,  2344.0,  2371.0,  2492.0,  2563.0,  2420.0,  2506.0,  2465.0,  2518.0,  2366.0,  2475.0,  2506.0,  2420.0,  2529.0,  2537.0,  2423.0,  2479.0,  2487.0,  2551.0,  22.0,  102.0,  434.0,  1003.0,  1543.0,  1911.0,  1987.0,  2191.0,  2381.0,  2322.0,  2540.0,  2278.0,  2514.0,  2337.0,  2431.0,  2486.0,  2464.0,  2382.0,  2360.0,  2519.0,  2399.0,  2509.0,  2464.0,  2412.0,  2443.0,  2555.0,  2511.0,  2589.0,  2497.0,  2436.0,  14.0,  127.0,  396.0,  863.0,  1483.0,  1874.0,  2039.0,  2209.0,  2274.0,  2303.0,  2385.0,  2323.0,  2464.0,  2375.0,  2436.0,  2367.0,  2432.0,  2375.0,  2512.0,  2412.0,  2351.0,  2433.0,  2387.0,  2441.0,  2557.0,  2476.0,  2408.0,  2527.0,  2372.0,  2493.0,  16.0,  131.0,  432.0,  816.0,  1434.0,  1893.0,  2001.0,  2071.0,  2140.0,  2290.0,  2393.0,  2340.0,  2496.0,  2308.0,  2290.0,  2360.0,  2477.0,  2381.0,  2398.0,  2368.0,  2417.0,  2452.0,  2408.0,  2442.0,  2452.0,  2385.0,  2424.0,  2433.0,  2332.0,  2537.0,  26.0,  100.0,  397.0,  743.0,  1319.0,  1694.0,  1971.0,  2169.0,  2246.0,  2277.0,  2358.0,  2373.0,  2350.0,  2402.0,  2365.0,  2450.0,  2437.0,  2430.0,  2405.0,  2431.0,  2408.0,  2522.0,  2429.0,  2442.0,  2490.0,  2535.0,  2445.0,  2524.0,  2505.0,  2532.0,  10.0,  143.0,  389.0,  778.0,  1189.0,  1639.0,  1960.0,  2125.0,  2416.0,  2399.0,  2403.0,  2285.0,  2363.0,  2413.0,  2421.0,  2414.0,  2454.0,  2432.0,  2492.0,  2344.0,  2409.0,  2518.0,  2459.0,  2488.0,  2459.0,  2375.0,  2561.0,  2452.0,  2395.0,  2334.0,  20.0,  159.0,  401.0,  740.0,  1178.0,  1578.0,  1845.0,  2023.0,  2195.0,  2356.0,  2244.0,  2388.0,  2382.0,  2495.0,  2390.0,  2354.0,  2392.0,  2478.0,  2488.0,  2357.0,  2361.0,  2455.0,  2458.0,  2378.0,  2456.0,  2407.0,  2496.0,  2445.0,  2428.0,  2462.0,  32.0,  158.0,  385.0,  703.0,  1184.0,  1549.0,  1994.0,  2102.0,  2134.0,  2307.0,  2294.0,  2297.0,  2424.0,  2400.0,  2381.0,  2443.0,  2474.0,  2473.0,  2497.0,  2446.0,  2636.0,  2574.0,  2341.0,  2379.0,  2382.0,  2444.0,  2475.0,  2575.0,  2456.0,  2469.0,  22.0,  133.0,  407.0,  741.0,  1163.0,  1522.0,  1830.0,  2046.0,  2167.0,  2329.0,  2328.0,  2376.0,  2379.0,  2401.0,  2374.0,  2501.0,  2337.0,  2506.0,  2405.0,  2558.0,  2433.0,  2541.0,  2663.0,  2488.0,  2350.0,  2527.0,  2487.0,  2516.0,  2498.0,  2478.0,  26.0,  137.0,  371.0,  759.0,  1138.0,  1522.0,  1885.0,  1993.0,  2130.0,  2228.0,  2246.0,  2232.0,  2373.0,  2416.0,  2471.0,  2514.0,  2438.0,  2556.0,  2453.0,  2481.0,  2338.0,  2411.0,  2532.0,  2430.0,  2512.0,  2447.0,  2388.0,  2424.0,  2550.0,  2367.0,  22.0,  183.0,  400.0,  696.0,  1105.0,  1516.0,  1834.0,  2099.0,  2129.0,  2256.0,  2264.0,  2290.0,  2432.0,  2400.0,  2473.0,  2431.0,  2453.0,  2415.0,  2472.0,  2419.0,  2398.0,  2551.0,  2448.0,  2472.0,  2468.0,  2623.0,  2423.0,  2454.0,  2521.0,  2412.0,  36.0,  136.0,  411.0,  830.0,  1150.0,  1502.0,  1828.0,  2157.0,  2215.0,  2197.0,  2247.0,  2058.0,  2319.0,  2468.0,  2457.0,  2311.0,  2439.0,  2353.0,  2445.0,  2403.0,  2509.0,  2427.0,  2484.0,  2472.0,  2440.0,  2481.0,  2516.0,  2534.0,  2595.0,  2440.0,  30.0,  178.0,  427.0,  743.0,  1068.0,  1596.0,  1758.0,  2040.0,  2202.0,  2294.0,  2322.0,  2319.0,  2372.0,  2329.0,  2374.0,  2424.0,  2511.0,  2409.0,  2609.0,  2410.0,  2442.0,  2555.0,  2500.0,  2465.0,  2436.0,  2489.0,  2487.0,  2442.0,  2497.0,  2414.0,  34.0,  146.0,  443.0,  761.0,  1119.0,  1548.0,  1883.0,  2164.0,  2238.0,  2237.0,  2186.0,  2276.0,  2276.0,  2272.0,  2422.0,  2356.0,  2521.0,  2490.0,  2397.0,  2561.0,  2414.0,  2452.0,  2456.0,  2423.0,  2390.0,  2360.0,  2399.0,  2412.0,  2442.0,  2459.0};
  vector <float> m_invpdf  = get_inverse_pdf(m_occ);
  vector <float> pT_invpdf = get_inverse_pdf(pT_occ);
  vector <float> invpdf    = get_inverse_pdf(occ);

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  if ( debug ) std::cout << " JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
  // Loop over jets
  for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
    reco::PFJetRef iJet( jets, iJ );
    if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
    if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
    bool passedGenSel = false;
    unsigned int iGenParticle = 0; 
    unsigned int NTau1Daughters = 0;
    unsigned int NTau2Daughters = 0;
    for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin(); iGen != genParticles->end(); ++iGen) {
      if ( abs(iGen->pdgId()) != 15 ) continue;
      if (iGen->numberOfMothers() != 1) continue;
      if (iGen->mother()->pdgId() != 25) continue;
      ++iGenParticle;
      if ( debug ) std::cout << "   GEN particle " << iGenParticle << " -> status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " nMoms: " <<iGen->numberOfMothers() << " | pt: "<< iGen->pt() << " eta: " <<iGen->eta() << " phi: " <<iGen->phi() << std::endl;
      float dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      if ( dR > 0.4 ) continue;

      float rand_sampler_pT  = rand() / float(RAND_MAX);
      float rand_sampler_m  = rand() / float(RAND_MAX);
      int pT_gen   = iGen->mother()->pt();
      float m_gen  = iGen->mother()->mass();
      float pT_wgt = lookup_pt_invpdf(pT_gen, pT_bins, pT_invpdf);
      if (debug) std::cout << " wgt pT " << pT_wgt  << " | rand_sampler_pT " << rand_sampler_pT << std::endl;
      float m_wgt = lookup_mass_invpdf(m_gen, m_bins, m_invpdf);
      if (debug) std::cout << " wgt m " << m_wgt  << " | rand_sampler_m "<< rand_sampler_m << std::endl;
      //if (rand_sampler_pT > pT_wgt) continue;
      //if (rand_sampler_m  > m_wgt) continue;

      passedGenSel = true;

      if (debug ) std::cout << "  >>>>>> Jet [" << iJ << "] matched particle [" << iGenParticle << "] -> pdgId: " << std::abs(iGen->pdgId()) << " | dR: " << dR << "| Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
      if (nMatchedGenParticles == 0) {
        if ( iGen->numberOfMothers() == 1 ) {
          aPdgId = std::abs(iGen->mother()->pdgId());
          if ( abs(iGen->mother()->pdgId()) == 25 && iGen->mother()->mass() < 15) {
            MatchedPseudoScalar = true;
            a_mass = iGen->mother()->mass();
            a_pt   = iGen->mother()->pt();
            dRa = reco::deltaR( iJet->eta(),iJet->phi(), iGen->mother()->eta(),iGen->mother()->phi() );
            if ( iGen->mother()->numberOfDaughters() == 2 ){ 
              if (abs(iGen->mother()->daughter(0)->pdgId()) == 15 && abs(iGen->mother()->daughter(1)->pdgId()) == 15){ 
                tausdR = reco::deltaR( iGen->mother()->daughter(0)->eta(),iGen->mother()->daughter(0)->phi(), iGen->mother()->daughter(1)->eta(),iGen->mother()->daughter(1)->phi() );
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
              }
            }
          }
        }
        if (debug ) std::cout << "   >>>>>> n1 dR = " << n1dR << " , n2 dR = " << n2dR << std::endl;
      }
      ++nMatchedGenParticles;
    } // primary gen particles
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
      v_tau_jetPdgIds_.push_back( aPdgId );
      v_jetadR.push_back( dRa );
      v_jetTaudR.push_back( tausdR );
      v_ma.push_back( a_mass );
      v_pta.push_back( a_pt );
      v_jetTau1dR.push_back( tau1dR );
      v_jetTau2dR.push_back( tau2dR );
      v_jetn1dR.push_back( n1dR );
      v_jetn2dR.push_back( n2dR );
      v_jetNGenTaus.push_back( nMatchedGenParticles );
      v_jetNrecoTaus.push_back( nMatchedRecoTaus );
      v_jetrecoTau1dR.push_back( recotau1dR );
      v_jetrecoTau2dR.push_back( recotau2dR );
      v_jetIsDiTau.push_back( MatchedPseudoScalar );

    }

  } // reco jets
  if ( debug ) std::cout << " Matched jets " << nMatchedJets << std::endl;

  // Check jet multiplicity
  if ( nMatchedJets < 2 ) return false;

  if ( debug ) std::cout << " >> has_jet_dijet_tau: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_jet_nJet->Fill( vJetIdxs.size() );

  v_tau_jet_pt_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jet_ma_.clear();
  v_tau_jet_pta_.clear();
  v_tau_jetIsDiTau_.clear();
  v_tau_jetadR_.clear();
  v_tau_jetTaudR_.clear();
  v_tau_jetTau1dR_.clear();
  v_tau_jetTau2dR_.clear();
  v_tau_jetNGenTaus_.clear();
  v_tau_jetNrecoTaus_.clear();
  v_tau_jetrecoTau1dR_.clear();
  v_tau_jetrecoTau2dR_.clear();
  v_tau_jetn1dR_.clear();
  v_tau_jetn2dR_.clear();
 
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_jet_eta->Fill( iJet->eta() );
    h_tau_jet_E->Fill( iJet->energy() );
    h_tau_jet_m0->Fill( iJet->mass() );
    h_tau_jet_isDiTau->Fill( v_jetIsDiTau[iJ] );
    h_tau_jet_dR->Fill( v_jetadR[iJ] );
    h_tau_jet_ma->Fill( v_ma[iJ] );
    h_tau_jet_pta->Fill( v_pta[iJ] );
    h_tau_jet_a_m_pt->Fill( v_ma[iJ], v_pta[iJ] );
    h_tau_jet_TaudR->Fill( v_jetTaudR[iJ] );
    h_tau_jet_Tau1dR->Fill( v_jetTau1dR[iJ] );
    h_tau_jet_Tau2dR->Fill( v_jetTau2dR[iJ] );
    h_tau_jet_NGenTaus->Fill( v_jetNGenTaus[iJ] );
    h_tau_jet_NrecoTaus->Fill( v_jetNrecoTaus[iJ] );
    h_tau_jet_recoTau1dR->Fill( v_jetrecoTau1dR[iJ] );
    h_tau_jet_recoTau2dR->Fill( v_jetrecoTau2dR[iJ] );
    h_tau_jet_n1dR->Fill( v_jetn1dR[iJ] );
    h_tau_jet_n2dR->Fill( v_jetn2dR[iJ] );

    // Fill branches 
    v_tau_jet_pt_.push_back( iJet->pt() );
    v_tau_jet_m0_.push_back( iJet->mass() );
    v_tau_jet_ma_.push_back( v_ma[iJ] );
    v_tau_jet_pta_.push_back( v_pta[iJ] );
    v_tau_jetIsDiTau_.push_back( v_jetIsDiTau[iJ] );
    v_tau_jetadR_.push_back( v_jetadR[iJ] );
    v_tau_jetTaudR_.push_back( v_jetTaudR[iJ] );
    v_tau_jetTau1dR_.push_back( v_jetTau1dR[iJ] );
    v_tau_jetTau2dR_.push_back( v_jetTau2dR[iJ] );
    v_tau_jetNGenTaus_.push_back( v_jetNGenTaus[iJ] );
    v_tau_jetNrecoTaus_.push_back( v_jetNrecoTaus[iJ] );
    v_tau_jetrecoTau1dR_.push_back( v_jetrecoTau1dR[iJ] );
    v_tau_jetrecoTau2dR_.push_back( v_jetrecoTau2dR[iJ] );
    v_tau_jetn1dR_.push_back( v_jetn1dR[iJ] );
    v_tau_jetn2dR_.push_back( v_jetn2dR[iJ] );

    // Gen jet constituents
    v_tau_subJetE_[iJ].clear();
    v_tau_subJetPx_[iJ].clear();
    v_tau_subJetPy_[iJ].clear();
    v_tau_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    if ( debug ) std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      if ( debug ) std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }

} // fillEvtSel_jet_dijet_tau()
