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
    unsigned int ibin = 0;
    unsigned int ibinx = 0;
    unsigned int ibiny = 0;
    unsigned int m1  = M_bins.size();
    unsigned int pt1 = pT_bins.size();
    unsigned int inv = invpdf.size();
    bool found_mass = false;
    bool found_end  = false;
    if (debug) std::cout << "Nbin = " << inv << " | Nbinx = " << m1 << " , Nbiny = " << pt1 << std::endl;
    for (unsigned int ibx = 0; ibx < m1; ibx++) {
        if (found_mass || found_end) { break; }
        for (unsigned int iby = 0; iby < pt1; iby++) {
            ibin = (ibx*pt1)+ iby;
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
    if (debug) std::cout << "mass gen   = " << Mgen  << " | bin = " << ibinx << " | mass bin = " << M_bins[ibinx]  << std::endl;
    if (debug) std::cout << "pt gen     = " << pTgen << " | bin = " << ibiny << " | pt bin   = " << pT_bins[ibiny] << std::endl;
    if (debug) std::cout << "Global bin = " << ibin <<  " | inv mass bin = " << invpdf[ibin] << std::endl;
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
  h_tau_jet_a_m_pt     = fs->make<TH2D>("h_a_m_pT"         , "m^{a} vs p_{T}^{a};m^{a} vs p_{T}^{a};Jets" ,  11, 3.7,  15.,  34, 30., 200.);
  h_tau_jet_ma         = fs->make<TH1D>("h_jet_ma"         , "m^{a};m^{a};Jets"                           ,  11, 3.7,  15.);
  h_tau_jet_pta        = fs->make<TH1D>("h_jet_pta"        , "p_{T}^{a};p_{T}^{a};Jets"                   ,  34, 30., 200.);
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
  //vector <int> pT_bins  = {26, 32, 38, 44, 50, 56, 62, 68, 74, 80, 86, 92, 98, 104, 110, 116, 122, 128, 134, 140, 146, 152, 158, 164, 170, 176, 182, 188, 194, 200};
  //vector <float> m_bins = {4.075, 4.55, 5.025, 5.5, 5.975, 6.45, 6.925, 7.4, 7.875, 8.35, 8.825, 9.3, 9.775, 10.25, 10.725, 11.2, 11.675, 12.15, 12.625, 13.1, 13.575, 14.05, 14.525, 15};
  //vector <int> occ      = {52, 380, 1020, 1257, 1665, 1836, 1981, 2130, 2144, 2112, 2103, 2163, 2183, 2160, 2153, 2213, 2139, 2076, 2251, 2308, 2147, 2184, 2241, 2039, 2148, 2271, 2189, 2148, 2249, 2089, 58, 478, 1116, 1434, 1878, 2093, 2005, 2277, 2370, 2332, 2490, 2433, 2505, 2354, 2428, 2488, 2555, 2512, 2544, 2498, 2545, 2571, 2371, 2592, 2462, 2425, 2458, 2522, 2512, 2505, 46, 443, 1109, 1468, 1921, 2060, 2246, 2354, 2316, 2265, 2417, 2383, 2463, 2459, 2541, 2499, 2454, 2403, 2450, 2602, 2285, 2503, 2564, 2435, 2554, 2554, 2409, 2393, 2629, 2383, 48, 458, 1059, 1583, 1724, 1944, 2114, 2356, 2288, 2296, 2387, 2319, 2309, 2518, 2425, 2395, 2376, 2390, 2619, 2599, 2349, 2421, 2474, 2552, 2502, 2469, 2382, 2462, 2559, 2425, 23, 324, 945, 1456, 1822, 2024, 2268, 2286, 2353, 2476, 2417, 2397, 2388, 2648, 2559, 2462, 2446, 2550, 2315, 2393, 2414, 2428, 2307, 2572, 2485, 2596, 2481, 2495, 2328, 2502, 33, 341, 933, 1443, 1753, 1837, 2139, 2239, 2193, 2381, 2348, 2362, 2452, 2368, 2590, 2514, 2422, 2486, 2517, 2407, 2426, 2517, 2481, 2512, 2543, 2534, 2514, 2539, 2371, 2459, 26, 290, 834, 1459, 1802, 1904, 2147, 2205, 2340, 2344, 2412, 2452, 2405, 2341, 2492, 2431, 2423, 2438, 2333, 2398, 2622, 2498, 2521, 2375, 2452, 2457, 2521, 2512, 2305, 2429, 22, 205, 853, 1416, 1781, 1948, 2092, 2245, 2187, 2399, 2419, 2357, 2460, 2443, 2477, 2516, 2519, 2432, 2492, 2425, 2466, 2557, 2453, 2471, 2485, 2499, 2476, 2464, 2460, 2460, 12, 176, 668, 1306, 1723, 1960, 2120, 2294, 2316, 2248, 2324, 2388, 2477, 2275, 2419, 2502, 2479, 2441, 2421, 2429, 2449, 2479, 2577, 2556, 2523, 2331, 2402, 2444, 2616, 2383, 10, 152, 562, 1180, 1667, 1946, 2132, 2101, 2193, 2332, 2451, 2449, 2330, 2409, 2481, 2499, 2480, 2349, 2439, 2421, 2437, 2536, 2413, 2550, 2418, 2480, 2379, 2528, 2449, 2594, 5, 123, 534, 1086, 1491, 1824, 2149, 2350, 2275, 2385, 2297, 2278, 2344, 2371, 2492, 2563, 2420, 2506, 2465, 2518, 2366, 2475, 2506, 2420, 2529, 2537, 2423, 2479, 2487, 2551, 22, 102, 434, 1003, 1543, 1911, 1987, 2191, 2381, 2322, 2540, 2278, 2514, 2337, 2431, 2486, 2464, 2382, 2360, 2519, 2399, 2509, 2464, 2412, 2443, 2555, 2511, 2589, 2497, 2436, 14, 127, 396, 863, 1483, 1874, 2039, 2209, 2274, 2303, 2385, 2323, 2464, 2375, 2436, 2367, 2432, 2375, 2512, 2412, 2351, 2433, 2387, 2441, 2557, 2476, 2408, 2527, 2372, 2493, 16, 131, 432, 816, 1434, 1893, 2001, 2071, 2140, 2290, 2393, 2340, 2496, 2308, 2290, 2360, 2477, 2381, 2398, 2368, 2417, 2452, 2408, 2442, 2452, 2385, 2424, 2433, 2332, 2537, 26, 100, 397, 743, 1319, 1694, 1971, 2169, 2246, 2277, 2358, 2373, 2350, 2402, 2365, 2450, 2437, 2430, 2405, 2431, 2408, 2522, 2429, 2442, 2490, 2535, 2445, 2524, 2505, 2532, 10, 143, 389, 778, 1189, 1639, 1960, 2125, 2416, 2399, 2403, 2285, 2363, 2413, 2421, 2414, 2454, 2432, 2492, 2344, 2409, 2518, 2459, 2488, 2459, 2375, 2561, 2452, 2395, 2334, 20, 159, 401, 740, 1178, 1578, 1845, 2023, 2195, 2356, 2244, 2388, 2382, 2495, 2390, 2354, 2392, 2478, 2488, 2357, 2361, 2455, 2458, 2378, 2456, 2407, 2496, 2445, 2428, 2462, 32, 158, 385, 703, 1184, 1549, 1994, 2102, 2134, 2307, 2294, 2297, 2424, 2400, 2381, 2443, 2474, 2473, 2497, 2446, 2636, 2574, 2341, 2379, 2382, 2444, 2475, 2575, 2456, 2469, 22, 133, 407, 741, 1163, 1522, 1830, 2046, 2167, 2329, 2328, 2376, 2379, 2401, 2374, 2501, 2337, 2506, 2405, 2558, 2433, 2541, 2663, 2488, 2350, 2527, 2487, 2516, 2498, 2478, 26, 137, 371, 759, 1138, 1522, 1885, 1993, 2130, 2228, 2246, 2232, 2373, 2416, 2471, 2514, 2438, 2556, 2453, 2481, 2338, 2411, 2532, 2430, 2512, 2447, 2388, 2424, 2550, 2367, 22, 183, 400, 696, 1105, 1516, 1834, 2099, 2129, 2256, 2264, 2290, 2432, 2400, 2473, 2431, 2453, 2415, 2472, 2419, 2398, 2551, 2448, 2472, 2468, 2623, 2423, 2454, 2521, 2412, 36, 136, 411, 830, 1150, 1502, 1828, 2157, 2215, 2197, 2247, 2058, 2319, 2468, 2457, 2311, 2439, 2353, 2445, 2403, 2509, 2427, 2484, 2472, 2440, 2481, 2516, 2534, 2595, 2440, 30, 178, 427, 743, 1068, 1596, 1758, 2040, 2202, 2294, 2322, 2319, 2372, 2329, 2374, 2424, 2511, 2409, 2609, 2410, 2442, 2555, 2500, 2465, 2436, 2489, 2487, 2442, 2497, 2414, 34, 146, 443, 761, 1119, 1548, 1883, 2164, 2238, 2237, 2186, 2276, 2276, 2272, 2422, 2356, 2521, 2490, 2397, 2561, 2414, 2452, 2456, 2423, 2390, 2360, 2399, 2412, 2442, 2459};
  vector <int> pT_bins  = {35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200};
  vector <float> m_bins = {4.72727272727, 5.75454545455, 6.78181818182, 7.80909090909, 8.83636363636, 9.86363636364, 10.8909090909, 11.9181818182, 12.9454545455, 13.9727272727, 15.0};
  vector <int> occ      = {2817, 4769, 5858, 6957, 7584, 7547, 8037, 8447, 8317, 8693, 8654, 8676, 8719, 8828, 8925, 8853, 8871, 9062, 8819, 9016, 8923, 8928, 8860, 9061, 8863, 8986, 9081, 9010, 8907, 8859, 8711, 8968, 8886, 9055, 2798, 4344, 5678, 6517, 7366, 7463, 7926, 8176, 8408, 8483, 8488, 8811, 8920, 8591, 8810, 8789, 8854, 9015, 9028, 9030, 8822, 8914, 8645, 9065, 8818, 9086, 9168, 9042, 9099, 8909, 8880, 8899, 8649, 8919, 2291, 4086, 5635, 6559, 7167, 7875, 7962, 8350, 8199, 8351, 8459, 8422, 8579, 8580, 8963, 8918, 8974, 8987, 8939, 8892, 9101, 9041, 8779, 8721, 9006, 8862, 8962, 8861, 8907, 8675, 9052, 8952, 9023, 8970, 1853, 3592, 5281, 6328, 7003, 7398, 7792, 8158, 8242, 8480, 8339, 8551, 8655, 8556, 8811, 8665, 8673, 8685, 8934, 8931, 8830, 8936, 8591, 9160, 8815, 8973, 8808, 8789, 8879, 8993, 8873, 8876, 8870, 8862, 1081, 2942, 4651, 5971, 6800, 7427, 7700, 8285, 8270, 8307, 8523, 8505, 8484, 8544, 8746, 8984, 8849, 9140, 8871, 8945, 8911, 8748, 8981, 8702, 8811, 8675, 8975, 8827, 8929, 8915, 8889, 8774, 9053, 8747, 908, 2330, 3789, 5370, 6568, 7293, 7739, 8015, 8031, 8273, 8664, 8505, 8629, 8750, 8671, 9071, 8797, 8762, 8871, 8832, 8825, 8831, 8936, 8734, 9036, 8785, 8862, 9035, 8779, 8996, 8914, 9145, 8516, 8696, 983, 1928, 3354, 4807, 6243, 6890, 7529, 7894, 8124, 8419, 8419, 8308, 8622, 8645, 8638, 8760, 8849, 8758, 8813, 8710, 8546, 8975, 8879, 8772, 8831, 8914, 8979, 8730, 8717, 8974, 8797, 9091, 8719, 8926, 937, 1909, 3135, 4541, 5682, 6645, 7538, 7608, 8009, 8318, 8574, 8377, 8544, 8574, 8562, 8806, 8685, 8479, 9018, 8627, 8650, 8808, 8700, 8742, 8636, 9143, 8662, 8872, 8666, 9113, 8972, 8914, 8986, 8831, 959, 1887, 2934, 4250, 5427, 6502, 7099, 7257, 7771, 8118, 8308, 8561, 8481, 8481, 8571, 8734, 8605, 8811, 8642, 8906, 8879, 8694, 8642, 8888, 8645, 8850, 8697, 8792, 8840, 9027, 8944, 8986, 8907, 9037, 984, 1981, 3103, 4232, 5227, 6360, 7246, 7670, 7916, 7989, 8437, 8580, 8369, 8385, 8566, 8573, 8761, 8936, 8886, 8832, 8904, 8831, 8741, 8739, 9076, 9049, 8761, 8821, 8836, 8922, 9077, 8790, 8733, 8703, 1059, 2025, 3094, 4232, 5424, 6284, 7221, 7860, 8187, 7985, 8357, 8053, 8649, 8541, 8539, 8847, 8690, 8608, 8668, 8878, 8720, 8765, 9002, 8681, 8834, 9108, 8681, 8961, 9051, 8935, 8732, 8986, 9065, 8837};
  //vector <float> m_invpdf  = get_inverse_pdf(m_occ);
  //vector <float> pT_invpdf = get_inverse_pdf(pT_occ);
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

      bool unbiasing = false;
      if ( unbiasing ) {
          // 2D
          float rand_sampler  = rand() / float(RAND_MAX);
          float pT_gen        = iGen->mother()->pt();
          float m_gen         = iGen->mother()->mass();
          float wgt           = lookup_invpdf(m_gen, m_bins, pT_gen, pT_bins, invpdf);
          if (debug) std::cout << " wgt " << wgt  << " | rand_sampler " << rand_sampler << std::endl;
          if (rand_sampler > wgt) continue;
          // 2 x 1D
          //float rand_sampler_pT  = rand() / float(RAND_MAX);
          //float rand_sampler_m   = rand() / float(RAND_MAX);
          //float pT_wgt           = lookup_pt_invpdf(pT_gen, pT_bins, pT_invpdf);
          //float m_wgt            = lookup_mass_invpdf(m_gen, m_bins, m_invpdf);
          //if (debug) std::cout << " wgt pT " << pT_wgt  << " | rand_sampler_pT " << rand_sampler_pT << std::endl;
          //if (debug) std::cout << " wgt m " << m_wgt  << " | rand_sampler_m "<< rand_sampler_m << std::endl;
          //if (rand_sampler_pT > pT_wgt) continue;
          //if (rand_sampler_m  > m_wgt) continue;
      }
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
