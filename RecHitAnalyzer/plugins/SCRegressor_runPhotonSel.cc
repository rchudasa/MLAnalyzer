#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesPhotonSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  hSC_pT = fs->make<TH1F>("SC_pT", "Pt", 27, 15., 150.);
  hMinDRgenRecoPho = fs->make<TH1F>("minDRgenRecoPho", "#DeltaR(#gamma_{gen},#gamma_{reco})_{min};#DeltaR;N", 100, 0., 25*0.0174);
  hMinDRrecoPtoGenPt = fs->make<TH1F>("minDRrecoPtoGenPt", "#DeltaR(#gamma_{gen},#gamma_{reco})_{min}, p_{T,reco}/p_{T,gen};p_{T,reco}/p_{T,gen};N", 60, -10., 10.);

  tree->Branch("SC_mass",   &vSC_mass_);
  tree->Branch("SC_DR",     &vSC_DR_);
  tree->Branch("SC_E",      &vSC_E_);
  tree->Branch("SC_pT",     &vSC_pT_);
  tree->Branch("SC_eta",    &vSC_eta_);
  tree->Branch("SC_phi",    &vSC_phi_);

}

// Define struct to handle mapping for gen pho<->matched reco photons<->matched presel photons
struct pho_map {
  unsigned int idx;
  std::vector<unsigned int> matchedRecoPhoIdxs;
  std::vector<unsigned int> matchedPreselPhoIdxs;
};
std::vector<pho_map> vPhos;

// Run event selection ___________________________________________________________________//
bool SCRegressor::runPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  ////////// Gen-level validation //////////

  // Identify particle gun gen phos from the event
  // (Assumed that event contains multiple pho guns)
  // WARNING: the MINIAODSIM prunedGenParticles collection clips particles with pt < 10 GeV from the truth table!
  float dR;
  std::vector<unsigned int> vGenPhoIdxs;
  for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );

    if ( iGen->pdgId() != 22 ) continue; //photon
    //if ( iGen->pdgId() != 11 ) continue; //electron
    if ( iGen->status() != 1 ) continue;

    vGenPhoIdxs.push_back( iG );

  } // genParticles
  if ( debug ) std::cout << " >> vGenPhoIdxs.size: " << vGenPhoIdxs.size() << std::endl;
  if ( vGenPhoIdxs.empty() ) return false;

  ////////// Build gen pho-reco photon mapping //////////

  float ptCut = 15., etaCut = 1.44;
  //float ptCut = 10., etaCut = 1.44;

  // Create mapping between gen pho<->matched reco photons<->matched presel photons
  // For each gen pho, find "reco" photons matched to it,
  // then check if that reco photon passes photon preselection criteria
  float minDR = 100.;
  float minDR_fpt = -10.;
  int minDR_idx = -1;
  vPhos.clear();
  // Loop over valid gen pho idxs
  for ( auto& iG : vGenPhoIdxs ) {

    reco::GenParticleRef iGenPho( genParticles, iG );
    if ( debug ) std::cout << " >> genPho[" << iG << "]" << " pt:" << iGenPho->pt() << " eta:" << iGenPho->eta() << std::endl;

    std::vector<unsigned int> vMatchedRecoPhoIdxs;
    std::vector<unsigned int> vMatchedPreselPhoIdxs;

    // Do dR match to closest reco photon
    minDR = 100.;
    minDR_fpt = -10;
    minDR_idx = -1;
    for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

      PhotonRef iRecoPho( photons, iP );

      // Definition of a "reco" photon--highly subject to interpretation
      //if ( iRecoPho->pt() < 5. ) continue;
      if ( iRecoPho->pt() < 10. ) continue; // pat/miniaod threshold is pt > 10 GeV

      dR = reco::deltaR( iRecoPho->eta(),iRecoPho->phi(), iGenPho->eta(),iGenPho->phi() );
      if ( dR > minDR ) continue;

      minDR = dR;
      minDR_idx = iP;
      minDR_fpt = iRecoPho->pt()/iGenPho->pt();
      if ( debug ) std::cout << "   >> minDR_idx:" << minDR_idx << " " << minDR << " pt:" << iRecoPho->pt() << " eta:" << iRecoPho->eta() << std::endl;

    } // reco photons
    hMinDRgenRecoPho->Fill( minDR );
    hMinDRrecoPtoGenPt->Fill( minDR_fpt );

    // Require minimum dR to declare match
    // Protects against matching to PU, although not a major issue since these will likely fail preselection
    // minDR only needs to be generous enough so that one of the gen photons match to a reco photon for analysis
    if ( minDR > 0.04 ) continue;

    // Declare reco photon matching to gen pho: only store unique reco idxs
    if ( std::find(vMatchedRecoPhoIdxs.begin(), vMatchedRecoPhoIdxs.end(), minDR_idx) != vMatchedRecoPhoIdxs.end() ) continue;
    vMatchedRecoPhoIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << "   >> !minDR_idx:" << minDR_idx << " f_pt(reco/gen):" << minDR_fpt << std::endl;

    // Check if matched reco photon passes preselection:
    PhotonRef iRecoPho( photons, minDR_idx );
    if ( std::abs(iRecoPho->pt()) <= ptCut ) continue;
    if ( std::abs(iRecoPho->eta()) >= etaCut ) continue;

    if ( iRecoPho->full5x5_r9() <= 0.5 ) continue;
    if ( iRecoPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iRecoPho->hasPixelSeed() == true ) continue;
    ///*
    //if ( iRecoPho->passElectronVeto() == true ) continue;
    //if ( iRecoPho->userFloat("phoChargedIsolation")/std::abs(iRecoPho->pt()) > 0.3 ) continue;

    ///*
    if ( iRecoPho->full5x5_r9() <= 0.85 ) {
      if ( iRecoPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
      //if ( iRecoPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      //if ( iRecoPho->photonIso() >= 4.0 ) continue;
      if ( iRecoPho->trkSumPtHollowConeDR03() >= 6. ) continue;
      //if ( iRecoPho->trackIso() >= 6. ) continue;
    }
    //*/
    vMatchedPreselPhoIdxs.push_back( minDR_idx );
    if ( debug ) std::cout << " >> presel photon: pT: " << iRecoPho->pt() << " eta: " << iRecoPho->eta() << std::endl;

    // Store this mapping
    pho_map iPho_obj = { iG, vMatchedRecoPhoIdxs, vMatchedPreselPhoIdxs };
    vPhos.push_back( iPho_obj );

  } // gen phos

  ////////// Apply selection criteria //////////

  // Ensure only 1 presel photon associated to each gen pho
  // Missing here: a check if no other reco photons (e.g. PU photons) are around the presel photon
  // NOTE: only unique reco idxs are stored
  vPreselPhoIdxs_.clear();
  for ( auto const& iPho : vPhos ) {

    if ( debug ) std::cout << " >> pho[" << iPho.idx
      << "], reco size:" << iPho.matchedRecoPhoIdxs.size()
      << ", presel size:"<< iPho.matchedPreselPhoIdxs.size()<< std::endl;

    // No presel OR >1 presel photon
    if ( iPho.matchedPreselPhoIdxs.empty() || iPho.matchedPreselPhoIdxs.size() > 1 ) continue;

    vPreselPhoIdxs_.push_back( iPho.matchedPreselPhoIdxs[0] );

  }
  if ( debug ) std::cout << " >> PreselPhos.size: " << vPreselPhoIdxs_.size() << std::endl;
  if ( vPreselPhoIdxs_.empty() ) return false;

  // Photon gets passed to cropping routine.
  if ( debug ) std::cout << " >> Passed selection. " << std::endl;
  return true;

} // runPhotonSel()

// Fill branches ___________________________________________________________________//
void SCRegressor::fillPhotonSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  ////////// Store gen-level pho kinematics //////////

  vSC_DR_.clear();
  vSC_mass_.clear();
  vSC_E_.clear();
  vSC_pT_.clear();
  vSC_eta_.clear();
  vSC_phi_.clear();
  for ( auto const& iPho : vPhos ) {

    // Skip phos which are not valid for regression
    if ( iPho.matchedPreselPhoIdxs.empty() || iPho.matchedPreselPhoIdxs.size() > 1 ) continue;
    /*std::vector<int>::iterator it;    
    //it =std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), iPho.matchedPreselPhoIdxs[0]) ;
    if(it!= vRegressPhoIdxs_.end()) std::cout <<"_______ Element " << it-vRegressPhoIdxs_.begin() << std::endl;*/

    //if matched preselcted photon does not pass zero suppression ( they are:regressed photon IDx) then continue
    if ( std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), iPho.matchedPreselPhoIdxs[0]) == vRegressPhoIdxs_.end() ) continue;

    reco::GenParticleRef iGen( genParticles, iPho.idx );

    vSC_DR_.push_back( 0. );
    vSC_mass_.push_back( 0. );
    vSC_E_.push_back( iGen->energy() );
    vSC_pT_.push_back( iGen->pt() );
    vSC_eta_.push_back( iGen->eta() );
    vSC_phi_.push_back( iGen->phi() );

    hSC_pT->Fill( iGen->pt() );

  } // gen phos

} // fillPhotonSel()
