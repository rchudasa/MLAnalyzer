#include "MLAnalyzer/RecHitAnalyzer/interface/SCRegressor.h"

// Initialize branches _____________________________________________________//
void SCRegressor::branchesNJetsSel ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("m0",        &m0_);
  //tree->Branch("FC_inputs", &vFC_inputs_);
  //tree->Branch("hltAccept", &hltAccept_);
  //tree->Branch("nRecoPho",  &nRecoPho_);
  //tree->Branch("minDR",     &vMinDR_);
  tree->Branch("jet_pt",      &vJet_pt_);
  tree->Branch("jet_eta",     &vJet_eta_);
  tree->Branch("jet_phi",     &vJet_phi_);
  tree->Branch("jet_energy",  &vJet_energy_);
  tree->Branch("jet_tightId", &vJet_tightId_);

  //hNpassed_kin      = fs->make<TH1F>("hNpassed_kin", "isPassed;isPassed;N", 2, 0., 2);
  //hNpassed_presel   = fs->make<TH1F>("hNpassed_presel", "isPassed;isPassed;N", 2, 0., 2);
  //hNpassed_mGG      = fs->make<TH1F>("hNpassed_mGG", "isPassed;isPassed;N", 2, 0., 2);
  //hNpassed_nRecoPho = fs->make<TH1F>("hNpassed_nRecoPho", "isPassed;isPassed;N", 2, 0., 2);
  //hNpassed_hlt      = fs->make<TH1F>("hNpassed_hlt", "isPassed;isPassed;N", 2, 0., 2);
}

// Run event selection ___________________________________________________________________//
bool SCRegressor::runNJetsSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<ElectronCollection> electrons;
  iEvent.getByToken(electronCollectionT_, electrons);

  edm::Handle<PhotonCollection> photons;
  iEvent.getByToken(photonCollectionT_, photons);

  edm::Handle<JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  ////////// Apply selection //////////

  //std::cout << " jet collection size:" << jets->size() << std::endl;

  bool looseJetID = false;
  bool tightJetID = false;
  std::vector<unsigned int> vRecoJetIdxs;
  for ( unsigned int iJ = 0; iJ < jets->size(); iJ++ ) {
    JetRef iJet( jets, iJ );
    std::cout << "!!!!!!!!! RUNNING JET SELECTION HERE!!!!!!!!!" << std::endl;
    if ( std::abs(iJet->pt()) < 20. ) continue;
    if ( std::abs(iJet->eta()) > 2.4 ) continue;
    //if ( std::abs(iJet->eta()) > 1.442 && std::abs(iJet->eta()) < 1.566 ) continue;
    // isLooseJet
    float NHF      = iJet->neutralHadronEnergyFraction();
    float NEMF     = iJet->neutralEmEnergyFraction();
    float NumConst = iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
    float CHF      = iJet->chargedHadronEnergyFraction();
    float CHM      = iJet->chargedMultiplicity();
    float CEMF     = iJet->chargedEmEnergyFraction();
    //float NNP      = iJet->neutralMultiplicity();
    looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99));
    if ( !looseJetID ) continue;
    //std::cout << " jet.pt:" << std::abs(iJet->pt()) << std::endl;
    //tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99));
    vRecoJetIdxs.push_back( iJ );
    //vDiJet += iJet->p4();
  }
  if ( vRecoJetIdxs.size() < 2 ) return false;
  //std::cout << " vRecoJetIdxs size:" << vRecoJetIdxs.size() << std::endl;

  bool isInHWindow = false;
  std::vector<unsigned int> vDiJetIdxs;
  for ( unsigned int iJ : vRecoJetIdxs ) {

    // Check ele passes tight ID
    JetRef iJet( jets, iJ );

    // Check is Z-mass compatible with other electrons
    isInHWindow = false;
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vDiJet;
    for ( unsigned int jJ : vRecoJetIdxs ) {

      if ( iJ == jJ ) continue;
      JetRef jJet( jets, jJ );

      vDiJet = iJet->p4() + jJet->p4();
      if ( vDiJet.mass() < 90. || vDiJet.mass() > 190. ) continue;
      //if ( vDiJet.mass() < 90. ) continue;

      m0_ = vDiJet.mass();
      //if ( std::abs(iJet->pt()) < m0_/3. ) continue;
      //if ( std::abs(jJet->pt()) < m0_/4. ) continue;
      vDiJetIdxs.push_back( iJ );
      vDiJetIdxs.push_back( jJ );
      isInHWindow = true;

      if ( isInHWindow ) break;
    }// j electrons

    // Stop at first compatible match
    if ( isInHWindow ) break;

  } // i electrons
  if ( !isInHWindow ) return false;
  if ( vDiJetIdxs.size() != 2 ) return false;
  //std::cout << " vDiJet.mass():" << m0_<< std::endl;

  // Apply pt/mGG cuts
  float ptOmCut[2] = {  3.,  4. };
  std::vector<unsigned int> vSelJetIdxs;
  //for ( unsigned int iJ : vDiJetIdxs ) {
  for ( unsigned int iJ = 0; iJ < vDiJetIdxs.size(); iJ++ ) {

    JetRef iJet( jets, vDiJetIdxs[iJ] );

    if ( std::abs(iJet->pt()) < m0_/ptOmCut[iJ] ) continue;

    vSelJetIdxs.push_back( vDiJetIdxs[iJ] );

  } // vDiJetIdxs
  if ( vSelJetIdxs.size() != 2 ) return false;
  //std::cout << " passed jet selection" << std::endl;

  /*
  //std::cout << " Pho collection size:" << photons->size() << std::endl;
  // Ensure jet passes photon presel
  std::vector<unsigned int> vRecoPhoIdxs;
  for ( unsigned int iP = 0; iP < photons->size(); iP++ ) {

    PhotonRef iPho( photons, iP );

    if ( std::abs(iPho->pt()) <= 25. ) continue;
    if ( std::abs(iPho->eta()) >= 1.442 ) continue;

    if ( iPho->full5x5_r9() <= 0.5 ) continue;
    if ( iPho->hadTowOverEm() >= 0.08 ) continue;
    if ( iPho->hasPixelSeed() == true ) continue;
    if ( iPho->full5x5_r9() <= 0.85 ) {
      if ( iPho->full5x5_sigmaIetaIeta() >= 0.015 ) continue;
      if ( iPho->userFloat("phoPhotonIsolation") >= 4.0 ) continue;
      //if ( iPho->photonIso() >= 4.0 ) continue;
      if ( iPho->trkSumPtHollowConeDR03() >= 6. ) continue;
      //if ( iPho->trackIso() >= 6. ) continue;
    }
    vRecoPhoIdxs.push_back( iP );

  } // photons
  if ( vRecoPhoIdxs.size() < 1 ) return false;
  //std::cout << " vRecoPhoIdxs.size():" << vRecoPhoIdxs.size() << std::endl;

  // Ensure presel photon doesnt overlap with eles
  float dR;
  bool isIsolated = true;
  for ( unsigned int iP : vRecoPhoIdxs ) {
    PhotonRef iPho( photons, iP );
    isIsolated = true;
    for ( unsigned int iE : vRecoEleIdxs ) {
      ElectronRef iEle( electrons, iE );
      dR = reco::deltaR( iPho->eta(),iPho->phi(), iEle->eta(),iEle->phi() );
      //std::cout << " dR to ele["<<iE<<"]:" << dR << std::endl;
      if ( dR < 0.3 ) {
        isIsolated = false;
        break;
      } // dR
    } // vRecoEleIdxs
    if ( !isIsolated ) continue;
    vPreselPhoIdxs_.push_back( iP );
  } // vRecoPhoIdxs
  //std::cout << " n passing NJets sel:" << vPreselPhoIdxs_.size() << std::endl;
  */

  vJet_pt_.clear();
  vJet_energy_.clear();
  vJet_eta_.clear();
  vJet_phi_.clear();
  vJet_tightId_.clear();
  for ( int iJ : vSelJetIdxs ) {

    JetRef iJet( jets, iJ );
    // Fill branch arrays
    vJet_pt_.push_back( iJet->pt() );
    vJet_energy_.push_back( iJet->energy() );
    vJet_eta_.push_back( iJet->eta() );
    vJet_phi_.push_back( iJet->phi() );

    float NHF      = iJet->neutralHadronEnergyFraction();
    float NEMF     = iJet->neutralEmEnergyFraction();
    float NumConst = iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
    float CHF      = iJet->chargedHadronEnergyFraction();
    float CHM      = iJet->chargedMultiplicity();
    float CEMF     = iJet->chargedEmEnergyFraction();
    //float NNP      = iJet->neutralMultiplicity();

    tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99));
    vJet_tightId_.push_back( tightJetID );
  }

  return true;
}

// Fill branches ___________________________________________________________________//
void SCRegressor::fillNJetsSel ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<JetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  /*
  // Fill kinematic variables
  float dphi[2] = { 0., 0. };
  vFC_inputs_.clear();
  for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
    PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
    vFC_inputs_.push_back( iPho->pt()/m0_ );
    vFC_inputs_.push_back( iPho->eta() );
    dphi[iP] = iPho->phi();
  }
  vFC_inputs_.push_back( TMath::Cos(reco::deltaPhi(dphi[0], dphi[1])) );

  // Get dR of closest reco photon to presel photon
  float minDR, dR;
  vMinDR_.clear();
  for ( unsigned int jP = 0; jP < vRegressPhoIdxs_.size(); jP++ ) {
    PhotonRef jPho( photons, vRegressPhoIdxs_[jP] );
    minDR = 100.;
    for ( unsigned int kP = 0; kP < photons->size(); kP++ ) {
      if ( std::find(vRegressPhoIdxs_.begin(), vRegressPhoIdxs_.end(), kP) != vRegressPhoIdxs_.end() ) continue;
      PhotonRef kPho( photons, kP );
      if ( std::abs(kPho->pt()) < 10. ) continue;
      dR = reco::deltaR( jPho->eta(),jPho->phi(), kPho->eta(),kPho->phi() );
      if ( dR < minDR ) minDR = dR;
    } //k
    vMinDR_.push_back( minDR );
  } //j
  */

}
