#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <algorithm>

using std::vector;
using namespace trigger;

struct gen_obj {
  unsigned int idx;
  double pt;
};
std::vector<gen_obj> vATaus;
std::vector<unsigned int> vGenATauIdxs;

std::vector<float> vA_E_;
std::vector<float> vA_pT_;
std::vector<float> vA_eta_;
std::vector<float> vA_phi_;
std::vector<float> vA_mass_;
std::vector<float> vA_DR_;
std::vector<float> vA_recoIdx_;
std::vector<float> vA_pdgId_;
std::vector<float> vA_mothPdgId_;
std::vector<float> vA_jetM_;
std::vector<float> vA_status_;
float mHgen_;
TH2F * hdPhidEta;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_ditau_h2aa4Tau ( TTree* tree, edm::Service<TFileService> &fs )
{
  tree->Branch("mHgen",     &mHgen_);
  //tree->Branch("FC_inputs", &vFC_inputs_);
  //tree->Branch("hltAccept", &hltAccept_);

  tree->Branch("A_mass",    &vA_mass_);
  tree->Branch("A_DR",      &vA_DR_);
  tree->Branch("A_E",       &vA_E_);
  tree->Branch("A_pT",      &vA_pT_);
  tree->Branch("A_eta",     &vA_eta_);
  tree->Branch("A_phi",     &vA_phi_);
  tree->Branch("A_recoIdx", &vA_recoIdx_);

  hdPhidEta = fs->make<TH2F>("dPhidEta_GG", "#Delta(#phi,#eta,m);#Delta#phi(#gamma,#gamma);#Delta#eta(#gamma,#gamma)",
              6, 0., 6.*0.0174, 6, 0., 6.*0.0174);
}

// Run event selection ___________________________________________________________________//
bool RecHitAnalyzer::runEvtSel_jet_dijet_ditau_h2aa4Tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  std::cout << "Coming in Gen loop" << std::endl;

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  vATaus.clear();
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vH;

for ( unsigned int iG = 0; iG < genParticles->size(); iG++ ) {

    reco::GenParticleRef iGen( genParticles, iG );
    if ( std::abs(iGen->pdgId()) != 25 ) continue;
    if ( iGen->numberOfDaughters() != 2 ) continue;
    if ( abs(iGen->daughter(0)->pdgId()) != 15 || abs(iGen->daughter(1)->pdgId()) != 15 ) continue;

    std::cout<<"*****************************************************"<< std::endl;
    std::cout<< "iG:" << iG << " ID:" << iGen->pdgId() << " A mass:" << iGen->mass() << std::endl;
    gen_obj Gen_obj = { iG, std::abs(iGen->pt()) };
    vATaus.push_back( Gen_obj );
    vH += iGen->p4();

  }

  std::cout<< "********** Gen Higgs mass " << vH.mass() << std::endl;
  
  if ( vATaus.size() != 2 ) return false;
  mHgen_ = vH.mass();
  std::cout << "*************************************** Higgs mass:" << mHgen_ << "**************************" << std::endl;
  // Sort As by pT, for abitrary N
  std::sort( vATaus.begin(), vATaus.end(), [](auto const &a, auto const &b) { return a.pt > b.pt; } );

  vGenATauIdxs.clear();
  for ( unsigned int iG = 0; iG < vATaus.size(); iG++ ) {
    reco::GenParticleRef iGen( genParticles, vATaus[iG].idx );
    if ( debug ) std::cout << " >> pT:" << iGen->pt() << " eta:" << iGen->eta() << " phi: " << iGen->phi() << " E:" << iGen->energy() << std::endl;
    vGenATauIdxs.push_back( vATaus[iG].idx );
  }


  return true;
}

// Fill branches ___________________________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_ditau_h2aa4Tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticleCollectionT_, genParticles);

  //edm::Handle<PhotonCollection> photons;
  //iEvent.getByToken(photonCollectionT_, photons);

  vA_E_.clear();
  vA_pT_.clear();
  vA_eta_.clear();
  vA_phi_.clear();
  vA_mass_.clear();
  vA_DR_.clear();
  vA_recoIdx_.clear();
  float dPhi, dEta;
  //float dPhi, dEta, dR, recoDR;
  int recoDR_idx;
  for ( unsigned int iG : vGenATauIdxs ) {

    reco::GenParticleRef iGen( genParticles, iG );

    vA_E_.push_back( std::abs(iGen->energy()) );
    vA_pT_.push_back( std::abs(iGen->pt()) );
    vA_eta_.push_back( iGen->eta() );
    vA_phi_.push_back( iGen->phi() );
    vA_mass_.push_back( iGen->mass() );
    vA_DR_.push_back( reco::deltaR(iGen->daughter(0)->eta(),iGen->daughter(0)->phi(), iGen->daughter(1)->eta(),iGen->daughter(1)->phi()) );

    dPhi = reco::deltaPhi( iGen->daughter(0)->phi(), iGen->daughter(1)->phi() );
    dEta = std::abs( iGen->daughter(0)->eta() - iGen->daughter(1)->eta() );
    hdPhidEta->Fill( dPhi, dEta );

    // Get index to dR-matched preselected photon
    //recoDR = 2*0.04;
    recoDR_idx = -1;
    // Want vA_recoIdx_ to store vector index in vRegressPhoIdxs_
    // i.e., vRegressPhoIdxs_[0]:leading reco pho, vRegressPhoIdxs_[1]:sub-leading reco pho
    // not position in original photon collection
    /*for ( unsigned int iP = 0; iP < vRegressPhoIdxs_.size(); iP++ ) {
      PhotonRef iPho( photons, vRegressPhoIdxs_[iP] );
      dR = reco::deltaR(iGen->eta(),iGen->phi(), iPho->eta(),iPho->phi());
      if ( dR < recoDR ) {
        recoDR = dR;
        recoDR_idx = iP;
      }
    }*/ // reco pho
    vA_recoIdx_.push_back( recoDR_idx );

  } // gen A

}
