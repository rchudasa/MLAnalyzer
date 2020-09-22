#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

using std::vector;

const unsigned nJets = 2; //TODO: use cfg level nJets_
TH1D *h_tau_jet_pT;
TH1D *h_tau_jet_E;
TH1D *h_tau_jet_eta;
TH1D *h_tau_jet_m0;
TH1D *h_tau_jet_nJet;
TH1D *h_nNonTau;
TH1D *h_nIsTau;
vector<float> v_tau_jet_m0_;
vector<float> v_tau_jet_pt_;
vector<float> v_tau_jetIsTau_;
vector<float> v_tau_jetPdgIds_;

vector<float> v_tau_subJetE_[nJets];
vector<float> v_tau_subJetPx_[nJets];
vector<float> v_tau_subJetPy_[nJets];
vector<float> v_tau_subJetPz_[nJets];

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEvtSel_jet_dijet_tau ( TTree* tree, edm::Service<TFileService> &fs ) {

  h_tau_jet_pT    = fs->make<TH1D>("h_jet_pT"  , "p_{T};p_{T};Particles", 100,  0., 500.);
  h_tau_jet_E     = fs->make<TH1D>("h_jet_E"   , "E;E;Particles"        , 100,  0., 800.);
  h_tau_jet_eta   = fs->make<TH1D>("h_jet_eta" , "#eta;#eta;Particles"  , 100, -5., 5.);
  h_tau_jet_nJet  = fs->make<TH1D>("h_jet_nJet", "nJet;nJet;Events"     ,  10,  0., 10.);
  h_tau_jet_m0    = fs->make<TH1D>("h_jet_m0"  , "m0;m0;Events"         , 100,  0., 100.);
  h_nNonTau        = fs->make<TH1D>("h_nNonTau"  , "nNonTau;nNonTau;Events" ,   3,  0.,   3.);
  h_nIsTau           = fs->make<TH1D>("h_nIsTau"     , "nIsTau;nIsTau;Events"       ,   3,  0.,   3.);

  tree->Branch("jetM",       &v_tau_jet_m0_);
  tree->Branch("jetPt",      &v_tau_jet_pt_);
  tree->Branch("jetIsTau",   &v_tau_jetIsTau_);
  tree->Branch("jetPdgIds",  &v_tau_jetPdgIds_);

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
  float dR;

  vJetIdxs.clear();
  v_tau_jetPdgIds_.clear();
  /*
  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetCollectionT_, genJets);
  std::vector<float> v_tau_jetFakePhoIdxs;
  */

  unsigned nNonTau = 0;
  unsigned nIsTau = 0;

  if ( debug ) std::cout << " >>>>>>>>>>>>>>>>>>>> evt:" << std::endl;
  for (reco::GenParticleCollection::const_iterator iGen = genParticles->begin();
      iGen != genParticles->end();
      ++iGen) {

    //std::cout << " DEBUG >> GEN particle status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " | pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;
    //if (!( (abs(iGen->pdgId()) == 15 && iGen->numberOfMothers() == 0) || ( iGen->status() == 23 || iGen->status() == 33 || iGen->status() == 43 || iGen->status() == 44 ))) continue;
    //if (!( (abs(iGen->pdgId()) == 15 && iGen->numberOfMothers() == 0) || ( iGen->status() == 23 || iGen->status() == 33 || iGen->status() == 44 ))) continue;
    if (!( (abs(iGen->pdgId()) == 15 && iGen->numberOfMothers() == 0) || ( iGen->status() == 23 || iGen->status() == 33 ))) continue;
    //if ( iGen->status() != 3 ) continue; // pythia6: 3 = hard scatter particle
    //if ( iGen->status() != 23 ) continue; // pythia8: 23 = outgoing particle from hard scatter
    std::cout << " DEBUG >> GEN particle status: " << iGen->status() << ", id: " << iGen->pdgId() << ", nDaught: " << iGen->numberOfDaughters() << " | pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;
    if ( debug ) std::cout << " >> id:" << iGen->pdgId() << " status:" << iGen->status() << " nDaught:" << iGen->numberOfDaughters() << " pt:"<< iGen->pt() << " eta:" <<iGen->eta() << " phi:" <<iGen->phi() << " nMoms:" <<iGen->numberOfMothers()<< std::endl;

    // Loop over jets
    std::cout << "  JETS IN THE EVENT = " << jets->size() << " | Selection requires minpT = " << minJetPt_ << " and maxEta = "<< maxJetEta_ << std::endl;
    for ( unsigned iJ(0); iJ != jets->size(); ++iJ ) {
      reco::PFJetRef iJet( jets, iJ );
      //std::cout << "Jet[" << iJ << "] : pt = " << iJet->pt() << " , eta = " << iJet->eta() << std::endl;
      if ( std::abs(iJet->pt())  < minJetPt_ ) continue;
      if ( std::abs(iJet->eta()) > maxJetEta_ ) continue;
      dR = reco::deltaR( iJet->eta(),iJet->phi(), iGen->eta(),iGen->phi() );
      //std::cout << "Jet[" << iJ << "] PASSED MIN PT and MAX ETA! |  dR = " << dR << std::endl;
      std::cout << " >>>>>> jet[" << iJ << "] Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi() << std::endl;
      if ( debug ) {
        std::cout << " >>>>>> jet[" << iJ << "] Pt:" << iJet->pt() << " jetEta:" << iJet->eta() << " jetPhi:" << iJet->phi()
        << " dR:" << dR << std::endl;
      }
      if ( dR > 0.4 ) {
        std::cout << "FAIL dR SELECTION, dR = " << dR << std::endl;
        continue;
      }
      std::cout << " >>>>>> DR matched: jet[" << iJ << "] pdgId: " << std::abs(iGen->pdgId()) << " | dR: " << dR << "| Pt: " << iJet->pt() << ", Eta: " << iJet->eta() << ", Phi: " << iJet->phi() << std::endl;
      vJetIdxs.push_back( iJ );
      v_tau_jetPdgIds_.push_back( std::abs(iGen->pdgId()) );
      if ( debug ) {
        std::cout << " >>>>>> DR matched: jet[" << iJ << "] pdgId:" << std::abs(iGen->pdgId()) << std::endl;
      }
      break;
    } // reco jets

  } // gen particles

  // Check jet multiplicity
  //std::cout << "    --- Jet multiplicity: vJetIdxs.size() = " << vJetIdxs.size() << " | nJets = " << nJets << std::endl;
  //if ( vJetIdxs.size() != nJets ) return false;
  if ( vJetIdxs.size() < 1 || vJetIdxs.size() > 2) return false; // Requiring at least one jet

  // Check jet identities
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {
    std::cout << "    >>>>>>>> Jet Pdg ID = " << v_tau_jetPdgIds_[iJ] << ", Idx = " << vJetIdxs[iJ] << std::endl;
    if ( abs(v_tau_jetPdgIds_[iJ]) == 15 ) nIsTau++;
    else nNonTau++;
  }
  //if ( vJetIdxs[0] == vJetIdxs[1] ) return false; // protect against double matching: only valid for nJets==2
  //if ( nIsTau+nNonTau != nJets ) return false; // require dijet
  //if ( nIsTau != nJets && nNonTau != nJets ) return false; // require gg or qq final state
  //if ( nIsTau != 1 && nNonTau != 1 ) return false; // require qg final state

  if ( debug ) std::cout << " >> has_jet_dijet_tau: passed" << std::endl;
  return true;

} // runEvtSel_jet_dijet_tau()

// Fill branches and histograms _____________________________________________________//
void RecHitAnalyzer::fillEvtSel_jet_dijet_tau ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetCollectionT_, jets);

  h_tau_jet_nJet->Fill( vJetIdxs.size() );

  unsigned nNonTau = 0;
  unsigned nIsTau = 0;
  v_tau_jet_pt_.clear();
  v_tau_jet_m0_.clear();
  v_tau_jetIsTau_.clear();
  for ( unsigned iJ(0); iJ != vJetIdxs.size(); ++iJ ) {

    reco::PFJetRef iJet( jets, vJetIdxs[iJ] );

    // Fill histograms 
    h_tau_jet_pT->Fill( std::abs(iJet->pt()) );
    h_tau_jet_eta->Fill( iJet->eta() );
    h_tau_jet_E->Fill( iJet->energy() );
    h_tau_jet_m0->Fill( iJet->mass() );

    // Fill branches 
    v_tau_jet_pt_.push_back( iJet->pt() );
    v_tau_jet_m0_.push_back( iJet->mass() );
    if ( abs(v_tau_jetPdgIds_[iJ]) == 15 ) {
      v_tau_jetIsTau_.push_back( true );
      nIsTau++;
    } else {
      v_tau_jetIsTau_.push_back( false );
      nNonTau++;
    }

    // Get jet constituents
    v_tau_subJetE_[iJ].clear();
    v_tau_subJetPx_[iJ].clear();
    v_tau_subJetPy_[iJ].clear();
    v_tau_subJetPz_[iJ].clear();
    //std::vector<reco::PFCandidatePtr> jetConstituents = iJet->getPFConstituents();
    unsigned int nConstituents = iJet->getPFConstituents().size();
    std::cout << " >>>>>> Jet " << iJ << " has "<< nConstituents << " constituents" << std::endl;
    for ( unsigned int j = 0; j < nConstituents; j++ ) {
      const reco::PFCandidatePtr subJet = iJet->getPFConstituent( j );
      std::cout << " >>>>>>>>>>>  Jet constituent " << j << "-> E:" << subJet->energy() << " px:" << subJet->px() << " py:" << subJet->py() << " pz:" << subJet->pz() << std::endl;
      v_tau_subJetE_[iJ].push_back( subJet->energy() );
      v_tau_subJetPx_[iJ].push_back( subJet->px() );
      v_tau_subJetPy_[iJ].push_back( subJet->py() );
      v_tau_subJetPz_[iJ].push_back( subJet->pz() );
    }
  }
  h_nNonTau->Fill( nNonTau );
  h_nIsTau->Fill( nIsTau );


} // fillEvtSel_jet_dijet_tau()
