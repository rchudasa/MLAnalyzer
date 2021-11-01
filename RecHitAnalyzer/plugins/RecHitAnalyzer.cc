// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
//
//

#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"

//
// constructors and destructor
//
RecHitAnalyzer::RecHitAnalyzer(const edm::ParameterSet& iConfig)
{
  // ECAL AND HCAL 
  //EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHitCollection"));
  EBRecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEBRecHitCollection"));
  //EBDigiCollectionT_      = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("selectedEBDigiCollection"));
  //EBDigiCollectionT_      = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBDigiCollection"));
  EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedEERecHitCollection"));
  //EERecHitCollectionT_    = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHitCollection"));
  HBHERecHitCollectionT_  = consumes<HBHERecHitCollection>(iConfig.getParameter<edm::InputTag>("reducedHBHERecHitCollection"));

  // VERTEXES AND TRACKS
  vertexCollectionT_      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
  trackCollectionT_       = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  TRKRecHitCollectionT_   = consumes<TrackingRecHitCollection>(iConfig.getParameter<edm::InputTag>("trackCollection"));
  transientTrackBuilderT_ = iConfig.getParameter<edm::ESInputTag>("transTrackBuilder");

  // PIXELS AND STRIPS
  siPixelRecHitCollectionT_        = consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("siPixelRecHitCollection"));
  siStripMatchedRecHitCollectionT_ = consumes<SiStripMatchedRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripMatchedRecHitCollection"));
  siStripRPhiRecHitCollectionT_    = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripRphiRecHits"));
  siStripStereoRecHitCollectionT_  = consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("siStripStereoRecHits"));


  //MET, JETS and PARTICLES
  //metCollectionT_         = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("metCollection"));
  //metPuppiCollectionT_    = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("slimmedMETsPuppi"));
  jetCollectionT_         = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("JetCollection"));

  genParticleCollectionT_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
  tauCollectionT_         = consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauCollection"));
  //eleCollectionT_         = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("eleCollection"));

  //OTHERS
  processName_            = iConfig.getUntrackedParameter<std::string>("processName","HLT");
  triggerResultsToken_    = consumes<edm::TriggerResults> (iConfig.getParameter<edm::InputTag>("triggerResultsTag"));
  rhoLabel_               = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"));
  //metSigAlgo_            = new metsig::METSignificance(iConfig);

  //jet/event configuration
  debug      = iConfig.getParameter<bool>("isDebug");
  isData_    = iConfig.getParameter<bool>("isData");
  isSignal_  = iConfig.getParameter<bool>("isSignal");
  isTrain_   = iConfig.getParameter<bool>("isTrain");
  isW_       = iConfig.getParameter<bool>("isW");
  mode_      = iConfig.getParameter<std::string>("mode");
  task_      = iConfig.getParameter<std::string>("task");
  minJetPt_  = iConfig.getParameter<double>("minJetPt");
  maxJetEta_ = iConfig.getParameter<double>("maxJetEta");
  z0PVCut_   = iConfig.getParameter<double>("z0PVCut");

  std::cout << " >> Mode set to " << mode_ << std::endl;
  if ( mode_ == "JetLevel" ) {
    doJets_ = true;
    nJets_ = iConfig.getParameter<int>("nJets");
    std::cout << "\t>> nJets set to " << nJets_ << std::endl;
  } else if ( mode_ == "EventLevel" ) {
    doJets_ = false;
  } else {
    std::cout << " >> Assuming EventLevel Config. " << std::endl;
    doJets_ = false;
  }



  // Initialize file writer
  // NOTE: initializing dynamic-memory histograms outside of TFileService
  // will cause memory leaks
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  h_sel = fs->make<TH1F>("h_sel", "isSelected;isSelected;Events", 2, 0., 2.);

  ///////////adjustable granularity stuff

  granularityMultiPhi[0]  = iConfig.getParameter<int>("granularityMultiPhi");
  granularityMultiEta[0]  = iConfig.getParameter<int>("granularityMultiEta");

  granularityMultiPhi[1] = 3;
  granularityMultiEta[1] = 3;

  for (unsigned int proj=0; proj<Nadjproj; proj++)
  {

    int totalMultiEta = granularityMultiEta[proj] * granularityMultiECAL;

    for (int i=0; i<eta_nbins_HBHE; i++)
    {
      double step=(eta_bins_HBHE[i+1]-eta_bins_HBHE[i])/totalMultiEta;
      for (int j=0; j<totalMultiEta; j++)
      {
        adjEtaBins[proj].push_back(eta_bins_HBHE[i]+step*j);
      }
    }
    adjEtaBins[proj].push_back(eta_bins_HBHE[eta_nbins_HBHE]);

    totalEtaBins[proj] = totalMultiEta*(eta_nbins_HBHE);
    totalPhiBins[proj] = granularityMultiPhi[proj] * granularityMultiECAL*HBHE_IPHI_NUM;

  }

  //////////// TTree //////////

  // These will be use to create the actual images
  RHTree = fs->make<TTree>("RHTree", "RecHit tree");
  if ( doJets_ ) {
    branchesEvtSel_jet( RHTree, fs );
  } 
  branchesEB           ( RHTree, fs );
  branchesEE           ( RHTree, fs );
  branchesHBHE         ( RHTree, fs );
  branchesECALatHCAL   ( RHTree, fs );
  branchesECALstitched ( RHTree, fs );
  branchesHCALatEBEE   ( RHTree, fs );
  branchesTracksAtEBEE(RHTree, fs);
  branchesTracksAtECALstitched( RHTree, fs);
  branchesTRKlayersAtECALstitched(RHTree, fs);

} // constructor
//
RecHitAnalyzer::~RecHitAnalyzer()
{
  //delete metSigAlgo_;  //FIXME
}
//
// member functions
// ------------ method called for each event  ------------
void
RecHitAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  nTotal++;
  using namespace edm;

  // ----- Apply event selection cuts ----- //

  bool passedSelection = false;
  if ( doJets_ ) {
    passedSelection = runEvtSel_jet( iEvent, iSetup );
  } 

  if ( !passedSelection ) {
    if ( debug ) std::cout << "!!!!!!!!!!! DID NOT PASS EVENT/JET SELECTION !!!!!!!!!!!" << std::endl;
    h_sel->Fill( 0. );;
    return;
  }

  fillEB( iEvent, iSetup );
  fillEE( iEvent, iSetup );
  fillHBHE( iEvent, iSetup );
  fillECALatHCAL( iEvent, iSetup );
  fillECALstitched( iEvent, iSetup );
  fillHCALatEBEE( iEvent, iSetup );
  fillTracksAtEBEE( iEvent, iSetup );
  for (unsigned int i=0;i<Nproj;i++)
  {
    fillTracksAtECALstitched( iEvent, iSetup, i );
  }
  for (unsigned int i=0;i<Nhitproj;i++)
  {
    fillTRKlayersAtECALstitched( iEvent, iSetup, i );
  }

  ////////////// 4-Momenta //////////
  //fillFC( iEvent, iSetup );

  // Fill RHTree
  RHTree->Fill();
  h_sel->Fill( 1. );
  nPassed++;

} // analyze()


// ------------ method called once each job just before starting event loop  ------------
void 
RecHitAnalyzer::beginJob()
{
  nTotal = 0;
  nPassed = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
RecHitAnalyzer::endJob() 
{
  std::cout << " selected: " << nPassed << "/" << nTotal << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecHitAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitAnalyzer);
