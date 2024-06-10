#ifndef RecHitAnalyzer_h
#define RecHitAnalyzer_h
// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
//
//
// Original Author:  Michael Andrews
//         Created:  Sat, 14 Jan 2017 17:45:54 GMT
//
//

// system include files
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
//#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
//#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
//#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
//#include "DataFormats/SiStripDetId/interface/TECDetId.h"
//#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
//#include "DataFormats/SiStripDetId/interface/TIDDetId.h"

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "DQM/SiPixelMonitorRecHit/interface/SiPixelRecHitModule.h"
//#include "DQM/SiPixelMonitorRecHit/interface/SiPixelClusterModule.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
//#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "DataFormats/HcalDetId/interface/HcalElectronicsId.h"
#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" // reco::PhotonCollection defined here
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/CommonMETData.h"
#include "RecoMET/METProducers/interface/METSignificanceProducer.h"
#include "RecoMET/METAlgorithms/interface/METSignificance.h"

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TSiStripRecHit2DLocalPos.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"
#include "TauAnalysis/ClassicSVfit/interface/FastMTT.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h" 

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Calibration/IsolatedParticles/interface/CaloPropagateTrack.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "FWCore/Utilities/interface/RegexMatch.h"

using namespace classic_svFit;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//static const unsigned int Nproj = 5;
static const unsigned int Nproj = 1;
static const unsigned int Nhitproj = 1;
//static const unsigned int Nhitproj = 2;
static const unsigned int Nadjproj = 2;

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
//class RecHitAnalyzer : public edm::EDAnalyzer  {
  public:
    explicit RecHitAnalyzer(const edm::ParameterSet&);
    ~RecHitAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    //switches 
    std::string mode_;  // EventLevel / JetLevel
    std::string task_;
    bool isMC_;
    bool isSignal_;
    bool isW_;
    bool isBoostedTop_;
    bool doJets_;


    // ----------member data ---------------------------
    // Tokens
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;
    edm::EDGetTokenT<TrackingRecHitCollection> TRKRecHitCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    edm::EDGetTokenT<reco::PhotonCollection> photonCollectionT_;
    edm::EDGetTokenT<reco::PFJetCollection> jetCollectionT_;
    edm::EDGetTokenT<edm::View<reco::Jet> > pfjetsToken_;
    edm::EDGetTokenT<reco::GenJetCollection> genJetCollectionT_;
    edm::EDGetTokenT<reco::TrackCollection> trackCollectionT_;
    edm::EDGetTokenT<reco::VertexCollection> vertexCollectionT_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secVertexCollectionT_;
    edm::ESInputTag transientTrackBuilderT_;
    edm::EDGetTokenT<edm::View<reco::Jet> > recoJetsT_;
    edm::EDGetTokenT<reco::JetTagCollection> jetTagCollectionT_;
    edm::EDGetTokenT<std::vector<reco::CandIPTagInfo> >    ipTagInfoCollectionT_;
    edm::EDGetTokenT<reco::PFMETCollection> metCollectionT_;
    edm::EDGetTokenT<reco::GsfElectronCollection> eleCollectionT_;

    edm::EDGetTokenT<reco::PFTauCollection> tauCollectionT_;
    edm::EDGetTokenT<reco::PFTauDiscriminator> tauDiscriminatorT_;
    edm::EDGetTokenT<reco::PFTauDiscriminator> tauDecayMode_;
    edm::EDGetTokenT<reco::PFTauDiscriminator> tauMVAIsolation_;
    edm::EDGetTokenT<reco::PFTauDiscriminator> tauMuonRejection_;
    edm::EDGetTokenT<reco::PFTauDiscriminator> tauElectronRejectionMVA6_;

    std::string   processName_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    HLTConfigProvider hltConfig_;

    edm::EDGetTokenT<double> rhoLabel_;
    std::string jetSFType_;      //to set
    std::string jetResPtType_;   //to set
    std::string jetResPhiType_;  //to set
 
    std::vector<edm::EDGetTokenT<edm::View<reco::Candidate> > > lepTokens_;
 
    typedef std::vector<reco::PFCandidate>  PFCollection;
    edm::EDGetTokenT<PFCollection> pfCollectionT_;

    edm::EDGetTokenT<edm::View<reco::Candidate> > pfCandidatesToken_;
    
    metsig::METSignificance* metSigAlgo_;

    edm::EDGetTokenT<SiPixelRecHitCollection> siPixelRecHitCollectionT_;
    edm::EDGetTokenT<SiStripMatchedRecHit2DCollection> siStripMatchedRecHitCollectionT_;
    edm::EDGetTokenT<SiStripRecHit2DCollection> siStripRPhiRecHitCollectionT_;
    edm::EDGetTokenT<SiStripRecHit2DCollection> siStripUnmatchedRPhiRecHitCollectionT_;
    edm::EDGetTokenT<SiStripRecHit2DCollection> siStripStereoRecHitCollectionT_;
    edm::EDGetTokenT<SiStripRecHit2DCollection> siStripUnmatchedStereoRecHitCollectionT_;

    //edm::InputTag trackTags_; //used to select what tracks to read from configuration file

    // Diagnostic histograms
    //TH2D * hEB_adc[EcalDataFrame::MAXSAMPLES]; 
    //TH1D * hHBHE_depth; 
    TH1F *h_sel;

    // Main TTree
    TTree* RHTree;

    // Objects used to fill RHTree branches
    //std::vector<float> vEB_adc_[EcalDataFrame::MAXSAMPLES];
    //std::vector<float> vFC_inputs_;
    //math::PtEtaPhiELorentzVectorD vPho_[2];
  
    // Selection and filling functions
    void branchesEvtSel         ( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet     ( TTree*, edm::Service<TFileService>& );
    void branchesEB             ( TTree*, edm::Service<TFileService>& );
    void branchesEE             ( TTree*, edm::Service<TFileService>& );
    void branchesHBHE           ( TTree*, edm::Service<TFileService>& );
    void branchesECALatHCAL     ( TTree*, edm::Service<TFileService>& );
    void branchesECALstitched   ( TTree*, edm::Service<TFileService>& );
    void branchesHCALatEBEE     ( TTree*, edm::Service<TFileService>& );
    void branchesTracksAtEBEE   ( TTree*, edm::Service<TFileService>& );
    void branchesTracksAtECALstitched   ( TTree*, edm::Service<TFileService>& );
    void branchesPFCandsAtEBEE   ( TTree*, edm::Service<TFileService>& );
    void branchesPFCandsAtECALstitched   ( TTree*, edm::Service<TFileService>& );
    void branchesTRKlayersAtEBEE( TTree*, edm::Service<TFileService>& );
    //void branchesTRKlayersAtECAL( TTree*, edm::Service<TFileService>& );
    void branchesTRKvolumeAtEBEE( TTree*, edm::Service<TFileService>& );
    //void branchesTRKvolumeAtECAL( TTree*, edm::Service<TFileService>& );
    void branchesJetInfoAtECALstitched   ( TTree*, edm::Service<TFileService>& );
    void branchesTRKlayersAtECALstitched( TTree*, edm::Service<TFileService>& );
    void branchesScalarInfo( TTree*, edm::Service<TFileService>& );

    bool runEvtSel          ( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet      ( const edm::Event&, const edm::EventSetup& );
    void fillEB             ( const edm::Event&, const edm::EventSetup& );
    void fillEE             ( const edm::Event&, const edm::EventSetup& );
    void fillHBHE           ( const edm::Event&, const edm::EventSetup& );
    void fillECALatHCAL     ( const edm::Event&, const edm::EventSetup& );
    void fillECALstitched   ( const edm::Event&, const edm::EventSetup& );
    void fillHCALatEBEE     ( const edm::Event&, const edm::EventSetup& );
    void fillTracksAtEBEE   ( const edm::Event&, const edm::EventSetup& );
    void fillTracksAtECALstitched   ( const edm::Event&, const edm::EventSetup&, unsigned int proj );
    //void fillTracksAtECALstitched   ( const edm::Event&, const edm::EventSetup& );
    void fillPFCandsAtEBEE   ( const edm::Event&, const edm::EventSetup& );
    void fillPFCandsAtECALstitched   ( const edm::Event&, const edm::EventSetup& );
    void fillTRKlayersAtEBEE( const edm::Event&, const edm::EventSetup& );
    //void fillTRKlayersAtECAL( const edm::Event&, const edm::EventSetup& );
    void fillTRKvolumeAtEBEE( const edm::Event&, const edm::EventSetup& );
    //void fillTRKvolumeAtECAL( const edm::Event&, const edm::EventSetup& );
    void fillJetInfoAtECALstitched   ( const edm::Event&, const edm::EventSetup& );
    //void fillTRKlayersAtECALstitched( TTree*, edm::Service<TFileService>& );
    void fillTRKlayersAtECALstitched( const edm::Event&, const edm::EventSetup&, unsigned int proj );
    void fillScalarInfo( const edm::Event&, const edm::EventSetup& );

    //bool debug;
    const reco::PFCandidate* getPFCand(edm::Handle<PFCollection> pfCands, float eta, float phi, float& minDr, bool debug = false);
    const reco::Track* getTrackCand(edm::Handle<reco::TrackCollection> trackCands, float eta, float phi, float& minDr, bool debug = false);
    int   getTruthLabel(const reco::PFJetRef& recJet, edm::Handle<reco::GenParticleCollection> genParticles, float dRMatch = 0.4, bool debug = false);
    float getBTaggingValue(const reco::PFJetRef& recJet, edm::Handle<edm::View<reco::Jet> >& recoJetCollection, edm::Handle<reco::JetTagCollection>& btagCollection, float dRMatch = 0.1, bool debug= false );

    unsigned int getLayer(const DetId& detid, const TrackerTopology* tTopo);

    // Jet level functions
    int  nJets_;
    double minJetPt_;
    double maxJetEta_;
    double z0PVCut_;
    std::vector<int> vJetIdxs;
    std::vector<int> passedJetIdxs;
    void branchesEvtSel_jet_dijet      ( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_tau( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_top( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_ditau( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_ditau_h2aa4Tau( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_tau_massregression( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_ele_massregression( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_ele_classification( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_background( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_dijet_gg_qq( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_qcd( TTree*, edm::Service<TFileService>& );
    void branchesEvtSel_jet_photonSel( TTree*, edm::Service<TFileService>& );
    bool runEvtSel_jet_dijet      ( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_tau( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_top( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_ditau( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_ditau_h2aa4Tau( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_tau_massregression( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_ele_massregression( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_ele_classification( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_background( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_dijet_gg_qq( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_qcd( const edm::Event&, const edm::EventSetup& );
    bool runEvtSel_jet_photonSel( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet      ( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_tau( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_top( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_ditau( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_ditau_h2aa4Tau( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_tau_massregression( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_ele_massregression( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_ele_classification( const edm::Event&, const edm::EventSetup& );
    //void fillEvtSel_jet_ele_classification( const edm::Event&, const edm::EventSetup&, std::vector<int>, std::vector<int> );
    void fillEvtSel_jet_background( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_gg_qq( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_qcd( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_photonSel( const edm::Event&, const edm::EventSetup& );

    //functions for secondary vertices
    inline float catchInfs(const float& in, float replace_value=0){
      if(std::isinf(in) || std::isnan(in))
	return replace_value;
      else if(in < -1e32 || in > 1e32)
	return replace_value;
      return in;
    }    
    static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv);
    static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv);
    static float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv);
    
    int nTotal, nPassed;
    //bool debug;

    std::map<uint32_t,SiPixelRecHitModule*> thePixelStructure;

    std::vector<int> findSubcrystal(const CaloGeometry* caloGeom, const float& eta, const float& phi, const int& granularityMultiEta, const int& granularityMultiPhi);
    void fillByBinNumber(TH2F * histo, const std::vector<int>& phi_eta, const float& value);
    void fillTRKlayerHelper (int layer_, unsigned int proj, TH2F *hSUBDET_ECAL[][Nadjproj], TH2F *hEvt_Adj_SUBDET[][Nadjproj], const CaloGeometry* caloGeom, const float& eta, const float& phi);
    unsigned int getLayer(const DetId& detid);
    
    unsigned int granularityMultiPhi[Nadjproj];
    unsigned int granularityMultiEta[Nadjproj];
    
    int totalEtaBins[Nadjproj];// = totalMultiEta*(eta_nbins_HBHE);
    int totalPhiBins[Nadjproj];// = granularityMultiPhi * granularityMultiECAL*HBHE_IPHI_NUM;
    std::vector<double> adjEtaBins[Nadjproj];
    //std::vector<double> adjPhiBins[Nadjproj];
    
}; // class RecHitAnalyzer

//
// constants, enums and typedefs
//
static const bool debug = true;
//static const bool debug = false;

static const int nEE = 2;
static const int nTOB = 6;
static const int nTEC = 9;
static const int nTIB = 4;
static const int nTID = 3;
static const int nBPIX = 4;
static const int nFPIX = 6;

static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
static const int EE_MIN_IX = EEDetId::IX_MIN;//1;
static const int EE_MIN_IY = EEDetId::IY_MIN;//1;
static const int EE_MAX_IX = EEDetId::IX_MAX;//100;
static const int EE_MAX_IY = EEDetId::IY_MAX;//100;
static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100
static const int HBHE_IETA_MAX_FINE = 20;
static const int HBHE_IETA_MAX_HB = hcaldqm::constants::IETA_MAX_HB;//16;
static const int HBHE_IETA_MIN_HB = hcaldqm::constants::IETA_MIN_HB;//1
static const int HBHE_IETA_MAX_HE = hcaldqm::constants::IETA_MAX_HE;//29;
static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17
static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;
static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;
static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;
static const int ECAL_IETA_MAX_EXT = 140;

static const float zs = 0.;

// EE-(phi,eta) projection eta edges
// These are generated by requiring 5 fictional crystals
// to uniformly span each HCAL tower in eta (as in EB).
static const double eta_bins_EEm[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                  {-3.    , -2.93  , -2.86  , -2.79  , -2.72  , -2.65  , -2.62  ,
                   -2.59  , -2.56  , -2.53  , -2.5   , -2.4644, -2.4288, -2.3932,
                   -2.3576, -2.322 , -2.292 , -2.262 , -2.232 , -2.202 , -2.172 ,
                   -2.1462, -2.1204, -2.0946, -2.0688, -2.043 , -2.0204, -1.9978,
                   -1.9752, -1.9526, -1.93  , -1.91  , -1.89  , -1.87  , -1.85  ,
                   -1.83  , -1.812 , -1.794 , -1.776 , -1.758 , -1.74  , -1.7226,
                   -1.7052, -1.6878, -1.6704, -1.653 , -1.6356, -1.6182, -1.6008,
                   -1.5834, -1.566 , -1.5486, -1.5312, -1.5138, -1.4964, -1.479 }; // 56
// EE+(phi,eta) projection eta edges
static const double eta_bins_EEp[5*(hcaldqm::constants::IETA_MAX_HE-1-HBHE_IETA_MAX_EB)+1] =
                   {1.479 ,  1.4964,  1.5138,  1.5312,  1.5486,  1.566 ,  1.5834,
                    1.6008,  1.6182,  1.6356,  1.653 ,  1.6704,  1.6878,  1.7052,
                    1.7226,  1.74  ,  1.758 ,  1.776 ,  1.794 ,  1.812 ,  1.83  ,
                    1.85  ,  1.87  ,  1.89  ,  1.91  ,  1.93  ,  1.9526,  1.9752,
                    1.9978,  2.0204,  2.043 ,  2.0688,  2.0946,  2.1204,  2.1462,
                    2.172 ,  2.202 ,  2.232 ,  2.262 ,  2.292 ,  2.322 ,  2.3576,
                    2.3932,  2.4288,  2.4644,  2.5   ,  2.53  ,  2.56  ,  2.59  ,
                    2.62  ,  2.65  ,  2.72  ,  2.79  ,  2.86  ,  2.93  ,  3.    }; // 56

// HBHE eta bin edges
static const double eta_bins_HBHE[2*(hcaldqm::constants::IETA_MAX_HE-1)+1] =
                  {-3.000, -2.650, -2.500, -2.322, -2.172, -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305,
                   -1.218, -1.131, -1.044, -0.957, -0.870, -0.783, -0.695, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087, 0.000,
                    0.087,  0.174,  0.261,  0.348,  0.435,  0.522,  0.609,  0.695,  0.783,  0.870,  0.957,  1.044,  1.131,  1.218,
                    1.305,  1.392,  1.479,  1.566,  1.653,  1.740,  1.830,  1.930,  2.043,  2.172,  2.322,  2.500,  2.650,  3.000}; // 57

// MGG80, pt/m0 cut
static const int runTotal[3] = {14907, 22323, 20195}; //57425
// MGG80
//static const int runTotal[3] = {21918, 29913, 32805}; //84636
//static const int runTotal[3] = {6575, 8974, 9842}; //25391
//static const int runTotal[3] = {28493, 38887, 42647}; //84636+25391
//static const int runTotal[3] = {45363, 67946, 61752}; //175061
// MGG90
//static const int runTotal[3] = {16308, 24538, 22206}; //63052
//static const int runTotal[3] = {4892, 7361, 6662}; //18915
//static const int runTotal[3] = {21200, 31899, 28868}; //63052+18915
//static const int runTotal[3] = {35141, 47885, 52576}; //135602


//static const std::string projections[Nproj] = {"_atECALfixIP", "", "_atECAL", "_atHCAL","_atECALfixIPfromPV", "_atECALtransientTrack"}; //57425
static const std::string projections[Nproj] = {"_atECALfixIP"}; //57425
//static const std::string hit_projections[Nhitproj] = {"_atPV", ""};
static const std::string hit_projections[Nhitproj] = {"_atPV"};
static const std::string adj_projections[Nadjproj] = {"_5x5", "_3x3"};
static const int eta_nbins_HBHE = 2*(HBHE_IETA_MAX_HE-1);
static const int granularityMultiECAL=5;

//
// static data member definitions
//

#endif
