#ifndef RecHitAnalyzer_h
#define RecHitAnalyzer_h
// -*- C++ -*-
//
// Package:    MLAnalyzer/RecHitAnalyzer
// Class:      RecHitAnalyzer
//
//

// system include files
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DQM/SiPixelMonitorRecHit/interface/SiPixelRecHitModule.h"

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
#include "DQM/HcalCommon/interface/Constants.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace classic_svFit;

//
// class declaration
//

static const unsigned int Nproj = 1;
static const unsigned int Nhitproj = 1;
static const unsigned int Nadjproj = 2;

class RecHitAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:
    explicit RecHitAnalyzer(const edm::ParameterSet&);
    ~RecHitAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    // Tokens
    edm::EDGetTokenT<EcalRecHitCollection> EBRecHitCollectionT_;
    edm::EDGetTokenT<EBDigiCollection>     EBDigiCollectionT_;
    edm::EDGetTokenT<EcalRecHitCollection> EERecHitCollectionT_;
    edm::EDGetTokenT<HBHERecHitCollection> HBHERecHitCollectionT_;

    edm::EDGetTokenT<TrackingRecHitCollection> TRKRecHitCollectionT_;
    edm::ESInputTag transientTrackBuilderT_;
    edm::EDGetTokenT<reco::TrackCollection> trackCollectionT_;
    edm::EDGetTokenT<reco::VertexCollection> vertexCollectionT_;

    edm::EDGetTokenT<SiPixelRecHitCollection> siPixelRecHitCollectionT_;
    edm::EDGetTokenT<SiStripMatchedRecHit2DCollection> siStripMatchedRecHitCollectionT_;
    edm::EDGetTokenT<SiStripRecHit2DCollection> siStripRPhiRecHitCollectionT_;
    edm::EDGetTokenT<SiStripRecHit2DCollection> siStripStereoRecHitCollectionT_;
 
    edm::EDGetTokenT<pat::JetCollection> jetCollectionT_;
    edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollectionT_;
    edm::EDGetTokenT<pat::METCollection> metCollectionT_;
    edm::EDGetTokenT<pat::METCollection> metPuppiCollectionT_;
    edm::EDGetTokenT<pat::ElectronCollection> eleCollectionT_;
    edm::EDGetTokenT<pat::TauCollection> tauCollectionT_;

    std::string   processName_;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    HLTConfigProvider hltConfig_;

    edm::EDGetTokenT<double> rhoLabel_;
 
    //metsig::METSignificance* metSigAlgo_;

    TH1F *h_sel;
    TTree* RHTree;

    // Selection and filling functions
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
    void branchesTRKvolumeAtEBEE( TTree*, edm::Service<TFileService>& );
    void branchesJetInfoAtECALstitched   ( TTree*, edm::Service<TFileService>& );
    void branchesTRKlayersAtECALstitched( TTree*, edm::Service<TFileService>& );

    bool runEvtSel_jet      ( const edm::Event&, const edm::EventSetup& );
    void fillEB             ( const edm::Event&, const edm::EventSetup& );
    void fillEE             ( const edm::Event&, const edm::EventSetup& );
    void fillHBHE           ( const edm::Event&, const edm::EventSetup& );
    void fillECALatHCAL     ( const edm::Event&, const edm::EventSetup& );
    void fillECALstitched   ( const edm::Event&, const edm::EventSetup& );
    void fillHCALatEBEE     ( const edm::Event&, const edm::EventSetup& );
    void fillTracksAtEBEE   ( const edm::Event&, const edm::EventSetup& );
    void fillTracksAtECALstitched   ( const edm::Event&, const edm::EventSetup&, unsigned int proj );
    void fillPFCandsAtEBEE   ( const edm::Event&, const edm::EventSetup& );
    void fillPFCandsAtECALstitched   ( const edm::Event&, const edm::EventSetup& );
    void fillTRKlayersAtEBEE( const edm::Event&, const edm::EventSetup& );
    void fillTRKvolumeAtEBEE( const edm::Event&, const edm::EventSetup& );
    void fillJetInfoAtECALstitched   ( const edm::Event&, const edm::EventSetup& );
    void fillTRKlayersAtECALstitched( const edm::Event&, const edm::EventSetup&, unsigned int proj );

    unsigned int getLayer(const DetId& detid, const TrackerTopology* tTopo);

    // Jet level functions
    std::string mode_;  // EventLevel / JetLevel
    std::string task_;  
    bool debug;  
    bool isData_;  
    bool isSignal_;  
    bool isTrain_;  
    bool isW_;  
    bool doJets_;
    int  nJets_;
    double minJetPt_;
    double maxJetEta_;
    double z0PVCut_;
    std::vector<int> vJetIdxs;
    void branchesEvtSel_jet_dijet_tau( TTree*, edm::Service<TFileService>& );
    bool runEvtSel_jet_dijet_tau( const edm::Event&, const edm::EventSetup& );
    void fillEvtSel_jet_dijet_tau( const edm::Event&, const edm::EventSetup& );

    int nTotal, nPassed;

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
