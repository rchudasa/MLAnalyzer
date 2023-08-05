#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
//#include <iostream>
//#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"


// Fill TRK rec hits ////////////////////////////////
// by layer at ECAL stitched

TH2F *hTOB_ECAL[nTOB][Nhitproj];
std::vector<float> vTOB_ECAL_[nTOB][Nhitproj];
TH2F *hEvt_EE_TOB[nTOB][nEE];

TH2F *hTEC_ECAL[nTEC][Nhitproj];
std::vector<float> vTEC_ECAL_[nTEC][Nhitproj];
TH2F *hEvt_EE_TEC[nTEC][nEE];

TH2F *hTIB_ECAL[nTIB][Nhitproj];
std::vector<float> vTIB_ECAL_[nTIB][Nhitproj];
TH2F *hEvt_EE_TIB[nTIB][nEE];

TH2F *hTID_ECAL[nTID][Nhitproj];
std::vector<float> vTID_ECAL_[nTID][Nhitproj];
TH2F *hEvt_EE_TID[nTID][nEE];

TH2F *hBPIX_ECAL[nBPIX][Nhitproj];
std::vector<float> vBPIX_ECAL_[nBPIX][Nhitproj];
std::vector<float> vBPIX_X_[nBPIX];
std::vector<float> vBPIX_Y_[nBPIX];
std::vector<float> vBPIX_Z_[nBPIX];
std::vector<float> vBPIX_Eta_[nBPIX];
std::vector<float> vBPIX_Phi_[nBPIX];
TH2F *hEvt_EE_BPIX[nBPIX][nEE];

TH2F *hFPIX_ECAL[nFPIX][Nhitproj];
std::vector<float> vFPIX_ECAL_[nFPIX][Nhitproj];
std::vector<float> vFPIX_X_[nFPIX];
std::vector<float> vFPIX_Y_[nFPIX];
std::vector<float> vFPIX_Z_[nFPIX];
std::vector<float> vFPIX_Eta_[nFPIX];
std::vector<float> vFPIX_Phi_[nFPIX];
TH2F *hEvt_EE_FPIX[nFPIX][nEE];

// Initialize branches ____________________________________________________________//
void RecHitAnalyzer::branchesTRKlayersAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images

  int layer;
  char hname[50], htitle[50];
  char bpix_X[50], bpix_Y[50], bpix_Z[50], bpix_eta[50], bpix_phi[50];
  char fpix_X[50], fpix_Y[50], fpix_Z[50], fpix_eta[50], fpix_phi[50];
  const double * eta_bins_EE[2] = {eta_bins_EEm,eta_bins_EEp};

  for ( unsigned int proj=0; proj < Nhitproj; proj++ ) {

    //TOB 
    for ( int iL(0); iL < nTOB; iL++ ) {
    // Branches for images
    layer = iL + 1;
    sprintf(hname, "TOB_layer%d_ECAL%s",layer,hit_projections[proj].c_str());
    tree->Branch(hname,        &vTOB_ECAL_[iL][proj]);
    // Histograms for monitoring
    sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
    hTOB_ECAL[iL][proj] = fs->make<TH2F>(hname, htitle,
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
    if (proj==0)
      for ( int iz(0); iz < nEE; iz++ ) {
      	const char *zside = (iz > 0) ? "p" : "m";
      	sprintf(hname, "evt_TOB_layer%d_EE%s",layer,zside);
      	sprintf(htitle,"N(ix,iy);ix;iy");
      	hEvt_EE_TOB[iL][iz] = new TH2F(hname, htitle,
      	EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      	5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
      } // iz
    } // iL
  

    //TEC
    for ( int iL(0); iL < nTEC; iL++ ) {
      // Branches for images
      layer = iL + 1;
      sprintf(hname, "TEC_layer%d_ECAL%s",layer,hit_projections[proj].c_str());
      tree->Branch(hname,        &vTEC_ECAL_[iL][proj]);
      // Histograms for monitoring
      sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
      hTEC_ECAL[iL][proj] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
      if (proj==0)
        for ( int iz(0); iz < nEE; iz++ ) {
          const char *zside = (iz > 0) ? "p" : "m";
          sprintf(hname, "evt_TEC_layer%d_EE%s",layer,zside);
          sprintf(htitle,"N(ix,iy);ix;iy");
          hEvt_EE_TEC[iL][iz] = new TH2F(hname, htitle,
          EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
          5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
        } // iz
    } // iL

    //TIB
    for ( int iL(0); iL < nTIB; iL++ ) {
      // Branches for images
      layer = iL + 1;
      sprintf(hname, "TIB_layer%d_ECAL%s",layer,hit_projections[proj].c_str());
      tree->Branch(hname,        &vTIB_ECAL_[iL][proj]);
      // Histograms for monitoring
      sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
      hTIB_ECAL[iL][proj] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
      if (proj==0)
        for ( int iz(0); iz < nEE; iz++ ) {
          const char *zside = (iz > 0) ? "p" : "m";
          sprintf(hname, "evt_TIB_layer%d_EE%s",layer,zside);
          sprintf(htitle,"N(ix,iy);ix;iy");
          hEvt_EE_TIB[iL][iz] = new TH2F(hname, htitle,
          EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
          5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
        } // iz
    } // iL

    //TID
    for ( int iL(0); iL < nTID; iL++ ) {
      // Branches for images
      layer = iL + 1;
      sprintf(hname, "TID_layer%d_ECAL%s",layer,hit_projections[proj].c_str());
      tree->Branch(hname,        &vTID_ECAL_[iL][proj]);
      // Histograms for monitoring
      sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
      hTID_ECAL[iL][proj] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
      if (proj==0)
        for ( int iz(0); iz < nEE; iz++ ) {
          const char *zside = (iz > 0) ? "p" : "m";
          sprintf(hname, "evt_TID_layer%d_EE%s",layer,zside);
          sprintf(htitle,"N(ix,iy);ix;iy");
          hEvt_EE_TID[iL][iz] = new TH2F(hname, htitle,
          EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
          5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
        } // iz
    } // iL
    
    //BPIX
    for ( int iL(0); iL < nBPIX; iL++ ) {
      // Branches for images
      layer = iL + 1;
      sprintf(hname, "BPIX_layer%d_ECAL%s",layer,hit_projections[proj].c_str());
      tree->Branch(hname,        &vBPIX_ECAL_[iL][proj]);
      
      sprintf(bpix_X, "BPIX_layer_X_%d",layer);
      tree->Branch(bpix_X,        &vBPIX_X_[iL]);
      sprintf(bpix_Y, "BPIX_layer_Y_%d",layer);
      tree->Branch(bpix_Y,        &vBPIX_Y_[iL]);
      sprintf(bpix_Z, "BPIX_layer_Z_%d",layer);
      tree->Branch(bpix_Z,        &vBPIX_Z_[iL]);
      sprintf(bpix_eta, "BPIX_layer_Eta_%d",layer);
      tree->Branch(bpix_eta,        &vBPIX_Eta_[iL]);
      sprintf(bpix_phi, "BPIX_layer_Phi_%d",layer);
      tree->Branch(bpix_phi,        &vBPIX_Phi_[iL]);


      // Histograms for monitoring
      sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
      hBPIX_ECAL[iL][proj] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
      if (proj==0)
        for ( int iz(0); iz < nEE; iz++ ) {
          const char *zside = (iz > 0) ? "p" : "m";
          sprintf(hname, "evt_BPIX_layer%d_EE%s",layer,zside);
          sprintf(htitle,"N(ix,iy);ix;iy");
          hEvt_EE_BPIX[iL][iz] = new TH2F(hname, htitle,
          EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
          5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
        } // iz
    } // iL

    //FPIX
    for ( int iL(0); iL < nFPIX; iL++ ) {
      // Branches for images
      layer = iL + 1;
      sprintf(hname, "FPIX_layer%d_ECAL%s",layer,hit_projections[proj].c_str());
      tree->Branch(hname,        &vFPIX_ECAL_[iL][proj]);

      sprintf(fpix_X, "FPIX_layer_X_%d",layer);
      tree->Branch(fpix_X,        &vFPIX_X_[iL]);
      sprintf(fpix_Y, "FPIX_layer_Y_%d",layer);
      tree->Branch(fpix_Y,        &vFPIX_Y_[iL]);
      sprintf(fpix_Z, "FPIX_layer_Z_%d",layer);
      tree->Branch(fpix_Z,        &vFPIX_Z_[iL]);
      sprintf(fpix_eta, "FPIX_layer_Eta_%d",layer);
      tree->Branch(fpix_eta,        &vFPIX_Eta_[iL]);
      sprintf(fpix_phi, "FPIX_layer_Phi_%d",layer);
      tree->Branch(fpix_phi,        &vFPIX_Phi_[iL]);

      // Histograms for monitoring
      sprintf(htitle,"N(i#phi,i#eta);i#phi;i#eta");
      hFPIX_ECAL[iL][proj] = fs->make<TH2F>(hname, htitle,
        EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT,-ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
      if (proj==0)
        for ( int iz(0); iz < nEE; iz++ ) {
          const char *zside = (iz > 0) ? "p" : "m";
          sprintf(hname, "evt_FPIX_layer%d_EE%s",layer,zside);
          sprintf(htitle,"N(ix,iy);ix;iy");
          hEvt_EE_FPIX[iL][iz] = new TH2F(hname, htitle,
          EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
          5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EE[iz] );
        } // iz
    } // iL

  }//proj
} // branchesEB()


// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________// 

void fillTRKLayerAtECAL_with_EEproj( TH2F *hEvt_EE_SUBDET, std::vector<float> & vSUBDET_ECAL_, TH2F *hSUBDET_ECAL, int ieta_global_offset, int ieta_signed_offset ){
  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float nEntries_=0.;
  for (int ieta = 1; ieta < hEvt_EE_SUBDET->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_SUBDET->GetNbinsX()+1; iphi++) {
      nEntries_ = hEvt_EE_SUBDET->GetBinContent( iphi, ieta );
      if ( (nEntries_ == 0.) ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
      // Fill vector for image
      vSUBDET_ECAL_[idx_] = nEntries_;
      // Fill histogram for monitoring
      hSUBDET_ECAL->Fill( iphi_, ieta_signed_, nEntries_ );
    } // iphi_
  } // ieta_
} // fillTracksAtECAL_with_EEproj


void fillTRKLayerAtECAL_with_EEproj( TH2F *hEvt_EE_SUBDET[][nEE], std::vector<float> vSUBDET_ECAL_[][Nhitproj], TH2F *hSUBDET_ECAL[][Nhitproj], int nSUBDET, unsigned int proj ){
  int ieta_global_offset,ieta_signed_offset;
  for(int nLayer=0; nLayer<nSUBDET; nLayer++){

    // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
    ieta_global_offset = 0;
    ieta_signed_offset = -ECAL_IETA_MAX_EXT;
    fillTRKLayerAtECAL_with_EEproj(hEvt_EE_SUBDET[nLayer][0], vSUBDET_ECAL_[nLayer][proj], hSUBDET_ECAL[nLayer][proj], ieta_global_offset, ieta_signed_offset);

    // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
    ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
    ieta_signed_offset = EB_IETA_MAX;
    fillTRKLayerAtECAL_with_EEproj(hEvt_EE_SUBDET[nLayer][1], vSUBDET_ECAL_[nLayer][proj], hSUBDET_ECAL[nLayer][proj], ieta_global_offset, ieta_signed_offset);
  }
}

void fillTRKLayerAtEB (DetId id, int layer_, unsigned int proj, TH2F *hSUBDET_ECAL[][Nhitproj], std::vector<float> vSUBDET_ECAL_[][Nhitproj] ) {
  int ieta_global_offset = 55;
  EBDetId ebId( id );
  int iphi_ = ebId.iphi() - 1;
  int ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
  int ieta_signed = ieta_;
  int ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
  int idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
  vSUBDET_ECAL_[layer_-1][proj][idx_] += 1.0;
  hSUBDET_ECAL[layer_-1][proj]->Fill( iphi_, ieta_signed, 1. );
}

void fillHelperAtEE ( float phi_, float eta_, int layer_, TH2F *hEvt_EE_SUBDET[][nEE]) {
  int iz_ = (eta_ > 0.) ? 1 : 0;
  hEvt_EE_SUBDET[layer_-1][iz_]->Fill( phi_, eta_);
}

unsigned int RecHitAnalyzer::getLayer(const DetId& detid, const TrackerTopology* tTopo) {
 //                                            uint16_t
 // +------------+---------------+---------------------------+-----------------+----------------+
 // |  tk/mu/mtd | sub-structure |     sub-sub-structure     |     stereo      |    hit type    |
 // +------------+---------------+---------------------------+-----------------+----------------+
 // |    11-10   | 9   8    7    |  6     5     4     3      |        2        |    1        0  |  bit
 // +------------+---------------+---------------------------+-----------------+----------------|
 // | tk  = 1    |    PXB = 1    | layer = 1-3               |                 | hit type = 0-3 |
 // | tk  = 1    |    PXF = 2    | disk  = 1-2               |                 | hit type = 0-3 |
 // | tk  = 1    |    TIB = 3    | layer = 1-4               | 0=rphi,1=stereo | hit type = 0-3 |
 // | tk  = 1    |    TID = 4    | wheel = 1-3               | 0=rphi,1=stereo | hit type = 0-3 |
 // | tk  = 1    |    TOB = 5    | layer = 1-6               | 0=rphi,1=stereo | hit type = 0-3 |
 // | tk  = 1    |    TEC = 6    | wheel = 1-9               | 0=rphi,1=stereo | hit type = 0-3 |
 // | mu  = 0    |    DT  = 1    | 4*(stat-1)+superlayer     |                 | hit type = 0-3 |
 // | mu  = 0    |    CSC = 2    | 4*(stat-1)+(ring-1)       |                 | hit type = 0-3 |
 // | mu  = 0    |    RPC = 3    | 4*(stat-1)+2*layer+region |                 | hit type = 0-3 |
 // | mu  = 0    |    GEM = 4    | 2*(stat-1)+2*(layer-1)    |                 | hit type = 0-3 |
 // | mu  = 0    |    ME0 = 5    | roll                      |                 | hit type = 0-3 |
 // | mtd = 2    |    BTL = 1    | moduleType = 1-3          |                 | hit type = 0-3 | 
 // | mtd = 2    |    ETL = 2    | ring = 1-12               |                 | hit type = 0-3 | 
 // +------------+---------------+---------------------------+-----------------+----------------+
  unsigned int subid=detid.subdetId();
  switch(subid){

    case PixelSubdetector::PixelBarrel:{//BPIX
      PXBDetId pdetId = PXBDetId(detid.rawId());
      return tTopo->pxbLayer(pdetId);
      //return pdetId.layer();
    }break;

    case PixelSubdetector::PixelEndcap:{//FPIX
      PXFDetId pdetId = PXFDetId(detid.rawId());
      return tTopo->pxfDisk(pdetId);
      //return pdetId.disk();
    }break;

    case SiStripDetId::TIB:{//TIB
      SiStripDetId pdetId = SiStripDetId(detid.rawId());
      return tTopo->tibLayer(pdetId);
      //return tTopo->layer(SiStripDetId(detid));
      //return tTopo->tibLayer(detid.rawId());
      //return tTopo->tibLayer(detid.rawId());
    }break;

    case SiStripDetId::TID:{//TID
      SiStripDetId pdetId = SiStripDetId(detid.rawId());
      return tTopo->tidWheel(pdetId);
      //return tTopo->layer(SiStripDetId(detid));
      //return tTopo->tidWheel(SiStripDetId(detid));
    }break;

    case SiStripDetId::TOB:{//TOB
      SiStripDetId pdetId = SiStripDetId(detid.rawId());
      return tTopo->tobLayer(pdetId);
      //return tTopo->layer(SiStripDetId(detid));
      //return tTopo->tobLayer(detid.rawId());
    }break;

    case SiStripDetId::TEC:{//TEC
      SiStripDetId pdetId = SiStripDetId(detid.rawId());
      return tTopo->tecWheel(pdetId);
      //return tTopo->layer(SiStripDetId(detid));
      //return tTopo->tecWheel(detid.rawId());
    }break;

  }
  return 999;
}



// Fill TRK rechits at ECAL stitched ______________________________________________________________//
void RecHitAnalyzer::fillTRKlayersAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int proj ) {
  
  //conf_=iConfig;

  float eta, phi;
  GlobalPoint pos;

  for ( int iL(0); iL < nTOB; iL++ ) {
    vTOB_ECAL_[iL][proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TOB[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nTEC; iL++ ) {
    vTEC_ECAL_[iL][proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TEC[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nTIB; iL++ ) {
    vTIB_ECAL_[iL][proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TIB[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nTID; iL++ ) {
    vTID_ECAL_[iL][proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_TID[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nBPIX; iL++ ) {
    vBPIX_X_[iL].clear();
    vBPIX_Y_[iL].clear();
    vBPIX_Z_[iL].clear();
    vBPIX_Eta_[iL].clear();
    vBPIX_Phi_[iL].clear();
    vBPIX_ECAL_[iL][proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_BPIX[iL][iz]->Reset();
  }
  for ( int iL(0); iL < nFPIX; iL++ ) {
    vFPIX_X_[iL].clear();
    vFPIX_Y_[iL].clear();
    vFPIX_Z_[iL].clear();
    vFPIX_Eta_[iL].clear();
    vFPIX_Phi_[iL].clear();
    vFPIX_ECAL_[iL][proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
    for ( int iz(0); iz < nEE; ++iz ) hEvt_EE_FPIX[iL][iz]->Reset();
  }

  //edm::Handle<TrackingRecHitCollection> TRKRecHitsH_;
  //iEvent.getByToken( TRKRecHitCollectionT_, TRKRecHitsH_ );
  // Provides access to global cell position
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  bool isPVgood=false;
  edm::Handle<reco::VertexCollection> vertexInfo;
  iEvent.getByToken(vertexCollectionT_, vertexInfo);
  //const reco::VertexCollection& vtxs = *vertexInfo;
  isPVgood = vertexInfo.product()->size()>0;
  reco::Vertex the_PV;
  if (isPVgood) the_PV = vertexInfo.product()->at(0);
  TVector3 pv_v(the_PV.x(),the_PV.y(),the_PV.z());

  //sipixel
  edm::Handle<SiPixelRecHitCollection>  recHitColl;
  iEvent.getByToken(siPixelRecHitCollectionT_, recHitColl);

  edm::Handle<SiStripMatchedRecHit2DCollection>  stripMatchedRecHitColl;
  iEvent.getByToken(siStripMatchedRecHitCollectionT_, stripMatchedRecHitColl);

  edm::Handle<SiStripRecHit2DCollection>  stripRPhiRecHitColl;
  iEvent.getByToken(siStripRPhiRecHitCollectionT_, stripRPhiRecHitColl);

  edm::Handle<SiStripRecHit2DCollection>  stripStereoRecHitColl;
  iEvent.getByToken(siStripStereoRecHitCollectionT_, stripStereoRecHitColl);

  edm::ESHandle<TrackerGeometry> geom;
  iSetup.get<TrackerDigiGeometryRecord>().get( geom );
  const TrackerGeometry& theTracker(*geom);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* const tTopo = tTopoHandle.product();

  //std::cout <<" FOUND "<<(recHitColl.product())->dataSize()<<" Pixel Hits" << std::endl;

  SiPixelRecHitCollection::const_iterator recHitIdIterator    = (recHitColl.product())->begin();
  SiPixelRecHitCollection::const_iterator recHitIdIteratorEnd = (recHitColl.product())->end();

  for ( ; recHitIdIterator != recHitIdIteratorEnd; recHitIdIterator++)
  {
    SiPixelRecHitCollection::DetSet detset = *recHitIdIterator;
    DetId detId = DetId(detset.detId()); // Get the Detid object
    //uint32_t TheID = detset.first;
    //DetId detId = DetId(TheID); // Get the Detid object
    unsigned int detType=detId.det(); // det type, tracker=1
    unsigned int subid=detId.subdetId(); //subdetector type, barrel=1, fpix=2
    if( detset.empty() ) {
      std::cout << "detset is empty" << std::endl;
      continue;
    }
    unsigned int layer = getLayer(detId, tTopo);
    //std::cout<<"Pixel Id = : "<<detId.rawId()<<" "<<detId.null()<<" , type = "<<detType<<" - subId ( 1->bpix | 2->fpix ) = "<< subid << " - Layer = " << layer << std::endl;
    const PixelGeomDetUnit* theGeomDet  = dynamic_cast<const PixelGeomDetUnit*> (theTracker.idToDetUnit(detId) );
    unsigned int iPixelHit = 0; 
    SiPixelRecHitCollection::DetSet::const_iterator pixeliter=detset.begin();
    SiPixelRecHitCollection::DetSet::const_iterator rechitRangeIteratorEnd   = detset.end();
    for(;pixeliter!=rechitRangeIteratorEnd;++pixeliter)
    {//loop on the rechit
      if (pixeliter->isValid())
      {
        iPixelHit++;
        LocalPoint lp = pixeliter->localPosition();
        GlobalPoint GP = theGeomDet->surface().toGlobal(Local3DPoint(lp));
        phi=0.;
        eta=0.;
        //std::cout << " " << iPixelHit << " | global position: x = " << GP.x() << " , y = "<< GP.y() << " , z = " << GP.z() <<std::endl;
        
        TVector3 GP_v(GP.x(),GP.y(),GP.z());
        //std::cout << " " << iPixelHit << " | P.V correction global position: x = " << GP_v.x() << " , y = "<< GP_v.y() << " , z = " << GP_v.z() <<std::endl;
        GP_v=GP_v-pv_v;
        phi=GP_v.Phi();
        eta=GP_v.Eta();

        if(subid==1){
	vBPIX_X_[layer].push_back(GP_v.x());
	vBPIX_Y_[layer].push_back(GP_v.y());
	vBPIX_Z_[layer].push_back(GP_v.z());
	vBPIX_Eta_[layer].push_back(GP_v.Eta());
	vBPIX_Phi_[layer].push_back(GP_v.Phi());
 	}

        if(subid==2){
        vFPIX_X_[layer].push_back(GP_v.x());
        vFPIX_Y_[layer].push_back(GP_v.y());
        vFPIX_Z_[layer].push_back(GP_v.z());
        vFPIX_Eta_[layer].push_back(GP_v.Eta());
        vFPIX_Phi_[layer].push_back(GP_v.Phi());

        }

        /*switch (proj)
        {
          case 1:
          {
            phi = GP.phi();
            eta = GP.eta();
            break;
          }
          case 0:
          {
            TVector3 GP_v(GP.x(),GP.y(),GP.z());
            GP_v=GP_v-pv_v;
            phi=GP_v.Phi();
            eta=GP_v.Eta();
            break;
          }
          default:
          {
            phi=0.;
            eta=0.;
            break;
          }
        }*/
        //if ( std::abs(eta) > 3. ) continue;
        DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
        if ( subid == PixelSubdetector::PixelBarrel ){
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hBPIX_ECAL, vBPIX_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_BPIX );
          }
        }
        else if ( subid == PixelSubdetector::PixelEndcap )
        {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hFPIX_ECAL, vFPIX_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_FPIX);
          }
        }
      } else std::cout << "!!!!!!!!!!!!!! NO PIXEL HITS ARE VALID !!!!!!!!!!!!!!" << std::endl;
    } //std::cout << "End of PixelRecHit " << iPixelHit << std::endl;
  }



  // --  siSTRIP --

  // MATCHED REC HIT COLLECTION
  for ( SiStripMatchedRecHit2DCollection::const_iterator detunit_iterator = stripMatchedRecHitColl->begin(), detunit_end = stripMatchedRecHitColl->end(); detunit_iterator != detunit_end; ++detunit_iterator) {
    SiStripMatchedRecHit2DCollection::DetSet rechitRange = *detunit_iterator;
    DetId detId = DetId(detunit_iterator->detId());
    unsigned int id = detunit_iterator->detId();
    unsigned int subid=detId.subdetId();
    unsigned int layer = getLayer(id, tTopo);
    //std::cout << "Strip Id = " << id << " - subId ( 3->TIB | 4->TID | 5->TOB | 6->TEC ) = " << subid << " - Layer = " << layer <<  std::endl;
    const StripGeomDetUnit* stripDet = (const StripGeomDetUnit*)theTracker.idToDet(detId);
    if(stripDet==0) {
      std::cout << "SiStripRecHitConverter: Detid=" << id << " not found, trying next one" << std::endl;
      continue;
    }
    const StripTopology * stripTopol = (StripTopology*)(&stripDet->topology()); ;
    SiStripMatchedRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorBegin = rechitRange.begin();
    SiStripMatchedRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorEnd   = rechitRange.end();
    SiStripMatchedRecHit2DCollection::DetSet::const_iterator stripiter=rechitRangeIteratorBegin;
    unsigned int iRecHit = 0;
    for(stripiter=rechitRangeIteratorBegin;stripiter!=rechitRangeIteratorEnd;++stripiter){//loop on the rechit
      if (stripiter->isValid()){
        iRecHit++;
        SiStripMatchedRecHit2D const rechit=*stripiter;
        const GeomDet* stripdet=rechit.det();
        //DetId stripid=rechit.geographicalId();
        //std::vector<const SiStripCluster*> clust=rechit.cluster();
        LocalPoint lp = rechit.localPosition();
        GlobalPoint GP = stripDet->surface().toGlobal(Local3DPoint(lp));
        //std::cout << " " << iRecHit << " | global position: x = " << GP.x() << " , y = "<< GP.y() << " , z = " << GP.z() <<std::endl;
        TVector3 GP_v(GP.x(),GP.y(),GP.z());
        GP_v=GP_v-pv_v;
        phi=GP_v.Phi();
        eta=GP_v.Eta();
        /*switch (proj)
        {
          case 0:
          {
            phi = GP.phi();
            eta = GP.eta();
            std::cout << "hits not corrected:" << std::endl;
            break;
          }
          case 1:
          {
            TVector3 GP_v(GP.x(),GP.y(),GP.z());
            GP_v=GP_v-pv_v;
            phi=GP_v.Phi();
            eta=GP_v.Eta();
            std::cout << "hit position corrected:" << std::endl;
            break;
          }
          default:
          {
            phi=0.;
            eta=0.;
            std::cout << "hit position default:" << std::endl;
            break;
          }
        }*/
        //if ( std::abs(eta) > 3. ) continue;
        DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
        if ( subid == StripSubdetector::TOB ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTOB_ECAL, vTOB_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TOB);
          }
        }
        else if ( subid == StripSubdetector::TIB ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTIB_ECAL, vTIB_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TIB);
          }
        }
        else if ( subid == StripSubdetector::TEC ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTEC_ECAL, vTEC_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TEC);
          }
        }
        else if ( subid == StripSubdetector::TID ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTID_ECAL, vTID_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TID);
          }
        }
      } else std::cout << "!!!!!!!!!!!!!! NO MATCHED STRIP HITS ARE VALID !!!!!!!!!!!!!!" << std::endl;
    } //std::cout << "End of StripReCHit " << iRecHit << std::endl;
  } // end loop over detectors 


  // RPHI REC HIT COLLECTION
  for ( SiStripRecHit2DCollection::const_iterator detunit_iterator = stripRPhiRecHitColl->begin(), detunit_end = stripRPhiRecHitColl->end(); detunit_iterator != detunit_end; ++detunit_iterator) {
    SiStripRecHit2DCollection::DetSet rechitRange = *detunit_iterator;
    DetId detId = DetId(detunit_iterator->detId());
    unsigned int id = detunit_iterator->detId();
    unsigned int subid=detId.subdetId();
    unsigned int layer = getLayer(id, tTopo);
    //std::cout << "Strip Id = " << id << " - subId ( 3->TIB | 4->TID | 5->TOB | 6->TEC ) = " << subid << " - Layer = " << layer <<  std::endl;
    const StripGeomDetUnit* stripDet = (const StripGeomDetUnit*)theTracker.idToDet(detId);
    if(stripDet==0) {
      std::cout << "SiStripRecHitConverter: Detid=" << id << " not found, trying next one" << std::endl;
      continue;
    }
    const StripTopology * stripTopol = (StripTopology*)(&stripDet->topology()); ;
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorBegin = rechitRange.begin();
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorEnd   = rechitRange.end();
    SiStripRecHit2DCollection::DetSet::const_iterator stripiter=rechitRangeIteratorBegin;
    unsigned int iRecHit = 0;
    for(stripiter=rechitRangeIteratorBegin;stripiter!=rechitRangeIteratorEnd;++stripiter){//loop on the rechit
      if (stripiter->isValid()){
        iRecHit++;
        SiStripRecHit2D const rechit=*stripiter;
        const GeomDet* stripdet=rechit.det();
        //DetId stripid=rechit.geographicalId();
        //std::vector<const SiStripCluster*> clust=rechit.cluster();
        LocalPoint lp = rechit.localPosition();
        GlobalPoint GP = stripDet->surface().toGlobal(Local3DPoint(lp));
        //std::cout << " " << iRecHit << " | global position: x = " << GP.x() << " , y = "<< GP.y() << " , z = " << GP.z() <<std::endl;
        TVector3 GP_v(GP.x(),GP.y(),GP.z());
        GP_v=GP_v-pv_v;
        phi=GP_v.Phi();
        eta=GP_v.Eta();
        /*switch (proj)
        {
          case 0:
          {
            phi = GP.phi();
            eta = GP.eta();
            break;
          }
          case 1:
          {
            TVector3 GP_v(GP.x(),GP.y(),GP.z());
            GP_v=GP_v-pv_v;
            phi=GP_v.Phi();
            eta=GP_v.Eta();
            break;
          }
          default:
          {
            phi=0.;
            eta=0.;
            break;
          }
        }*/
        //if ( std::abs(eta) > 3. ) continue;
        DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
        if ( subid == StripSubdetector::TOB ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTOB_ECAL, vTOB_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TOB);
          }
        }
        else if ( subid == StripSubdetector::TIB ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTIB_ECAL, vTIB_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TIB);
          }
        }
        else if ( subid == StripSubdetector::TEC ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTEC_ECAL, vTEC_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TEC);
          }
        }
        else if ( subid == StripSubdetector::TID ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTID_ECAL, vTID_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TID);
          }
        }
      } else std::cout << "!!!!!!!!!!!!!! NO RPHI STRIP HITS ARE VALID !!!!!!!!!!!!!!" << std::endl;
    } //std::cout << "End of StripReCHit " << iRecHit << std::endl;
  } // end loop over detectors 

  // STEREO REC HIT COLLECTION

  for ( SiStripRecHit2DCollection::const_iterator detunit_iterator = stripStereoRecHitColl->begin(), detunit_end = stripStereoRecHitColl->end(); detunit_iterator != detunit_end; ++detunit_iterator) {
    SiStripRecHit2DCollection::DetSet rechitRange = *detunit_iterator;
    DetId detId = DetId(detunit_iterator->detId());
    unsigned int id = detunit_iterator->detId();
    unsigned int subid=detId.subdetId();
    unsigned int layer = getLayer(id, tTopo);
    //std::cout << "Strip Id = " << id << " - subId ( 3->TIB | 4->TID | 5->TOB | 6->TEC ) = " << subid << " - Layer = " << layer <<  std::endl;
    const StripGeomDetUnit* stripDet = (const StripGeomDetUnit*)theTracker.idToDet(detId);
    if(stripDet==0) {
      std::cout << "SiStripRecHitConverter: Detid=" << id << " not found, trying next one" << std::endl;
      continue;
    }
    const StripTopology * stripTopol = (StripTopology*)(&stripDet->topology()); ;
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorBegin = rechitRange.begin();
    SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorEnd   = rechitRange.end();
    SiStripRecHit2DCollection::DetSet::const_iterator stripiter=rechitRangeIteratorBegin;
    unsigned int iRecHit = 0;
    for(stripiter=rechitRangeIteratorBegin;stripiter!=rechitRangeIteratorEnd;++stripiter){//loop on the rechit
      if (stripiter->isValid()){
        iRecHit++;
        SiStripRecHit2D const rechit=*stripiter;
        const GeomDet* stripdet=rechit.det();
        //DetId stripid=rechit.geographicalId();
        //std::vector<const SiStripCluster*> clust=rechit.cluster();
        LocalPoint lp = rechit.localPosition();
        GlobalPoint GP = stripDet->surface().toGlobal(Local3DPoint(lp));
        //std::cout << " " << iRecHit << " | global position: x = " << GP.x() << " , y = "<< GP.y() << " , z = " << GP.z() <<std::endl;
        TVector3 GP_v(GP.x(),GP.y(),GP.z());
        GP_v=GP_v-pv_v;
        phi=GP_v.Phi();
        eta=GP_v.Eta();
        
        /*switch (proj)
        {
          case 0:
          {
            phi = GP.phi();
            eta = GP.eta();
            break;
          }
          case 1:
          {
            TVector3 GP_v(GP.x(),GP.y(),GP.z());
            GP_v=GP_v-pv_v;
            phi=GP_v.Phi();
            eta=GP_v.Eta();
            break;
          }
          default:
          {
            phi=0.;
            eta=0.;
            break;
          }
        }*/
        //if ( std::abs(eta) > 3. ) continue;
        DetId ecalId( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
        if ( subid == StripSubdetector::TOB ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTOB_ECAL, vTOB_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TOB);
          }
        }
        else if ( subid == StripSubdetector::TIB ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTIB_ECAL, vTIB_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TIB);
          }
        }
        else if ( subid == StripSubdetector::TEC ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTEC_ECAL, vTEC_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TEC);
          }
        }
        else if ( subid == StripSubdetector::TID ) {
          if ( ecalId.subdetId() == EcalBarrel ){
            fillTRKLayerAtEB ( ecalId, layer, proj, hTID_ECAL, vTID_ECAL_ );
          }
          else if ( ecalId.subdetId() == EcalEndcap ){
            fillHelperAtEE ( phi, eta, layer, hEvt_EE_TID);
          }
        }
      } else std::cout << "!!!!!!!!!!!!!! NO STEREO STRIP HITS ARE VALID !!!!!!!!!!!!!!" << std::endl;
    } //std::cout << "End of StripReCHit " << iRecHit << std::endl;
  } // end loop over detectors 



  fillTRKLayerAtECAL_with_EEproj( hEvt_EE_BPIX, vBPIX_ECAL_, hBPIX_ECAL, nBPIX, proj);
  fillTRKLayerAtECAL_with_EEproj( hEvt_EE_FPIX, vFPIX_ECAL_, hFPIX_ECAL, nFPIX, proj);
  fillTRKLayerAtECAL_with_EEproj( hEvt_EE_TOB, vTOB_ECAL_, hTOB_ECAL, nTOB, proj);
  fillTRKLayerAtECAL_with_EEproj( hEvt_EE_TIB, vTIB_ECAL_, hTIB_ECAL, nTIB, proj);
  fillTRKLayerAtECAL_with_EEproj( hEvt_EE_TEC, vTEC_ECAL_, hTEC_ECAL, nTEC, proj);
  fillTRKLayerAtECAL_with_EEproj( hEvt_EE_TID, vTID_ECAL_, hTID_ECAL, nTID, proj);

} // fillEB()

