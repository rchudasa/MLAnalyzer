#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
//#include "Calibration/IsolatedParticles/interface/DetIdFromEtaPhi.h"

// Fill stitched EEm_EB_EEp image /////////////////////............/
// Store all ECAL event rechits into a stitched EEm_EB_EEp image 
// segmented in iphi,ieta spannning the full -3 < eta < 3. 
// Use EB-like granularity giving an extended range ieta=[-140,140].
//
// For endcaps, project EE hits into a helper histogram binned by 
// phi,eta before filling the full extended ECAL(iphi,eta) image.
// For barrel, fill EB hits directly since geometries are 1:1. 
//
// 'ieta_global' keeps track of the global ieta index count used
// for filling the extended image vector vECAL_tracksPt.
// 'ieta_signed' keeps track of the position along [-140,140] used
// for filling the monitoring histogram hECAL_tracksPt.

TH2F *hEvt_EE_tracks[nEE];
TH2F *hEvt_EE_tracksPt[nEE];
TH2F *hEvt_EE_tracksQPt[nEE];
TH2F *hEvt_EE_tracksD0[nEE];
TH2F *hEvt_EE_tracksDz[nEE];
TH2F *hEvt_EE_tracksPt_max[nEE];
TH2F *hEvt_EE_tracksD0_max[nEE];
TH2F *hEvt_EE_tracksDz_max[nEE];


TProfile2D *hECAL_tracks[Nproj];
TProfile2D *hECAL_tracksPt[Nproj];
TProfile2D *hECAL_tracksQPt[Nproj];
TProfile2D *hECAL_tracksD0[Nproj];
TProfile2D *hECAL_tracksDz[Nproj];

std::vector<float> vECAL_tracksPt_[Nproj];
std::vector<float> vECAL_tracksQPt_[Nproj];
std::vector<float> vECAL_tracksD0_[Nproj];
std::vector<float> vECAL_tracksDz_[Nproj];
std::vector<float> vECAL_tracks_[Nproj];
std::vector<float> vECAL_tracksPt_max_[Nproj];
std::vector<float> vECAL_tracksD0_max_[Nproj];
std::vector<float> vECAL_tracksDz_max_[Nproj];
// //at ECAL
// std::vector<float> vECAL_tracksPt_atECAL_;
// std::vector<float> vECAL_tracksQPt_atECAL_;
// std::vector<float> vECAL_tracksD0_atECAL_;
// std::vector<float> vECAL_tracksDz_atECAL_;
// std::vector<float> vECAL_tracks_atECAL_;
// std::vector<float> vECAL_tracksPt_max_atECAL_;
// std::vector<float> vECAL_tracksD0_max_atECAL_;
// std::vector<float> vECAL_tracksDz_max_atECAL_;
// //at HCAL
// std::vector<float> vECAL_tracksPt_atHCAL_;
// std::vector<float> vECAL_tracksQPt_atHCAL_;
// std::vector<float> vECAL_tracksD0_atHCAL_;
// std::vector<float> vECAL_tracksDz_atHCAL_;
// std::vector<float> vECAL_tracks_atHCAL_;
// std::vector<float> vECAL_tracksPt_max_atHCAL_;
// std::vector<float> vECAL_tracksD0_max_atHCAL_;
// std::vector<float> vECAL_tracksDz_max_atHCAL_;

// TODO Take the D0 and Dz values from the highest pT track for endcap regions

// Initialize branches _______________________________________________________________//
void RecHitAnalyzer::branchesTracksAtECALstitched ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Intermediate helper histogram (single event only)
  // Track Intermediate Hists
  hEvt_EE_tracks[0] = new TH2F("evt_EEm_tracks", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracks[1] = new TH2F("evt_EEp_tracks", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  // Track Parameter Intermediate Hists
  hEvt_EE_tracksPt[0] = new TH2F("evt_EEm_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPt[1] = new TH2F("evt_EEp_tracksPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksQPt[0] = new TH2F("evt_EEm_tracksQPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksQPt[1] = new TH2F("evt_EEp_tracksQPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksD0[0] = new TH2F("evt_EEm_tracksD0", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksD0[1] = new TH2F("evt_EEp_tracksD0", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksDz[0] = new TH2F("evt_EEm_tracksDz", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksDz[1] = new TH2F("evt_EEp_tracksDz", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  // Track Max Parameter Intermediate Hists
  hEvt_EE_tracksPt_max[0] = new TH2F("evt_EEm_tracksPt_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksPt_max[1] = new TH2F("evt_EEp_tracksPt_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksD0_max[0] = new TH2F("evt_EEm_tracksD0_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksD0_max[1] = new TH2F("evt_EEp_tracksD0_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );
  hEvt_EE_tracksDz_max[0] = new TH2F("evt_EEm_tracksDz_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEm );
  hEvt_EE_tracksDz_max[1] = new TH2F("evt_EEp_tracksDz_maxPt", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX, -TMath::Pi(), TMath::Pi(),
      5*(HBHE_IETA_MAX_HE-1-HBHE_IETA_MAX_EB), eta_bins_EEp );

  
  for (unsigned int proj=0; proj<Nproj; proj++)
  {
    // Branches for images
    tree->Branch((std::string("ECAL_tracks")+projections[proj]).c_str(),      &(vECAL_tracks_[proj]));
    tree->Branch((std::string("ECAL_tracksPt")+projections[proj]).c_str(),    &(vECAL_tracksPt_[proj]));
    tree->Branch((std::string("ECAL_tracksQPt")+projections[proj]).c_str(),    &(vECAL_tracksPt_[proj]));
    tree->Branch((std::string("ECAL_tracksD0")+projections[proj]).c_str(),    &(vECAL_tracksD0_[proj]));
    tree->Branch((std::string("ECAL_tracksDz")+projections[proj]).c_str(),    &(vECAL_tracksDz_[proj]));
    tree->Branch((std::string("ECAL_tracksPt_maxPt")+projections[proj]).c_str(),    &(vECAL_tracksPt_max_[proj]));
    tree->Branch((std::string("ECAL_tracksD0_maxPt")+projections[proj]).c_str(),    &(vECAL_tracksD0_max_[proj]));
    tree->Branch((std::string("ECAL_tracksDz_maxPt")+projections[proj]).c_str(),    &(vECAL_tracksDz_max_[proj]));

    // Histograms for monitoring
    hECAL_tracks[proj] = fs->make<TProfile2D>((std::string("ECAL_tracks")+projections[proj]).c_str(), "E(i#phi,i#eta);i#phi;i#eta",
        EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
    hECAL_tracksPt[proj] = fs->make<TProfile2D>((std::string("ECAL_tracksPt")+projections[proj]).c_str(), "E(i#phi,i#eta);i#phi;i#eta",
        EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
    hECAL_tracksQPt[proj] = fs->make<TProfile2D>((std::string("ECAL_tracksQPt")+projections[proj]).c_str(), "E(i#phi,i#eta);i#phi;i#eta",
        EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
    hECAL_tracksD0[proj] = fs->make<TProfile2D>((std::string("ECAL_tracksD0")+projections[proj]).c_str(), "E(i#phi,i#eta);i#phi;i#eta",
        EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
    hECAL_tracksDz[proj] = fs->make<TProfile2D>((std::string("ECAL_tracksDz")+projections[proj]).c_str(), "E(i#phi,i#eta);i#phi;i#eta",
        EB_IPHI_MAX,    EB_IPHI_MIN-1, EB_IPHI_MAX,
        2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
  }

} // branchesTracksAtECALstitched()

// Function to map EE(phi,eta) histograms to ECAL(iphi,ieta) vector _______________________________//
//void fillTracksAtECAL_with_EEproj ( TH2F *hEvt_EE_tracks_, TH2F *hEvt_EE_tracksPt_, TH2F *hEvt_EE_tracksQPt_, TH2F *hEvt_EE_tracksD0_, TH2F *hEvt_EE_tracksDz_,  TH2F *hEvt_EE_tracksPt_max_, TH2F *hEvt_EE_tracksD0_max_, TH2F *hEvt_EE_tracksDz_max_, int ieta_global_offset, int ieta_signed_offset ) {
void fillTracksAtECAL_with_EEproj ( int side, int ieta_global_offset, int ieta_signed_offset, int proj ) {

  int ieta_global_, ieta_signed_;
  int ieta_, iphi_, idx_;
  float track_;
  float trackPt_;
  float trackQPt_;
  float trackD0_; 
  float trackDz_;
  float trackPt_max_;
  float trackD0_max_;
  float trackDz_max_;
  for (int ieta = 1; ieta < hEvt_EE_tracksPt[side]->GetNbinsY()+1; ieta++) {
    ieta_ = ieta - 1;
    ieta_global_ = ieta_ + ieta_global_offset;
    ieta_signed_ = ieta_ + ieta_signed_offset;
    for (int iphi = 1; iphi < hEvt_EE_tracksPt[side]->GetNbinsX()+1; iphi++) {

      track_  = hEvt_EE_tracks[side]->GetBinContent( iphi, ieta );
      trackPt_ = hEvt_EE_tracksPt[side]->GetBinContent( iphi, ieta );
      trackQPt_ = hEvt_EE_tracksQPt[side]->GetBinContent( iphi, ieta );
      trackD0_ = hEvt_EE_tracksD0[side]->GetBinContent( iphi, ieta );
      trackDz_ = hEvt_EE_tracksDz[side]->GetBinContent( iphi, ieta );
      trackPt_max_ = hEvt_EE_tracksPt_max[side]->GetBinContent( iphi, ieta );
      trackD0_max_ = hEvt_EE_tracksD0_max[side]->GetBinContent( iphi, ieta );
      trackDz_max_ = hEvt_EE_tracksDz_max[side]->GetBinContent( iphi, ieta );
      if ( trackPt_ <= zs ) continue;
      // NOTE: EB iphi = 1 does not correspond to physical phi = -pi so need to shift!
      iphi_ = iphi  + 5*38; // shift
      iphi_ = iphi_ > EB_IPHI_MAX ? iphi_-EB_IPHI_MAX : iphi_; // wrap-around
      iphi_ = iphi_ - 1;
      idx_  = ieta_global_*EB_IPHI_MAX + iphi_;
      // Fill vector for image
      vECAL_tracks_[proj][idx_] = track_;
      vECAL_tracksPt_[proj][idx_] = trackPt_;
      vECAL_tracksQPt_[proj][idx_] = trackQPt_;
      vECAL_tracksD0_[proj][idx_] = trackD0_;
      vECAL_tracksDz_[proj][idx_] = trackDz_;
      vECAL_tracksPt_max_[proj][idx_] = trackPt_max_;
      vECAL_tracksD0_max_[proj][idx_] = trackD0_max_;
      vECAL_tracksDz_max_[proj][idx_] = trackDz_max_;
      // Fill histogram for monitoring
      hECAL_tracks[proj]->Fill( iphi_, ieta_signed_, track_ );
      hECAL_tracksPt[proj]->Fill( iphi_, ieta_signed_, trackPt_ );
      hECAL_tracksQPt[proj]->Fill( iphi_, ieta_signed_, trackQPt_ );
      hECAL_tracksD0[proj]->Fill( iphi_, ieta_signed_, trackD0_ );
      hECAL_tracksDz[proj]->Fill( iphi_, ieta_signed_, trackDz_ );

    } // iphi_
  } // ieta_

} // fillTracksAtECAL_with_EEproj

// Fill stitched EE-, EB, EE+ rechits ________________________________________________________//
void RecHitAnalyzer::fillTracksAtECALstitched ( const edm::Event& iEvent, const edm::EventSetup& iSetup, unsigned int proj ) {
  int iphi_, ieta_, iz_, idx_;
  int ieta_global, ieta_signed;
  int ieta_global_offset, ieta_signed_offset;
  float eta, phi, trackPt_, trackQPt_,trackD0_, trackDz_;
  GlobalPoint pos;

  edm::ESHandle<MagneticField> magfield;
  iSetup.get<IdealMagneticFieldRecord>().get(magfield);

  vECAL_tracks_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksPt_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksQPt_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksD0_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksDz_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksPt_max_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksD0_max_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  vECAL_tracksDz_max_[proj].assign( 2*ECAL_IETA_MAX_EXT*EB_IPHI_MAX, 0. );
  for ( int iz(0); iz < nEE; ++iz ) {
    hEvt_EE_tracks[iz]->Reset();
    hEvt_EE_tracksPt[iz]->Reset();
    hEvt_EE_tracksQPt[iz]->Reset();
    hEvt_EE_tracksD0[iz]->Reset();
    hEvt_EE_tracksDz[iz]->Reset();
    hEvt_EE_tracksPt_max[iz]->Reset();
    hEvt_EE_tracksD0_max[iz]->Reset();
    hEvt_EE_tracksDz_max[iz]->Reset();
  }

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_ );
  edm::Handle<EcalRecHitCollection> EERecHitsH_;
  iEvent.getByToken( EERecHitCollectionT_, EERecHitsH_ );
  edm::ESHandle<CaloGeometry> caloGeomH_;
  iSetup.get<CaloGeometryRecord>().get( caloGeomH_ );
  const CaloGeometry* caloGeom = caloGeomH_.product();

  edm::Handle<reco::TrackCollection> tracksH_;
  iEvent.getByToken( trackCollectionT_, tracksH_ );

  bool isPVgood=false;
  edm::Handle<reco::VertexCollection> pvColl;
  iEvent.getByToken(vertexCollectionT_, pvColl);
 
  isPVgood = pvColl.product()->size()>0;
  reco::Vertex the_PV;
  if (isPVgood) the_PV = pvColl.product()->at(0);


  reco::Track::TrackQuality tkQt_ = reco::Track::qualityByName("highPurity");

  int bin;
  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) {

    if(iTk->pt()<=0.5)continue; 
    if(iTk->charge()==0) continue;// NO neutral objects
    if ( !(iTk->quality(tkQt_)) ) continue;


    if (proj==4)
    {
      bool pv_match=false;
      for ( reco::Vertex::trackRef_iterator iTkPV = the_PV.tracks_begin();
            iTkPV != the_PV.tracks_end(); ++iTkPV ) {
        if (fabs(iTk->pt()-iTkPV->get()->pt())<0.001 &&
            fabs(iTk->eta()-iTkPV->get()->eta())<0.001 &&
            fabs(iTk->phi()-iTkPV->get()->phi())<0.001)
        {
          pv_match=true;
          break;
        }
      }
      if (!pv_match)continue;
    }    


    bool isPropagationOk=false;
    eta = 0.;
    phi = 0.;
    switch (proj)
    {
      case 0:
      {
        eta = iTk->eta();
        phi = iTk->phi();
        isPropagationOk=true;
      }
      break;

      case 1: case 3: case 4:
      {
        auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
        isPropagationOk=propagatedECALTrack.ok;
        if (propagatedECALTrack.ok)
        {
          eta = propagatedECALTrack.direction.eta();
          phi = propagatedECALTrack.direction.phi();
        }
      }
      break;

      case 2:
      {
        auto propagatedHCALTrack = spr::propagateTrackToHCAL(&*iTk, magfield.product());
        isPropagationOk=propagatedHCALTrack.ok;
        if (propagatedHCALTrack.ok)
        {
          eta = propagatedHCALTrack.direction.eta();
          phi = propagatedHCALTrack.direction.phi();
        }
      }
      break;

      default:
      {
        isPropagationOk=false;
      }
      break;
    }
    
    if ( std::abs(eta) > 3. || !isPropagationOk) continue;
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalBarrel ) continue;
    if ( id.subdetId() == EcalEndcap ) {
      iz_ = (eta > 0.) ? 1 : 0;
      switch(proj)
      {
        case 0: case 1: case 2:
        {
          trackD0_ = iTk->d0();
          trackDz_ = iTk->dz();
          break;
        }
        case 3: case 4: default:
        {
          trackD0_ = -iTk->dxy(the_PV.position());
          trackDz_ = iTk->dz(the_PV.position());
          break;
        }
      }
      // Fill intermediate helper histogram by eta,phi
      hEvt_EE_tracks[iz_]->Fill( phi, eta );
      hEvt_EE_tracksPt[iz_]->Fill( phi, eta, iTk->pt() );
      hEvt_EE_tracksQPt[iz_]->Fill( phi, eta, iTk->pt()*iTk->charge() );
      hEvt_EE_tracksD0[iz_]->Fill( phi, eta, trackD0_ );
      hEvt_EE_tracksDz[iz_]->Fill( phi, eta, trackDz_ );
      bin = hEvt_EE_tracks[iz_]->FindBin( phi, eta );
      if ( iTk->pt() > hEvt_EE_tracksPt_max[iz_]->GetBinContent( bin ) ) {
        hEvt_EE_tracksPt_max[iz_]->SetBinContent( bin, iTk->pt() );
        hEvt_EE_tracksD0_max[iz_]->SetBinContent( bin, trackD0_ );
        hEvt_EE_tracksDz_max[iz_]->SetBinContent( bin, trackDz_ );
      }
    }
  } // tracks


  // Map EE-(phi,eta) to bottom part of ECAL(iphi,ieta)
  ieta_global_offset = 0;
  ieta_signed_offset = -ECAL_IETA_MAX_EXT;
  fillTracksAtECAL_with_EEproj( 0, ieta_global_offset, ieta_signed_offset, proj );

  // Fill middle part of ECAL(iphi,ieta) with the EB rechits.
  ieta_global_offset = 55;

  for ( reco::TrackCollection::const_iterator iTk = tracksH_->begin();
        iTk != tracksH_->end(); ++iTk ) { 
    if ( !(iTk->quality(tkQt_)) ) continue;
    if(iTk->pt()<=0.5)continue; 
    if(iTk->charge()==0) continue;


    if (proj==4)
    {
      bool pv_match=false;
      for ( reco::Vertex::trackRef_iterator iTkPV = the_PV.tracks_begin();
            iTkPV != the_PV.tracks_end(); ++iTkPV ) {
        if (fabs(iTk->pt()-iTkPV->get()->pt())<0.001 &&
            fabs(iTk->eta()-iTkPV->get()->eta())<0.001 &&
            fabs(iTk->phi()-iTkPV->get()->phi())<0.001)
        {
          pv_match=true;
          break;
        }
      }
      if (!pv_match)continue;
    }


    bool isPropagationOk=false;
    eta = 0.;
    phi = 0.;
    switch (proj)
    {
      case 0:
      {
        eta = iTk->eta();
        phi = iTk->phi();
        isPropagationOk=true;
      }
      break;

      case 1: case 3: case 4:
      {
        auto propagatedECALTrack = spr::propagateTrackToECAL(&*iTk, magfield.product());
        isPropagationOk=propagatedECALTrack.ok;
        if (propagatedECALTrack.ok)
        {
          eta = propagatedECALTrack.direction.eta();
          phi = propagatedECALTrack.direction.phi();
        }
      }
      break;

      case 2:
      {
        auto propagatedHCALTrack = spr::propagateTrackToHCAL(&*iTk, magfield.product());
        isPropagationOk=propagatedHCALTrack.ok;
        if (propagatedHCALTrack.ok)
        {
          eta = propagatedHCALTrack.direction.eta();
          phi = propagatedHCALTrack.direction.phi();
        }
      }
      break;

      default:
      {
        isPropagationOk=false;
      }
      break;
    }
    
    if ( std::abs(eta) > 3. || !isPropagationOk) continue;

    trackPt_ = iTk->pt();
    trackQPt_ = iTk->pt()*iTk->charge();
    switch(proj)
    {
      case 0: case 1: case 2:
      {
        trackD0_ = iTk->d0();
        trackDz_ = iTk->dz();
        break;
      }
      case 3: case 4: default:
      {
        trackD0_ = -iTk->dxy(the_PV.position());
        trackDz_ = iTk->dz(the_PV.position());
        break;
      }
    }
    
    DetId id( spr::findDetIdECAL( caloGeom, eta, phi, false ) );
    if ( id.subdetId() == EcalEndcap ) continue;
    if ( id.subdetId() == EcalBarrel ) { 
      EBDetId ebId( id );
      iphi_ = ebId.iphi() - 1;
      ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
      if ( trackPt_ <= zs ) continue;
      // Fill vector for image
      ieta_signed = ieta_;
      ieta_global = ieta_ + EB_IETA_MAX + ieta_global_offset;
      idx_ = ieta_global*EB_IPHI_MAX + iphi_; 
      vECAL_tracks_[proj][idx_] += 1.;
      vECAL_tracksPt_[proj][idx_] += trackPt_;
      vECAL_tracksQPt_[proj][idx_] += trackQPt_;
      vECAL_tracksD0_[proj][idx_] += trackD0_;
      vECAL_tracksDz_[proj][idx_] += trackDz_;
      if ( trackPt_ > vECAL_tracksPt_max_[proj][idx_] ) {
        vECAL_tracksPt_max_[proj][idx_] = trackPt_;
        vECAL_tracksD0_max_[proj][idx_] = trackD0_;
        vECAL_tracksDz_max_[proj][idx_] = trackDz_;
      }
      // Fill histogram for monitoring
      hECAL_tracks[proj]->Fill( iphi_, ieta_signed, 1.0 );
      hECAL_tracksPt[proj]->Fill( iphi_, ieta_signed, trackPt_ );
      hECAL_tracksQPt[proj]->Fill( iphi_, ieta_signed, trackQPt_ );
      hECAL_tracksD0[proj]->Fill( iphi_, ieta_signed, trackD0_ );
      hECAL_tracksDz[proj]->Fill( iphi_, ieta_signed, trackDz_ );
    }

  } // EB



  // Map EE+(phi,eta) to upper part of ECAL(iphi,ieta)
  ieta_global_offset = ECAL_IETA_MAX_EXT + EB_IETA_MAX;
  ieta_signed_offset = EB_IETA_MAX;
  fillTracksAtECAL_with_EEproj( 1, ieta_global_offset, ieta_signed_offset, proj );

  // Get average D0 and Dz for each position
  for (unsigned int idx_=0;idx_<vECAL_tracks_[proj].size();idx_++) {
    if (vECAL_tracks_[proj][idx_] != 0) {
      vECAL_tracksD0_[proj][idx_] = vECAL_tracksD0_[proj][idx_] / vECAL_tracks_[proj][idx_];
      vECAL_tracksDz_[proj][idx_] = vECAL_tracksDz_[proj][idx_] / vECAL_tracks_[proj][idx_];
    }
  }

} // fillTracksAtECALstitched()
