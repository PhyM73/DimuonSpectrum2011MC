// -*- C++ -*-
//
// Original Package:    DimuonSpectrum2011
// Original Class:      DimuonSpectrum2011
//
/**\class DimuonSpectrum2011MC DimuonSpectrum2011MC.cc Demo/DimuonSpectrum2011/src/DimuonSpectrum2011MC.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:
//         Created:  Mon May  4 15:24:13 CEST 2015
//         Finalized: February 24, 2016  by   A. Geiser
//                    with contributions from I. Dutta,
//                                            H. Hirvonsalo
//                                            B. Sheeran
// $Id$
// ..
//
// ***************************************************************************
// version of DEMO setup provided by CMS open data access team               *
// expanded/upgraded to contain a pedagocigal analysis example for the       *
// dimuon mass spectrum (MUO-10-004)                                         *
//                                                                           *
// Note that the published spectrum is reproduced approximately, but not     *
// exactly, since the data sets only partially overlap, and, for reasons of  *
// simplicity, there is no trigger selection beyond the one implicit in the  *
// Mu data set, only global muons are used, and only part of the muon        *
// quality cuts are applied                                                  *
// ***************************************************************************

// system include files
#include <memory>

// user include files, general
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//------ EXTRA HEADER FILES--------------------//
#include "math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

// for Root histogramming
#include "TH1.h"

// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

// for Moun information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
// #include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// class declaration
//

class DimuonSpectrum2011MC: public edm::EDAnalyzer {
public:
        explicit DimuonSpectrum2011MC(const edm::ParameterSet&);
        ~DimuonSpectrum2011MC();

private:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();
        bool providesGoodLumisection(const edm::Event& iEvent);
        bool htl138active(int);
        bool eta21pt1510(const reco::MuonCollection::const_iterator, 
                         const reco::MuonCollection::const_iterator, double);
        bool istight(const reco::MuonCollection::const_iterator,const math::XYZPoint);
        bool isolation(const reco::MuonCollection::const_iterator, 
                       const reco::MuonCollection::const_iterator);

// ----------member data ---------------------------

// declare Root histograms
// for a description of their content see below
TH1D *h1;
TH1D *h2;
TH1D *h3;
TH1D *h4;

TH1D *h5;
TH1D *h6;
TH1D *h66;
TH1D *h661;
TH1D *h662;

TH1D *h10;

TH1D *h53;
TH1D *h54;
TH1D *h55;

TH1D *h60;
TH1D *h61;

TH1D *h7;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

DimuonSpectrum2011MC::DimuonSpectrum2011MC(const edm::ParameterSet& iConfig) {

// ***************************************************************************
// This is the main analysis routine
// The goal is to approximately reproduce the dimuon mass spectrum from 
// MUO-10-004 and reproduce the validation test in 10.1103/PhysRevD.100.015021
// ***************************************************************************

//now do what ever initialization is needed
edm::Service<TFileService> fs;

// ************************************
// book histograms and set axis labels
// (called once for initialization)
// ************************************

// monitoring histograms for muons,
// intended for muons from Mu sample

// muon multiplicity
h10 = fs->make<TH1D>("Mmultiplicity", "Mmultiplicity", 8, 0, 8);
h10->GetXaxis()->SetTitle("Number of Muons");
h10->GetYaxis()->SetTitle("Number of Events");

// muon momentum
h1 = fs->make<TH1D>("GMmomentum", "GM_Momentum", 240, 0., 120.);
h1->GetXaxis()->SetTitle("Global Muon Momentum (in GeV/c)");
h1->GetYaxis()->SetTitle("Number of Events");

// muon Transverse_momentum
h2 = fs->make<TH1D>("GM_Transverse_momentum", "TransverseMomentum", 240, 0., 120.);
h2->GetXaxis()->SetTitle("Transverse Momentum of global muons (in GeV/c)");
h2->GetYaxis()->SetTitle("Number of Events");

// muon pseudorapity
h3 = fs->make<TH1D>("GM_eta", "GM_Eta", 140, -3.5, 3.5);
h3->GetXaxis()->SetTitle("Eta of global muons (in radians)");
h3->GetYaxis()->SetTitle("Number of Events");

// muon azimuth angle
h4 = fs->make<TH1D>("GM_phi", "GM_phi", 314, -3.15, 3.15);
h4->GetXaxis()->SetTitle("Phi");
h4->GetYaxis()->SetTitle("Number of Events");

// dimuon mass spectrum, 8 to 12 GeV (upsilon)
h5 = fs->make<TH1D>("GMmass_up" , "GMmass_up" , 40, 8. , 12. );
h5->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h5->GetYaxis()->SetTitle("Number of Events");

// dimuon mass spectrum up to 120 GeV (high mass range: upsilon, Z)
h6 = fs->make<TH1D>("GMmass" , "GMmass" , 150, 0. , 150. );
h6->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h6->GetYaxis()->SetTitle("Number of Events");

// muon track chi2
h53 = fs->make<TH1D>("GM_chi2", "GM_Chi2", 300, 0, 150);
h53->GetXaxis()->SetTitle("Chi2 values");
h53->GetYaxis()->SetTitle("Number of Events");

// muon track number of degrees of freedom
h54 = fs->make<TH1D>("GM_ndof", "GM_ndof", 100, 0, 100);
h54->GetXaxis()->SetTitle("Ndof values");
h54->GetYaxis()->SetTitle("Number of Events");

// muon track chi2 normalized to number of degrees of freedom
h55 = fs->make<TH1D>("GM_normalizedchi2", "GM_normalizedChi2", 200, 0, 20);
h55->GetXaxis()->SetTitle("NormalizedChi2 values");
h55->GetYaxis()->SetTitle("Number of Events");

// muon track, number of valid hits
h60 = fs->make<TH1D>("GM_validhits", "GM_ValidHits", 100, 0., 100);
h60->GetXaxis()->SetTitle("Number of valid hits");
h60->GetYaxis()->SetTitle("Number of Events");

// muon track, number of pixel hits
h61 = fs->make<TH1D>("GM_pixelhits", "GM_pixelhits", 14, 0., 14);
h61->GetXaxis()->SetTitle("Munber of pixel hits");
h61->GetYaxis()->SetTitle("Number of Events");


// dimuon mass spectrum up to 120 GeV after impose bound
h66 = fs->make<TH1D>("GM_mass_cut", "GM mass Cut", 70, 10., 150.);
h66->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h66->GetYaxis()->SetTitle("Number of Events");

// dimuon mass spectrum up to 120 GeV after impose IP bound, 
h661 = fs->make<TH1D>("GM_mass_cut_IP", "GM mass Cut IP", 70, 10., 150.);
h661->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h661->GetYaxis()->SetTitle("Number of Events");

// dimuon mass spectrum up to 120 GeV after impose IP&ISo bound
h662 = fs->make<TH1D>("GM_mass_cut_IP_IS", "GM mass Cut IP IS", 70, 10., 150.);
h662->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h662->GetYaxis()->SetTitle("Number of Events");

// cut flow for the analysis of xsec_Zmumu
h7 = fs->make<TH1D>("Mu13Mu8", "Mu13Mu8", 12, 0, 12);
h7->GetYaxis()->SetTitle("Number of Events");

}


DimuonSpectrum2011MC::~DimuonSpectrum2011MC() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

// member functions

// ------------ method called for each event  ------------
void DimuonSpectrum2011MC::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

// **********************************************
// here each relevant event will get analyzed
// **********************************************

using namespace edm;
using namespace reco;
using namespace std;


#ifdef THIS_IS_AN_EVENT_EXAMPLE
        Handle<ExampleData> pIn;
        iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
        ESHandle<SetupData> pSetup;
        iSetup.get<SetupRecord>().get(pSetup);
#endif


// Event is to be analyzed
  // LogInfo("Demo")
  // << "Starting to analyze \n"
  // << "Event number: " << (iEvent.id()).event()
  // << ", Run number: " << iEvent.run()
  // << ", Lumisection: " << iEvent.luminosityBlock();

//------------------Load (relevant) Event information------------------------//
// INFO: Getting Data From an Event
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter4#GetData
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEDMGetDataFromEvent#get_ByLabel
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAodDataTable
// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideRecoDataTable


// INFO: Muons
// NB: note that when using keyword "Muons" getByLabel-function returns reco::MuonCollection
  Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons", muons);

  Handle<reco::VertexCollection> primvtxHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primvtxHandle);
  reco::VertexCollection primvtx;
  if (primvtxHandle.isValid()) {
      primvtx = *primvtxHandle;
  } 
  else{
      LogInfo("Demo")<< "No primary vertex available from EventSetup \n";
  }
  math::XYZPoint point(primvtx[0].position());

  if (htl138active(iEvent.run())){
    h7->Fill(0);
  
    bool bsac = false;
    bool tight = false; 
    // bool opps = false;
    bool zreg = false;
    bool pt20 = false;
    bool iso = false;
//------------------analysing Muons (muons-TrackCollection)----------//

// WHAT: Fill histogram of the number of Muon-Tracks in the current Event.
// WHY:  for monitoring purposes
    h10->Fill(muons->size());

// WHAT: Declare variables used later
    double sqm1, s1, s2, s;

// WHAT: Set square of muon mass
// WHY:  Needed in later calculations
    sqm1 = (0.105658) * (0.105658);

// WHAT: Loop over all the Muons of current Event
// WHY:  to select good candidates to be used in invariant mass calculation
    for (reco::MuonCollection::const_iterator it = muons->begin();
      it != muons->end(); it++) {
//       if (it->isGlobalMuon() && (it->globalTrack()).isNonnull()){
// // WHAT: Fill histograms for the following attributes from the current Muon-Track:
// // - p (momentum vector magnitude)
// // - pt (track transverse momentum)
// // - eta (pseudorapidity of momentum vector)
// // - chi-square
// // - ndof (number of degrees of freedom of the fit)
// // - normalizedChi2 (normalized chi-square == chi-squared divided by ndof
// //                   OR chi-squared * 1e6 if ndof is zero)
//         h1->Fill(it->p());
//         h2->Fill(it->pt());
//         h3->Fill(it->eta());
//         h4->Fill(it->phi());
//         h53->Fill(it->globalTrack()->chi2());
//         h54->Fill(it->globalTrack()->ndof());
//         h55->Fill(it->globalTrack()->normalizedChi2());

// // WHAT: Get HitPattern-object for Track of current Muon
// // WHY:  in order to count the number of hits on the track
//         const reco::HitPattern& p = it->globalTrack()->hitPattern();

// // WHAT: Fill number of ValidHits and PixelHits in current globalMuon-Track
// //       into histogram
// // WHY:  to check distribution before cuts
//         h60->Fill(p.numberOfValidHits());
//         h61->Fill(p.numberOfValidPixelHits());
//       }

  
//-----------------prepare variables to determine quality cuts---------------//

// loop over globalMuon-Tracks satisfying quality cuts //


// WHAT: If current globalMuon-Track satisfies quality-cut-criteria, it is
//       compared to other globalMuon-Tracks that come after this current one.
//       (succeeding globalMuon-Tracks that are in the muons-MuonCollection)
//       need at least two candidates to calculate dimuon mass
      if (muons->size() >= 2){
      //     && p.numberOfValidHits() >= 12
      //     && p.numberOfValidPixelHits() >= 2
      //     && it->globalTrack()->normalizedChi2() < 10.0) {
// NTS: Stores iterator for current globalMuon-Track and advances it by one.
//      In other words, the needed preparation to be able to compare all the
//      other globalMuon-Tracks after
//      the current one to the current globalMuon-Track with iterator it.
        reco::MuonCollection::const_iterator i = it;
        i++;

// loop over 2nd muon candidate
        for (; i != muons->end(); i++) {
          s1 = sqrt(((it->p())*(it->p()) + sqm1) * ((i->p())*(i->p()) + sqm1));
          s2 = it->px()*i->px() + it->py()*i->py() + it->pz()*i->pz();
          s = sqrt(2.0 * (sqm1 + (s1 - s2)));
          
          if (eta21pt1510(it,i,s)){
          bsac = true;
    
            if (istight(it,point) && istight(i,point)){
              tight = true;

// WHAT: Compare electric charges of the current two globalMuon-Tracks
//       (Iterators "it" and "i")
// WHY:  Need to find out if the charges of the current two globalMuons-Tracks
//      are like or unlike charge, since the decaying parents are neutral
              if (it->charge() == -(i->charge()) ){// unlike charges
// and cut on quality of 2nd muon candidate

              // os = true;

//----------Calculate invariant mass-----------------//
// WHAT: Calculate invariant mass of globalMuon-Tracks under comparison
// (Iterators "it" and "i")
// WHY:  in order to fill the mass histogram

// WHAT: Store the invariant mass of two muons with unlike sign charges in
//       linear scale
// WHY:  in order to see the various mass peaks on linear scale
                h66->Fill(s);
                if (isolation(it,i)) {
                  h661->Fill(s);
                }
              } // end of unlike charge if

              if (s >= 60. && s <= 120.) {
                zreg = true;
                if (it->pt()>20. && i->pt()>20.) {
                  pt20 = true;
                  if (isolation(it,i)) {
                    iso = true;
                  }
                }
              }
// WHAT: Store the invariant mass of two muons with unlike charges in log scale
// WHY: Reproduce the "Invariant mass spectrum of dimuons in events"-plot
//      from MUO-10-004
           
            } // end of if(istight)
          } // end of if(eta21pt15pt10)
        } //end of for(;i!=muons....)
      } //end of if(muons->size >=2 )
    } //end of reco ::MuonCollection loop
  if (bsac == true){
    h7->Fill(1);
  }
  if (tight == true){
    h7->Fill(2);
  }
  if (zreg == true){
    h7->Fill(3);
  }
  if (pt20 == true){
    h7->Fill(4);
  }
  if (iso == true){
    h7->Fill(5);
  }
  // if (accept == true){
  //   h73->Fill(0);
  // }
  } // end of HLT mu13mu8 active
} //DimuonSpectrum2011MC: analyze ends


// ------------ method called once each job just before starting event loop  ------------
void DimuonSpectrum2011MC::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DimuonSpectrum2011MC::endJob() {
}

bool DimuonSpectrum2011MC::htl138active (int run){
  if ((run >= 165088 && run <= 167043) ||
      (run >= 167078 && run <= 167913) ||
      (run >= 170249 && run <= 173198) ||
      (run >= 173236 && run <= 178380) ||
      (run >= 178420 && run <= 179889) ||
      (run >= 179959 && run <= 180252)){ // http://opendata.cern.ch/record/2698 
    return true;
  }
  return false;
}

bool DimuonSpectrum2011MC::eta21pt1510 (const reco::MuonCollection::const_iterator m1, 
       const reco::MuonCollection::const_iterator m2, double s){
  double pt = sqrt( pow(m1->px()+m2->px(), 2.0) + pow(m1->py()+m2->py(), 2.0) );

  if ((fabs(m1->eta()) < 2.1 && fabs(m2->eta()) < 2.1)
      && (m1->pt() > 10. && m2->pt() > 10.)
      && (m1->pt() > 15. || m2->pt() > 15.)
      && (pt < s)){
    return true;
  } // import bounds and selections in 10.1103/PhysRevD.100.015021
  return false;
}

bool DimuonSpectrum2011MC::istight (const reco::MuonCollection::const_iterator muon, math::XYZPoint point){
  if (muon->isGlobalMuon()){
    if ( muon->globalTrack()->normalizedChi2() < 10. 
      && muon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0 
      && muon->numberOfMatchedStations() > 1 
      && fabs(muon->innerTrack()->dxy(point)) < 0.2 
      && fabs(muon->innerTrack()->dz(point)) < 1.0 
      && muon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0 
      && muon->innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 ){
        return true;
    }
  } 
  // IP requirement, the reconstructed muon tracks must intersect the primary vertex
  // within d < 2mm in the x–y plane and z < 10mm in the z direction.
  return false;
}


bool DimuonSpectrum2011MC::isolation (const reco::MuonCollection::const_iterator m1, 
       const reco::MuonCollection::const_iterator m2){
  if (m1->isIsolationValid() && m2->isIsolationValid()) {
    double iso1=(m1->isolationR03().hadEt+m1->isolationR03().emEt+m1->isolationR03().sumPt)/m1->pt();
    double iso2=(m2->isolationR03().hadEt+m2->isolationR03().emEt+m2->isolationR03().sumPt)/m2->pt();
    if (iso1 < 0.15 && iso2 < 0.15){
      return true;
    }
  }
  return false;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DimuonSpectrum2011MC);                    
