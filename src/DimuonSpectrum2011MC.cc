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
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonMETCorrectionData.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
// #include "DataFormats/MuonReco/interface/MuonSelectors.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for trigger information 
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

//for generator information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

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
        bool eta21pt1510(const reco::Muon&, const reco::Muon&);
        bool istight(const reco::Muon&, const math::XYZPoint);
        bool isolation(const reco::Muon&);

// ----------member data ---------------------------

// declare Root histograms
// for a description of their content see below

TH1D *h66;
TH1D *h661;

TH1D *h10;

TH1D *h7;

// triggerExpression::Data triggerCache;
// std::unique_ptr<triggerExpression::Evaluator> triggerSelector;

};

//
// constants, enums and typedefs
//
  const double sqmums = (0.105658)*(0.105658); // square of muon mass (in GeV^2/c^4)

//
// static data member definitions
//

//
// constructors and destructor
//

DimuonSpectrum2011MC::DimuonSpectrum2011MC(const edm::ParameterSet& iConfig){

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


// dimuon mass spectrum up to 120 GeV after impose bound
h66 = fs->make<TH1D>("GM_mass_tight", "GTM mass ", 70, 10., 150.);
h66->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h66->GetYaxis()->SetTitle("Number of Events");

// dimuon mass spectrum up to 120 GeV after impose IP bound, 
h661 = fs->make<TH1D>("GM_mass_tight_iso", "GTM mass Iso", 70, 10., 150.);
h661->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h661->GetYaxis()->SetTitle("Number of Events");

// cut flow for the analysis of xsec_Zmumu
h7 = fs->make<TH1D>("Cut_Flow", "Cut Flow", 12, 0, 12);
h7->GetYaxis()->SetTitle("Number of Events");

}


DimuonSpectrum2011MC::~DimuonSpectrum2011MC() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// member functions

bool DimuonSpectrum2011MC::eta21pt1510 (const reco::Muon& m1, const reco::Muon& m2){

  if ((fabs(m1.eta()) < 2.1 && fabs(m2.eta()) < 2.1)
      && (m1.pt() > 10. && m2.pt() > 10.)
      && (m1.pt() > 15. || m2.pt() > 15.)){
    return true;
  } // baseline acceptance in 10.1103/PhysRevD.100.015021
  return false;
}


bool DimuonSpectrum2011MC::istight (const reco::Muon& muon, math::XYZPoint point){
  // Global muon with additional muon quality reqirements.
  // Starting from 50X release this set of selection is into an omni-comprehensive selector 
  // in DataFormats/MuonReco/interface/MuonSelectors.h
  // See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Tight_Muon_selection
  if (muon.isGlobalMuon()){
    if ( muon.globalTrack()->normalizedChi2() < 10. 
      && muon.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 
      && muon.numberOfMatchedStations() > 1 
      && fabs(muon.innerTrack()->dxy(point)) < 0.2 
      && fabs(muon.innerTrack()->dz(point)) < 1.0 
      && muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
      && muon.innerTrack()->hitPattern().numberOfValidTrackerHits() > 10 ){
        return true;
    }
  } 
  return false;
}


bool DimuonSpectrum2011MC::isolation (const reco::Muon& muon){
  if (muon.isIsolationValid()) {
    double iso=(muon.isolationR03().hadEt + muon.isolationR03().emEt + muon.isolationR03().sumPt)/muon.pt();
    if (iso < 0.15) return true;
  }
  return false;
}


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

  // INFO: Primary Vertex
  Handle<reco::VertexCollection> primvtxHandle;
  iEvent.getByLabel("offlinePrimaryVertices", primvtxHandle);
  reco::VertexCollection primvtx;
  if (primvtxHandle.isValid()) {
      primvtx = *primvtxHandle;
  } 
  else{ LogInfo("Demo")<< "No primary vertex available from EventSetup \n"; return; }

  math::XYZPoint point(primvtx[0].position());

  // INFO: GenParticle
  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);

  // WHAT: Fill histogram of the number of Muon-Tracks in the current Event.
  // WHY:  for monitoring purposes
  h10->Fill(muons->size());

  for(size_t i = 0; i < genParticles->size(); ++i) {
    const reco::GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    int st = p.status();  
    // const reco::Candidate * mom = p.mother();
    // double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
    // double vx = p.vx(), vy = p.vy(), vz = p.vz();
    // int charge = p.charge();
    // int n = p.numberOfDaughters();
    cout<<"id: "<<id<<" status: "<<st<<endl;
    // cout<<"id: "<<id<<" status: "<<st<<" mom: "<<mom.pdgId()
  }


  if (muons->size() >= 2) {

    h7->Fill(0);
  
    bool bsac = false;
    bool tight= false; 
    bool opps = false;
    bool zreg = false;
    bool pt20 = false;
    bool iso  = false;
//------------------analysing Muons (muons-TrackCollection)----------//


  // WHAT: Declare variables used later
    double s1, s2, s;

  // WHAT: Loop over all the Muons of current Event
  // WHY:  to select good candidates to be used in invariant mass calculation
    for (reco::MuonCollection::const_iterator it = muons->begin();
      it != muons->end(); it++) {

      // NTS: Stores iterator for current globalMuon-Track and advances it by one.
      //      In other words, the needed preparation to be able to compare all the
      //      other globalMuon-Tracks after
      //      the current one to the current globalMuon-Track with iterator it.
      reco::MuonCollection::const_iterator i = it;
      i++;

      // loop over 2nd muon candidate
      for (; i != muons->end(); i++) {

//-------------------------Calculate invariant mass-----------------------------//
// WHAT: Calculate invariant mass of globalMuon-Tracks under comparison
// (Iterators "it" and "i")
// WHY:  in order to fill the mass histogram
        s1 = sqrt(((it->p())*(it->p()) + sqmums) * ((i->p())*(i->p()) + sqmums));
        s2 = it->px()*i->px() + it->py()*i->py() + it->pz()*i->pz();
        s = sqrt(2.0 * (sqmums + (s1 - s2)));

//--------------------determine quality cuts----------------------//

// WHAT: If current globalMuon-Track satisfies quality-cut-criteria, it is
//       compared to other globalMuon-Tracks that come after this current one.
//       (succeeding globalMuon-Tracks that are in the muons-MuonCollection)
//       need at least two candidates to calculate dimuon mass

        if (eta21pt1510(*it,*i)) { bsac = true;
    
          if (istight(*it,point) && istight(*i,point)) { tight = true;

// WHAT: Compare electric charges of the current two globalMuon-Tracks
//       (Iterators "it" and "i")
// WHY:  Need to find out if the charges of the current two globalMuons-Tracks
//       are like or unlike charge, since the decaying parents are neutral
            if (it->charge() == -(i->charge()) ){ opps = true;
              // unlike charges

// WHAT: Store the invariant mass of two muons with unlike sign charges in
//       linear scale
// WHY:  in order to see the various mass peaks on linear scale
// WHAT: Store the invariant mass of two muons with unlike charges
// WHY: Reproduce the "Invariant mass spectrum of dimuons in events"-plot
//      from MUO-10-004
              double pt = sqrt( pow(it->px()+i->px(), 2.0) + pow(it->py()+i->py(), 2.0) );
              if (pt<s){
                h66->Fill(s);
                if (isolation(*it) && isolation(*i)) {
                  h661->Fill(s);
                }
              }
            } // end of unlike charge if

            if (s >= 60. && s <= 120.) {
              zreg = true;
              if (it->pt()>20. && i->pt()>20.) {
                pt20 = true;
                if (isolation(*it) && isolation(*i)) {
                  iso = true;
                }
              }
            }
          } // end of if(istight)
        } // end of if(eta21pt15pt10)
      } //end of for(;i!=muons....)
    } //end of reco ::MuonCollection loop
    if (bsac == true){
      h7->Fill(1);
    }
    if (tight == true){
      h7->Fill(2);
      if (opps == true) h7->Fill(3);
      else {h7->Fill(8);}
    }
    if (zreg == true){
      if (opps == true) {h7->Fill(4);}
      else {h7->Fill(9);}
    }
    if (pt20 == true){
      if (opps == true) {h7->Fill(5);}
      else {h7->Fill(10);}
    }
    if (iso == true){
      if (opps == true) {h7->Fill(6);}
      else {h7->Fill(11);}
    }
  } //end of if (size() >=2 )
} //DimuonSpectrum2011MC: analyze ends


// ------------ method called once each job just before starting event loop  ------------
void DimuonSpectrum2011MC::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DimuonSpectrum2011MC::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(DimuonSpectrum2011MC);   