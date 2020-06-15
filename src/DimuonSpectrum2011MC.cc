// -*- C++ -*-
//
// Package:    DimuonSpectrum2011MC
// Class:      DimuonSpectrum2011MC
//
/**\class DimuonSpectrum2011MC DimuonSpectrum2011MC.cc Demo/DimuonSpectrum2011/src/DimuonSpectrum2011MC.cc

 Description: [one line class summary]

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:
//         Created:  Mon  May 4, 15:24:13 CEST 2015
//         Finalized: February 24, 2016  by   A. Geiser
//                    with contributions from I. Dutta,
//                                            H. Hirvonsalo
//                                            B. Sheeran

// Author:
//         Created:  Mon  Apr 6, 2020
//         Finalized: Jun 12, 2020  by   F.Q. Meng (still in development)
// $Id$
// ..
//

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
// #include "DataFormats/MuonReco/interface/MuonSelectors.h"

// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// for trigger information 
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

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
        bool isolation15(const reco::Muon&);
        bool search(const reco::Muon&, const math::XYZPoint);


// ----------member data ---------------------------

// declare Root histograms
// for a description of their content see below

TH1D *h10;

TH1D *h6;
TH1D *h66[6];

TH1D *h7;

// declare the trigger selector
triggerExpression::Data triggerCache;
std::unique_ptr<triggerExpression::Evaluator> triggerSelector;

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

DimuonSpectrum2011MC::DimuonSpectrum2011MC(const edm::ParameterSet& iConfig):
    triggerCache(triggerExpression::Data(edm::InputTag("TriggerResults","","HLT"), 
                 edm::InputTag(""), 1, false, false, false)) ,
    triggerSelector(triggerExpression::parse("HLT_Mu13_Mu8*")) {

// ***************************************************************************
// This is the main analysis routine
// The goal is to reproduce the validation test and the resonance search
// in 10.1103/PhysRevD.100.015021  
// ***************************************************************************

//now do what ever initialization is needed
edm::Service<TFileService> fs;

// ************************************
// book histograms and set axis labels
// (called once for initialization)
// ************************************

// monitoring histograms for muons, intended for muons from Mu sample
// muon multiplicity
h10 = fs->make<TH1D>("Mmultiplicity", "Mmultiplicity", 8, 0, 8);
h10->GetXaxis()->SetTitle("Number of Muons");
h10->GetYaxis()->SetTitle("Number of Events");


// dimuon mass spectrum up to 150 GeV for tight muons after impose Isolaiton requires 
// Perform the comparison between the CMS2011a data set and Monte Carlo Samples
h6 = fs->make<TH1D>("GM_mass_tight_iso", "GTM mass Iso", 70, 10., 150.);
h6->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
h6->GetYaxis()->SetTitle("Number of Events");


const char* name[6] = {"GM_mass_iso_0", "GM_mass_iso_25", "GM_mass_iso_60", 
                       "GM_mass_pro_0", "GM_mass_pro_25", "GM_mass_pro_60"};
// dimuon mass spectrum in 2 GeV bins. 
for (int i=0; i<6; i++){
  h66[i] = fs->make<TH1D>(name[i], "GM mass", 37, 10.5, 84.5);
  h66[i]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  h66[i]->GetYaxis()->SetTitle("Number of Events");
}


// cut flow for the analysis of xsec_Zmumu and the search work
h7 = fs->make<TH1D>("Cut_Flow", "Cut Flow", 14, 0, 14);
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
  // Perform the Tight Muon selection explicitly
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


bool DimuonSpectrum2011MC::isolation15 (const reco::Muon& m){
  if (m.isIsolationValid()) {
    double iso=(m.isolationR03().hadEt+m.isolationR03().emEt+m.isolationR03().sumPt)/m.pt();
    if (iso < 0.15){ return true; }
  }
  return false;
}


bool DimuonSpectrum2011MC::search (const reco::Muon& muon, math::XYZPoint point){
  // the tight IP cuts for resonance search
  if (fabs(muon.innerTrack()->dxy(point)) < 0.025 
      && fabs(muon.innerTrack()->dz(point)) < 0.2) return true; 
  return false;
}

// ------------ method called for each event  ------------//
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

  // WHAT: Fill histogram of the number of Muon-Tracks in the current Event.
  // WHY:  for monitoring purposes
  h10->Fill(muons->size());

  // INFO: Use the trigger result as a evnet selector, see the link below: 
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TriggerResultsFilter#Use_as_a_Selector_AN1
  // Pass the Event and EventSetup to the cache object
  if (triggerSelector and triggerCache.setEvent(iEvent, iSetup)){
    // if the L1 or HLT configurations have changed, (re)initialize the filters 
    // (including during the first event)
    if (triggerCache.configurationUpdated())
      triggerSelector ->init(triggerCache);
  }

  bool trigger_result = (*triggerSelector)(triggerCache);
  if (trigger_result && (muons->size() >= 2)) {

    h7->Fill(0);
  
    // sign for cut flow 
    bool bsac = false;  // baseline acceptance
    bool tight= false;  // tight muon
    bool opps = false;  // opposite sign
    bool zreg = false;  // z-mass region
    bool pt20 = false;  // pt > 20 Gev/c
    bool iso  = false;  // Iso < 0.15
    bool sea  = false;  // search region

    //------------------analysing Muons (muons-TrackCollection)--------------------//

    // WHAT: Declare variables used later
    double s1, s2, s;

    // Loop over all the Muons of current Event
    for (reco::MuonCollection::const_iterator it = muons->begin();
      it != muons->end(); it++) {

      reco::MuonCollection::const_iterator i = it;
      i++;

      // Loop over 2nd muon candidate
      for (; i != muons->end(); i++) {

        //-------------------------Calculate invariant mass----------------------//
        // WHAT: Calculate invariant mass of globalMuon-Tracks under comparison
        // (Iterators "it" and "i")
        // WHY:  in order to fill the mass histogram
        s1 = sqrt(((it->p())*(it->p()) + sqmums) * ((i->p())*(i->p()) + sqmums));
        s2 = it->px()*i->px() + it->py()*i->py() + it->pz()*i->pz();
        s = sqrt(2.0 * (sqmums + (s1 - s2)));

        //--------------------determine quality cuts-----------------------------//

        // WHAT: If these Muon-Tracks satisfy the quality-cut-criteria, the cut flow 
        //       is recorded and their invariant mass is collected.

        if (eta21pt1510(*it,*i)) { bsac = true;
    
          if (istight(*it,point) && istight(*i,point)) { tight = true;

            // WHAT: Compare electric charges of the current two globalMuon-Tracks
            //       (Iterators "it" and "i")
            if (it->charge() == -(i->charge()) ){ opps = true;

              double pt = sqrt( pow(it->px()+i->px(), 2.0) + pow(it->py()+i->py(), 2.0) );
              if (pt<s && isolation15(*it) && isolation15(*i)){
                h6->Fill(s);
                // WHAT: Store the invariant mass of two muons with unlike sign charges
              }

              if (search(*it,point) && search(*i,point) && s >= 11. && s <= 83. ) { 
                sea = true;
                if (isolation15(*it,point) && isolation15(*i,point)){ 
                  // isolated sample
                  h66[0]->Fill(s);
                  if (pt>25.) h66[1]->Fill(s);
                  if (pt>60.) h66[2]->Fill(s);
                }
                if (fabs(it->innerTrack()->dxy(point)) < 0.01 && fabs(i->innerTrack()->dxy(point)) < 0.01){
                  // prompt sample
                  h66[3]->Fill(s);
                  if (pt>25.) h66[4]->Fill(s);
                  if (pt>60.) h66[5]->Fill(s);
                }
              }

            } // end of unlike charge if

            if (s >= 60. && s <= 120.) { zreg = true;
              if (it->pt()>20. && i->pt()>20.) { pt20 = true;
                if (isolation15(*it) && isolation15(*i)) { iso = true; }
              }
            }

          } // end of if(istight)
        } // end of if(eta21pt15pt10)
      } //end of for(;i!=muons->end();...)
    } //end of reco::MuonCollection loop
    // Fill the histo of the cut flow
    if (bsac == true) h7->Fill(1);
    if (tight == true){ h7->Fill(2);
      if (opps == true) h7->Fill(3); else h7->Fill(8); 
      }
    if (zreg == true){ if (opps == true) h7->Fill(4); else h7->Fill(9);  }
    if (pt20 == true){ if (opps == true) h7->Fill(5); else h7->Fill(10); }
    if (iso == true) { if (opps == true) h7->Fill(6); else h7->Fill(11); }
    if (sea == true) { h7->Fill(13); }
  } // end of trigger_result
} //DimuonSpectrum2011MC: analyze ends


// ------------ method called once each job just before starting event loop  ------------
void DimuonSpectrum2011MC::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DimuonSpectrum2011MC::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(DimuonSpectrum2011MC);   