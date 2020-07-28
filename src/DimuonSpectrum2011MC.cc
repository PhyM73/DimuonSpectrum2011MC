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
//         Created:   May  4 15:24:13 CEST 2015
//         Finalized: February 24, 2016  by   A. Geiser
//                    with contributions from I. Dutta,
//                                            H. Hirvonsalo
//                                            B. Sheeran
//
// Author:
//         Created:   Apr 6,  2020
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
// #include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
// #include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
// #include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"

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
        bool isolation15(const reco::Muon&);
        double invmass(const reco::Candidate&, const reco::Candidate&);
        const reco::Candidate* daughter_afsr(const reco::Candidate* );
// ----------member data ---------------------------

// declare Root histograms
// for a description of their content see below

TH1D *h10;

TH1D *h6;

TH1D *h8;

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


// the numerator and the denominator of the acceptance for the analysis of xsec_Zmumu
h8 = fs->make<TH1D>("Accept", "Acceptance", 2, 0, 2);
h8->GetYaxis()->SetTitle("Number of Events");

}


DimuonSpectrum2011MC::~DimuonSpectrum2011MC() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}


// member functions

bool DimuonSpectrum2011MC::eta21pt1510 (const reco::Muon& m1, const reco::Muon& m2){

  if ((fabs(m1.eta()) < 2.1 && fabs(m2.eta()) < 2.1) && m1.pt() > 15. && m2.pt() > 10.){
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


bool DimuonSpectrum2011MC::isolation15 (const reco::Muon& muon){
  if (muon.isIsolationValid()) {
    double iso=(muon.isolationR03().hadEt + muon.isolationR03().emEt + muon.isolationR03().sumPt)/muon.pt();
    if (iso < 0.15) return true;
  }
  return false;
}


double DimuonSpectrum2011MC::invmass (const reco::Candidate& p1, const reco::Candidate& p2){
  double  s1 = sqrt(((p1.p())*(p1.p()) + sqmums) * ((p2.p())*(p2.p()) + sqmums));
  double  s2 = p1.px()*p2.px() + p1.py()*p2.py() + p1.pz()*p2.pz();
  double  s = sqrt(2.0 * (sqmums + (s1 - s2)));
  return s;
}


const reco::Candidate* DimuonSpectrum2011MC::daughter_afsr(const reco::Candidate* p ){
  for(size_t i = 0; i < p->numberOfDaughters();++i){
    if (p->daughter(i)->pdgId() == p->pdgId()){
      return daughter_fsr(p->daughter(i));
    }
  }
  return p;
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

  // Fill histogram of the number of Muon-Tracks in the current Event for monitoring purposes
  h10->Fill(muons->size());

  if (muons->size() >= 2) {

    //------------------analysing Muons (muons-TrackCollection)----------//

    // select the leading muon and subleading muon
    reco::MuonCollection::const_iterator it = muons->begin();
    reco::Muon muon1 = *it;
    it++;
    reco::Muon muon2 = *it;
    it++;
    if (muon2.pt() > muon1.pt()){
      reco::Muon m = muon2;
      muon2 = muon1;
      muon1 = m;
    }
    for (; it != muons->end(); it++) {
      // Loop over all the remain Muons (if any) of current Event
      if (it->pt()>muon1.pt()){
        muon2 = muon1;
        muon1 = *it;
      } else if (it->pt()>muon2.pt()){
        muon2 = *it;
      }
    }  

    //-------------------------Calculate invariant mass-----------------------------//
    double s1, s2, s;
    s1 = sqrt(((muon1.p())*(muon1.p()) + sqmums) * ((muon2.p())*(muon2.p()) + sqmums));
    s2 = muon1.px()*muon2.px() + muon1.py()*muon2.py() + muon1.pz()*muon2.pz();
    s = sqrt(2.0 * (sqmums + (s1 - s2)));

    //--------------------determine quality cuts----------------------//

    // If these Muon-Tracks satisfy the quality-cut-criteria, their invariant mass is collected
    if (eta21pt1510(muon1,muon2)) {
      if (istight(muon1,point) && istight(muon2,point)) {
        if (muon1.charge() == -(muon2.charge()) ){ 

          double pt = sqrt( pow(muon1.px()+muon2.px(), 2.0) + pow(muon1.py()+muon2.py(), 2.0) );
          if (pt<s && isolation15(muon1) && isolation15(muon2)) {
            h6->Fill(s);
          }
        } // end of unlike charge if
      } // end of if(istight)
    } // end of if(eta21pt15pt10)
  } //end of if (size() >=2 )

  //------------------------------evaluate Acceptance-----------------------------------//

  // INFO: GenParticle
  Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel("genParticles", genParticles);

  reco::GenParticle muonbeforeFSR1, muonbeforeFSR2;
  const reco::Candidate *muonafterFSR1, *muonafterFSR2;
  int count = 0;

  for(reco::GenParticleCollection::const_iterator itp = genParticles->begin();
      itp != genParticles->end() && itp->status() == 3 ; itp++) {

    if(abs(itp->pdgId()) == 13 && itp->mother()->pdgId() == 23){
      if (count == 0) {
        muonbeforeFSR1 = *itp;
        for(size_t i = 0; i < itp->numberOfDaughters();i++){
          if (itp->daughter(i)->pdgId()==itp->pdgId()) muonafterFSR1 = daughter_afsr(itp->daughter(i));
        }
        count++;
      }
      else { 
        muonbeforeFSR2 = *itp;
        for(size_t i = 0; i < itp->numberOfDaughters();i++){
          if (itp->daughter(i)->pdgId()==itp->pdgId()) muonafterFSR2 = daughter_afsr(itp->daughter(i));
        }
        count++;
      }    
    }
  }    

  if (count == 2) {
    double mass = invmass(muonbeforeFSR1, muonbeforeFSR2);
    if (mass > 60. && mass < 120.) h8->Fill(0); //the denominator of the acceptance
    if (muonafterFSR1->pt() > 20 && muonafterFSR2->pt() > 20
        && fabs(muonafterFSR1->eta()) < 2.1 && fabs(muonafterFSR2->eta()) < 2.1 ){
      double m = invmass(*muonafterFSR1, *muonafterFSR2);
      if (m > 60. && m < 120.) h8->Fill(1); //the nominator of the acceptance
    }
    
  }

} //DimuonSpectrum2011MC: analyze ends


// ------------ method called once each job just before starting event loop  ------------
void DimuonSpectrum2011MC::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DimuonSpectrum2011MC::endJob() {
}


//define this as a plug-in
DEFINE_FWK_MODULE(DimuonSpectrum2011MC);   