// -*- C++ -*-
//
// Package:    MuonAnalyzer
// Class:      MuonAnalyzer
//
/**\class MuonAnalyzer MuonAnalyzer.cc MuonAnalyzer/MuonAnalyzer/src/MuonAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Luca MARTINI
//         Created:  Thu Jan 10 16:14:20 CET 2013
// $Id: MuonAnalyzer.cc,v 1.6 2013/02/18 16:13:37 lmartini Exp $
//
//


// system include files
#include <memory>

#include <math.h>
#include <iterator>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <DataFormats/TrackReco/interface/TrackBase.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;

struct SimpleTrack
{
  SimpleTrack(const double & inP, const double & inPt, const double & inEta, const double & inPhi,
              const int inNdof, const double & inDoca,
              const double & inDocaSignificance,
              const int inVertexId, const bool inHighPurity,
              const double & inDeltaRMu1, const double & inDeltaRMu2, const double inDeltaRCand) :
      p(inP), pt(inPt), eta(inEta), phi(inPhi),
      ndof(inNdof), doca(inDoca),
      docaSignificance(inDocaSignificance),
      vertexId(inVertexId), highPurity(inHighPurity),
      deltaRMu1(inDeltaRMu1), deltaRMu2(inDeltaRMu2), deltaRCand(inDeltaRCand)
  {}

  bool operator<(const SimpleTrack & compTrack) const {
    return( doca < compTrack.doca );
  }

  double p;
  double pt;
  double eta;
  double phi;
  int ndof;
  double doca;
  double docaSignificance;
  int vertexId;
  bool highPurity;
  double deltaRMu1;
  double deltaRMu2;
  double deltaRCand;
};

//
// class declaration
//

class MuonAnalyzer : public edm::EDAnalyzer {
public:
  explicit MuonAnalyzer(const edm::ParameterSet&);
  ~MuonAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  void initialize_tree_vars();
  void fillGen(const edm::Event& iEvent);
  void fillTriggerHLTpaths(const edm::Event& iEvent);
  void fillTriggerL3(const edm::Event& iEvent, const edm::EventSetup &iSetup);
  void fillRAWTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void fillReco(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  int findVertexId(const reco::Vertex &theOriginalPV,
                   const std::vector<reco::Vertex> &priVtxs,
                   const reco::Track & track, const double & maxDeltaR );

  vector<int> is_MC_matched(const edm::Event& iEvent, const Muon * mu);
  bool is_trigger_matched(const edm::Event& iEvent, const TLorentzVector mu);

  bool IsTightMuon(const reco::Muon& muon);

  // ----------member data ---------------------------
  string outputname_;
  Bool_t doMC_;
  Bool_t doTrigger_;
  Bool_t doRAWTrigger_;
  Bool_t doReco_;
  vector< Int_t> pdgId_;
  vector <string> HLTPaths_;

  map <string, bool> mHLTPaths_;

  string HLTString_;

  InputTag onlineBS_;
  InputTag offlineBS_;
  InputTag PVs_;
  InputTag GenParticle_;

  Int_t ngen;
  Int_t nsv_hlt;
  Int_t nL3;
  Int_t nmuons;
  Int_t ndimuons;
  Int_t npv;

  TFile * output_;
  TTree * tree_;
  TClonesArray * gen_particle_4mom;
  TClonesArray * L3_particle_4mom;
  TClonesArray * reco_muon_4mom;
  TClonesArray * reco_dimuon_4mom;
  TClonesArray * hlt_vtx_4mom;

  Int_t charge[100];
  Int_t dicharge[100];
  Int_t pdgId[100];
  Int_t pdgIdMom[100];
  Int_t pdgIdGrandma[100];

  Int_t gen_pdgId[100];
  Int_t gen_pdgIdMom[100];
  Int_t gen_status[100];

  Int_t L3_particle_id[100];

  Float_t normChi2[100];
  Double_t vtxProb[100];
  Float_t cosAlpha[100];
  Float_t LxySignificance[100];
  Float_t Lxy[100];
  Float_t Lxyerr[100];

  bool isTightMuon[100];
  bool isTriggerMatched[100];
  Double_t dimuon_iso[100];

  Int_t mu1Index[100];
  Int_t mu2Index[100];


};

//
// constants, enums and typedefs
//
const Double_t muon_mass = 0.105658;
const Double_t MuMass2(muon_mass * muon_mass);
//
// static data member definitions
//

//
// constructors and destructor
//
MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed

  outputname_   = iConfig.getParameter<string> ("OutputFileName");
  doReco_       = iConfig.getParameter<bool> ("doReco");
  doMC_         = iConfig.getParameter<bool> ("doMC");
  doTrigger_    = iConfig.getParameter<bool> ("doTrigger");
  doRAWTrigger_ = iConfig.getParameter<bool> ("doRAWTrigger");
  pdgId_        = iConfig.getParameter <vector <int> > ("pdgId");
  HLTString_    = iConfig.getParameter<string> ("HLTString");
  HLTPaths_     = iConfig.getParameter <vector <string> > ("HLTPaths");
  onlineBS_     = iConfig.getParameter <InputTag> ("OnlineBeamSpot");
  offlineBS_    = iConfig.getParameter<InputTag>("OfflineBeamSpot");
  PVs_          = iConfig.getParameter <InputTag> ("PrimaryVertexCollection");
  GenParticle_  = iConfig.getParameter <InputTag> ("GenParticleCollection");

  output_ = new TFile (outputname_.c_str(), "RECREATE" );
  tree_   = new TTree("tree","tree");

  reco_muon_4mom    = new TClonesArray("TLorentzVector", 100);
  reco_dimuon_4mom  = new TClonesArray("TLorentzVector", 100);
  gen_particle_4mom = new TClonesArray("TLorentzVector", 100);
  L3_particle_4mom  = new TClonesArray("TLorentzVector", 100);
  hlt_vtx_4mom      = new TClonesArray("TLorentzVector", 100);

  for (unsigned int i = 0; i < HLTPaths_.size(); i++) {
    mHLTPaths_[HLTPaths_[i]] = false;
  }

}


MuonAnalyzer::~MuonAnalyzer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  delete tree_;
  output_->Close();
  delete output_;

}


//
// member functions
//

// ------------ method called for each event  ------------
void MuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  initialize_tree_vars();
  if (doMC_) fillGen(iEvent);
  if (doTrigger_) {
    fillTriggerHLTpaths(iEvent);
    fillTriggerL3(iEvent, iSetup);
  }
  if (doReco_) fillReco(iEvent, iSetup);
  if (doRAWTrigger_) fillRAWTrigger(iEvent, iSetup);

  //////////////////
  tree_->Fill(); ///
  //////////////////
}

void MuonAnalyzer::fillReco(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  /// vertex
  Handle<VertexCollection> priVtxs;
  iEvent.getByLabel(PVs_, priVtxs);
  npv = priVtxs->size();
  vector< Vertex> priVtxs_v;
  for( auto PVit = priVtxs->begin(); PVit != priVtxs->end(); ++PVit ) {
    const Vertex * v = &(*PVit);
    priVtxs_v.push_back(*v);
  }

  Handle<TrackCollection> tkColl;
  iEvent.getByLabel("generalTracks", tkColl);
  //get the transient track builder:
  ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);

  Handle<BeamSpot> theBeamSpot;
  iEvent.getByLabel(offlineBS_, theBeamSpot);
  BeamSpot bs = *theBeamSpot;

  string theMuonLabel = "muons";

  //	EDGetTokenT<MuonCollection> theMuonLabel;
  //	theMuonLabel = consumes<MuonCollection>(InputTag("muons"));

  Handle<MuonCollection> muons;
  //	iEvent.getByToken(theMuonLabel, muons);
  iEvent.getByLabel(theMuonLabel, muons);

  if (!muons.isValid()) {
    cout << "theMuonLabel is not valid" << endl;
    return;
  }


  // first muon
  MuonCollection::const_iterator muon_i;
  for (muon_i = muons->begin(); muon_i!=muons->end(); ++muon_i) { // loop over all muons
    if (nmuons >= 100) {
      cout << "array too small = 100" << endl;
      return;
    }
    const Muon * mu = &(*muon_i);
    if (!mu->isGlobalMuon() || !mu->isTrackerMuon()) continue;
    reco::TrackBase::TrackQuality trackQualityhighPur = reco::TrackBase::qualityByName("highPurity");
    if (!mu->bestTrack()->quality(trackQualityhighPur)) continue;

    const Track * muon_track = mu->bestTrack();
    TLorentzVector muon_4mom(0., 0., 0., 0.);
    if (muon_track->pt() < 2.5) continue;
    muon_4mom.SetPtEtaPhiM(muon_track->pt(), muon_track->eta(), muon_track->phi(), muon_mass);

    new((*reco_muon_4mom)[nmuons]) TLorentzVector(muon_4mom);
    charge[nmuons] = (*muon_i).charge();

    if (doMC_) {
      vector <int> family(3, -1);
      family = is_MC_matched(iEvent, mu);
      pdgId[nmuons] = family[0];
      pdgIdMom[nmuons] = family[1];
      pdgIdGrandma[nmuons] = family[2];
    }

    if (doTrigger_) {
      isTriggerMatched[nmuons] = is_trigger_matched(iEvent, muon_4mom);
    }

    // ID
    isTightMuon[nmuons] = IsTightMuon(*mu);
    nmuons++;

    // second muon and make jpsi
    MuonCollection::const_iterator muon_j;
    for (muon_j = muon_i+1; muon_j!=muons->end(); ++muon_j) { // loop over all muons
      const Muon * mu_j = &(*muon_j);
      if (!mu_j->isGlobalMuon() || !mu_j->isTrackerMuon()) continue;
      if (!mu_j->bestTrack()->quality(trackQualityhighPur)) continue;
      const Track * muon_track_j = mu_j->bestTrack();
      TLorentzVector muon_4mom_j(0., 0., 0., 0.);
      if (muon_track_j->pt() < 2.5) continue;
      muon_4mom_j.SetPtEtaPhiM(muon_track_j->pt(), muon_track_j->eta(), muon_track_j->phi(), muon_mass);

      TLorentzVector dimuon_4mom(0., 0., 0., 0.);
      dimuon_4mom = muon_4mom + muon_4mom_j;
      new((*reco_dimuon_4mom)[ndimuons]) TLorentzVector(dimuon_4mom);
      dicharge[ndimuons] = (*muon_i).charge() + (*muon_j).charge();

      unsigned int index1 = std::distance( muons->begin(), muon_i);
      unsigned int index2 = std::distance( muons->begin(), muon_j);
      mu1Index[ndimuons] = -1;
      mu2Index[ndimuons] = -1;
      if (index1 < muons->size()) mu1Index[ndimuons] = index1;
      if (index2 < muons->size()) mu2Index[ndimuons] = index2;

      // SV
      vector<TransientTrack> t_tks;
      t_tks.push_back(theTTBuilder->build(mu->bestTrack()));
      t_tks.push_back(theTTBuilder->build(mu_j->bestTrack()));
      KalmanVertexFitter kvf;
      TransientVertex tSV = kvf.vertex(t_tks);
      if (!tSV.isValid()) continue;
      Vertex SV = tSV;

      // associate PV, closest in extrapolated Z

      float minDz = 999999.;

      TwoTrackMinimumDistance ttmd;
      bool status = ttmd.calculate( GlobalTrajectoryParameters(GlobalPoint(SV.position().x(), SV.position().y(), SV.position().z()), GlobalVector(dimuon_4mom.Px(), dimuon_4mom.Py(), dimuon_4mom.Pz()), TrackCharge(dicharge[ndimuons]), &(*magneticField)), GlobalTrajectoryParameters(GlobalPoint(bs.position().x(), bs.position().y(), bs.position().z()), GlobalVector(bs.dxdz(), bs.dydz(), 1.), TrackCharge(0), &(*magneticField)) );
      float extrapZ = -9E20;
      if (status) extrapZ = ttmd.points().first.z();

      Vertex thePrimaryV;
      // Select the closest PV in z
      for( auto tPV = priVtxs->begin(); tPV != priVtxs->end(); ++tPV ) {
        float deltaZ = fabs(extrapZ - tPV->position().z()) ;
        if ( deltaZ < minDz ) {
          minDz = deltaZ;
          // thePrimaryV = Vertex(*tPV);
          thePrimaryV = *tPV;
        }
      }

 // count the number of tracks with pT > 900 MeV attached to the chosen vertex
      double sumPTPV = 0., sumPTNoVtx = 0.;
      double sumPTPV20 = 0., sumPTNoVtx20 = 0.;
      double sumPTkMu1 = 0., sumPTkMu2 = 0.;
      // double sumPTkMu1_20 = 0., sumPTkMu2_20 = 0.;
      // int countTksOfPV = 0, countTksOfNoVtx = 0;
      int Ntrk = 0, Ntrkhp = 0, Ntrk20 = 0, Ntrkhp20 = 0;
      int Ntrk1sigma = 0, Ntrk1sigmahp = 0, Ntrk1sigma20 = 0, Ntrk1sigmahp20 = 0;
      int Ntrk2sigma = 0, Ntrk2sigmahp = 0, Ntrk2sigma20 = 0, Ntrk2sigmahp20 = 0;
      int Ntrk3sigma = 0, Ntrk3sigmahp = 0, Ntrk3sigma20 = 0, Ntrk3sigmahp20 = 0;
      double minDca = 9999.;
      double sumNdofPV = 0., sumNdofNoVtx = 0.;
      // double PVndof = thePrimaryV.ndof();

      // TrajectoryStateClosestToPoint mu1TS = t_tks[0].impactPointTSCP();
      // TrajectoryStateClosestToPoint mu2TS = t_tks[1].impactPointTSCP();

      // Fill a collection of all tracks and:
      // - if they belong to the PV associated to the candidate set vertexId to 1
      // - if they belong to another PV set the vertexId to 2
      // - if they belong to no PV set the vertexId to 0
      double maxDeltaR = 1.e-5;
      std::vector<SimpleTrack> simpleTracks;
      for ( auto track = tkColl->begin(); track != tkColl->end(); ++track ) {
        if( track->extra().key() == mu->track().key() ) continue;
        if( track->extra().key() == mu_j->track().key() ) continue;

        // Set defaults
        double doca = 999.;
        double docaSig = 999.;
        int vertexId = 999;
        bool highPurity = track->quality(reco::TrackBase::highPurity);
        // Compute doca and doca significance
        TransientTrack tt = theTTBuilder->build(*track);
        // Looking in the code it seems absoluteImpactParameter3D uses only the position of the vertex

        pair<bool, Measurement1D> tkPVdist = IPTools::absoluteImpactParameter3D(tt, SV);
        if( tkPVdist.first ) {
          doca = tkPVdist.second.value();
          if(tkPVdist.second.error() != 0) docaSig = doca/tkPVdist.second.error();
        }
        vertexId = findVertexId( thePrimaryV, priVtxs_v, *track, maxDeltaR );

        // Muon isolation
        // --------------
        double deltaRMu1 = 1.;
        double deltaRMu2 = 1.;
        double deltaRCand = 1.;
        if( vertexId != 2 ) {
          deltaRCand = deltaR(track->eta(), track->phi(), dimuon_4mom.Eta(), dimuon_4mom.Phi());
          deltaRMu1 = deltaR(track->eta(),track->phi(), mu->eta(), mu->phi());
          deltaRMu2 = deltaR(track->eta(),track->phi(), mu_j->eta(), mu_j->phi());
        }
        simpleTracks.push_back(SimpleTrack(track->p(), track->pt(), track->eta(), track->phi(), track->ndof(),                                   doca, docaSig, vertexId, highPurity,
                                                 deltaRMu1, deltaRMu2, deltaRCand));
      }

      // Sort the simpleTracks. Smallest doca first.
      std::sort(simpleTracks.begin(), simpleTracks.end());


      // Count the number of tracks from the candidate PV or no PV and with pt > 0.5 GeV, doca < 300 microns and doca significance > 1.
      // Using only the first 20 tracks
      double maxDoca = 0.03;
      int tkCount = 0;
      for( auto tk = simpleTracks.begin(); tk != simpleTracks.end(); ++tk, ++tkCount ) {

        // From here on use only tracks from the candidate PV or no PV
        if( tk->vertexId == 2 ) continue;

        // if( tk->vertexId == 1 ) ++countTksOfPV;
        // else if( tk->vertexId == 0 ) ++countTksOfNoVtx;

        double deltaRMuIsoCut = 0.5;

        // Consistent
        if( tk->vertexId == 1 && tk->p > 0.5 ) {
          if( tk->deltaRMu1 < deltaRMuIsoCut ) sumPTkMu1 += tk->p;
          if( tk->deltaRMu2 < deltaRMuIsoCut ) sumPTkMu2 += tk->p;
        }
        if( tk->vertexId == 0 && tk->doca < 0.1 && tk->p > 0.5 && tkCount < 20 ) {
          if( tk->deltaRCand < deltaRMuIsoCut ) sumPTkMu1 += tk->p;
          if( tk->deltaRCand < deltaRMuIsoCut ) sumPTkMu2 += tk->p;
        }

        // Take the minimum doca
        if( minDca > tk->doca ) minDca = tk->doca;

        // Number of close tracks
        if( tk->pt > 0.5 && tk->doca < maxDoca ) {
          ++Ntrk;
          if( tk->highPurity ) ++Ntrkhp;
          if( tkCount < 20 ) {
            ++Ntrk20;
            if( tk->highPurity ) ++Ntrkhp20;
          }

          // Number of close tracks with doca > 1 sigma
          if( tk->docaSignificance > 1. ) {
            ++Ntrk1sigma;
            if( tk->highPurity ) ++Ntrk1sigmahp;
            if( tkCount < 20 ) {
              ++Ntrk1sigma20;
              if( tk->highPurity ) ++Ntrk1sigmahp20;
            }

            // Number of close tracks with doca > 2 sigma
            if( tk->docaSignificance > 2. ) {
              ++Ntrk2sigma;
              if( tk->highPurity ) ++Ntrk2sigmahp;
              if( tkCount < 20 ) {
                ++Ntrk2sigma20;
                if( tk->highPurity ) ++Ntrk2sigmahp20;
              }

              // Number of close tracks with doca > 3 sigma
              if( tk->docaSignificance > 3. ) {
                ++Ntrk3sigma;
                if( tk->highPurity ) ++Ntrk3sigmahp;
                if( tkCount < 20 ) {
                  ++Ntrk3sigma20;
                  if( tk->highPurity ) ++Ntrk3sigmahp20;
                }
              }
            }
          }
        }

        // SumPT of tracks from candidate PV or no PV
        if( tk->pt < 0.9 ) continue;
        if( deltaR(tk->eta ,tk->phi, dimuon_4mom.Eta(), dimuon_4mom.Phi()) > 0.7 ) continue;
        if( tk->vertexId == 1 ) {
          sumNdofPV += tk->ndof;
          // std::cout << "track(" << tkCount << ") from PV pt = " << tk->pt << std::endl;
          sumPTPV += tk->pt;
          if( tkCount < 20 ) sumPTPV20 += tk->pt;
        }
        // In this case the doca is also required to be less than 500 microns
        else if( tk->vertexId == 0 && tk->doca < 0.05 ) {
          sumNdofNoVtx += tk->ndof;
          // std::cout << "track(" << tkCount << ") from no PV pt = " << tk->pt << std::endl;
          sumPTNoVtx += tk->pt;
          if( tkCount < 20 ) sumPTNoVtx20 += tk->pt;
        }
      }
      double Iso = dimuon_4mom.Pt()/(dimuon_4mom.Pt()+sumPTNoVtx+sumPTPV);
      //      double Iso20 = myCand.pt()/(myCand.pt()+sumPTNoVtx20+sumPTPV);
      dimuon_iso[ndimuons] = Iso;

      ndimuons++;
    }
  }



}

vector <int> MuonAnalyzer::is_MC_matched(const edm::Event& iEvent, const Muon * mu) {
  string theGenLabel = "genParticles";

  //	EDGetTokenT<GenParticleCollection>  theGenLabel;
  //	theGenLabel = consumes<GenParticleCollection>(InputTag("genParticles"));

  edm::Handle<GenParticleCollection> genParticles;
  //	iEvent.getByToken(theGenLabel, genParticles);
  iEvent.getByLabel(GenParticle_, genParticles);

  vector <int> family_tree(3, -1);
  for(size_t i = 0; i < genParticles->size(); ++i) {
    const GenParticle & p = (*genParticles)[i];
    if (p.status() == 1) {
      TLorentzVector p_4mom(0., 0., 0., 0.);
      TLorentzVector muon_4mom(0., 0., 0., 0.);
      if (p.pt() < 2.5 || abs(p.eta()) > 2.4) continue;
      p_4mom.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      const Track * muon_track = mu->bestTrack();
      muon_4mom.SetPtEtaPhiM(muon_track->pt(), muon_track->eta(), muon_track->phi(), muon_mass);
      Double_t deltaR = p_4mom.DeltaR(muon_4mom);
      if (deltaR < 0.1 && p.charge() == mu->charge()) {
        int pdgid = p.pdgId();
        const Candidate * mom = p.mother();
        const Candidate * grandma = mom->mother();
        while (grandma->pdgId() > 80 && grandma->pdgId() < 101) {
          grandma = grandma->mother();
        }
        family_tree[0] = pdgid;
        family_tree[1] = mom->pdgId();
        family_tree[2] = grandma->pdgId();
        break;
        //        cout << "deltaR = " << deltaR << "  pdgId = " << pdgid << "  mother = " << mom->pdgId() << "  grandma = " << grandma->pdgId() << endl;

      }
    }
  }
  return family_tree;
}

void MuonAnalyzer::initialize_tree_vars(){

  reco_muon_4mom->Clear();
  reco_dimuon_4mom->Clear();
  gen_particle_4mom->Clear();
  L3_particle_4mom->Clear();
  hlt_vtx_4mom->Clear();

  ngen = 0;
  nL3 = 0;
  nmuons = 0;
  ndimuons = 0;
  nsv_hlt = 0;
  npv = 0;

  for (unsigned int i = 0; i < HLTPaths_.size(); i++) mHLTPaths_[HLTPaths_[i]] = false;
}

// ------------ method called once each job just before starting event loop  ------------
void MuonAnalyzer::beginJob() {
  output_->cd();

  if (doReco_) {
    tree_->Branch("nmuons",                  &nmuons,       "nmuons/I");
    tree_->Branch("reco_muon_charge",        charge,        "charge[nmuons]/I");

    tree_->Branch("ndimuons",                  &ndimuons,      "ndimuons/I");
    tree_->Branch("reco_dimuon_charge",        dicharge,        "charge[ndimuons]/I");
    tree_->Branch("reco_muon_4mom",          "TClonesArray", &reco_muon_4mom, 32000, 0);
    tree_->Branch("reco_muon_isTightMuon",   isTightMuon,      "isTightMuon[nmuons]/O");
    tree_->Branch("reco_dimuon_4mom",          "TClonesArray", &reco_dimuon_4mom, 32000, 0);
    tree_->Branch("reco_dimuon_iso",        dimuon_iso,        "reco_dimuon_iso[ndimuons]/D");
    tree_->Branch("reco_dimuon_mu1Index", mu1Index, "reco_dimuon_mu1Index[ndimuons]/I");
    tree_->Branch("reco_dimuon_mu2Index", mu2Index, "reco_dimuon_mu2Index[ndimuons]/I");

    tree_->Branch("npv", &npv, "npv/I");

    if (doMC_) {
      tree_->Branch("reco_muon_pdgId",         pdgId,         "pdgId[nmuons]/I");
      tree_->Branch("reco_muon_pdgId_mother",  pdgIdMom,      "pdgIdMom[nmuons]/I");
      tree_->Branch("reco_muon_pdgId_grandma", pdgIdGrandma,  "pdgIdGrandma[nmuons]/I");
    }

    if (doTrigger_) {
      tree_->Branch("reco_muon_isTriggerMatched",   isTriggerMatched,      "isTriggerMatched[nmuons]/O");
    }
  }

  if (doMC_) {
    tree_->Branch("ngen",          &ngen,  "ngen/I");
    tree_->Branch("gen_particle_4mom",       "TClonesArray",   &gen_particle_4mom, 32000, 0);
    tree_->Branch("gen_particle_pdgId",      gen_pdgId,        "gen_pdgId[ngen]/I");
    tree_->Branch("gen_particle_pdgIdMom",   gen_pdgIdMom,   "gen_pdgIdMom[ngen]/I");
    tree_->Branch("gen_particle_status",     gen_status,       "gen_status[ngen]/I");
  }

  if (doTrigger_) {
    tree_->Branch("nL3",          &nL3,  "nL3/I");
    tree_->Branch("L3_particle_4mom",       "TClonesArray",   &L3_particle_4mom, 32000, 0);
    tree_->Branch("L3_particle_id",      L3_particle_id,        "L3_particle_id[nL3]/I");
    for (unsigned int i = 0; i < HLTPaths_.size(); i++) {
      tree_->Branch(HLTPaths_[i].c_str(),   &mHLTPaths_[HLTPaths_[i]],      Form("%s/O", HLTPaths_[i].c_str()));
    }
  }
  if (doRAWTrigger_) {
    tree_->Branch("nsv_hlt",                  &nsv_hlt,       "nsv_hlt/I");
    tree_->Branch("hlt_vtx_normChi2",        normChi2,        "hlt_vtx_normChi2[nsv_hlt]/F");
    tree_->Branch("hlt_vtx_Prob",        vtxProb,        "hlt_vtx_Prob[nsv_hlt]/D");
    tree_->Branch("hlt_vtx_cosAlpha",        cosAlpha,        "hlt_vtx_cosAlpha[nsv_hlt]/F");
    tree_->Branch("hlt_vtx_LxyS",        LxySignificance,        "hlt_vtx_LxyS[nsv_hlt]/F");
    tree_->Branch("hlt_vtx_Lxy",        Lxy,        "hlt_vtx_Lxy[nsv_hlt]/F");
    tree_->Branch("hlt_vtx_Lxyerr",        Lxyerr,        "hlt_vtx_Lxyerr[nsv_hlt]/F");
    tree_->Branch("hlt_vtx_4mom",       "TClonesArray",   &hlt_vtx_4mom, 32000, 0);
  }
}

void MuonAnalyzer::fillGen(const edm::Event& iEvent) {
  string theGenLabel = "genParticles";
  Handle< GenParticleCollection > genParticles;

  //	EDGetTokenT<GenParticleCollection>  theGenLabel;
  //	theGenLabel = consumes<GenParticleCollection>(InputTag("genParticles"));
  //	edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC

  //	iEvent.getByToken(theGenLabel, genParticles);
  iEvent.getByLabel(GenParticle_, genParticles);

  for (vector<int>::iterator it = pdgId_.begin(); it != pdgId_.end(); ++it) {
    int pdgId_i = *it;
    for(size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & p = (*genParticles)[i];
      const Candidate * mom = p.mother();
      if (abs(p.pdgId()) == pdgId_i) {
      	gen_status[ngen] = p.status();
      	gen_pdgId[ngen] = p.pdgId();
        gen_pdgIdMom[ngen] = mom->pdgId();
      	TLorentzVector moment(0., 0., 0., 0.);
      	moment.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      	new((*gen_particle_4mom)[ngen]) TLorentzVector(moment);
      	ngen++;
      }
    }
  }
}

void MuonAnalyzer::fillTriggerHLTpaths(const edm::Event& iEvent) {
  //// HLT bits
  Handle<TriggerResults> hltresults;
  InputTag HLTbits("TriggerResults", "", HLTString_.c_str());
  iEvent.getByLabel(HLTbits, hltresults);

  if (hltresults.isValid()) {
    int ntrigs=hltresults->size();
    if (ntrigs==0) cout << "%HLTInfo -- No trigger name given in TriggerResults of the input " << endl;
    edm::TriggerNames triggerNames_ = iEvent.triggerNames(* hltresults);

    for (unsigned int i = 0; i < HLTPaths_.size(); i++) {

      for ( int j=0; j!=ntrigs; j++) {

        string trigName = triggerNames_.triggerName(j);
        size_t found = trigName.find(HLTPaths_[i]);
        if (found != std::string::npos) {
          mHLTPaths_[HLTPaths_[i]] = hltresults->accept(j);
          break;
        }
      }
    }
  }
  else cout << "hltresults is not valid" << endl;

}

void MuonAnalyzer::fillTriggerL3(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //// L3 objects
  //	EDGetTokenT<trigger::TriggerEvent>  trigEventToken;
  //	trigEventToken = consumes<trigger::TriggerEvent>(trigEventTag);
  edm::InputTag trigEventTag("hltTriggerSummaryAOD", "", HLTString_.c_str()); //make sure have correct process on MC

  Handle<trigger::TriggerEvent> trigEvent;
  //	iEvent.getByToken(trigEventToken, trigEvent);
  iEvent.getByLabel(trigEventTag, trigEvent);

  if (!trigEvent.isValid()) {
    cout << "Trigger summary product not found! Collection returns false always" << endl;
    return;
  }

  InputTag L3NameCollection("hltL3MuonCandidates", "", trigEventTag.process());

  trigger::size_type Index(0);
  Index = trigEvent->collectionIndex(L3NameCollection);
  if (Index < trigEvent->sizeCollections()) {
    const trigger::Keys& Keys(trigEvent->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( trigEvent->getObjects().at(i) );

      if (abs(obj.id()) == 13) {
        TLorentzVector L3muon_4mom(0., 0., 0., 0.);
        L3muon_4mom.SetPtEtaPhiM(obj.pt(), obj.eta(), obj.phi(), obj.mass());
        new((*L3_particle_4mom)[nL3]) TLorentzVector(L3muon_4mom);
        L3_particle_id[nL3] = obj.id();
        nL3++;
        for (trigger::size_type j = n0+1; j != n1; ++j) {
          const trigger::TriggerObject& obj_j( trigEvent->getObjects().at(j) );
          if (abs(obj_j.id()) == 13) {
            TLorentzVector L3muon_j_4mom(0., 0., 0., 0.);
            L3muon_j_4mom.SetPtEtaPhiM(obj_j.pt(), obj_j.eta(), obj_j.phi(), obj_j.mass());
          }
        }
      }
    }
  }
}

void MuonAnalyzer::fillRAWTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  //  InputTag L3NameCollection("hltL3MuonCandidates", "", "HLT");
  InputTag L3NameCollection("hltL3MuonCandidates", "", HLTString_.c_str());

  Handle<reco::RecoChargedCandidateCollection> mucands;
  iEvent.getByLabel(L3NameCollection, mucands);

  if (!mucands.isValid()) return;

  //get the transient track builder:
  edm::ESHandle<TransientTrackBuilder> theB;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  std::auto_ptr<VertexCollection> vertexCollection(new VertexCollection());

  // look at all mucands, check cuts and make vertices
  double e1,e2;
  Particle::LorentzVector p,p1,p2;

  RecoChargedCandidateCollection::const_iterator cand1;
  RecoChargedCandidateCollection::const_iterator cand2;

  for (cand1=mucands->begin(); cand1!=mucands->end(); cand1++) {
    TrackRef tk1 = cand1->get<TrackRef>();
    cand2 = cand1; cand2++;
    for (; cand2!=mucands->end(); cand2++) {
      TrackRef tk2 = cand2->get<TrackRef>();
      if (cand1->charge()*cand2->charge()>0) continue;

      // Combined dimuon system
      e1 = sqrt(cand1->momentum().Mag2()+MuMass2);
      e2 = sqrt(cand2->momentum().Mag2()+MuMass2);
      p1 = Particle::LorentzVector(cand1->px(),cand1->py(),cand1->pz(),e1);
      p2 = Particle::LorentzVector(cand2->px(),cand2->py(),cand2->pz(),e2);
      p = p1+p2;

      // do the vertex fit
      vector<TransientTrack> t_tks;
      TransientTrack ttkp1 = (*theB).build(&tk1);
      TransientTrack ttkp2 = (*theB).build(&tk2);
      t_tks.push_back(ttkp1);
      t_tks.push_back(ttkp2);

      KalmanVertexFitter kvf;
      TransientVertex tv = kvf.vertex(t_tks);
      if (!tv.isValid()) continue;

      Vertex vertex = tv;
      //      cout << vertex.x() << " " << vertex.y() << " " << vertex.z() << endl;
      // put vertex in the event
      vertexCollection->push_back(vertex);
    }
  }

  // get beam spot
  reco::BeamSpot vertexBeamSpot;
  edm::Handle<reco::BeamSpot> onlineBeamSpotHande;
  iEvent.getByLabel("hltOnlineBeamSpot",onlineBeamSpotHande);
  vertexBeamSpot = *onlineBeamSpotHande;


  for(reco::VertexCollection::iterator it = vertexCollection->begin(); it!= vertexCollection->end(); it++) {
    reco::Vertex displacedVertex = *it;
    if( (displacedVertex.chi2()>=0.0) && (displacedVertex.ndof()>0) ) vtxProb[nsv_hlt] = TMath::Prob(displacedVertex.chi2(), displacedVertex.ndof() );

    normChi2[nsv_hlt] = displacedVertex.normalizedChi2();

    // get the two muons from the vertex
    reco::Vertex::trackRef_iterator trackIt = displacedVertex.tracks_begin();
    reco::TrackRef vertextkRef1 = (*trackIt).castTo<reco::TrackRef>() ;
    // the second one
    trackIt++;
    reco::TrackRef vertextkRef2 = (*trackIt).castTo<reco::TrackRef>();

    // first find these two tracks in the muon collection
    reco::RecoChargedCandidateCollection::const_iterator cand1;
    reco::RecoChargedCandidateCollection::const_iterator cand2;

    int iFoundRefs = 0;
    for (reco::RecoChargedCandidateCollection::const_iterator cand=mucands->begin(); cand!=mucands->end(); cand++) {
      reco::TrackRef tkRef = cand->get<reco::TrackRef>();
      if(tkRef == vertextkRef1) {cand1 = cand; iFoundRefs++;}
      if(tkRef == vertextkRef2) {cand2 = cand; iFoundRefs++;}
    }
    if(iFoundRefs != 2) throw cms::Exception("BadLogic") << "HLTDisplacedmumuFilter: ERROR: the Jpsi vertex must have exactly two muons by definition." << std::endl;

    TLorentzVector cand1_4mom(0., 0., 0., 0.);
    TLorentzVector cand2_4mom(0., 0., 0., 0.);
    cand1_4mom.SetPtEtaPhiM(cand1->pt(), cand1->eta(), cand1->phi(), muon_mass);
    cand2_4mom.SetPtEtaPhiM(cand2->pt(), cand2->eta(), cand2->phi(), muon_mass);
    TLorentzVector vtx_4mom(0., 0., 0., 0.);
    vtx_4mom = cand1_4mom + cand2_4mom;
    new((*hlt_vtx_4mom)[nsv_hlt]) TLorentzVector(vtx_4mom);

    // calculate two-track transverse momentum
    math::XYZVector pperp(cand1->px() + cand2->px(),
                          cand1->py() + cand2->py(),
                          0.);


    reco::Vertex::Point vpoint=displacedVertex.position();
    //translate to global point, should be improved
    GlobalPoint secondaryVertex (vpoint.x(), vpoint.y(), vpoint.z());

    reco::Vertex::Error verr = displacedVertex.error();
    // translate to global error, should be improved
    GlobalError err(verr.At(0,0), verr.At(1,0), verr.At(1,1), verr.At(2,0), verr.At(2,1), verr.At(2,2) );

    GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) + (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()),
                                          -1*((vertexBeamSpot.y0() - secondaryVertex.y())+ (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 0);

    Lxy[nsv_hlt] = displacementFromBeamspot.perp();
    Lxyerr[nsv_hlt] = sqrt(err.rerr(displacementFromBeamspot));
    LxySignificance[nsv_hlt] = Lxy[nsv_hlt]/Lxyerr[nsv_hlt];

    //calculate the angle between the decay length and the mumu momentum
    reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);

    cosAlpha[nsv_hlt] = vperp.Dot(pperp)/(vperp.R()*pperp.R());

    nsv_hlt++;
  }
}

bool MuonAnalyzer::IsTightMuon(const reco::Muon & muon) {
  //if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;
  if(!muon.isGlobalMuon() ) return false;
  bool muID = isGoodMuon(muon, muon::GlobalMuonPromptTight) && (muon.numberOfMatchedStations() > 1);
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
  //  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  return muID && hits/* && ip*/;
}

bool MuonAnalyzer::is_trigger_matched(const edm::Event& iEvent, const TLorentzVector mu) {


  //	EDGetTokenT<trigger::TriggerEvent>  trigEventToken;
  //	trigEventToken = consumes<trigger::TriggerEvent>(trigEventTag);
  //	edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","", HLTString_.c_str()); //make sure have correct process on MC

  Handle<trigger::TriggerEvent> trigEvent;
  //	iEvent.getByToken(trigEventToken, trigEvent);
  iEvent.getByLabel(trigEventTag, trigEvent);

  if (!trigEvent.isValid()) {
    cout << "Trigger summary product not found! Collection returns false always";
    return false;
  }

  //std::string filterName("hltVertexmumuFilterBs47"); //hltVertexmumuFilterBs345
  string L3NameCollection("hltL3MuonCandidates");

  trigger::size_type Index(0);
  Index = trigEvent->collectionIndex(edm::InputTag(L3NameCollection, "", trigEventTag.process()));
  if (Index < trigEvent->sizeCollections()) {
    const trigger::Keys& Keys(trigEvent->collectionKeys());
    const trigger::size_type n0 (Index == 0? 0 : Keys.at(Index-1));
    const trigger::size_type n1 (Keys.at(Index));
    for (trigger::size_type i = n0; i != n1; ++i) {
      const trigger::TriggerObject& obj( trigEvent->getObjects().at(i) );
      if (abs(obj.id()) == 13) {
        TLorentzVector L3muon_4mom(0., 0., 0., 0.);
        L3muon_4mom.SetPtEtaPhiM(obj.pt(), obj.eta(), obj.phi(), obj.mass());
        if (mu.DeltaR(L3muon_4mom) < 0.1) return true;
      }
    }
  }
  return false;
}

// ------------ method called once each job just after ending the event loop  ------------
void MuonAnalyzer::endJob() {
  output_->Write();
}

// ------------ method called when starting to processes a run  ------------
void
MuonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void
MuonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void
MuonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
MuonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int MuonAnalyzer::findVertexId(const reco::Vertex & theOriginalPV,
                               const std::vector<reco::Vertex > &priVtxs,
                               const reco::Track & track, const double & maxDeltaR)
{
  // Check if the track belongs to the PV associated to the candidate
  for(auto PVtk = theOriginalPV.tracks_begin(); PVtk != theOriginalPV.tracks_end(); ++PVtk) {
    if( PVtk->isNonnull() ) {
      if( PVtk->key() == track.extra().key() ) {
        // if( deltaR(track.eta(), track.phi(), (*PVtk)->eta(), (*PVtk)->phi()) < maxDeltaR ) return 1;
        return 1;
      }
    }
  }
  // Check if the track is associated to any other PV
  for(auto itv = priVtxs.begin(), itvend = priVtxs.end(); itv != itvend; ++itv){
    for( auto vtxTk = itv->tracks_begin(); vtxTk != itv->tracks_end(); ++vtxTk) {
      if( vtxTk->isNonnull() ) {
        if( vtxTk->key() == track.extra().key() ) {
          return 2;
        }
        // if( deltaR(track.eta(), track.phi(), (*vtxTk)->eta(), (*vtxTk)->phi()) < maxDeltaR ) return 2;
      }
    }
  }
  // Track is not associated to any PV
  return 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzer);
