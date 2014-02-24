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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
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

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

using namespace std;
using namespace edm;
using namespace reco;

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

      virtual void initialize_tree_vars();
      virtual void fillGen(const edm::Event& iEvent);
      virtual void fillReco(const edm::Event& iEvent);
      
      
      virtual vector<int> is_MC_matched(const edm::Event& iEvent, const Muon * mu);
      virtual bool is_trigger_matched(const edm::Event& iEvent, const TLorentzVector mu);
      
      virtual bool IsTightMuon(const reco::Muon& muon);
      
      // ----------member data ---------------------------
      string outputname_;
      Bool_t doMC_;
      Bool_t doReco_;
      vector< Int_t> pdgId_;

      Int_t ngen;
      TFile * output_;
      TTree * tree_;
      TClonesArray * gen_particle_4mom;
      TClonesArray * reco_muon_4mom;
      TClonesArray * reco_dimuon_4mom;
      Int_t nmuons;
      Int_t ndimuons;
      Int_t charge[100];
      Int_t dicharge[100];
      Int_t pdgId[100];
      Int_t pdgIdMom[100];
      Int_t pdgIdGrandma[100];
      
      Int_t gen_pdgId[100];
      Int_t gen_status[100];

      bool isTightMuon[100];
      bool isTriggerMatched[100];

};

//
// constants, enums and typedefs
//
const Double_t muon_mass = 0.105658;
//
// static data member definitions
//

//
// constructors and destructor
//
MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  outputname_ = iConfig.getParameter<string> ("OutputFileName");
  doMC_ = iConfig.getParameter<bool> ("doMC");
  doReco_ = iConfig.getParameter<bool> ("doReco");
  pdgId_ = iConfig.getParameter <vector <int> > ("pdgId");

  output_ = new TFile (outputname_.c_str(), "RECREATE" );
  tree_ = new TTree("tree","tree");

  reco_muon_4mom = new TClonesArray("TLorentzVector", 100);
  reco_dimuon_4mom = new TClonesArray("TLorentzVector", 100);
  gen_particle_4mom = new TClonesArray("TLorentzVector", 100);
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
  if (doReco_) fillReco(iEvent);
  ////////////////
  tree_->Fill();//
  ////////////////
}

void MuonAnalyzer::fillReco(const edm::Event& iEvent) {
	
  string theMuonLabel = "muons";
  
//	EDGetTokenT<MuonCollection> theMuonLabel;
//	theMuonLabel = consumes<MuonCollection>(InputTag("muons"));
	
	edm::Handle<MuonCollection> muons;
//	iEvent.getByToken(theMuonLabel, muons);
	iEvent.getByLabel(theMuonLabel, muons);
  
  if (!muons.isValid()) {
    cout << "theMuonLabel is not valid" << endl;
    return;
  }
  
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
    isTriggerMatched[nmuons] = is_trigger_matched(iEvent, muon_4mom);

    if (doMC_) {
      vector <int> family(3, -1);
      family = is_MC_matched(iEvent, mu);
      pdgId[nmuons] = family[0];
      pdgIdMom[nmuons] = family[1];
      pdgIdGrandma[nmuons] = family[2];
    }

    // ID
    isTightMuon[nmuons] = IsTightMuon(*mu);
    nmuons++;
    
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
	iEvent.getByLabel(theGenLabel, genParticles);
  
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

  ngen = 0;
  nmuons = 0;
  ndimuons = 0;
}

// ------------ method called once each job just before starting event loop  ------------
void MuonAnalyzer::beginJob() {
  output_->cd();

  tree_->Branch("nmuons",                  &nmuons,       "nmuons/I");
  tree_->Branch("reco_muon_charge",        charge,        "charge[nmuons]/I");
  
  tree_->Branch("ndimuons",                  &ndimuons,      "ndimuons/I");
  tree_->Branch("reco_dimuon_charge",        dicharge,        "charge[ndimuons]/I");

  if (doMC_) {
  	tree_->Branch("ngen",          &ngen,  "ngen/I");
  	tree_->Branch("gen_particle_4mom",       "TClonesArray",   &gen_particle_4mom, 32000, 0);
  	tree_->Branch("gen_particle_pdgId",      gen_pdgId,        "pdgId[ngen]/I");
  	tree_->Branch("gen_particle_status",     gen_status,       "gen_status[ngen]/I");
 
  	tree_->Branch("reco_muon_pdgId",         pdgId,         "pdgId[nmuons]/I");
  	tree_->Branch("reco_muon_pdgId_mother",  pdgIdMom,      "pdgIdMom[nmuons]/I");
  	tree_->Branch("reco_muon_pdgId_grandma", pdgIdGrandma,  "pdgIdGrandma[nmuons]/I");
  }

  tree_->Branch("reco_muon_4mom",          "TClonesArray", &reco_muon_4mom, 32000, 0);
  tree_->Branch("reco_muon_isTightMuon",   isTightMuon,      "isTightMuon[nmuons]/O");
  tree_->Branch("reco_muon_isTriggerMatched",   isTriggerMatched,      "isTriggerMatched[nmuons]/O");
  
  tree_->Branch("reco_dimuon_4mom",          "TClonesArray", &reco_dimuon_4mom, 32000, 0);

}

void MuonAnalyzer::fillGen(const edm::Event& iEvent) {
  string theGenLabel = "genParticles";
  Handle< GenParticleCollection > genParticles;
  
//	EDGetTokenT<GenParticleCollection>  theGenLabel;
//	theGenLabel = consumes<GenParticleCollection>(InputTag("genParticles"));
//	edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
	
  //	iEvent.getByToken(theGenLabel, genParticles);
	iEvent.getByLabel(theGenLabel, genParticles);
  
  for (vector<int>::iterator it = pdgId_.begin(); it != pdgId_.end(); ++it) {
    int pdgId_i = *it;
    for(size_t i = 0; i < genParticles->size(); ++i) {
      const GenParticle & p = (*genParticles)[i];
//      cout << p.pdgId()  << "  " << p.status() << endl;
      if (abs(p.pdgId()) == pdgId_i) {
      	gen_status[ngen] = p.status();
      	gen_pdgId[ngen] = p.pdgId();
      	TLorentzVector moment(0., 0., 0., 0.);
      	moment.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), p.mass());
      	new((*gen_particle_4mom)[ngen]) TLorentzVector(moment);
      	ngen++;
      }
    }
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
//	edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
//	trigEventToken = consumes<trigger::TriggerEvent>(trigEventTag);
	edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
	
	Handle<trigger::TriggerEvent> trigEvent;
//	iEvent.getByToken(trigEventToken, trigEvent);
	iEvent.getByLabel(trigEventTag, trigEvent);

	if (!trigEvent.isValid()) {
		//cout << "Trigger summary product not found! Collection returns false always";
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

//define this as a plug-in
DEFINE_FWK_MODULE(MuonAnalyzer);
