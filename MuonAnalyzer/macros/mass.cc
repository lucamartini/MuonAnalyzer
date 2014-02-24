#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TAxis.h"

#include <iostream> 

using namespace std;

void plot(TH1D * h_5314, TH1D * h_620SLHC7) {
	TCanvas * c = new TCanvas("c", "c", 600, 600);
	h_5314->SetMaximum(1.1*max(h_5314->GetMaximum(), h_620SLHC7->GetMaximum()));
	h_5314->GetYaxis()->SetTitleOffset(1.4);
	h_5314->SetLineColor(kBlue);
	h_5314->SetFillColor(kBlue);
	h_5314->SetFillStyle(3004);
	h_5314->Draw();
	h_620SLHC7->SetLineColor(kRed);
	h_620SLHC7->SetFillColor(kRed);
	h_620SLHC7->SetFillStyle(3005);
	h_620SLHC7->Draw("sames");
	c->Update();
	TPaveStats *st_5314 = (TPaveStats*)h_5314->FindObject("stats");
	st_5314->SetOptStat(1101);
	st_5314->SetLineColor(kBlue);
	TPaveStats *st_620SLHC7 = (TPaveStats*)h_620SLHC7->FindObject("stats");
	st_620SLHC7->SetOptStat(1101);
	st_620SLHC7->SetLineColor(kRed);
	st_620SLHC7->SetY1NDC(st_620SLHC7->GetY1NDC()-0.2);
	st_620SLHC7->SetY2NDC(st_620SLHC7->GetY2NDC()-0.2);
	string name(h_5314->GetName());
	name.erase(name.end()-7, name.end());
	c->Print(Form("%s.pdf", name.c_str()));
	delete c;
}

int main() {
	
	double eta = 1.4;
	
	// Declaration of leaf types
	Int_t           nmuons;
	Int_t           reco_muon_charge[100];   //[nmuons]
	Int_t           ndimuons;
	Int_t           reco_dimuon_charge[100];   //[ndimuons]
	Int_t           ngen;
	TClonesArray    *gen_particle_4mom;
	Int_t           gen_particle_pdgId[100];   //[ngen]
	Int_t           gen_particle_status[100];   //[ngen]
	Int_t           reco_muon_pdgId[100];   //[nmuons]
	Int_t           reco_muon_pdgId_mother[100];   //[nmuons]
	Int_t           reco_muon_pdgId_grandma[100];   //[nmuons]
	TClonesArray    *reco_muon_4mom;
	Bool_t          reco_muon_isTightMuon[100];   //[nmuons]
	Bool_t          reco_muon_isTriggerMatched[100];   //[nmuons]
	TClonesArray    *reco_dimuon_4mom;
	
	TBranch        *b_nmuons;   //!
	TBranch        *b_reco_muon_charge;   //!
	TBranch        *b_ndimuons;   //!
	TBranch        *b_reco_dimuon_charge;   //!
	TBranch        *b_ngen;   //!
	TBranch        *b_gen_particle_4mom;   //!
	TBranch        *b_gen_particle_pdgId;   //!
	TBranch        *b_gen_particle_status;   //!
	TBranch        *b_reco_muon_pdgId;   //!
	TBranch        *b_reco_muon_pdgId_mother;   //!
	TBranch        *b_reco_muon_pdgId_grandma;   //!
	TBranch        *b_reco_muon_4mom;   //!
	TBranch        *b_reco_muon_isTightMuon;   //!
	TBranch        *b_reco_muon_isTriggerMatched;   //!
	TBranch        *b_reco_dimuon_4mom;   //!
	
	gen_particle_4mom = 0;
	reco_muon_4mom = 0;
	reco_dimuon_4mom = 0;
	// Set branch addresses and branch pointers
	
	TFile* input_5314 = new TFile("muontree_5314_v2.root");
	TTree* tree_5314 = (TTree*)input_5314->Get("tree");
	
	TH1D * mass_barrel_h_5314 = new TH1D("mass_barrel_h_5314", Form("barrel mass (|#eta_{#mu^{1}}| < %.1f AND |#eta_{#mu^{2}}| < %.1f);mass (GeV); Normalized", eta, eta), 100, 4.9, 5.9);
	TH1D * mass_endcap_h_5314 = new TH1D("mass_endcap_h_5314", Form("endcap mass (|#eta_{#mu^{1}}| > %.1f OR |#eta_{#mu^{2}}| > %.1f);mass (GeV); Normalized", eta, eta), 100, 4.9, 5.9);
	TH1D * pt_res_barrel_h_5314 = new TH1D("pt_res_barrel_h_5314", "barrel p_{T} resolution (|#eta_{#mu}| < 1.2);p_{T}^{gen}-p_{T}^{reco} (GeV); Normalized", 100, -1., 1.);
	TH1D * pt_res_intermediate_h_5314 = new TH1D("pt_res_intermediate_h_5314", "intermediate p_{T} resolution (1.2 < |#eta_{#mu}| < 1.6);p_{T}^{gen}-p_{T}^{reco} (GeV); Normalized", 100, -1., 1.);
	TH1D * pt_res_endcap_h_5314 = new TH1D("pt_res_endcap_h_5314", "endcap p_{T} resolution (1.6 < |#eta_{#mu}| < 2.4);p_{T}^{gen}-p_{T}^{reco} (GeV); Normalized", 100, -1., 1.);
	
	tree_5314->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
	tree_5314->SetBranchAddress("reco_muon_charge", reco_muon_charge, &b_reco_muon_charge);
	tree_5314->SetBranchAddress("ndimuons", &ndimuons, &b_ndimuons);
	tree_5314->SetBranchAddress("reco_dimuon_charge", reco_dimuon_charge, &b_reco_dimuon_charge);
	tree_5314->SetBranchAddress("ngen", &ngen, &b_ngen);
	tree_5314->SetBranchAddress("gen_particle_4mom", &gen_particle_4mom, &b_gen_particle_4mom);
	tree_5314->SetBranchAddress("gen_particle_pdgId", gen_particle_pdgId, &b_gen_particle_pdgId);
	tree_5314->SetBranchAddress("gen_particle_status", gen_particle_status, &b_gen_particle_status);
	tree_5314->SetBranchAddress("reco_muon_pdgId", reco_muon_pdgId, &b_reco_muon_pdgId);
	tree_5314->SetBranchAddress("reco_muon_pdgId_mother", reco_muon_pdgId_mother, &b_reco_muon_pdgId_mother);
	tree_5314->SetBranchAddress("reco_muon_pdgId_grandma", reco_muon_pdgId_grandma, &b_reco_muon_pdgId_grandma);
	tree_5314->SetBranchAddress("reco_muon_4mom", &reco_muon_4mom, &b_reco_muon_4mom);
	tree_5314->SetBranchAddress("reco_muon_isTightMuon", reco_muon_isTightMuon, &b_reco_muon_isTightMuon);
	tree_5314->SetBranchAddress("reco_muon_isTriggerMatched", reco_muon_isTriggerMatched, &b_reco_muon_isTriggerMatched);
	tree_5314->SetBranchAddress("reco_dimuon_4mom", &reco_dimuon_4mom, &b_reco_dimuon_4mom);
	
	Long64_t nentries_5314 = tree_5314->GetEntries();
	
	for (Long64_t i = 0; i < nentries_5314; i++) {
		tree_5314->GetEntry(i);
		for (int j = 0; j < nmuons; j++) {
			TLorentzVector * mu_j = (TLorentzVector*)reco_muon_4mom->At(j);
			if (abs(reco_muon_pdgId[j]) == 13 && abs(reco_muon_pdgId_mother[j]) == 531 && reco_muon_isTightMuon[j] == true && mu_j->Pt() > 4.) {
				double deltar_min = 100.;
				double pt_gen = -1000;
				double pt_reco = mu_j->Pt();
				double eta_reco = mu_j->Eta();
				for (int jj = 0; jj < ngen; jj++) {
					if (abs(gen_particle_pdgId[jj]) == 13 && gen_particle_status[jj] == 1) {
						TLorentzVector * mu_jj = (TLorentzVector*)gen_particle_4mom->At(jj);
						Double_t deltar = mu_j->DeltaR(*mu_jj);
						if (deltar < deltar_min) {
							deltar_min = deltar;
							pt_gen = mu_jj->Pt();
						}
					}
				}
				if (deltar_min < 0.1) {
					double delta_pt = pt_gen - pt_reco;
					if (fabs(eta_reco) < 1.2) pt_res_barrel_h_5314->Fill(delta_pt);
					else if (fabs(eta_reco) < 1.6) pt_res_intermediate_h_5314->Fill(delta_pt);
					else if (fabs(eta_reco) < 2.4) pt_res_endcap_h_5314->Fill(delta_pt);
				}
				for (int k = j+1; k < nmuons; k++) {
					TLorentzVector * mu_k = (TLorentzVector*)reco_muon_4mom->At(k);
					if (abs(reco_muon_pdgId[k]) == 13 && abs(reco_muon_pdgId_mother[k]) == 531 && reco_muon_isTightMuon[k] == true && mu_k->Pt() > 4.) {
						if (reco_muon_charge[j] + reco_muon_charge[k] == 0){
							TLorentzVector dimuon(0.,0.,0.,0.);
							dimuon = *mu_j + *mu_k;
							if ((abs(mu_j->Eta()) < eta && abs(mu_k->Eta()) < eta)) {
								mass_barrel_h_5314->Fill(dimuon.M());
							}
							else {
								TLorentzVector dimuon(0.,0.,0.,0.);
								dimuon = *mu_j + *mu_k;
								mass_endcap_h_5314->Fill(dimuon.M());
							}
						}
					}
				}
			}
		}
	}	
//	input_5314->Close();
	mass_barrel_h_5314->Scale(1./mass_barrel_h_5314->GetEntries());
	mass_endcap_h_5314->Scale(1./mass_endcap_h_5314->GetEntries());
	
	pt_res_barrel_h_5314->Scale(1./pt_res_barrel_h_5314->GetEntries());
	pt_res_intermediate_h_5314->Scale(1./pt_res_intermediate_h_5314->GetEntries());
	pt_res_endcap_h_5314->Scale(1./pt_res_endcap_h_5314->GetEntries());
	
	TFile* input_620SLHC7 = new TFile("muontree_620SLHC7_v2.root");
	TTree* tree_620SLHC7 = (TTree*)input_620SLHC7->Get("tree");
	
	TH1D * mass_barrel_h_620SLHC7 = new TH1D("mass_barrel_h_620SLHC7", Form("barrel mass (|#eta_{#mu^{1}}| < %.1f AND |#eta_{#mu^{2}}| < %.1f);mass (GeV); Normalized", eta, eta), 100, 4.9, 5.9);
	TH1D * mass_endcap_h_620SLHC7 = new TH1D("mass_endcap_h_620SLHC7", Form("endcap mass (|#eta_{#mu^{1}}| > %.1f OR |#eta_{#mu^{2}}| > %.1f);mass (GeV); Normalized", eta, eta), 100, 4.9, 5.9);
	TH1D * pt_res_barrel_h_620SLHC7 = new TH1D("pt_res_barrel_h_620SLHC7", "barrel p_{T} resolution (|#eta_{#mu}| < 1.2);p_{T}^{gen}-p_{T}^{reco} (GeV); Normalized", 100, -1., 1.);
	TH1D * pt_res_intermediate_h_620SLHC7 = new TH1D("pt_res_intermediate_h_620SLHC7", "intermediate p_{T} resolution (1.2 < |#eta_{#mu}| < 1.6);p_{T}^{gen}-p_{T}^{reco} (GeV); Normalized", 100, -1., 1.);
	TH1D * pt_res_endcap_h_620SLHC7 = new TH1D("pt_res_endcap_h_620SLHC7", "endcap p_{T} resolution (1.6 < |#eta_{#mu}| < 2.4);p_{T}^{gen}-p_{T}^{reco} (GeV); Normalized", 100, -1., 1.);
	
	tree_620SLHC7->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
	tree_620SLHC7->SetBranchAddress("reco_muon_charge", reco_muon_charge, &b_reco_muon_charge);
	tree_620SLHC7->SetBranchAddress("ndimuons", &ndimuons, &b_ndimuons);
	tree_620SLHC7->SetBranchAddress("reco_dimuon_charge", reco_dimuon_charge, &b_reco_dimuon_charge);
	tree_620SLHC7->SetBranchAddress("ngen", &ngen, &b_ngen);
	tree_620SLHC7->SetBranchAddress("gen_particle_4mom", &gen_particle_4mom, &b_gen_particle_4mom);
	tree_620SLHC7->SetBranchAddress("gen_particle_pdgId", gen_particle_pdgId, &b_gen_particle_pdgId);
	tree_620SLHC7->SetBranchAddress("gen_particle_status", gen_particle_status, &b_gen_particle_status);
	tree_620SLHC7->SetBranchAddress("reco_muon_pdgId", reco_muon_pdgId, &b_reco_muon_pdgId);
	tree_620SLHC7->SetBranchAddress("reco_muon_pdgId_mother", reco_muon_pdgId_mother, &b_reco_muon_pdgId_mother);
	tree_620SLHC7->SetBranchAddress("reco_muon_pdgId_grandma", reco_muon_pdgId_grandma, &b_reco_muon_pdgId_grandma);
	tree_620SLHC7->SetBranchAddress("reco_muon_4mom", &reco_muon_4mom, &b_reco_muon_4mom);
	tree_620SLHC7->SetBranchAddress("reco_muon_isTightMuon", reco_muon_isTightMuon, &b_reco_muon_isTightMuon);
	tree_620SLHC7->SetBranchAddress("reco_muon_isTriggerMatched", reco_muon_isTriggerMatched, &b_reco_muon_isTriggerMatched);
	tree_620SLHC7->SetBranchAddress("reco_dimuon_4mom", &reco_dimuon_4mom, &b_reco_dimuon_4mom);
	
	Long64_t nentries_620SLHC7 = tree_620SLHC7->GetEntries();
	
	for (Long64_t i = 0; i < nentries_620SLHC7; i++) {
		tree_620SLHC7->GetEntry(i);
		for (int j = 0; j < nmuons; j++) {
			TLorentzVector * mu_j = (TLorentzVector*)reco_muon_4mom->At(j);
			if (abs(reco_muon_pdgId[j]) == 13 && abs(reco_muon_pdgId_mother[j]) == 531 && reco_muon_isTightMuon[j] == true && mu_j->Pt() > 4.) {
				
				double deltar_min = 100.;
				double pt_gen = -1000;
				double pt_reco = mu_j->Pt();
				double eta_reco = mu_j->Eta();
				for (int jj = 0; jj < ngen; jj++) {
					if (abs(gen_particle_pdgId[jj]) == 13 && gen_particle_status[jj] == 1) {
						TLorentzVector * mu_jj = (TLorentzVector*)gen_particle_4mom->At(jj);
						Double_t deltar = mu_j->DeltaR(*mu_jj);
						if (deltar < deltar_min) {
							deltar_min = deltar;
							pt_gen = mu_jj->Pt();
						}
					}
				}
				if (deltar_min < 0.1) {
					double delta_pt = pt_gen - pt_reco;
					if (fabs(eta_reco) < 1.2) pt_res_barrel_h_620SLHC7->Fill(delta_pt);
					else if (fabs(eta_reco) < 1.6) pt_res_intermediate_h_620SLHC7->Fill(delta_pt);
					else if (fabs(eta_reco) < 2.4) pt_res_endcap_h_620SLHC7->Fill(delta_pt);
				}
				for (int k = j+1; k < nmuons; k++) {
					TLorentzVector * mu_k = (TLorentzVector*)reco_muon_4mom->At(k);
					if (abs(reco_muon_pdgId[k]) == 13 && abs(reco_muon_pdgId_mother[k]) == 531 && reco_muon_isTightMuon[k] == true && mu_k->Pt() > 4.) {
						if (reco_muon_charge[j] + reco_muon_charge[k] == 0){
							TLorentzVector dimuon(0.,0.,0.,0.);
							dimuon = *mu_j + *mu_k;
							if ((abs(mu_j->Eta()) < eta && abs(mu_k->Eta()) < eta)) {
								mass_barrel_h_620SLHC7->Fill(dimuon.M());
							}
							else {
								TLorentzVector dimuon(0.,0.,0.,0.);
								dimuon = *mu_j + *mu_k;
								mass_endcap_h_620SLHC7->Fill(dimuon.M());
							}
						}
					}
				}
			}
		}
	}	
//	input_620SLHC7->Close();
	mass_barrel_h_620SLHC7->Scale(1./mass_barrel_h_620SLHC7->GetEntries());
	mass_endcap_h_620SLHC7->Scale(1./mass_endcap_h_620SLHC7->GetEntries());
		
	pt_res_barrel_h_620SLHC7->Scale(1./pt_res_barrel_h_620SLHC7->GetEntries());
	pt_res_intermediate_h_620SLHC7->Scale(1./pt_res_intermediate_h_620SLHC7->GetEntries());
	pt_res_endcap_h_620SLHC7->Scale(1./pt_res_endcap_h_620SLHC7->GetEntries());
		
		
	plot(mass_barrel_h_5314, mass_barrel_h_620SLHC7);
	plot(mass_endcap_h_5314, mass_endcap_h_620SLHC7);
	plot(pt_res_barrel_h_5314, pt_res_barrel_h_620SLHC7);
	plot(pt_res_intermediate_h_5314, pt_res_intermediate_h_620SLHC7);
	plot(pt_res_endcap_h_5314, pt_res_endcap_h_620SLHC7);
	
	return 1;
}