#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuEfficiencyHists::LQToTopMuEfficiencyHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1D>("N_ele", "N_{e}", 10, -0.5, 9.5); 
  book<TH1D>("h_N_gen_ele", "N_{e} generated", 10, -0.5, 9.5); 

  book<TH1D>("N_mu", "N_{#mu}", 10, -0.5, 9.5); 
  book<TH1D>("h_N_gen_mu", "N_{#mu} generated", 10, -0.5, 9.5); 


  // general
  book<TH1D>("H_T", "H_{T}", 48, 0, 4200);
  book<TH1D>("Parton_H_T", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_2gen_muons", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_2gen_ele", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_2matched_muons", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_2matched_ele", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_2rec_muons", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_2rec_ele", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_not2matched_muons", "H_{T}", 48, 0, 4200);
  book<TH1D>("H_T_not2matched_ele", "H_{T}", 48, 0, 4200);
  book<TH1F>("jets_faking_ele_pt", "p_{T} of jets faking an electron [GeV]", 50, 20, 1500);
  book<TH1F>("jets_faking_mu_pt", "p_{T} of jets faking a muon [GeV]", 50, 20, 1500);
  book<TH1F>("jets_2rec_ele_pt", "p_{T} of jets with 2 reconstructed electrons [GeV]", 50, 20, 1500);
  book<TH1F>("jets_2rec_mu_pt", "p_{T} of jets with 2 reconstructed muons [GeV]", 50, 20, 1500);
  book<TH1F>("ele_2rec_ele_pt", "p_{T} of reconstucted electrons with #geq 2 electrons [GeV]", 50, 20, 1500);
  book<TH1F>("mu_2rec_mu_pt", "p_{T} of reconstructed muons with #geq 2 muons [GeV]", 50, 20, 1500);

  book<TH1D>("H_T_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("Parton_H_T_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_2gen_muons_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_2gen_ele_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_2matched_muons_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_2matched_ele_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_2rec_muons_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_2rec_ele_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_not2matched_muons_rebin", "H_{T}", 20, bins_from350);
  book<TH1D>("H_T_not2matched_ele_rebin", "H_{T}", 20, bins_from350);

  book<TH2D>("H_T_rec_gen_2rec_muons", "H_{T};H_{T}^{rec};H_{T}^{gen}", 48, 0, 4200, 48, 0, 4200);
  book<TH2D>("H_T_rec_gen_2rec_ele", "H_{T};H_{T}^{rec};H_{T}^{gen}", 48, 0, 4200, 48, 0, 4200);
  book<TH2D>("H_T_rec_gen_2rec_muons_rebin", "H_{T};H_{T}^{rec};H_{T}^{gen}", 20, bins_from350, 20, bins_from350);
  book<TH2D>("H_T_rec_gen_2rec_ele_rebin", "H_{T};H_{T}^{rec};H_{T}^{gen}", 20, bins_from350, 20, bins_from350);

  is_mc = ctx.get("dataset_type") == "MC";
}


void LQToTopMuEfficiencyHists::fill(const Event & event){
  if(is_mc){
    // Don't forget to always use the weight when filling.
    double weight = event.weight;


    //calculate HT
    auto met = event.met->pt();
    double ht = 0.0;
    double ht_jets = 0.0;
    double ht_lep = 0.0;
    for(const auto & jet : *event.jets){
      ht_jets += jet.pt();
    }
    for(const auto & electron : *event.electrons){
      ht_lep += electron.pt();
    }
    for(const auto & muon : *event.muons){
      ht_lep += muon.pt();
    }

    ht = ht_lep + ht_jets + met;

    //calcutate parton HT
    double partonHT = 0;
    constexpr const int invalid_daughter = (unsigned short)(-1);
    //cout << endl << endl << "---------- CALCULATING PARTON HT ----------" << endl;
    for(const auto & gp : *event.genparticles){
      if(gp.daughter1() != invalid_daughter || gp.daughter2() != invalid_daughter) continue;
      // if we are here, it means we have a final state particle.
      // Add to HT in case it is a parton (quark -- including b but not top as tops are never final state particles -- or gluon -- or ele/mu -- or its respective neutrino).
      // Note that the exact HT definition depends on the madgraph configuration, but this
      // should cover the most common case.
      int id = abs(gp.pdgId());
      if((id >= 1 && id <= 5) || (id == 21) || (id>=11 && id <= 14)){
	partonHT += gp.pt();
	//cout << "pdgId: " << id << ", pt: " << gp.pt() << ", current HT: " << partonHT << endl;
      }
    }

    //ht = partonHT;

    /*
    //find right bin for rec ht re-binned
    int idx = 0;
    while(ht >= bins_from350[idx]) idx++;
    if(!(partonHT >= bins_from350[idx] && partonHT < bins_from350[idx+1])) return;
    */

    hist("H_T")->Fill(ht, weight);
    hist("H_T_rebin")->Fill(ht, weight);
    hist("Parton_H_T")->Fill(partonHT, weight);
    hist("Parton_H_T_rebin")->Fill(partonHT, weight);

    unsigned int Nmuons = event.muons->size();
    hist("N_mu")->Fill(Nmuons, weight);

    unsigned int Nele = event.electrons->size();
    hist("N_ele")->Fill(Nele, weight);

    int N_gen_ele = 0, N_gen_mu = 0;
    for(const auto & gp : *event.genparticles){
      int id = abs(gp.pdgId());
      if(id == 11) N_gen_ele++;
      if(id == 13) N_gen_mu++;
    }
    hist("h_N_gen_ele")->Fill(N_gen_ele,weight);
    hist("h_N_gen_mu")->Fill(N_gen_mu,weight);
  

    //fill HT for each event that contains at least 2 gen muons --> denominator
    if(N_gen_mu >= 2){
      hist("H_T_2gen_muons")->Fill(ht, weight);
      hist("H_T_2gen_muons_rebin")->Fill(ht, weight);
    }
    if(N_gen_ele >= 2){
      hist("H_T_2gen_ele")->Fill(ht,weight);
      hist("H_T_2gen_ele_rebin")->Fill(ht,weight);
    }

    //cases, where at least 2 muons have been reconstucted
    if(Nmuons >= 2){
      //search for at least two matched muons
      vector<double> v_dr_min;
      vector<bool> v_matched;
      int N_matched = 0;
      for(unsigned int i=0; i<Nmuons; i++){
	double dr_min = 999999;
	bool matched = false;
	for(unsigned int j=0; j<event.genparticles->size(); j++){
	  auto genp = event.genparticles->at(j);
	  if(abs(genp.pdgId()) == 13){
	    double dr = deltaR(event.muons->at(i),genp);
	    if(dr<dr_min){
	      dr_min = dr;
	    }
	  }
	}
	if(dr_min <= 0.1){
	  matched = true;
	  N_matched++;
	}
	else{
	  //if lepton is a fake, try to match it to a jet to fill faking-jet pt
	  //search for jets within 0.4
	  int idx_matching_jet = -1;
	  for(unsigned int j=0; j<event.jets->size(); j++){
	    double dr = deltaR(event.muons->at(i),event.jets->at(j));
	    if(dr < 0.4){
	      idx_matching_jet = j;
	    }
	  }

	  //-1 is filled in case no jet could be matched to the fake lepton
	  double faking_pt = -1;
	  if(idx_matching_jet > -1) faking_pt = event.jets->at(idx_matching_jet).pt();
	  hist("jets_faking_mu_pt")->Fill(faking_pt,weight);
	}
	v_dr_min.push_back(dr_min);
	v_matched.push_back(matched);
      }
      if(v_matched.size() != Nmuons) throw runtime_error("In LQToTopMuEfficiencyHists.cxx: v_matched does not contain as many entries as Nmuons");

      if(N_matched >= 2){
	hist("H_T_2matched_muons")->Fill(ht,weight);
	hist("H_T_2matched_muons_rebin")->Fill(ht,weight);
      }
      else{
	hist("H_T_not2matched_muons")->Fill(ht,weight);
	hist("H_T_not2matched_muons_rebin")->Fill(ht,weight);
      }
    }

    //cases, where at least 2 electrons have been reconstucted
    if(Nele >= 2){
      //search for at least two matched muons
      vector<double> v_dr_min;
      vector<bool> v_matched;
      int N_matched = 0;
      for(unsigned int i=0; i<Nele; i++){
	double dr_min = 999999;
	bool matched = false;
	for(unsigned int j=0; j<event.genparticles->size(); j++){
	  auto genp = event.genparticles->at(j);
	  if(abs(genp.pdgId()) == 11){
	    double dr = deltaR(event.electrons->at(i),genp);
	    if(dr<dr_min){
	      dr_min = dr;
	    }
	  }
	}
	if(dr_min <= 0.1){
	  matched = true;
	  N_matched++;
	}
	else{
	  //if lepton is a fake, try to match it to a jet to fill faking-jet pt
	  //search for jets within 0.4
	  int idx_matching_jet = -1;
	  for(unsigned int j=0; j<event.jets->size(); j++){
	    double dr = deltaR(event.electrons->at(i),event.jets->at(j));
	    if(dr < 0.4){
	      idx_matching_jet = j;
	    }
	  }

	  //-1 is filled in case no jet could be matched to the fake lepton
	  double faking_pt = -1;
	  if(idx_matching_jet > -1) faking_pt = event.jets->at(idx_matching_jet).pt();
	  hist("jets_faking_ele_pt")->Fill(faking_pt,weight);
	}
	v_dr_min.push_back(dr_min);
	v_matched.push_back(matched);
      }
      if(v_matched.size() != Nele) throw runtime_error("In LQToTopMuEfficiencyHists.cxx: v_matched does not contain as many entries as Nele");

      if(N_matched >= 2){
	hist("H_T_2matched_ele")->Fill(ht,weight);
	hist("H_T_2matched_ele_rebin")->Fill(ht,weight);
      }
      else{
	hist("H_T_not2matched_ele")->Fill(ht,weight);
	hist("H_T_not2matched_ele_rebin")->Fill(ht,weight);
      }
    }

    //simply fill ht, pt_leptons for N_lep >= 2 to be 100% consistent and not rely on main selections with triggers etc.
    if(Nmuons >= 2){
    //if(N_gen_mu >= 2){
      hist("H_T_2rec_muons")->Fill(ht,weight);
      hist("H_T_2rec_muons_rebin")->Fill(ht,weight);
      ((TH2D*)hist("H_T_rec_gen_2rec_muons"))->Fill(ht, partonHT ,weight);
      ((TH2D*)hist("H_T_rec_gen_2rec_muons_rebin"))->Fill(ht, partonHT ,weight);
      for(const auto & jet : *event.jets){
	double pt = jet.pt();
	hist("jets_2rec_mu_pt")->Fill(pt,weight);
      }
      for(const auto & mu : *event.muons){
	double pt = mu.pt();
	hist("mu_2rec_mu_pt")->Fill(pt,weight);
      }
    }
    if(Nele >= 2){
    //if(N_gen_ele >= 2){
     hist("H_T_2rec_ele")->Fill(ht,weight);
     hist("H_T_2rec_ele_rebin")->Fill(ht,weight);
      ((TH2D*)hist("H_T_rec_gen_2rec_ele"))->Fill(ht, partonHT ,weight);
      ((TH2D*)hist("H_T_rec_gen_2rec_ele_rebin"))->Fill(ht, partonHT ,weight);
     for(const auto & jet : *event.jets){
       double pt = jet.pt();
       hist("jets_2rec_ele_pt")->Fill(pt,weight);
      }
      for(const auto & ele : *event.electrons){
	double pt = ele.pt();
	hist("ele_2rec_ele_pt")->Fill(pt,weight);
      }
    }


  } //is_mc
} //Methode


LQToTopMuEfficiencyHists::~LQToTopMuEfficiencyHists(){}
