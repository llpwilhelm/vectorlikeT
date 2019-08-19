#include "../include/LQToTopMuHFLeptonHists.h"
#include "../include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"

#include <math.h>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuHFLeptonHists::LQToTopMuHFLeptonHists(Context & ctx, const string & dirname, ElectronId & EleIdIso_, MuonId & MuIdIso_): Hists(ctx, dirname), EleIdIso(EleIdIso_), MuIdIso(MuIdIso_){
  book<TH1F>("flavor", "flavor of non-isolated lepton; -1: #mu, +1: e", 3, -1.5, 1.5);  
  book<TH1F>("pt", "p_{T} of non-isolated lepton", 30, 0, 600); 
  book<TH1F>("iso_ele", "isolation of non-isolated electron", 50, 0, 0.5);
  book<TH1F>("iso_mu", "isolation of non-isolated muon", 50, 0, 0.5);
  book<TH1F>("drmin_lep_genb", "#Delta R_{min} (non-iso lep, gen b)", 100, 0, 5); 
  book<TH1F>("pt_b", "p_{T} of non-isolated lepton, matched to gen b", 30, 0, 600);   
  book<TH1F>("drmin_lep_genc", "#Delta R_{min} (non-iso lep, gen c)", 100, 0, 5); 
  book<TH1F>("pt_c", "p_{T} of non-isolated lepton, matched to gen c", 30, 0, 600);   
  book<TH1F>("drmin_lep_gentau", "#Delta R_{min} (non-iso lep, gen #tau)", 100, 0, 5); 
  book<TH1F>("pt_tau", "p_{T} of non-isolated lepton, matched to gen #tau", 30, 0, 600);  
  book<TH1F>("drmin_lep_genlep", "#Delta R_{min} (non-iso lep, gen lep)", 100, 0, 5); 
  book<TH1F>("pt_lep", "p_{T} of non-isolated lepton, matched to true lepton", 30, 0, 600);  
  book<TH1F>("pt_fake", "p_{T} of non-isolated lepton, unmatched", 30, 0, 600);
  book<TH1F>("pt_noprompt_notau", "p_{T} of non-isolated lepton, non-prompt, not from #tau decay", 30, 0, 600);
  book<TH1F>("pt_prompt_tau", "p_{T} of non-isolated lepton, prompt or from #tau decay", 30, 0, 600);
  
  
  is_mc = ctx.get("dataset_type") == "MC";

}


void LQToTopMuHFLeptonHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.

  double weight = event.weight;


  //Find non-iso lepton
  if(event.muons->size() > 1 && event.electrons->size() > 1) throw runtime_error("In LQToTopMuHFLeptonHists.cxx: There is more than one ele AND more than one muon. This should not be the case, please fix.");
  if(event.muons->size() < 1 && event.electrons->size() < 1) throw runtime_error("In LQToTopMuHFLeptonHists.cxx: There are no electrons or muons. This should not be the case, please fix.");
  if(!(event.muons->size() == 2 && event.electrons->size() ==1) && !(event.muons->size()==1 && event.electrons->size() == 2)) throw runtime_error("In LQToTopMuHFLeptonHists.cxx: Not 1mu&2ele or 2mu&1ele. Please fix.");
  bool ismu = event.muons->size() > 1;

  Particle lepton;
  double isolation = -1.;
  if(ismu){
    for(const auto & mu : *event.muons){
      if(MuIdIso(mu, event)) continue;
      lepton = mu;
      isolation = mu.relIso();
    }
  }
  else{
    for(const auto & ele : *event.electrons){
      if(EleIdIso(ele, event)) continue;
      lepton = ele;
      isolation = ele.relIsorho(event.rho);
    }
  }

  //Now 'lepton' is the non-isolated additional lepton

  if(ismu) hist("flavor")->Fill(-1., weight);
  else     hist("flavor")->Fill(1., weight);

  hist("pt")->Fill(lepton.pt(), weight);

  if(isolation == -1.) throw runtime_error("In LQToTopMuHFLeptonHists.cxx: isolation has not been filled, how can this be?");
  if(ismu) hist("iso_mu")->Fill(isolation,weight);
  else     hist("iso_ele")->Fill(isolation,weight);

  if(is_mc){
    double drmin_b = 999, drmin_c = 999, drmin_tau = 999, drmin_lep = 999;
    for(const auto & gp : *event.genparticles){
      double dr = deltaR(gp,lepton);
      int id = fabs(gp.pdgId());
      if(id == 5){
	if(dr < drmin_b) drmin_b = dr;
      }
      else if(id == 4){
	if(dr < drmin_c) drmin_c = dr;
      }
      else if(id == 15){
	if(dr < drmin_tau) drmin_tau = dr;
      }
      else if(id == 13 && ismu){
	if(dr < drmin_lep) drmin_lep = dr;
      }
      else if(id == 11 && !ismu){
	if(dr < drmin_lep) drmin_lep = dr;
      }
    }
    hist("drmin_lep_genb")->Fill(drmin_b, weight);
    hist("drmin_lep_genc")->Fill(drmin_c, weight);
    hist("drmin_lep_gentau")->Fill(drmin_tau, weight);
    hist("drmin_lep_genlep")->Fill(drmin_lep, weight);

    bool is_b = drmin_b < drmin_c && drmin_b < drmin_tau && drmin_b < drmin_lep && drmin_b < 0.2;
    bool is_c = drmin_c < drmin_b && drmin_c < drmin_tau && drmin_c < drmin_lep && drmin_c < 0.2;
    bool is_tau = drmin_tau < drmin_b && drmin_tau < drmin_c && drmin_tau < drmin_lep && drmin_tau < 0.2;
    bool is_lep = drmin_lep < drmin_b && drmin_lep < drmin_c && drmin_lep < drmin_tau && drmin_lep < 0.1;

    int sum = is_b + is_c + is_tau + is_lep;
    if(sum>1) throw runtime_error("In LQToTopMuHFLeptonHists.cxx: Non-iso lepton matched to more than one gen particle. Please fix.");

    if(is_b) hist("pt_b")->Fill(lepton.pt(), weight);
    else if(is_c) hist("pt_c")->Fill(lepton.pt(), weight);
    else if(is_tau) hist("pt_tau")->Fill(lepton.pt(), weight);
    else if(is_lep) hist("pt_lep")->Fill(lepton.pt(), weight);
    else hist("pt_fake")->Fill(lepton.pt(), weight);

    if(!is_lep && !is_tau) hist("pt_noprompt_notau")->Fill(lepton.pt(), weight);
    else hist("pt_prompt_tau")->Fill(lepton.pt(), weight);
  }

} //Methode



LQToTopMuHFLeptonHists::~LQToTopMuHFLeptonHists(){}
