#include "../include/LQToTopMuFakeRateHists.h"
#include "../include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuFakeRateHists::LQToTopMuFakeRateHists(Context & ctx, const string & dirname, bool is_ele_): Hists(ctx, dirname), is_ele(is_ele_){
  // book all histograms here
  // jets



  book<TH1F>("pt_alljets", "p_{T}^{jets} [GeV]", 45, 0, 1500);   
  book<TH1F>("pt_matchingjets", "p_{T}^{jets} [GeV]", 45, 0, 1500);  
  book<TH1F>("pt_matchingjets_fakeele", "p_{T}^{jets} [GeV]", 45, 0, 1500);   
  book<TH1F>("pt_matchingjets_realele", "p_{T}^{jets} [GeV]", 45, 0, 1500); 
  book<TH1F>("pt_matchingjets_promptele", "p_{T}^{jets} [GeV]", 45, 0, 1500); 
  book<TH1F>("pt_matchingjets_realnonpromptele", "p_{T}^{jets} [GeV]", 45, 0, 1500); 
  double bins_jetpt[4] = {0,100,200,800};
  book<TH1F>("pt_alljets_rebin", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin", "p_{T}^{jets} [GeV]", 3, bins_jetpt);  
  book<TH1F>("pt_matchingjets_rebin_fakeele", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin_realele", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin_promptele", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin_realnonpromptele", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 

  book<TH1F>("ele_type", "Electron types", 5, -2.5, 2.5); 
  book<TH1F>("N_faking_jets", "Weighted number of faking jets", 1, -0.5, 0.5);   
  book<TH1F>("N_fake_ele", "Weighted number of fake electrons", 1, -0.5, 0.5);   




  book<TH1F>("pt_matchingjets_mu", "p_{T}^{jets} [GeV]", 45, 0, 1500);  
  book<TH1F>("pt_matchingjets_fakemu", "p_{T}^{jets} [GeV]", 45, 0, 1500);   
  book<TH1F>("pt_matchingjets_realmu", "p_{T}^{jets} [GeV]", 45, 0, 1500); 
  book<TH1F>("pt_matchingjets_promptmu", "p_{T}^{jets} [GeV]", 45, 0, 1500); 
  book<TH1F>("pt_matchingjets_realnonpromptmu", "p_{T}^{jets} [GeV]", 45, 0, 1500); 

  book<TH1F>("pt_matchingjets_mu_rebin", "p_{T}^{jets} [GeV]", 3, bins_jetpt);  
  book<TH1F>("pt_matchingjets_rebin_fakemu", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin_realmu", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin_promptmu", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 
  book<TH1F>("pt_matchingjets_rebin_realnonpromptmu", "p_{T}^{jets} [GeV]", 3, bins_jetpt); 

  book<TH1F>("mu_type", "Muon types", 5, -2.5, 2.5); 
  book<TH1F>("N_faking_jets_mu", "Weighted number of faking jets", 1, -0.5, 0.5);   
  book<TH1F>("N_fake_mu", "Weighted number of fake muons", 1, -0.5, 0.5); 
  book<TH1F>("dR_min_fakemu_gen_b", "dR_{min}(fake mu, gen b)", 100,0,5);

  book<TH1F>("Int_alljets", "Weighted number of jets", 1, -0.5, 0.5);  
  book<TH1F>("Int_matchingjets", "Weighted number of jets matching a #mu", 1, -0.5, 0.5);  
  book<TH1F>("Int_matchingjets_fakemu", "Weighted number of jets matching a faked #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_matchingjets_realmu", "Weighted number of jets matching a real #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_matchingjets_promptmu", "Weighted number of jets matching a prompt #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_matchingjets_realnonpromptmu", "Weighted number of jets matching a real, non-prompt #mu", 1, -0.5, 0.5); 

  book<TH1F>("Int_events_denom", "Weighted number of events with 0 or 1 #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_events_1mu", "Weighted number of events with == 1 #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_events_1mu_real", "Weighted number of events with == 1 real #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_events_1mu_fake", "Weighted number of events with == 1 fake #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_events_1mu_prompt", "Weighted number of events with == 1 prompt #mu", 1, -0.5, 0.5); 
  book<TH1F>("Int_events_1mu_realnonprompt", "Weighted number of events with == 1 real, non-prompt #mu", 1, -0.5, 0.5); 
 

  is_mc = ctx.get("dataset_type") == "MC";
  ZEE_finder.reset(new ZEEFinder());

}


void LQToTopMuFakeRateHists::fill(const Event & event){
  double weight = event.weight;


  if(event.electrons->size() >= 2){
    pair<unsigned int,unsigned int> best_ele = ZEE_finder->search(event);

    //contains booleans in the same order as event.electrons, specifying if the electron is part of the 'best combination'
    vector<bool> is_best_ele;
    for(unsigned int i=0; i<event.electrons->size(); i++){
      bool is_best = false;
      if(i == best_ele.first || i == best_ele.second) is_best = true;
      is_best_ele.push_back(is_best);
    }

    //contains booleans in the same order as event.electrons, specifying if the electron can be matched to a gen-electron, for MC only
    vector<bool> is_real_ele, is_prompt_ele;
    /*
    for(unsigned int i=0; i<event.electrons->size(); i++){
      bool is_real = false;
      if(is_mc){
	for(unsigned int j=0; j<event.genparticles->size(); j++){
	  //consider only electrons
	  if(fabs(event.genparticles->at(j).pdgId()) != 11) continue;
	  
	  //try to match
	  if(deltaR(event.genparticles->at(j), event.electrons->at(i)) < 0.1) is_real = true;
	}
      }
    */
      unsigned int n_genele = 0;
      if(is_mc){
	for(const auto & gp : *event.genparticles){
	  if(fabs(gp.pdgId()) != 11) continue;
	  n_genele++;
	}
      }
      else{
	//for data, the vector also needs to be filled (un-clever nesting of conditions later on)
	for(unsigned int i=0; i<event.electrons->size(); i++){
	  is_real_ele.push_back(false);
	  is_prompt_ele.push_back(false);
	}
      }

      if(n_genele == event.electrons->size()){
	for(unsigned int i=0; i<event.electrons->size(); i++){
	  is_real_ele.push_back(true);
	  is_prompt_ele.push_back(true);
	}
      }
      else{
	if(is_mc){
	  //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
	  unsigned int n_matched_to_ele = 0, n_matched_to_taus = 0, n_matched_to_b = 0, n_matched_to_c = 0;
	  int i=0;
	  for(const auto & ele : *event.electrons){
	    bool is_matched = false;
	    bool is_prompt = false;
	    for(const auto & gp : *event.genparticles){
	      if(fabs(gp.pdgId()) == 11){
		if(deltaR(gp,ele) < 0.1 && !is_matched){
		  is_matched = true;
		  is_prompt = true;
		  n_matched_to_ele++;
		}
	      }
	      else if(fabs(gp.pdgId()) == 15){ 
		if(deltaR(gp,ele) < 0.2 && !is_matched){
		  is_matched = true;
		  n_matched_to_taus++;
		}
	      }
	      else if(fabs(gp.pdgId()) == 5){ 
		if(deltaR(gp,ele) < 0.2 && !is_matched){
		  is_matched = true;
		  n_matched_to_b++;
		}
	      }
	      else if(fabs(gp.pdgId()) == 4){ 
		if(deltaR(gp,ele) < 0.2 && !is_matched){
		  is_matched = true;
		  n_matched_to_c++;
		}
	      }
	    }
	    //cout << "This muon is matched: " << is_matched << endl;
	    is_real_ele.push_back(is_matched);
	    is_prompt_ele.push_back(is_prompt);
	    if(is_mc && !is_matched && !is_best_ele[i]) hist("N_fake_ele")->Fill(0.,weight);
	    i++;
	  }
	}
      }
      //}

    
    //The same for muons
      vector<bool> is_real_muon, is_prompt_muon;
    if(event.muons->size() >= 1){  
      unsigned int n_genmu = 0;
      if(is_mc){
	for(const auto & gp : *event.genparticles){
	  if(fabs(gp.pdgId()) != 13) continue;
	  n_genmu++;
	}
      }
      else{
	//for data, the vector also needs to be filled (un-clever nesting of conditions later on)
	for(unsigned int i=0; i<event.muons->size(); i++){
	  is_real_muon.push_back(false);
	  is_prompt_muon.push_back(false);
	}
      }

      if(n_genmu == event.muons->size()){
	for(unsigned int i=0; i<event.muons->size(); i++){
	  is_real_muon.push_back(true);
	  is_prompt_muon.push_back(true);
	}
      }
      else{
	if(is_mc){
	  //if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
	  unsigned int n_matched_to_muons = 0, n_matched_to_taus = 0, n_matched_to_b = 0, n_matched_to_c = 0;
	  for(const auto & mu : *event.muons){
	    bool is_matched = false;
	    bool is_prompt = false;
	    for(const auto & gp : *event.genparticles){
	      if(fabs(gp.pdgId()) == 13){
		if(deltaR(gp,mu) < 0.1 && !is_matched){
		  is_matched = true;
		  is_prompt = true;
		  n_matched_to_muons++;
		}
	      }
	      else if(fabs(gp.pdgId()) == 15){ 
		if(deltaR(gp,mu) < 0.2 && !is_matched){
		  is_matched = true;
		  n_matched_to_taus++;
		}
	      }
	      else if(fabs(gp.pdgId()) == 5){ 
		if(deltaR(gp,mu) < 0.2 && !is_matched){
		  is_matched = true;
		  n_matched_to_b++;
		}
	      }
	      else if(fabs(gp.pdgId()) == 4){ 
		if(deltaR(gp,mu) < 0.2 && !is_matched){
		  is_matched = true;
		  n_matched_to_c++;
		}
	      }
	    }
	    //cout << "This muon is matched: " << is_matched << endl;
	    is_real_muon.push_back(is_matched);
	    is_prompt_muon.push_back(is_prompt);
	  }
	}
      }
    }






    if(is_mc){
      for(unsigned int i=0; i<event.muons->size(); i++){
	if(!is_real_muon[i]){
	  double dr_min = 99999;
	  for(const auto & gp : *event.genparticles){
	    if(fabs(gp.pdgId()) != 5) continue;
	    double dr = deltaR(event.muons->at(i),gp);
	    if(dr<dr_min) dr_min = dr;
	  }
	  hist("dR_min_fakemu_gen_b")->Fill(dr_min,weight);	
	}
      }
    }

    //fill ele-type histogram
    if(is_mc){
      for(unsigned int i=0; i<event.electrons->size(); i++){
	if(is_real_ele[i] && is_best_ele[i])       hist("ele_type")->Fill(-2,weight);
	else if(is_real_ele[i] && !is_best_ele[i]) hist("ele_type")->Fill(-1,weight);
	else if(!is_real_ele[i] && is_best_ele[i]) hist("ele_type")->Fill(1,weight);
	else if(!is_real_ele[i] && !is_best_ele[i])hist("ele_type")->Fill(2,weight);
	else throw runtime_error("In LQToTopMuFakeRateHists.cxx: Uncovered combination of real and best electrons.");
      }
    }

    //fill muon-type histogram
    if(is_mc){
      if(event.muons->size() >= 1){
	for(unsigned int i=0; i<event.muons->size(); i++){
	  if(is_real_muon[i])  hist("mu_type")->Fill(-1,weight);
	  else if(!is_real_muon[i]) hist("mu_type")->Fill(2,weight);
	}
      }
    }


    int idx_matched_jet = -1;
    int idx_matched_jet_fake_ele = -1;
    int idx_matched_jet_real_ele = -1;
    int idx_matched_jet_prompt_ele = -1;
    int idx_matched_jet_realnonprompt_ele = -1;
    int idx_matched_jet_mu = -1;
    int idx_matched_jet_fake_mu = -1;
    int idx_matched_jet_real_mu = -1;
    int idx_matched_jet_prompt_mu = -1;
    int idx_matched_jet_realnonprompt_mu = -1;
    if(is_ele){
      //make sure to have a max. of 3 electrons in the event. Otherwise this part will not work as intended.
      //find the ele that is not in the best pair
      double dr_min = 999;
      for(unsigned int i=0; i<event.electrons->size(); i++){
	if(is_best_ele[i]) continue;
      
	//find jet matching this ele within 0.4. If ambiguous, choose closest jet.
	for(unsigned int j=0; j<event.jets->size(); j++){
	  double dr = deltaR(event.electrons->at(i), event.jets->at(j));
	  if(dr <= 0.4){
	    if(event.electrons->size() > 3) throw runtime_error("In LQToTopMuFakeRateHists.cxx: More than 3 electrons in the event when identifying THE ONE faking jet. In case of >1 fake electrons, there should also be >1 faking jet. Modifiy the procedure or cut on Nele <= 3.");
	    if(dr < dr_min){
	      idx_matched_jet = j;
	      dr_min = dr;
	    
	      //if that ele is truly a fake (on gen-lvl), fill separate histogram to subtract from data later on (MC only)
	      if(!is_real_ele[i] && is_mc)     idx_matched_jet_fake_ele = j;
	      else if(is_real_ele[i] && is_mc) idx_matched_jet_real_ele = j;

	      if(is_prompt_ele[i] && is_mc)                         idx_matched_jet_prompt_ele = j;
	      else if(!is_prompt_ele[i] && is_real_ele[i] && is_mc) idx_matched_jet_realnonprompt_ele = j;
	    }
	  }
	}
      }
    }
    else{
      //make sure to have a max of 1 muon in the event, otherwise this part will not work as intended
      //same for muons
      if(event.muons->size() >= 1){
	double dr_min_mu = 999;
	for(unsigned int i=0; i<event.muons->size(); i++){
      
	  //find jet matching this ele within 0.4. If ambiguous, choose closest jet.
	  for(unsigned int j=0; j<event.jets->size(); j++){
	    double dr = deltaR(event.muons->at(i), event.jets->at(j));
	    if(dr <= 0.4){
	      if(event.muons->size() > 1) throw runtime_error("In LQToTopMuFakeRateHists.cxx: More than 1 muon in the event when identifying THE ONE faking jet. Procedure only works in cases of == 1 muon.");
	      if(dr < dr_min_mu){
		idx_matched_jet_mu = j;
		dr_min_mu = dr;
	    
		//if that mu is truly a fake (on gen-lvl), fill separate histogram to subtract from data later on (MC only)
		if(!is_real_muon[i] && is_mc)     idx_matched_jet_fake_mu = j;
		else if(is_real_muon[i] && is_mc) idx_matched_jet_real_mu = j;

		if(is_prompt_muon[i] && is_mc)                          idx_matched_jet_prompt_mu = j;
		else if(!is_prompt_muon[i] && is_real_muon[i] && is_mc) idx_matched_jet_realnonprompt_mu = j;
	      }
	    }
	  }
	}
      }

      if(idx_matched_jet_mu != -1){
	hist("pt_matchingjets_mu")->Fill(event.jets->at(idx_matched_jet_mu).pt(),weight);
	hist("pt_matchingjets_mu_rebin")->Fill(event.jets->at(idx_matched_jet_mu).pt(),weight); 
	hist("Int_matchingjets")->Fill(0.,weight);             

	if(idx_matched_jet_fake_mu != -1){
	  hist("pt_matchingjets_rebin_fakemu")->Fill(event.jets->at(idx_matched_jet_fake_mu).pt(),weight);
	  hist("pt_matchingjets_fakemu")->Fill(event.jets->at(idx_matched_jet_fake_mu).pt(),weight);
	  hist("Int_matchingjets_fakemu")->Fill(0.,weight);
	  hist("N_faking_jets_mu")->Fill(0.,weight);
	}
	else if(idx_matched_jet_real_mu != -1){
	  hist("pt_matchingjets_rebin_realmu")->Fill(event.jets->at(idx_matched_jet_real_mu).pt(),weight);
	  hist("pt_matchingjets_realmu")->Fill(event.jets->at(idx_matched_jet_real_mu).pt(),weight);
	  hist("Int_matchingjets_realmu")->Fill(0.,weight);
	}
	else if(is_mc) throw runtime_error("In LQToTopMuFakeRateHists.cxx: This MC-mu seems to be neither real nor fake. What is it?");

	if(idx_matched_jet_prompt_mu != -1){
	  hist("pt_matchingjets_rebin_promptmu")->Fill(event.jets->at(idx_matched_jet_prompt_mu).pt(),weight);
	  hist("pt_matchingjets_promptmu")->Fill(event.jets->at(idx_matched_jet_prompt_mu).pt(),weight);
	  hist("Int_matchingjets_promptmu")->Fill(0.,weight);
	}
	else if(idx_matched_jet_realnonprompt_mu != -1){
	  hist("pt_matchingjets_rebin_realnonpromptmu")->Fill(event.jets->at(idx_matched_jet_realnonprompt_mu).pt(),weight);
	  hist("pt_matchingjets_realnonpromptmu")->Fill(event.jets->at(idx_matched_jet_realnonprompt_mu).pt(),weight);
	  hist("Int_matchingjets_realnonpromptmu")->Fill(0.,weight);
	}
      }
    }

    
    if(idx_matched_jet != -1){
      hist("pt_matchingjets")->Fill(event.jets->at(idx_matched_jet).pt(),weight);
      hist("pt_matchingjets_rebin")->Fill(event.jets->at(idx_matched_jet).pt(),weight);    

      if(idx_matched_jet_fake_ele != -1){
	hist("pt_matchingjets_rebin_fakeele")->Fill(event.jets->at(idx_matched_jet_fake_ele).pt(),weight);
	hist("pt_matchingjets_fakeele")->Fill(event.jets->at(idx_matched_jet_fake_ele).pt(),weight);
	hist("N_faking_jets")->Fill(0.,weight);
      }
      else if(idx_matched_jet_real_ele != -1){
	hist("pt_matchingjets_rebin_realele")->Fill(event.jets->at(idx_matched_jet_real_ele).pt(),weight);
	hist("pt_matchingjets_realele")->Fill(event.jets->at(idx_matched_jet_real_ele).pt(),weight);
      }
      else if(is_mc) throw runtime_error("In LQToTopMuFakeRateHists.cxx: This MC-ele seems to be neither real nor fake. What is it?");

      if(idx_matched_jet_prompt_ele != -1){
	hist("pt_matchingjets_rebin_promptele")->Fill(event.jets->at(idx_matched_jet_prompt_ele).pt(),weight);
	hist("pt_matchingjets_promptele")->Fill(event.jets->at(idx_matched_jet_prompt_ele).pt(),weight);
      }
      else if(idx_matched_jet_realnonprompt_ele != -1){
	hist("pt_matchingjets_rebin_realnonpromptele")->Fill(event.jets->at(idx_matched_jet_realnonprompt_ele).pt(),weight);
	hist("pt_matchingjets_realnonpromptele")->Fill(event.jets->at(idx_matched_jet_realnonprompt_ele).pt(),weight);
      }
    }




    //Fill pt of all jets for the denominator, again only for at least 2 electrons
    for(auto & jet : *event.jets){
      hist("pt_alljets")->Fill(jet.pt(),weight);
      hist("pt_alljets_rebin")->Fill(jet.pt(),weight);
      hist("Int_alljets")->Fill(0.,weight);
    }

    if(!is_ele){
      hist("Int_events_denom")->Fill(0.,weight);

      if(event.muons->size() == 1) hist("Int_events_1mu")->Fill(0.,weight);
      else if(event.muons->size() != 0) throw runtime_error("In LQToTopMuFakeRateHists.cxx: There is more than 1 muon in the event. Fix this.");

      if(is_mc){
	if(event.muons->size() == 1){
	  
	  if(!is_real_muon[0]) hist("Int_events_1mu_fake")->Fill(0.,weight);
	  else                 hist("Int_events_1mu_real")->Fill(0.,weight);

	  if(is_prompt_muon[0])                          hist("Int_events_1mu_prompt")->Fill(0.,weight);
	  else if(!is_prompt_muon[0] && is_real_muon[0]) hist("Int_events_1mu_realnonprompt")->Fill(0.,weight);

	}
      }
    }

  }
} //Methode



LQToTopMuFakeRateHists::~LQToTopMuFakeRateHists(){}
