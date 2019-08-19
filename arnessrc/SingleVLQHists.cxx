#include "../include/SingleVLQHists.h"
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

SingleVLQHists::SingleVLQHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  // jets
  book<TH1F>("N_jets", "N_{jets}", 21, -0.5, 20.5);
  book<TH1F>("pt_jet", "p_{T}^{jets} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3} [GeV]", 50, 0, 1500);
  book<TH1F>("eta_jet", "#eta^{jets}", 50, -2.5, 2.5);
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 50, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 50, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 50, -2.5, 2.5);
  book<TH1F>("phi_jet", "#phi^{jets}", 35, -3.5, 3.5);
  book<TH1F>("phi_jet1", "#phi^{jet 1}", 35, -3.5, 3.5);
  book<TH1F>("phi_jet2", "#phi^{jet 2}", 35, -3.5, 3.5);
  book<TH1F>("phi_jet3", "#phi^{jet 3}", 35, -3.5, 3.5);
  book<TH1F>("m_jet", "m^{jets}", 50, 0, 500);
  book<TH1F>("m_jet1", "m^{jet 1}", 50, 0, 500);
  book<TH1F>("m_jet2", "m^{jet 2}", 50, 0, 500);
  book<TH1F>("m_jet3", "m^{jet 3}", 50, 0, 500);
  book<TH1F>("csv_jet", "CSV^{jets}", 20, 0, 1);
  book<TH1F>("csv_jet1", "CSV^{jet 1}", 20, 0, 1);
  book<TH1F>("csv_jet2", "CSV^{jet 2}", 20, 0, 1);
  book<TH1F>("csv_jet3", "CSV^{jet 3}", 20, 0, 1);
  book<TH1F>("N_bJets_loose", "N_{jets}^{CSV loose}", 11, -0.5, 10.5);
  book<TH1F>("N_bJets_med", "N_{jets}^{CSV medium}", 11, -0.5, 10.5);
  book<TH1F>("N_bJets_tight", "N_{jets}^{CSV tight}", 11, -0.5, 10.5);


  book<TH1F>("N_mu", "N^{#mu}", 11, -0.5, 10.5);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_mu1", "p_{T}^{#mu 1} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_mu2", "p_{T}^{#mu 2} [GeV]", 50, 0, 1500);
  book<TH1F>("eta_mu", "#eta^{#mu}", 50, -2.5, 2.5);
  book<TH1F>("eta_mu1", "#eta^{#mu 1}", 50, -2.5, 2.5);
  book<TH1F>("eta_mu2", "#eta^{#mu 2}", 50, -2.5, 2.5);
  book<TH1F>("phi_mu", "#phi^{#mu}", 35, -3.5, 3.5);
  book<TH1F>("phi_mu1", "#phi^{#mu 1}", 35, 3.5, 3.5);
  book<TH1F>("phi_mu2", "#phi^{#mu 2}", 35, 3.5, 3.5);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_mu1", "#mu 1 rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_mu2", "#mu 2 rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_mu_rebin", "#mu rel. Iso ", 400, 0, 5);
  book<TH1F>("reliso_mu1_rebin", "#mu 1 rel. Iso ", 400, 0, 5);
  book<TH1F>("reliso_mu2_rebin", "#mu 2 rel. Iso ", 400, 0, 5);
  book<TH1F>("N_ele", "N^{e}", 11, -0.5, 10.5);
  book<TH1F>("pt_ele", "p_{T}^{e} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_ele1", "p_{T}^{e 1} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_ele2", "p_{T}^{e 2} [GeV]", 50, 0, 1500);
  book<TH1F>("eta_ele", "#eta^{e}", 50, -2.5, 2.5);
  book<TH1F>("eta_ele1", "#eta^{ele 1}", 50, -2.5, 2.5);
  book<TH1F>("eta_ele2", "#eta^{ele 2}", 50, -2.5, 2.5);
  book<TH1F>("phi_ele", "#phi^{e}", 35, -3.5, 3.5);
  book<TH1F>("phi_ele1", "#phi^{e 1}", 35, -3.5, 3.5);
  book<TH1F>("phi_ele2", "#phi^{e 2}", 35, -3.5, 3.5);
  book<TH1F>("reliso_ele", "e rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_ele1", "e 1 rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_ele2", "e 2 rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_ele_rebin", "e rel. Iso", 400, 0, 5);
  book<TH1F>("reliso_ele1_rebin", "e 1 rel. Iso", 400, 0, 5);
  book<TH1F>("reliso_ele2_rebin", "e 2 rel. Iso", 400, 0, 5);
  book<TH1F>("M_mumu", "M_{#mu#mu} [GeV]",75 , 0, 500);
  book<TH1F>("M_ee", "M_{ee} [GeV]",75 , 0, 500);


  book<TH1F>("NPV", "number of primary vertices", 91, -0.50, 90.5);
  book<TH1F>("MET", "missing E_{T} [GeV]", 50, 0, 7000);
  book<TH1F>("MET_rebin", "missing E_{T} [GeV]", 50, 0, 1500);
  book<TH1F>("MET_rebin2", "missing E_{T} [GeV]", 30, 0, 1500);
  book<TH1F>("MET_rebin3", "missing E_{T} [GeV]", 15, 0, 1500);
  book<TH1F>("ST", "S_{T} [GeV]", 50, 0, 7000);
  book<TH1F>("ST_rebin", "S_{T} [GeV]", 200, 0, 5000);
  book<TH1F>("ST_rebin2", "S_{T} [GeV]", 100, 0, 5000);
  book<TH1F>("ST_rebin3", "S_{T} [GeV]", 50, 0, 5000);
  book<TH1F>("STjets", "S_{T}^{jets} [GeV]", 50, 0, 7000);
  book<TH1F>("STjets_rebin", "S_{T}^{jets} [GeV]", 200, 0, 5000);
  book<TH1F>("STjets_rebin2", "S_{T}^{jets} [GeV]", 100, 0, 5000);
  book<TH1F>("STjets_rebin3", "S_{T}^{jets} [GeV]", 50, 0, 5000);
  book<TH1F>("STlep", "S_{T}^{lep} [GeV]", 50, 0, 7000);
  book<TH1F>("STlep_rebin", "S_{T}^{lep} [GeV]", 50, 0, 1500);
  book<TH1F>("STlep_rebin2", "S_{T}^{lep} [GeV]", 30, 0, 1500);
  book<TH1F>("STlep_rebin3", "S_{T}^{lep} [GeV]", 15, 0, 1500);

  book<TH1F>("M_Tprime", "M_{T} [GeV]", 300, 0, 3000);
  book<TH1F>("M_Tprime_rebin", "M_{T} [GeV]", 150, 0, 3000);
  book<TH1F>("M_Tprime_rebin2", "M_{T} [GeV]", 60, 0, 3000);
  vector<double> bins_limits = {0, 300, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1200, 3000};
  book<TH1F>("M_Tprime_rebinlimit", "M_{T} [GeV]", bins_limits.size()-1, &bins_limits[0]);
  book<TH1F>("chi2", "#chi^{2}", 100, 0, 200);
  book<TH1F>("chi2_rebin", "#chi^{2}", 50, 0, 200);
  book<TH1F>("chi2_rebin2", "#chi^{2}", 20, 0, 200);
  book<TH1F>("chi2_rebin3", "#chi^{2}", 50, 0, 50);
  book<TH1F>("chi2_rebin4", "#chi^{2}", 25, 0, 50);
  book<TH1F>("chi2_rebin5", "#chi^{2}", 10, 0, 50);
  book<TH1F>("chi2_rebin6", "#chi^{2}", 50, 0, 10);
  book<TH1F>("chi2_rebin7", "#chi^{2}", 25, 0, 10);
  book<TH1F>("chi2_rebin8", "#chi^{2}", 10, 0, 10);
  book<TH1F>("chi2_top", "#chi^{2}_{t}", 100, 0, 200);
  book<TH1F>("chi2_top_rebin", "#chi^{2}_{t}", 50, 0, 200);
  book<TH1F>("chi2_top_rebin2", "#chi^{2}_{t}", 20, 0, 200);
  book<TH1F>("chi2_top_rebin3", "#chi^{2}_{t}", 50, 0, 50);
  book<TH1F>("chi2_top_rebin4", "#chi^{2}_{t}", 25, 0, 50);
  book<TH1F>("chi2_top_rebin5", "#chi^{2}_{t}", 10, 0, 50);
  book<TH1F>("chi2_top_rebin6", "#chi^{2}_{t}", 50, 0, 10);
  book<TH1F>("chi2_top_rebin7", "#chi^{2}_{t}", 25, 0, 10);
  book<TH1F>("chi2_top_rebin8", "#chi^{2}_{t}", 10, 0, 10);
  book<TH1F>("chi2_higgs", "#chi^{2}_{H}", 100, 0, 200);
  book<TH1F>("chi2_higgs_rebin", "#chi^{2}_{H}", 50, 0, 200);
  book<TH1F>("chi2_higgs_rebin2", "#chi^{2}_{H}", 20, 0, 200);
  book<TH1F>("chi2_higgs_rebin3", "#chi^{2}_{H}", 50, 0, 50);
  book<TH1F>("chi2_higgs_rebin4", "#chi^{2}_{H}", 25, 0, 50);
  book<TH1F>("chi2_higgs_rebin5", "#chi^{2}_{H}", 10, 0, 50);
  book<TH1F>("chi2_higgs_rebin6", "#chi^{2}_{H}", 50, 0, 10);
  book<TH1F>("chi2_higgs_rebin7", "#chi^{2}_{H}", 25, 0, 10);
  book<TH1F>("chi2_higgs_rebin8", "#chi^{2}_{H}", 10, 0, 10);
  book<TH1F>("N_jets_higgs","number of jets in the higgs hypothesis",6,-0.5,5.5);
  book<TH1F>("M_higgs", "M_{H} [GeV]", 50, 0, 500);
  book<TH1F>("N_jets_top","number of jets in the top hypothesis",6,-0.5,5.5);
  book<TH1F>("M_top", "M_{t} [GeV]", 50, 0, 500);
  book<TH1F>("dR_bH_bH", "#Delta R (b_{H}, b_{H})", 60, 0, 3);
  book<TH1F>("dR_bt_Wt", "#Delta R (b_{t}, W_{t})", 60, 0, 3);
  book<TH1F>("dR_t_H", "#Delta R (t, H)", 60, 0, 3);
  book<TH1F>("dPhi_bH_bH", "#Delta #Phi (b_{H}, b_{H})", 60, 0, 3);
  book<TH1F>("dPhi_bt_Wt", "#Delta #Phi (b_{t}, W_{t})", 60, 0, 3);
  book<TH1F>("dPhi_t_H", "#Delta #Phi (t, H)", 60, 0, 3);

  book<TH1F>("sum_event_weights", "BinContent = sum(eventweights)", 1, 0.5, 1.5);



  //For MLQ reconstruction
  h_hyps           = ctx.get_handle<std::vector<SingleVLQReconstructionHypothesis>>("TprimeHypotheses");
  h_is_tprime_reco = ctx.get_handle<bool>("is_tprime_reco");
  is_mc = ctx.get("dataset_type") == "MC";

}


void SingleVLQHists::fill(const Event & event){

  double weight = event.weight;
  hist("sum_event_weights")->Fill(1., weight);


  /*
  █      ██ ███████ ████████ ███████
  █      ██ ██         ██    ██
  █      ██ █████      ██    ███████
  █ ██   ██ ██         ██         ██
  █  █████  ███████    ██    ███████
  */


  vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);

  for(unsigned int i=0; i<jets->size(); i++){
    hist("pt_jet")->Fill(jets->at(i).pt(),weight);
    hist("eta_jet")->Fill(jets->at(i).eta(),weight);
    hist("phi_jet")->Fill(jets->at(i).phi(),weight);
    hist("m_jet")->Fill(jets->at(i).v4().M(),weight);
    hist("csv_jet")->Fill(jets->at(i).btag_combinedSecondaryVertex(), weight);
    if(i==0){
      hist("pt_jet1")->Fill(jets->at(i).pt(),weight);
      hist("eta_jet1")->Fill(jets->at(i).eta(),weight);
      hist("phi_jet1")->Fill(jets->at(i).phi(),weight);
      hist("m_jet1")->Fill(jets->at(i).v4().M(),weight);
      hist("csv_jet1")->Fill(jets->at(i).btag_combinedSecondaryVertex(), weight);
    }
    else if(i==1){
      hist("pt_jet2")->Fill(jets->at(i).pt(),weight);
      hist("eta_jet2")->Fill(jets->at(i).eta(),weight);
      hist("phi_jet2")->Fill(jets->at(i).phi(),weight);
      hist("m_jet2")->Fill(jets->at(i).v4().M(),weight);
      hist("csv_jet2")->Fill(jets->at(i).btag_combinedSecondaryVertex(), weight);
    }
    else if(i==2){
      hist("pt_jet3")->Fill(jets->at(i).pt(),weight);
      hist("eta_jet3")->Fill(jets->at(i).eta(),weight);
      hist("phi_jet3")->Fill(jets->at(i).phi(),weight);
      hist("m_jet3")->Fill(jets->at(i).v4().M(),weight);
      hist("csv_jet3")->Fill(jets->at(i).btag_combinedSecondaryVertex(), weight);
    }
  }

  int Nbjets_loose = 0, Nbjets_medium = 0, Nbjets_tight = 0;
  CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  CSVBTag Btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  CSVBTag Btag_tight = CSVBTag(CSVBTag::WP_TIGHT);

  for (unsigned int i =0; i<jets->size(); i++) {
    if(Btag_loose(jets->at(i),event))  Nbjets_loose++;
    if(Btag_medium(jets->at(i),event)) Nbjets_medium++;
    if(Btag_tight(jets->at(i),event))  Nbjets_tight++;
  }

  hist("N_bJets_loose")->Fill(Nbjets_loose,weight);
  hist("N_bJets_med")->Fill(Nbjets_medium,weight);
  hist("N_bJets_tight")->Fill(Nbjets_tight,weight);


  /*
  ███    ███ ██    ██  ██████  ███    ██ ███████
  ████  ████ ██    ██ ██    ██ ████   ██ ██
  ██ ████ ██ ██    ██ ██    ██ ██ ██  ██ ███████
  ██  ██  ██ ██    ██ ██    ██ ██  ██ ██      ██
  ██      ██  ██████   ██████  ██   ████ ███████
  */


  vector<Muon>* muons = event.muons;
  int Nmuons = muons->size();
  hist("N_mu")->Fill(Nmuons, weight);

  for(int i=0; i<Nmuons; i++){

    hist("pt_mu")->Fill(muons->at(i).pt(),weight);
    hist("eta_mu")->Fill(muons->at(i).eta(),weight);
    hist("phi_mu")->Fill(muons->at(i).phi(),weight);
    hist("reliso_mu")->Fill(muons->at(i).relIso(),weight);
    hist("reliso_mu_rebin")->Fill(muons->at(i).relIso(),weight);
    if(i==0){
      hist("pt_mu1")->Fill(muons->at(i).pt(),weight);
      hist("eta_mu1")->Fill(muons->at(i).eta(),weight);
      hist("phi_mu1")->Fill(muons->at(i).phi(),weight);
      hist("reliso_mu1")->Fill(muons->at(i).relIso(),weight);
      hist("reliso_mu1_rebin")->Fill(muons->at(i).relIso(),weight);
    }
    else if(i==1){
      hist("pt_mu2")->Fill(muons->at(i).pt(),weight);
      hist("eta_mu2")->Fill(muons->at(i).eta(),weight);
      hist("phi_mu2")->Fill(muons->at(i).phi(),weight);
      hist("reliso_mu2")->Fill(muons->at(i).relIso(),weight);
      hist("reliso_mu2_rebin")->Fill(muons->at(i).relIso(),weight);
    }
  }

  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j <= i) continue;
      hist("M_mumu")->Fill((muons->at(i).v4() + muons->at(j).v4()).M() ,weight);
    }
  }

  /*
  ███████ ██      ███████  ██████ ████████ ██████   ██████  ███    ██ ███████
  ██      ██      ██      ██         ██    ██   ██ ██    ██ ████   ██ ██
  █████   ██      █████   ██         ██    ██████  ██    ██ ██ ██  ██ ███████
  ██      ██      ██      ██         ██    ██   ██ ██    ██ ██  ██ ██      ██
  ███████ ███████ ███████  ██████    ██    ██   ██  ██████  ██   ████ ███████
  */


  vector<Electron>* electrons = event.electrons;
  int Nelectrons = electrons->size();
  hist("N_ele")->Fill(Nelectrons, weight);

  for(int i=0; i<Nelectrons; i++){
    hist("pt_ele")->Fill(electrons->at(i).pt(),weight);
    hist("eta_ele")->Fill(electrons->at(i).eta(),weight);
    hist("phi_ele")->Fill(electrons->at(i).phi(),weight);
    hist("reliso_ele")->Fill(electrons->at(i).relIso(),weight);
    hist("reliso_ele_rebin")->Fill(electrons->at(i).relIso(),weight);
    if(i==0){
      hist("pt_ele1")->Fill(electrons->at(i).pt(),weight);
      hist("eta_ele1")->Fill(electrons->at(i).eta(),weight);
      hist("phi_ele1")->Fill(electrons->at(i).phi(),weight);
      hist("reliso_ele1")->Fill(electrons->at(i).relIso(),weight);
      hist("reliso_ele1_rebin")->Fill(electrons->at(i).relIso(),weight);
    }
    else if(i==1){
      hist("pt_ele2")->Fill(electrons->at(i).pt(),weight);
      hist("eta_ele2")->Fill(electrons->at(i).eta(),weight);
      hist("phi_ele2")->Fill(electrons->at(i).phi(),weight);
      hist("reliso_ele2")->Fill(electrons->at(i).relIso(),weight);
      hist("reliso_ele2_rebin")->Fill(electrons->at(i).relIso(),weight);
    }
  }

  for(int i=0; i<Nelectrons; i++){
    for(int j=0; j<Nelectrons; j++){
      if(j <= i) continue;
      hist("M_ee")->Fill((electrons->at(i).v4() + electrons->at(j).v4()).M() ,weight);
    }
  }



  /*
  ██████  ███████ ███    ██ ███████ ██████   █████  ██
  ██      ██      ████   ██ ██      ██   ██ ██   ██ ██
  ██  ███ █████   ██ ██  ██ █████   ██████  ███████ ██
  ██   ██ ██      ██  ██ ██ ██      ██   ██ ██   ██ ██
  ██████  ███████ ██   ████ ███████ ██   ██ ██   ██ ███████
  */



  int Npvs = event.pvs->size();
  hist("NPV")->Fill(Npvs, weight);

  double met = event.met->pt();
  double st = 0., st_jets = 0., st_lep = 0.;
  for(const auto & jet : *jets)           st_jets += jet.pt();
  for(const auto & electron : *electrons) st_lep += electron.pt();
  for(const auto & muon : *muons)         st_lep += muon.pt();
  st = st_jets + st_lep + met;

  hist("MET")->Fill(met, weight);
  hist("MET_rebin")->Fill(met, weight);
  hist("MET_rebin2")->Fill(met, weight);
  hist("MET_rebin3")->Fill(met, weight);
  hist("ST")->Fill(st, weight);
  hist("ST_rebin")->Fill(st, weight);
  hist("ST_rebin2")->Fill(st, weight);
  hist("ST_rebin3")->Fill(st, weight);
  hist("STjets")->Fill(st_jets, weight);
  hist("STjets_rebin")->Fill(st_jets, weight);
  hist("STjets_rebin2")->Fill(st_jets, weight);
  hist("STjets_rebin3")->Fill(st_jets, weight);
  hist("STlep")->Fill(st_lep, weight);
  hist("STlep_rebin")->Fill(st_lep, weight);
  hist("STlep_rebin2")->Fill(st_lep, weight);
  hist("STlep_rebin3")->Fill(st_lep, weight);



  /*
  █ ████████ ██████  ██████  ██ ███    ███ ███████
  █    ██    ██   ██ ██   ██ ██ ████  ████ ██
  █    ██    ██████  ██████  ██ ██ ████ ██ █████
  █    ██    ██      ██   ██ ██ ██  ██  ██ ██
  █    ██    ██      ██   ██ ██ ██      ██ ███████
  */





  bool is_tprime_reco = event.get(h_is_tprime_reco);
  if(is_tprime_reco){
    std::vector<SingleVLQReconstructionHypothesis> hyps = event.get(h_hyps);
    const SingleVLQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
    double chi2 = hyp->discriminator("Chi2");
    double chi2_top = hyp->discriminator("Chi2_top");
    double chi2_higgs = hyp->discriminator("Chi2_higgs");
    hist("chi2")->Fill(chi2,weight);
    hist("chi2_rebin")->Fill(chi2,weight);
    hist("chi2_rebin2")->Fill(chi2,weight);
    hist("chi2_rebin3")->Fill(chi2,weight);
    hist("chi2_rebin4")->Fill(chi2,weight);
    hist("chi2_rebin5")->Fill(chi2,weight);
    hist("chi2_rebin6")->Fill(chi2,weight);
    hist("chi2_rebin7")->Fill(chi2,weight);
    hist("chi2_rebin8")->Fill(chi2,weight);
    hist("chi2_top")->Fill(chi2_top,weight);
    hist("chi2_top_rebin")->Fill(chi2_top,weight);
    hist("chi2_top_rebin2")->Fill(chi2_top,weight);
    hist("chi2_top_rebin3")->Fill(chi2_top,weight);
    hist("chi2_top_rebin4")->Fill(chi2_top,weight);
    hist("chi2_top_rebin5")->Fill(chi2_top,weight);
    hist("chi2_top_rebin6")->Fill(chi2_top,weight);
    hist("chi2_top_rebin7")->Fill(chi2_top,weight);
    hist("chi2_top_rebin8")->Fill(chi2_top,weight);
    hist("chi2_higgs")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin2")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin3")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin4")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin5")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin6")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin7")->Fill(chi2_higgs,weight);
    hist("chi2_higgs_rebin8")->Fill(chi2_higgs,weight);

    double m_tprime = 0.;
    double m_top = 0.;
    double m_higgs = 0.;

    if( (hyp->Tprime_v4()).isTimelike() ) {m_tprime = (hyp->Tprime_v4()).M();}
    else {m_tprime = sqrt( -(hyp->Tprime_v4()).mass2());}

    if( (hyp->toplep_v4()).isTimelike() ) {m_top = (hyp->toplep_v4()).M();}
    else {m_top = sqrt( -(hyp->toplep_v4()).mass2());}

    if( (hyp->higgs_v4()).isTimelike() ) {m_higgs = (hyp->higgs_v4()).M();}
    else {m_higgs = sqrt( -(hyp->higgs_v4()).mass2());}

    hist("M_Tprime")->Fill(m_tprime,weight);
    hist("M_Tprime_rebin")->Fill(m_tprime,weight);
    hist("M_Tprime_rebin2")->Fill(m_tprime,weight);
    hist("M_Tprime_rebinlimit")->Fill(m_tprime,weight);
    hist("M_higgs")->Fill(m_higgs,weight);
    hist("M_top")->Fill(m_top,weight);

    double n_jets_higgs = hyp->higgs_jets().size();
    double n_jets_top = hyp->toplep_jets().size();
    hist("N_jets_higgs")->Fill(n_jets_higgs,weight);
    hist("N_jets_top")->Fill(n_jets_top,weight);

    if(n_jets_higgs != 2 || n_jets_top != 1) throw runtime_error("Incorrect number of b-tagged jets in either H or t hypothesis (or both)");
    LorentzVector W = (hyp->lepton_v4() + hyp->neutrino_v4());
    float dR_bH_bH = deltaR(hyp->higgs_jets().at(0),hyp->higgs_jets().at(1));
    float dR_bt_Wt = deltaR(W, hyp->toplep_jets().at(0));
    float dR_t_H = deltaR(hyp->toplep_v4(), hyp->higgs_v4());
    float dPhi_bH_bH = deltaPhi(hyp->higgs_jets().at(0),hyp->higgs_jets().at(1));
    float dPhi_bt_Wt = deltaPhi(W, hyp->toplep_jets().at(0));
    float dPhi_t_H = deltaPhi(hyp->toplep_v4(), hyp->higgs_v4());
    hist("dR_bH_bH")->Fill(dR_bH_bH, weight);
    hist("dR_bt_Wt")->Fill(dR_bt_Wt, weight);
    hist("dR_t_H")->Fill(dR_t_H, weight);
    hist("dPhi_bH_bH")->Fill(dPhi_bH_bH, weight);
    hist("dPhi_bt_Wt")->Fill(dPhi_bt_Wt, weight);
    hist("dPhi_t_H")->Fill(dPhi_t_H, weight);
  }

} //Methode



SingleVLQHists::~SingleVLQHists(){}
