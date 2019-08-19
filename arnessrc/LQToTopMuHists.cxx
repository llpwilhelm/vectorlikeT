#include "../include/LQToTopMuHists.h"
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

LQToTopMuHists::LQToTopMuHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);  
  book<TH1F>("eta_jets", "#eta^{jets}", 25, -2.5, 2.5);
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);
  book<TH1F>("pt_jets", "p_{T}^{jets} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet1", "p_{T}^{jet 1} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet2", "p_{T}^{jet 2} [GeV]", 50, 20, 1500);
  book<TH1F>("pt_jet3", "p_{T}^{jet 3} [GeV]", 50, 20, 1500);
  book<TH1F>("N_bJets_loose", "N_{Bjets}^{loose}", 11, -0.5, 10.5);
  book<TH1F>("N_bJets_med", "N_{Bjets}^{medium}", 11, -0.5, 10.5);
  book<TH1F>("N_bJets_tight", "N_{Bjets}^{tight}", 11, -0.5, 10.5);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_mu_zoom", "p_{T}^{#mu} [GeV]", 34, 0, 1020);
  book<TH1F>("eta_mu", "#eta^{#mu}", 25, -2.5, 2.5);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);
  book<TH1F>("reliso_mu_rebin", "#mu rel. Iso ", 400, 0, 5);
  book<TH1F>("relisorho_ele", "e rel. Iso (eff. area corrected)", 40, 0, 0.5);
  book<TH1F>("relisorho_ele_rebin", "e rel. Iso (eff. area corrected)", 400, 0, 5);
  book<TH1F>("M_mumu", "M_#mu#mu [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("M_ee", "M_{ee} [GeV^{2}]",50 , 0, 1000);
  book<TH1F>("Pt_mu_sum", "#Sum p_{T}^{#mu} [GeV]", 50, 0, 7000);
  double bins_HTlept_low[6] = {0, 300, 600, 900, 1200, 7000};
  book<TH1F>("Pt_mu_sum_rebin", "#Sum p_{T}^{#mu} [GeV]", 5, bins_HTlept_low);
  book<TH1F>("Pt_lept1", "leading lepton p_{T [GeV]}", 75, 0, 1500);
  book<TH1F>("Pt_lept2", "subleading lepton p_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_lept12", "leading lepton p_{T} + subleading lepton p_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_lept12_rebin", "leading lepton p_{T} + subleading lepton p_{T} [GeV]", 5,bins_HTlept_low);
  book<TH1F>("Pt_mu1", "p_{T}^{leading #mu} [GeV]", 75, 0, 1500);
  double bins_pt_low[26] = {0,30,60,90,120,150,180,210,240,270,300,350,400,450,500,550,600,650,700,750,800,900,1000,1100,1300,1500};
  book<TH1F>("Pt_mu1_rebin", "P_{T}^{leading #mu} [GeV]", 25, bins_pt_low);
  book<TH1F>("Pt_mu1_NoEle", "p_{T}^{leading #mu}, no Ele [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_mu1_NoEle_rebin", "P_{T}^{leading #mu}, no Ele [GeV]", 25, bins_pt_low);

  book<TH1F>("pt_ele", "p_{T}^{ele} [GeV]", 50, 0, 1500);
  book<TH1F>("pt_ele_zoom", "p_{T}^{ele} [GeV]", 34, 0, 1020);


  // general
  book<TH1F>("N_pv", "number of primary vertices", 51, -0.50, 50.5);
  book<TH1F>("N_pv_zoom", "number of primary vertices", 31, -0.50, 30.5);
  book<TH1F>("E_Tmiss", "missing E_{T} [GeV]", 75, 0, 1500);
  book<TH1F>("E_Tmiss_0Ele2Mu", "missing E_{T} [GeV] for N_{e}=0, N_{#mu}=2", 75, 0,1500);
  book<TH1F>("H_T", "H_{T} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_from350", "H_{T} [GeV]", 24, 0, 4200);
  double bins_from350[21] = {0,175,350,525,700,875,1050,1225,1400,1575,1750,1925,2100,2275,2450,2625,2800,2975,3325,3675,4200}; //same binning as _from350 up to 2975, then two double-size and one triple-size bin
  book<TH1F>("H_T_from350_all_filled", "H_{T} [GeV]", 20, bins_from350);
  double bins_from350_allfilled[16] = {0,175,350,525,700,875,1050,1225,1400,1575,1750,1925,2100,2450,2800,3000}; //same binning as _from350_all_filled up to 2100, then larger bins
  book<TH1F>("H_T_from350_all_filled_rebin", "H_{T} [GeV]", 15, bins_from350_allfilled);
  book<TH1F>("H_T_from350_rebin", "H_{T} [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_from350_rebin2", "H_{T} [GeV]", 12, 0, 4200);
  book<TH1F>("Parton_H_T", "H_{T} [GeV] on parton level", 80,0,7000);
  double bins_low_1Ele[12] = {0,350,500,700,900,1100,1300,1500,1750,2000,2500,7000};
  double bins_low_NoEle[23] = {0,200,350,500,650,800,950,1100,1250,1400,1550,1700,1850,2000,2150,2300,2450,2600,2750,2900,3050,3200,7000};
  double bins_low_NoEle2[11] = {0,350,500,650,800,950,1100,1250,1450,1750,2050};
  book<TH1F>("H_T_rebin", "H_{T} [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_rebin2", "H_{T} [GeV]", 100, 0, 7000);
  book<TH1F>("H_T_rebin3", "H_{T} [GeV]", 10,bins_low_NoEle2);
  book<TH1F>("H_T_jets", "H_{T}^{jets} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_jets_from350_rebin", "H_{T}^{jets} [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_jets_from350_all_filled_rebin", "H_{T}^{jets} [GeV]", 15, bins_from350_allfilled);
  book<TH1F>("H_T_jets_rebin", "H_{T}^{jets} rebinned [GeV]", 5, bins_HTlept_low);
  book<TH1F>("H_T_lept", "H_{T}^{leptons} [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_lept_from350_rebin", "H_{T}^{leptons} [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_lept_from350_all_filled_rebin", "H_{T}^{leptons} [GeV]", 15, bins_from350_allfilled);
  book<TH1F>("H_T_lept_zoom", "H_{T}^{leptons} [GeV]", 40, 0, 4000);
  book<TH1F>("H_T_lept_rebin", "H_{T}^{leptons} rebinned [GeV]", 5, bins_HTlept_low);
  book<TH1F>("H_T_comb_NoEle", "H_{T}, no Ele [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_NoEle_from350", "H_{T}, no Ele [GeV]", 24, 0, 4200);
  book<TH1F>("H_T_comb_NoEle_from350_rebin", "H_{T}, no Ele [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_comb_NoEle_from350_rebin2", "H_{T}, no Ele [GeV]", 12, 0, 4200);
  book<TH1F>("H_T_comb_NoEle_from350_all_filled_rebin", "H_{T}, no Ele [GeV]", 15, bins_from350_allfilled);
  book<TH1F>("H_T_comb_NoEle_rebin", "H_{T}, no Ele [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_comb_NoEle_rebin2", "H_{T}, no Ele [GeV]", 10, bins_low_NoEle2);
  book<TH1F>("Integral_NoEle", "BinContent = sum(eventweights), NoEle", 1, 0.5, 1.5);
  book<TH1F>("H_T_comb_1Ele", "H_{T}, N_{Ele} #geq 1 [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_1Ele_from350", "H_{T}, N_{Ele} #geq 1 [GeV]", 24, 0, 4200);
  book<TH1F>("H_T_comb_1Ele_from350_rebin", "H_{T}, N_{Ele} #geq 1 [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_comb_1Ele_from350_rebin2", "H_{T}, N_{Ele} #geq 1 [GeV]", 12, 0, 4200);
  book<TH1F>("H_T_comb_1Ele_rebin", "H_{T}, N_{Ele} #geq 1 [GeV]", 11, bins_low_1Ele);
  book<TH1F>("H_T_comb_1Ele_rebin2", "H_{T}, N_{Ele} #geq 1, same binning as for N_{Ele} = 0 [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_comb_1Ele_rebin3", "H_{T}, N_{Ele} #geq 1, same binning as for N_{Ele} = 0 [GeV]", 10, bins_low_NoEle2);
  book<TH1F>("Integral_1Ele", "BinContent = sum(eventweights), 1Ele", 1, 0.5, 1.5);
  book<TH1F>("H_T_2mu_from350_all_filled_rebin", "H_{T}, 2 #mu [GeV]", 15, bins_from350_allfilled);
  book<TH1F>("M_LQ_comb", "M_{LQ,mean} [GeV]", 60, 0, 3000);
  double bins_mlq_low[17] = {100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,1000,2000};
  book<TH1F>("M_LQ_comb_rebin", "M_{LQ,mean} [GeV]", 16, bins_mlq_low);
  double bins_mlq_low2[6] = {0,200,400,600,800,1000};
  book<TH1F>("M_LQ_comb_rebin2", "M_{LQ,mean} [GeV]", 5, bins_mlq_low2);
  double bins_mlq_low3[7] = {0,250,350,450,600,750,1000};
  book<TH1F>("M_LQ_comb_all_filled", "M_{LQ,mean} [GeV]", 6, bins_mlq_low3);
  book<TH1F>("chi2", "#chi^{2}", 100, 0,200);
  book<TH1F>("M_LQ_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV]", 50, -500, 500);
  book<TH1F>("M_LQ_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep})/M_{LQ,mean} [GeV]", 50, -0.5, 0.5);
  book<TH1F>("M_LQLQ", "M_{LQLQ} [GeV]", 100, 0, 5000);
  double bins_mlqlq_low[12] = {0,400,500,600,700,800,900,1000,1250,1500,2000,5000};
  book<TH1F>("M_LQLQ_rebin", "M_{LQLQ} [GeV]", 11, bins_mlqlq_low);
  book<TH1F>("N_jets_had","number of jets in the hadronic top hypothesis",6,-0.5,5.5);
  book<TH1F>("N_jets_lep","number of jets in the leptonic top hypothesis",6,-0.5,5.5);

  book<TH1F>("dR_toplep_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_tophad_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_mulep","#Delta R(t^{lep},#mu^{lep})",50,0,5);
  book<TH1F>("dR_toplep_muhad","#Delta R(t^{had},#mu^{had})",50,0,5);
  book<TH1F>("dR_tophad_muX", "1: t^{had} closer to #mu^{had}, -1: closer to #mu^{lep}", 3,-1.5, 1.5);
  //book<TH1F>("dummy","dummy",50,0,5);

  book<TH1F>("M_LQ_Muonic_comb", "M_{LQ,mean} [GeV]", 60, 0, 3000);
  book<TH1F>("M_LQ_Muonic_comb_rebin", "M_{LQ,mean} [GeV]", 16, bins_mlq_low);
  book<TH1F>("M_LQ_Muonic_comb_rebin2", "M_{LQ,mean} [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_Muonic_comb_all_filled", "M_{LQ,mean} [GeV]", 6, bins_mlq_low3);
  book<TH1F>("chi2_Muonic", "#chi^{2}", 100, 0,200);
  book<TH1F>("M_LQ_Muonic_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV]", 50, -500, 500);
  book<TH1F>("M_LQ_Muonic_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep})/M_{LQ,mean} [GeV]", 50, -0.5, 0.5);
  book<TH1F>("M_LQLQ_Muonic", "M_{LQLQ} [GeV]", 100, 0, 5000);
  book<TH1F>("M_LQLQ_Muonic_rebin", "M_{LQLQ} [GeV]", 11, bins_mlqlq_low);

  book<TH1F>("M_LQ_MuEle_comb", "M_{LQ,mean} [GeV]", 60, 0, 3000);
  book<TH1F>("M_LQ_MuEle_comb_rebin", "M_{LQ,mean} [GeV]", 16, bins_mlq_low);
  book<TH1F>("M_LQ_MuEle_comb_rebin2", "M_{LQ,mean} [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_MuEle_comb_all_filled", "M_{LQ,mean} [GeV]", 6, bins_mlq_low3);
  book<TH1F>("chi2_MuEle", "#chi^{2}", 100, 0,200);
  book<TH1F>("chi2_MuEle_rebin", "#chi^{2}", 40, 0,200);
  book<TH1F>("chi2_MuEle_rebin2", "#chi^{2}", 20, 0,200);
  book<TH1F>("chi2_MuEle_rebin3", "#chi^{2}", 10, 0,200);
  book<TH2F>("chi2_MuEle_nothad_vs_total", "#chi^{2};#chi^2 (no thad);#chi^{2} (total)", 100, 0, 200, 100, 0, 200);
  book<TH2F>("chi2_MuEle_nothad_vs_total_rebin", "#chi^{2};#chi^2 (no thad);#chi^{2} (total)", 40, 0, 200, 40, 0, 200);
  book<TH2F>("chi2_MuEle_nothad_vs_total_rebin2", "#chi^{2};#chi^2 (no thad);#chi^{2} (total)", 20, 0, 200, 20, 0, 200);
  book<TH2F>("chi2_MuEle_nothad_vs_total_rebin3", "#chi^{2};#chi^2 (no thad);#chi^{2} (total)", 10, 0, 200, 10, 0, 200);
  book<TH1F>("M_LQ_MuEle_comb_1hadjet", "M_{LQ,mean}, 1 jet in had. hypo. [GeV]", 60, 0, 3000);
  book<TH1F>("M_LQ_MuEle_comb_rebin_1hadjet", "M_{LQ,mean}, 1 jet in had. hypo.  [GeV]", 16, bins_mlq_low);
  book<TH1F>("M_LQ_MuEle_comb_rebin2_1hadjet", "M_{LQ,mean}, 1 jet in had. hypo.  [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_MuEle_comb_all_filled_1hadjet", "M_{LQ,mean}, 1 jet in had. hypo.  [GeV]", 6, bins_mlq_low3);
  book<TH1F>("M_LQ_MuEle_comb_2hadjet", "M_{LQ,mean}, 2 jets in had. hypo. [GeV]", 60, 0, 3000);
  book<TH1F>("M_LQ_MuEle_comb_rebin_2hadjet", "M_{LQ,mean}, 2 jets in had. hypo.  [GeV]", 16, bins_mlq_low);
  book<TH1F>("M_LQ_MuEle_comb_rebin2_2hadjet", "M_{LQ,mean}, 2 jets in had. hypo.  [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_MuEle_comb_all_filled_2hadjet", "M_{LQ,mean}, 2 jets in had. hypo.  [GeV]", 6, bins_mlq_low3);
  book<TH1F>("M_LQ_MuEle_comb_3hadjet", "M_{LQ,mean}, 3 jets in had. hypo. [GeV]", 60, 0, 3000);
  book<TH1F>("M_LQ_MuEle_comb_rebin_3hadjet", "M_{LQ,mean}, 3 jets in had. hypo.  [GeV]", 16, bins_mlq_low);
  book<TH1F>("M_LQ_MuEle_comb_rebin2_3hadjet", "M_{LQ,mean}, 3 jets in had. hypo.  [GeV]", 5, bins_mlq_low2);
  book<TH1F>("M_LQ_MuEle_comb_all_filled_3hadjet", "M_{LQ,mean}, 3 jets in had. hypo.  [GeV]", 6, bins_mlq_low3);
  book<TH1F>("M_tophad_rec_MuEle_1hadjet", "M_{top}^{had}, 1 jet in had. hypo [GeV]", 30, 0, 300);
  book<TH1F>("M_tophad_rec_MuEle_1hadjet_rebin", "M_{top}^{had}, 1 jet in had. hypo [GeV]", 15, 0, 300);
  book<TH1F>("M_tophad_rec_MuEle_1hadjet_rebin2", "M_{top}^{had}, 1 jet in had. hypo [GeV]", 10, 0, 300);
  book<TH1F>("M_tophad_rec_MuEle_1hadjet_rebin3", "M_{top}^{had}, 1 jet in had. hypo [GeV]", 5, 0, 300);
  book<TH1F>("M_tophad_rec_MuEle_1hadjet_rebin4", "M_{top}^{had}, 1 jet in had. hypo [GeV]", 20, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 70, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 150, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin2", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 30, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin3", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 35, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin4", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin5", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 15, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin6", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 50, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin7", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 100, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin8", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 20, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin9", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 20, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin10", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 40, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_1hadjet_rebin11", "p_{T}^{jet} in 1-jet had. hypo  [GeV]", 10, 0, 300);
  book<TH1F>("M_tophad_rec_MuEle_2hadjet", "M_{top}^{had}, 2 jets in had. hypo [GeV]", 70, 0, 700);
  book<TH1F>("M_tophad_rec_MuEle_2hadjet_rebin", "M_{top}^{had}, 2 jets in had. hypo [GeV]", 35, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 70, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 150, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin2", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 30, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin3", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 35, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin4", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin5", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 15, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin6", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 50, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin7", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 100, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin8", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 20, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin9", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 20, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin10", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 40, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_2hadjet_rebin11", "p_{T}^{jet} in 2-jet had. hypo  [GeV]", 10, 0, 300);
  book<TH1F>("M_tophad_rec_MuEle_3hadjet", "M_{top}^{had}, 3 jets in had. hypo [GeV]", 70, 0, 700);
  book<TH1F>("M_tophad_rec_MuEle_3hadjet_rebin", "M_{top}^{had}, 3 jets in had. hypo [GeV]", 35, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 70, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 150, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin2", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 30, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin3", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 35, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin4", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 75, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin5", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 15, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin6", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 50, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin7", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 100, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin8", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 20, 0, 300);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin9", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 20, 0, 700);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin10", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 40, 0, 1500);
  book<TH1F>("Pt_tophad_rec_MuEle_3hadjet_rebin11", "p_{T}^{jet} in 3-jet had. hypo  [GeV]", 10, 0, 300);
  book<TH1F>("M_LQ_MuEle_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV]", 50, -500, 500);
  book<TH1F>("M_LQ_MuEle_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep})/M_{LQ,mean} [GeV]", 50, -0.5, 0.5);
  book<TH1F>("M_LQLQ_MuEle", "M_{LQLQ} [GeV]", 100, 0, 5000);
  book<TH1F>("M_LQLQ_MuEle_rebin", "M_{LQLQ} [GeV]", 11, bins_mlqlq_low);

 book<TH1F>("MT_Lep_Met", "M_{T} (lepton, E_{T}^{miss} [GeV])", 16, 0, 800);
 book<TH1F>("dPhi_Lep_Met", "#Delta #Phi (lepton, E_{T}^{miss} [GeV])", 16, 0, 3.2);

  book<TH1F>("H_T_comb_NoMLQ", "H_{T}, no MLQ [GeV]", 50, 0, 7000);
  book<TH1F>("H_T_comb_NoMLQ_from350", "H_{T}, no MLQ [GeV]", 24, 0, 4200);
  book<TH1F>("H_T_comb_NoMLQ_from350_rebin", "H_{T}, no MLQ [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_comb_NoMLQ_from350_rebin_noweights", "H_{T}, no MLQ [GeV]", 48, 0, 4200);
  book<TH1F>("H_T_comb_NoMLQ_from350_rebin2", "H_{T}, no MLQ [GeV]", 12, 0, 4200);
  book<TH1F>("H_T_comb_NoMLQ_from350_all_filled_rebin", "H_{T}, no MLQ [GeV]", 15, bins_from350_allfilled);
  book<TH1F>("H_T_comb_NoMLQ_rebin", "H_{T}, no MLQ [GeV]", 22, bins_low_NoEle);
  book<TH1F>("H_T_comb_NoMLQ_rebin2", "H_{T}, no MLQ [GeV]", 10, bins_low_NoEle2);
  book<TH1F>("Integral_NoMLQ", "BinContent = sum(eventweights), No MLQ", 1, 0.5, 1.5);

  book<TH2F>("dR_ee_HT","#Delta R(e,e) vs. H_{T};H_{T};#Delta R(e,e)",20,bins_from350,50,0,5);
  book<TH2F>("dR_mumu_HT","#Delta R(#mu,#mu) vs. H_{T};H_{T};#Delta R(#mu,#mu)",20,bins_from350,50,0,5);

 
  //electron and muon fakes
  book<TH1F>("ele_type", "0 real ele, 1 fake ele", 2,-0.5,1.5);
  book<TH1F>("mu_type", "0 real #mu, 1 fake #mu", 2,-0.5,1.5);

  book<TH1F>("ele_type_mlq_reco", "MLQ is reconstructed, 0 real ele, 1 fake ele", 2,-0.5,1.5);
  book<TH1F>("mu_type_mlq_reco", "MLQ is reconstructed, 0 real #mu, 1 fake #mu", 2,-0.5,1.5);

  book<TH1F>("more_gen_ele_mlq_reco", "MLQ is reconstructed, 0 <, 1 #geq", 2,-0.5,1.5);
  book<TH1F>("more_gen_mu_mlq_reco", "MLQ is reconstructed, 0 <, 1 #geq", 2,-0.5,1.5);

  book<TH1F>("dr_min_ele_genele", "#DeltaR(e, closest gen-e)", 100,0,5);
  book<TH1F>("dr_min_mu_genmu", "#DeltaR(#mu, closest gen-#mu)", 100,0,5);

  book<TH1F>("dr_min_ele_genele_mlq_reco", "MLQ is reconstructed, #DeltaR(e, closest gen-e)", 100,0,5);
  book<TH1F>("dr_min_mu_genmu_mlq_reco", "MLQ is reconstructed, #DeltaR(#mu, closest gen-#mu)", 100,0,5);

  book<TH1F>("dr_min_ele_geneletau", "#DeltaR(e, closest gen-e/#tau)", 100,0,5);
  book<TH1F>("dr_min_mu_genmutau", "#DeltaR(#mu, closest gen-#mu/#tau)", 100,0,5);

  book<TH1F>("dr_min_ele_geneletau_mlq_reco", "MLQ is reconstructed, #DeltaR(e, closest gen-e/#tau)", 100,0,5);
  book<TH1F>("dr_min_mu_genmutau_mlq_reco", "MLQ is reconstructed, #DeltaR(#mu, closest gen-#mu/#tau)", 100,0,5);

  book<TH1F>("jets_faking_ele_pt", "p_{T} of jets faking an electron [GeV]", 50, 20, 1500);
  book<TH1F>("jets_faking_mu_pt", "p_{T} of jets faking a muon [GeV]", 50, 20, 1500);

  book<TH1F>("dr_min_mu_gen_b", "dR min (muon, gen b-quark) for each muon", 100,0,5);
  book<TH1F>("dr_min_mu_gen_b_MLQ", "dR min (muon, gen b-quark) for each muon, MLQ reco.", 100,0,5);
  book<TH1F>("dr_min_mu_gen_b_NoMLQ", "dR min (muon, gen b-quark) for each muon, MLQ not reco.", 100,0,5);

  book<TH1F>("n_fake_mu", "number of fake muons in all events", 1, -0.5, 0.5);
  book<TH1F>("n_fake_mu_MLQ", "number of fake muons, MLQ reconstructed", 1, -0.5, 0.5);
  book<TH1F>("n_fake_mu_NoMLQ", "number of fake muons, MLQ not reconstructed", 1, -0.5, 0.5);

  book<TH1F>("n_fake_mu_HF", "number of fake muons in all events, not from HF", 1, -0.5, 0.5);
  book<TH1F>("n_fake_mu_MLQ_HF", "number of fake muons, MLQ reconstructed, not from HF", 1, -0.5, 0.5);
  book<TH1F>("n_fake_mu_NoMLQ_HF", "number of fake muons, MLQ not reconstructed, not from HF", 1, -0.5, 0.5);

  book<TH2F>("n_jets_hadhyp_vs_chi2",";number of jets in hadronic top hypothesis; #chi^{2} of best LQ hypothesis", 6, 0, 6, 50, 0, 500);

  //event weights: sum
  book<TH1F>("sum_event_weights", "BinContent = sum(eventweights)", 1, 0.5, 1.5);
  book<TH1F>("sum_event_weights_3lep", "cut and count, #geq 3 leptons", 1, 0.5, 1.5);
  book<TH1F>("sum_event_weights_MLQ_MuEle", "cut and count, MLQ reconstructed", 1, 0.5, 1.5);

 

  //For MLQ reconstruction
  h_hyps        = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
  h_muonic_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassMuonicLQReconstruction");
  //h_hadr_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassHadronicLQReconstruction");
  m_discriminator_name ="Chi2";
  //m_discriminator_name ="CorrectMatch";

  is_mc = ctx.get("dataset_type") == "MC";

}


void LQToTopMuHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.

  double weight = event.weight;

  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);

  for(unsigned int i=0; i<event.jets->size(); i++){
    hist("pt_jets")->Fill(jets->at(i).pt(),weight);
    hist("eta_jets")->Fill(jets->at(i).eta(),weight);
  }
  
  if(Njets>=1){
    hist("eta_jet1")->Fill(jets->at(0).eta(), weight);
  }
  if(Njets>=2){
    hist("eta_jet2")->Fill(jets->at(1).eta(), weight);
  }
  if(Njets>=3){
    hist("eta_jet3")->Fill(jets->at(2).eta(), weight);
  }
  if(Njets>=4){
    hist("eta_jet4")->Fill(jets->at(3).eta(), weight);
  }
  if(Njets>=1){
    hist("pt_jet1")->Fill(jets->at(0).pt(), weight);
  }
  if(Njets>=2){
    hist("pt_jet2")->Fill(jets->at(1).pt(), weight);
  }
  if(Njets>=3){
    hist("pt_jet3")->Fill(jets->at(2).pt(), weight);
  }

  //# b-jets
  std::vector<Jet> bjets_loose, bjets_med, bjets_tight;
  CSVBTag Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
  CSVBTag Btag_medium = CSVBTag(CSVBTag::WP_MEDIUM);
  CSVBTag Btag_tight = CSVBTag(CSVBTag::WP_TIGHT);


  for (unsigned int i =0; i<jets->size(); ++i) {
    if(Btag_loose(jets->at(i),event)) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_loose.push_back(jets->at(i));
    }
    if(Btag_medium(jets->at(i),event)) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_med.push_back(jets->at(i));
    }
    if(Btag_tight(jets->at(i),event)) { //loose: >0.605, medium: >0.890, tight: >0.970
      bjets_tight.push_back(jets->at(i));
    }
  }

  int NbJets_loose = bjets_loose.size();
  hist("N_bJets_loose")->Fill(NbJets_loose,weight);
  int NbJets_med = bjets_med.size();
  hist("N_bJets_med")->Fill(NbJets_med,weight);
  int NbJets_tight = bjets_tight.size();
  hist("N_bJets_tight")->Fill(NbJets_tight,weight);


  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  double sum_mu_pt = 0;
  for (const Muon & thismu : *event.muons){
    hist("pt_mu")->Fill(thismu.pt(), weight);
    hist("pt_mu_zoom")->Fill(thismu.pt(), weight);
    hist("eta_mu")->Fill(thismu.eta(), weight);
    hist("reliso_mu")->Fill(thismu.relIso(), weight);
    hist("reliso_mu_rebin")->Fill(thismu.relIso(), weight);
    sum_mu_pt += thismu.pt();
  }
  hist("Pt_mu_sum")->Fill(sum_mu_pt,weight);
  hist("Pt_mu_sum_rebin")->Fill(sum_mu_pt,weight);
  if(Nmuons>0){
    hist("Pt_mu1")->Fill(event.muons->at(0).pt(), weight);
    hist("Pt_mu1_rebin")->Fill(event.muons->at(0).pt(), weight);
  }

  for(const auto & e : *event.electrons){
    hist("relisorho_ele")->Fill(e.relIsorho(event.rho), weight);
    hist("relisorho_ele_rebin")->Fill(e.relIsorho(event.rho), weight);
  }

  double pt_lept1 = 0, pt_lept2 = 0;
  if(event.electrons->size() > 1){
    if(event.muons->size() > 0){
      if(event.muons->at(0).pt()>event.electrons->at(0).pt()) {
	pt_lept1 = event.muons->at(0).pt();
	if(event.muons->size() > 1){
	  if(event.muons->at(1).pt()>event.electrons->at(0).pt()) pt_lept2 = event.muons->at(1).pt();
	  else pt_lept2 = event.electrons->at(0).pt();
	}
	else pt_lept2 = event.electrons->at(0).pt();
      }
      else{
	pt_lept1 = event.electrons->at(0).pt();
	if(event.electrons->at(1).pt()>=event.muons->at(0).pt()) pt_lept2 = event.electrons->at(1).pt();
	else pt_lept2 = event.muons->at(0).pt();
      }
    }
    else {
      pt_lept1 = event.electrons->at(0).pt();
      pt_lept2 = event.electrons->at(1).pt();
    }
  }
  if(event.electrons->size() == 1 && event.muons->size() > 0){
    if(event.muons->at(0).pt()>event.electrons->at(0).pt()){
      pt_lept1 = event.muons->at(0).pt();
      if(event.muons->size() > 1){
	if(event.muons->at(1).pt()>event.electrons->at(0).pt()) pt_lept2 = event.muons->at(1).pt();
	else pt_lept2 = event.electrons->at(0).pt();
      }
      else pt_lept2 = event.electrons->at(0).pt();
    }
    else{
      pt_lept1 = event.electrons->at(0).pt();
      pt_lept2 = event.muons->at(0).pt();
    }
  }
  else if(event.electrons->size() == 1 && event.muons->size() == 0){
    pt_lept1 = event.electrons->at(0).pt();
    pt_lept2 = 0;
  }
  else if(event.electrons->size() == 0){
    if(event.muons->size() > 0){
      pt_lept1 = event.muons->at(0).pt();
      if(event.muons->size() > 1) pt_lept2 = event.muons->at(1).pt();
    }
  }



  hist("Pt_lept1")->Fill(pt_lept1,weight);
  if(event.electrons->size()+event.muons->size() > 1) hist("Pt_lept2")->Fill(pt_lept2,weight);
  hist("Pt_lept12")->Fill(pt_lept1+pt_lept2,weight);
  hist("Pt_lept12_rebin")->Fill(pt_lept1+pt_lept2,weight);

  for (const Electron & thisele : *event.electrons){
    hist("pt_ele")->Fill(thisele.pt(), weight);
    hist("pt_ele_zoom")->Fill(thisele.pt(), weight);

  }

  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);
  hist("N_pv_zoom")->Fill(Npvs, weight);

  //HT
  auto met = event.met->pt();
  hist("E_Tmiss")->Fill(met, weight);
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
  hist("H_T_jets")->Fill(ht_jets,weight);
  hist("H_T_jets_from350_rebin")->Fill(ht_jets,weight);
  if(ht_jets <= 2900) hist("H_T_jets_from350_all_filled_rebin")->Fill(ht_jets,weight);
  else                hist("H_T_jets_from350_all_filled_rebin")->Fill(2900,weight);
  hist("H_T_lept")->Fill(ht_lep,weight);
  hist("H_T_lept_from350_rebin")->Fill(ht_lep,weight);
  if(ht_lep <= 2900) hist("H_T_lept_from350_all_filled_rebin")->Fill(ht_lep,weight);
  else               hist("H_T_lept_from350_all_filled_rebin")->Fill(2900,weight);
  hist("H_T_lept_zoom")->Fill(ht_lep,weight);
  hist("H_T_jets_rebin")->Fill(ht_jets,weight);
  hist("H_T_lept_rebin")->Fill(ht_lep,weight);
  hist("H_T")->Fill(ht, weight);
  hist("H_T_from350")->Fill(ht, weight);
  if(ht <= 4000) hist("H_T_from350_all_filled")->Fill(ht, weight);
  else hist("H_T_from350_all_filled")->Fill(4000, weight);
  if(ht <= 2900) hist("H_T_from350_all_filled_rebin")->Fill(ht, weight);
  else hist("H_T_from350_all_filled_rebin")->Fill(2900, weight);
  hist("H_T_from350_rebin")->Fill(ht, weight);
  hist("H_T_from350_rebin2")->Fill(ht, weight);
  hist("H_T_rebin")->Fill(ht, weight);
  hist("H_T_rebin2")->Fill(ht,weight);
  if(ht <= 2000) hist("H_T_rebin3")->Fill(ht,weight);
  else hist("H_T_rebin3")->Fill(2000,weight);

  //partonlvl HT:
  if(is_mc){
    double partonHT = 0;
    constexpr const int invalid_daughter = (unsigned short)(-1);
    for(const auto & gp : *event.genparticles){
      if(gp.daughter1() != invalid_daughter || gp.daughter2() != invalid_daughter) continue;
      // if we are here, it means we have a final state particle.
      // Add to HT in cas it is a parton (quark -- including b but not top as tops are never final state particles -- or gluon -- or ele/mu -- or its respective neutrino).
      // Note that the exact HT definition depends on the madgraph configuration, but this
      // should cover the most common case.
      int id = abs(gp.pdgId());
      if((id >= 1 && id <= 5) || (id == 21) || (id>=11 && id <= 14)){
	partonHT += gp.pt();
      }
    }
    hist("Parton_H_T")->Fill(partonHT, weight);
  }

  // M_mumu Invariant Mass
  double M_mumu;
  LorentzVector muons[Nmuons];
  for(int i=0; i<Nmuons; i++){
    muons[i] = event.muons->at(i).v4();
  }
  for(int i=0; i<Nmuons; i++){
    for(int j=0; j<Nmuons; j++){
      if(j > i){
	M_mumu = (muons[i] + muons[j]).M();
	hist("M_mumu")->Fill(M_mumu, weight);	
      }
    }
  }

  //M_ee
  double M_ee;
  const int Nele = event.electrons->size();
  LorentzVector electrons[Nele];
  for(int i=0; i<Nele; i++){
    electrons[i] = event.electrons->at(i).v4();
  }
  for(int i=0; i<Nele; i++){
    for(int j=0; j<Nele; j++){
      if(j > i){
	M_ee = (electrons[i] + electrons[j]).M();
	hist("M_ee")->Fill(M_ee, weight);	
      }
    }
  }

  if(event.electrons->size() >= 2){
    for(unsigned int i=0; i<event.electrons->size(); i++){
      for(unsigned int j=0; j<event.electrons->size(); j++){
	if(j<=i) continue;

	double dr_ee = deltaR(event.electrons->at(i),event.electrons->at(j));
	((TH2F*)hist("dR_ee_HT"))->Fill(ht,dr_ee,weight);
      }
    }
  }

  if(event.muons->size() >= 2){
    for(unsigned int i=0; i<event.muons->size(); i++){
      for(unsigned int j=0; j<event.muons->size(); j++){
	if(j<=i) continue;

	double dr_mumu = deltaR(event.muons->at(i),event.muons->at(j));
	((TH2F*)hist("dR_mumu_HT"))->Fill(ht,dr_mumu,weight);
      }
    }
  }


  if(event.electrons->size() + event.muons->size() >= 3) hist("sum_event_weights_3lep")->Fill(1.,weight);
  else hist("H_T_2mu_from350_all_filled_rebin")->Fill(ht,weight);

  
  //check for at least 1 muon pair with opposite charge
  bool charge_opposite = false;
  for(unsigned int i=0; i<event.muons->size(); i++){
    for(unsigned int j=0; j<event.muons->size(); j++){
      if(j>i){
	if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
	  charge_opposite = true;
	}
      }
    }
  }

  bool reconstruct_mlq_ele = (Nele >= 1 && event.muons->size() >= 2 && charge_opposite && event.jets->size() >= 2);
  if(reconstruct_mlq_ele){   
    std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps); 
    const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
    double chi2 = hyp->discriminator(m_discriminator_name);
    double chi2_nothad = hyp->discriminator(m_discriminator_name+"_tlep_MQLdiff_rel");
    hist("chi2")->Fill(chi2,weight);
    hist("chi2_MuEle")->Fill(chi2,weight);
    hist("chi2_MuEle_rebin")->Fill(chi2,weight);
    hist("chi2_MuEle_rebin2")->Fill(chi2,weight);
    hist("chi2_MuEle_rebin3")->Fill(chi2,weight);
    ((TH2F*)hist("chi2_MuEle_nothad_vs_total"))->Fill(chi2_nothad,chi2,weight);
    ((TH2F*)hist("chi2_MuEle_nothad_vs_total_rebin"))->Fill(chi2_nothad,chi2,weight);
    ((TH2F*)hist("chi2_MuEle_nothad_vs_total_rebin2"))->Fill(chi2_nothad,chi2,weight);
    ((TH2F*)hist("chi2_MuEle_nothad_vs_total_rebin3"))->Fill(chi2_nothad,chi2,weight);

    double mLQlep_rec = 0;
    double mLQhad_rec = 0;
    double mLQmed_rec = 0;
    double mLQdiff = 0;
    double mLQdiff_rel = 0;
    double mLQLQ = 0;

    if( (hyp->LQlep_v4()).isTimelike() ) {mLQlep_rec = (hyp->LQlep_v4()).M();}
    else {mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());}
    if( (hyp->LQhad_v4()).isTimelike() ) {mLQhad_rec = (hyp->LQhad_v4()).M();}
    else {mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());}
    


    double n_jets_had = hyp->tophad_jets().size();
    double n_jets_lep = hyp->toplep_jets().size();
    hist("N_jets_had")->Fill(n_jets_had,weight);
    hist("N_jets_lep")->Fill(n_jets_lep,weight);
    ((TH2F*)hist("n_jets_hadhyp_vs_chi2"))->Fill(n_jets_had,chi2,weight);
    

    mLQdiff = mLQhad_rec - mLQlep_rec;
    mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;
    hist("M_LQ_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_comb_rebin")->Fill(mLQmed_rec, weight);
    if(mLQmed_rec < 900)   hist("M_LQ_comb_rebin2")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_comb_rebin2")->Fill(900., weight);
    if(mLQmed_rec < 900)   hist("M_LQ_comb_all_filled")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_comb_all_filled")->Fill(900., weight);
    hist("M_LQ_diff")->Fill(mLQdiff, weight);
    mLQdiff_rel = mLQdiff / mLQmed_rec;
    hist("M_LQ_diff_rel")->Fill(mLQdiff_rel,weight);
    mLQLQ = (hyp->LQlep_v4()+hyp->LQhad_v4()).M();
    hist("M_LQLQ")->Fill(mLQLQ,weight);
    hist("M_LQLQ_rebin")->Fill(mLQLQ,weight);

    hist("M_LQ_MuEle_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_MuEle_comb_rebin")->Fill(mLQmed_rec, weight);
    if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_MuEle_comb_rebin2")->Fill(900., weight);
    if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_MuEle_comb_all_filled")->Fill(900., weight);


    if(n_jets_had==1){
      hist("M_LQ_MuEle_comb_1hadjet")->Fill(mLQmed_rec, weight);
      hist("M_LQ_MuEle_comb_rebin_1hadjet")->Fill(mLQmed_rec, weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2_1hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_rebin2_1hadjet")->Fill(900., weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled_1hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_all_filled_1hadjet")->Fill(900., weight);
      hist("M_tophad_rec_MuEle_1hadjet")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin2")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin3")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin4")->Fill(hyp->tophad_v4().M(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin2")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin3")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin4")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin5")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin6")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin7")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin8")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin9")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin10")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin11")->Fill(hyp->tophad_v4().Pt(),weight);
    }
    else if(n_jets_had==2){
      hist("M_LQ_MuEle_comb_2hadjet")->Fill(mLQmed_rec, weight);
      hist("M_LQ_MuEle_comb_rebin_2hadjet")->Fill(mLQmed_rec, weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2_2hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_rebin2_2hadjet")->Fill(900., weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled_2hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_all_filled_2hadjet")->Fill(900., weight);
      hist("M_tophad_rec_MuEle_2hadjet")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_2hadjet_rebin")->Fill(hyp->tophad_v4().M(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin2")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin3")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin4")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin5")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin6")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin7")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin8")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin9")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin10")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin11")->Fill(hyp->tophad_v4().Pt(),weight);
    }
    else if(n_jets_had==3){
      hist("M_LQ_MuEle_comb_3hadjet")->Fill(mLQmed_rec, weight);
      hist("M_LQ_MuEle_comb_rebin_3hadjet")->Fill(mLQmed_rec, weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2_3hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_rebin2_3hadjet")->Fill(900., weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled_3hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_all_filled_3hadjet")->Fill(900., weight);
      hist("M_tophad_rec_MuEle_3hadjet")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_3hadjet_rebin")->Fill(hyp->tophad_v4().M(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin2")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin3")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin4")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin5")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin6")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin7")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin8")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin9")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin10")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin11")->Fill(hyp->tophad_v4().Pt(),weight);
    }



    hist("M_LQ_MuEle_diff")->Fill(mLQdiff, weight);
    hist("M_LQ_MuEle_diff_rel")->Fill(mLQdiff_rel,weight);
    hist("M_LQLQ_MuEle")->Fill(mLQLQ,weight);
    hist("M_LQLQ_MuEle_rebin")->Fill(mLQLQ,weight);
    
    //deltaR between tops and associated muons
    double dR_lep = deltaR(hyp->toplep_v4(),hyp->mu_lep_v4());
    double dR_had = deltaR(hyp->tophad_v4(),hyp->mu_had_v4());
    double dR_lephad = deltaR(hyp->toplep_v4(),hyp->mu_had_v4());
    double dR_hadlep = deltaR(hyp->tophad_v4(),hyp->mu_lep_v4());

    hist("dR_toplep_mulep")->Fill(dR_lep,weight);
    hist("dR_tophad_muhad")->Fill(dR_had,weight);
    hist("dR_toplep_muhad")->Fill(dR_lephad,weight);
    hist("dR_tophad_mulep")->Fill(dR_hadlep,weight);
    if(dR_had<dR_hadlep)  hist("dR_tophad_muX")->Fill(1,weight);
    else                  hist("dR_tophad_muX")->Fill(-1,weight);

    //mT between ele and met
    Particle ele = hyp->electron();
    double mt = sqrt(2*ele.pt()*event.met->pt()*(1-cos(deltaPhi(*event.met,ele))));
    hist("MT_Lep_Met")->Fill(mt,weight);
    hist("dPhi_Lep_Met")->Fill(deltaPhi(*event.met,ele),weight);

  }
  else if(Nele == 0){   //Fill HT, if Nele = 0, else
                        //reconstruct MLQ and fill MLQmean
    hist("H_T_comb_NoEle")->Fill(ht, weight);
    hist("H_T_comb_NoEle_from350")->Fill(ht, weight);
    hist("H_T_comb_NoEle_from350_rebin")->Fill(ht, weight);
    hist("H_T_comb_NoEle_from350_rebin2")->Fill(ht, weight);
    hist("H_T_comb_NoEle_rebin")->Fill(ht, weight);
    if(ht <= 2000) hist("H_T_comb_NoEle_rebin2")->Fill(ht, weight);
    else hist("H_T_comb_NoEle_rebin2")->Fill(2000., weight);
    if(ht <= 2900) hist("H_T_comb_NoEle_from350_all_filled_rebin")->Fill(ht, weight);
    else hist("H_T_comb_NoEle_from350_all_filled_rebin")->Fill(2900, weight);
    if(Nmuons>0){
      hist("Pt_mu1_NoEle")->Fill(event.muons->at(0).pt(), weight);
      hist("Pt_mu1_NoEle_rebin")->Fill(event.muons->at(0).pt(), weight);
    }
    hist("Integral_NoEle")->Fill(1,weight);
  }

  if(Nele >= 1){
    hist("H_T_comb_1Ele")->Fill(ht, weight);
    hist("H_T_comb_1Ele_from350")->Fill(ht, weight);
    hist("H_T_comb_1Ele_from350_rebin")->Fill(ht, weight);
    hist("H_T_comb_1Ele_from350_rebin2")->Fill(ht, weight);
    hist("H_T_comb_1Ele_rebin")->Fill(ht, weight);
    hist("H_T_comb_1Ele_rebin2")->Fill(ht, weight);
    if(ht <= 2000) hist("H_T_comb_1Ele_rebin3")->Fill(ht, weight);
    else hist("H_T_comb_1Ele_rebin3")->Fill(2000, weight);
    hist("Integral_1Ele")->Fill(1,weight);
  }

  double sum_mu_charge = 0;
  for(const auto & mu : *event.muons) sum_mu_charge += mu.charge();

  bool reconstruct_mlq_mu = (event.electrons->size() == 0 && event.muons->size() == 3 && fabs(sum_mu_charge) == 1 && event.jets->size() >= 2);
  if(reconstruct_mlq_mu){   
    std::vector<LQReconstructionHypothesis> muonic_hyps = event.get(h_muonic_hyps); 
    const LQReconstructionHypothesis* hyp = get_best_hypothesis( muonic_hyps, m_discriminator_name );
    double chi2 = hyp->discriminator(m_discriminator_name);
    hist("chi2_Muonic")->Fill(chi2,weight);
    hist("chi2_MuEle")->Fill(chi2,weight);
    hist("chi2_MuEle_rebin")->Fill(chi2,weight);
    hist("chi2_MuEle_rebin2")->Fill(chi2,weight);
    hist("chi2_MuEle_rebin3")->Fill(chi2,weight);

    double mLQlep_rec = 0;
    double mLQhad_rec = 0;
    double mLQmed_rec = 0;
    double mLQdiff = 0;
    double mLQdiff_rel = 0;
    double mLQLQ = 0;

    if( (hyp->LQlep_v4()).isTimelike() ) {mLQlep_rec = (hyp->LQlep_v4()).M();}
    else {mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());}
    if( (hyp->LQhad_v4()).isTimelike() ) {mLQhad_rec = (hyp->LQhad_v4()).M();}
    else {mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());}
    

    

    mLQdiff = mLQhad_rec - mLQlep_rec;
    mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;
    hist("M_LQ_Muonic_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_Muonic_comb_rebin")->Fill(mLQmed_rec, weight);
    if(mLQmed_rec < 900)   hist("M_LQ_Muonic_comb_rebin2")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_Muonic_comb_rebin2")->Fill(900., weight);
    if(mLQmed_rec < 900)   hist("M_LQ_Muonic_comb_all_filled")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_Muonic_comb_all_filled")->Fill(900., weight);
    hist("M_LQ_Muonic_diff")->Fill(mLQdiff, weight);
    mLQdiff_rel = mLQdiff / mLQmed_rec;
    hist("M_LQ_Muonic_diff_rel")->Fill(mLQdiff_rel,weight);
    mLQLQ = (hyp->LQlep_v4()+hyp->LQhad_v4()).M();
    hist("M_LQLQ_Muonic")->Fill(mLQLQ,weight);
    hist("M_LQLQ_Muonic_rebin")->Fill(mLQLQ,weight);

    double n_jets_had = hyp->tophad_jets().size();
    double n_jets_lep = hyp->toplep_jets().size();
    hist("N_jets_had")->Fill(n_jets_had,weight);
    hist("N_jets_lep")->Fill(n_jets_lep,weight);
    ((TH2F*)hist("n_jets_hadhyp_vs_chi2"))->Fill(n_jets_had,chi2,weight);

    if(n_jets_had==1){
      hist("M_LQ_MuEle_comb_1hadjet")->Fill(mLQmed_rec, weight);
      hist("M_LQ_MuEle_comb_rebin_1hadjet")->Fill(mLQmed_rec, weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2_1hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_rebin2_1hadjet")->Fill(900., weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled_1hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_all_filled_1hadjet")->Fill(900., weight);
      hist("M_tophad_rec_MuEle_1hadjet")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin2")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin3")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_1hadjet_rebin4")->Fill(hyp->tophad_v4().M(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin2")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin3")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin4")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin5")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin6")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin7")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin8")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin9")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin10")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_1hadjet_rebin11")->Fill(hyp->tophad_v4().Pt(),weight);
    }
    else if(n_jets_had==2){
      hist("M_LQ_MuEle_comb_2hadjet")->Fill(mLQmed_rec, weight);
      hist("M_LQ_MuEle_comb_rebin_2hadjet")->Fill(mLQmed_rec, weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2_2hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_rebin2_2hadjet")->Fill(900., weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled_2hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_all_filled_2hadjet")->Fill(900., weight);
      hist("M_tophad_rec_MuEle_2hadjet")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_2hadjet_rebin")->Fill(hyp->tophad_v4().M(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin2")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin3")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin4")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin5")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin6")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin7")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin8")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin9")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin10")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_2hadjet_rebin11")->Fill(hyp->tophad_v4().Pt(),weight);
    }
    else if(n_jets_had==3){
      hist("M_LQ_MuEle_comb_3hadjet")->Fill(mLQmed_rec, weight);
      hist("M_LQ_MuEle_comb_rebin_3hadjet")->Fill(mLQmed_rec, weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2_3hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_rebin2_3hadjet")->Fill(900., weight);
      if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled_3hadjet")->Fill(mLQmed_rec, weight);
      else                   hist("M_LQ_MuEle_comb_all_filled_3hadjet")->Fill(900., weight);
      hist("M_tophad_rec_MuEle_3hadjet")->Fill(hyp->tophad_v4().M(),weight);
      hist("M_tophad_rec_MuEle_3hadjet_rebin")->Fill(hyp->tophad_v4().M(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin2")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin3")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin4")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin5")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin6")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin7")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin8")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin9")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin10")->Fill(hyp->tophad_v4().Pt(),weight);
      hist("Pt_tophad_rec_MuEle_3hadjet_rebin11")->Fill(hyp->tophad_v4().Pt(),weight);
    }


    //combined histograms
    hist("M_LQ_MuEle_comb")->Fill(mLQmed_rec, weight);
    hist("M_LQ_MuEle_comb_rebin")->Fill(mLQmed_rec, weight);
    if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_rebin2")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_MuEle_comb_rebin2")->Fill(900., weight);
    if(mLQmed_rec < 900)   hist("M_LQ_MuEle_comb_all_filled")->Fill(mLQmed_rec, weight);
    else                   hist("M_LQ_MuEle_comb_all_filled")->Fill(900., weight);
    hist("M_LQ_MuEle_diff")->Fill(mLQdiff, weight);
    hist("M_LQ_MuEle_diff_rel")->Fill(mLQdiff_rel,weight);
    hist("M_LQLQ_MuEle")->Fill(mLQLQ,weight);
    hist("M_LQLQ_MuEle_rebin")->Fill(mLQLQ,weight);


    //mT between ele and met
    Particle mu = hyp->muon();
    double mt = sqrt(2*mu.pt()*event.met->pt()*(1-cos(deltaPhi(*event.met,mu))));
    hist("MT_Lep_Met")->Fill(mt,weight);
    hist("dPhi_Lep_Met")->Fill(deltaPhi(*event.met,mu),weight);

  }

  if(!reconstruct_mlq_ele && !reconstruct_mlq_mu){
    //fill HT
    hist("H_T_comb_NoMLQ")->Fill(ht, weight);
    hist("H_T_comb_NoMLQ_from350")->Fill(ht, weight);
    hist("H_T_comb_NoMLQ_from350_rebin")->Fill(ht, weight);
    hist("H_T_comb_NoMLQ_from350_rebin_noweights")->Fill(ht);
    hist("H_T_comb_NoMLQ_from350_rebin2")->Fill(ht, weight);
    hist("H_T_comb_NoMLQ_rebin")->Fill(ht, weight);
    if(ht <= 2000) hist("H_T_comb_NoMLQ_rebin2")->Fill(ht, weight);
    else hist("H_T_comb_NoMLQ_rebin2")->Fill(2000., weight);
    if(ht <= 2900) hist("H_T_comb_NoMLQ_from350_all_filled_rebin")->Fill(ht, weight);
    else hist("H_T_comb_NoMLQ_from350_all_filled_rebin")->Fill(2900, weight);
    hist("Integral_NoMLQ")->Fill(1,weight);
  }

  if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("sum_event_weights_MLQ_MuEle")->Fill(1.,weight);
 

    hist("sum_event_weights")->Fill(1, weight);
    if(is_mc){
      unsigned int n_genp_eletau = 0, n_genp_ele = 0;
      int idx_e = 0;
      for(const auto & ele : *event.electrons){
	double type = -1;
	double dr_min = 9999999;
	double dr_min_genele = 999;
	for(const auto & gp : *event.genparticles){
	  if(fabs(gp.pdgId()) == 11 || fabs(gp.pdgId()) == 15){
	    if(idx_e == 0) n_genp_eletau++;
	    if(idx_e == 0 && fabs(gp.pdgId()) == 11) n_genp_ele++;
	    double dr = deltaR(gp,ele);
	    if(dr < dr_min){
	      dr_min = dr;
	    }
	    if(dr < dr_min_genele && fabs(gp.pdgId()) == 11) dr_min_genele = dr;
	  }
	}
	if(dr_min <= 0.1) type = 0;
	else type = 1;
	hist("ele_type")->Fill(type,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("ele_type_mlq_reco")->Fill(type,weight);
	hist("dr_min_ele_genele")->Fill(dr_min_genele,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("dr_min_ele_genele_mlq_reco")->Fill(dr_min_genele,weight);
	hist("dr_min_ele_geneletau")->Fill(dr_min,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("dr_min_ele_geneletau_mlq_reco")->Fill(dr_min,weight);

	//consider only fake leptons
	if(type == 1){
	  //search for jets within 0.4
	  int idx_matching_jet = -1;
	  for(unsigned int i=0; i<event.jets->size(); i++){
	    double dr = deltaR(ele,event.jets->at(i));
	    if(dr < 0.4){
	      idx_matching_jet = i;
	    }
	  }

	  //-1 is filled in case no jet could be matched to the fake lepton
	  double faking_pt = -1;
	  if(idx_matching_jet > -1) faking_pt = event.jets->at(idx_matching_jet).pt();
	  hist("jets_faking_ele_pt")->Fill(faking_pt,weight);
	}
	idx_e++;
      }

      if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("more_gen_ele_mlq_reco")->Fill(n_genp_eletau >= event.electrons->size(),weight);

      unsigned int n_genp_mutau = 0, n_genp_mu = 0;
      int idx_mu=0;
      for(const auto & mu : *event.muons){
	double type = -1;
	double dr_min = 9999999;
	double dr_min_genmu = 999;
	for(const auto & gp : *event.genparticles){
	  if(fabs(gp.pdgId()) == 13 || fabs(gp.pdgId()) == 15){
	    if(idx_mu==0)                           n_genp_mutau++;
	    if(idx_mu==0 && fabs(gp.pdgId()) == 13) n_genp_mu++;
	    double dr = deltaR(gp,mu);
	    if(dr < dr_min){
	      dr_min = dr;
	    }
	    if(dr < dr_min_genmu && fabs(gp.pdgId()) == 13) dr_min_genmu = dr;
	  }
	}
	if(dr_min <= 0.1) type = 0;
	else type = 1;
	hist("mu_type")->Fill(type,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("mu_type_mlq_reco")->Fill(type,weight);
	hist("dr_min_mu_genmu")->Fill(dr_min_genmu,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("dr_min_mu_genmu_mlq_reco")->Fill(dr_min_genmu,weight);
	hist("dr_min_mu_genmutau")->Fill(dr_min,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("dr_min_mu_genmutau_mlq_reco")->Fill(dr_min,weight);

	//consider only fake leptons
	if(type == 1){
	  //search for jets within 0.4
	  int idx_matching_jet = -1;
	  for(unsigned int i=0; i<event.jets->size(); i++){
	    double dr = deltaR(mu,event.jets->at(i));
	    if(dr < 0.4){
	      idx_matching_jet = i;
	    }
	  }

	  //-1 is filled in case no jet could be matched to the fake lepton
	  double faking_pt = -1;
	  if(idx_matching_jet > -1) faking_pt = event.jets->at(idx_matching_jet).pt();
	  hist("jets_faking_mu_pt")->Fill(faking_pt,weight);
	}
	idx_mu++;
      }
      
      if(reconstruct_mlq_ele || reconstruct_mlq_mu) {
	hist("more_gen_mu_mlq_reco")->Fill(n_genp_mutau >= event.muons->size(),weight);
      }

      //look at muons from b-quark decays --> dr_min between a rec muon and a gen b-quark for all muons in the event
      for(const auto & mu : *event.muons){
	double dr_min= 99999;
	for(const auto & gp : *event.genparticles){
	  if(fabs(gp.pdgId()) == 5){
	    double dr = deltaR(mu,gp);
	    if(dr < dr_min) dr_min = dr;
	  }
	}
	hist("dr_min_mu_gen_b")->Fill(dr_min,weight);
	if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("dr_min_mu_gen_b_MLQ")->Fill(dr_min,weight);
	else                                          hist("dr_min_mu_gen_b_NoMLQ")->Fill(dr_min,weight);
      }

      vector<bool> is_fake_mu;
      unsigned int n_genmu = 0;
      for(const auto & gp : *event.genparticles){
	if(fabs(gp.pdgId()) != 13) continue;
	n_genmu++;
      }

      if(n_genmu != event.muons->size()){
	//if ngen and nreco are unequal, try to match muons to muons within 0.1 and to taus within 0.2
	unsigned int n_matched_to_muons = 0, n_matched_to_taus = 0, n_matched_to_b = 0, n_matched_to_c = 0;
	for(const auto & mu : *event.muons){
	  bool is_matched = false;
	  for(const auto & gp : *event.genparticles){
	    if(fabs(gp.pdgId()) == 13){
	      if(deltaR(gp,mu) < 0.1 && !is_matched){
		is_matched = true;
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
		n_matched_to_b++;
	      }
	    }
	  }
	  is_fake_mu.push_back(!is_matched);
	}
	if(n_matched_to_taus + n_matched_to_muons != event.muons->size()){
	  //find unmatched muons
	  for(unsigned int i=0; i<event.muons->size() - n_matched_to_taus - n_matched_to_muons; i++){
	    hist("n_fake_mu")->Fill(0.,weight);
	    if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("n_fake_mu_MLQ")->Fill(0.,weight);
	    else hist("n_fake_mu_NoMLQ")->Fill(0.,weight);
	  }
	}
	if(n_matched_to_taus + n_matched_to_muons + n_matched_to_b + n_matched_to_c != event.muons->size()){
	  //find unmatched muons
	  for(unsigned int i=0; i<event.muons->size() - n_matched_to_taus - n_matched_to_muons - n_matched_to_b - n_matched_to_c; i++){
	    hist("n_fake_mu_HF")->Fill(0.,weight);
	    if(reconstruct_mlq_ele || reconstruct_mlq_mu) hist("n_fake_mu_MLQ_HF")->Fill(0.,weight);
	    else hist("n_fake_mu_NoMLQ_HF")->Fill(0.,weight);
	  }
	}
      }
    }


} //Methode



LQToTopMuHists::~LQToTopMuHists(){}
