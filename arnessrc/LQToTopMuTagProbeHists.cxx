#include "../include/LQToTopMuTagProbeHists.h"
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

LQToTopMuTagProbeHists::LQToTopMuTagProbeHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  
  // book all histograms here
  book<TH1F>("pt_ele_fine", "p_{T}^{ele} [GeV]", 200, 0, 1000);  
  double bins[13] = {0,30,65,100,150,200,250,300,350,400,500,700,1000};
  book<TH1F>("pt_ele_binned", "p_{T}^{ele} [GeV]", 12, bins); 
  book<TH1F>("eta_ele", "#eta^{ele}", 24, -2.4, 2.4);  
  book<TH1F>("eta_ele_binned", "#eta^{ele}", 8, -2.4, 2.4);  
  book<TH1D>("dR_ele_mu", "#Delta R(ele, mu)", 50, 0, 5);
  book<TH1D>("dRmin_ele_jet", "#Delta R^{min}(ele, jet)", 50, 0, 5);
  book<TH1D>("dRmin_ele_obj", "#Delta R^{min}(ele, obj)", 50, 0, 5);
}


void LQToTopMuTagProbeHists::fill(const Event & event){
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  if(event.electrons->size() != 1) throw runtime_error("In LQToTopMuTagProbeHists.cxx: Not ==1 electron --> this is needed for Tag&Probe selection.");
  double pt = event.electrons->at(0).pt();
  double eta = event.electrons->at(0).eta();

  hist("pt_ele_fine")->Fill(pt, weight);
  hist("pt_ele_binned")->Fill(pt, weight);
  hist("eta_ele")->Fill(eta, weight);
  hist("eta_ele_binned")->Fill(eta, weight);

  double dR_e_mu = deltaR(event.electrons->at(0),event.muons->at(0));
  hist("dR_ele_mu")->Fill(dR_e_mu,weight);

  double dR_min_e_jet = 999999;
  for(unsigned int i=0; i<event.jets->size(); i++){
    double dr = deltaR(event.electrons->at(0),event.jets->at(i));
    if(dr < dR_min_e_jet) dR_min_e_jet = dr;
  }
  hist("dRmin_ele_jet")->Fill(dR_min_e_jet,weight);

  double dR_min_e_obj = min(dR_e_mu, dR_min_e_jet);
  hist("dRmin_ele_obj")->Fill(dR_min_e_obj,weight);

} //Methode



LQToTopMuTagProbeHists::~LQToTopMuTagProbeHists(){}
