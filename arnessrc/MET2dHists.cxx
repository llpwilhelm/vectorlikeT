#include "UHH2/LQToTopMu/include/MET2dHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

MET2dHists::MET2dHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  book<TH2F>("met_ptmu", "missing E_{T} vs p_{T}^{#mu}", 50, 0, 1000, 200, 0, 1000);  
  book<TH2F>("met_ptjet", "missing E_{T} vs p_{T}^{jet}", 50, 0, 1000, 40, 0, 1500);
  book<TH2F>("met_ptmu1", "missing E_{T} vs p_{T}^{leading #mu}", 50, 0, 1000, 200, 0, 1000);
  book<TH2F>("met_ptjet1", "missing E_{T} vs p_{T}^{leading jet}", 50, 0, 1000, 40, 0, 1500);
  book<TH2F>("met_etamu", "missing E_{T} vs #eta^{#mu}", 50, 0, 1000, 40, -2.5, 2.5);
  book<TH2F>("met_etajet", "missing E_{T} vs #eta^{jet}", 50, 0, 1000, 40, -2.5, 2.5);
  book<TH2F>("met_ht", "missing E_{T} vs H_{T}", 50, 0, 1000, 100, 0, 5000);
} // end of constructor


void MET2dHists::fill(const Event & event){

  double weight = event.weight;
  double met = event.met->pt();
  double ptmu1 = event.muons->at(0).pt();
  double ptjet1 = event.jets->at(0).pt();

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    ht += jet.pt();
  }
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
  }
  ht = ht_lep + ht_jets + met;

 

  ((TH2F*)hist("met_ptmu1"))->Fill(met, ptmu1, weight);
  ((TH2F*)hist("met_ptjet1"))->Fill(met, ptjet1, weight);
  ((TH2F*)hist("met_ht"))->Fill(met, ht, weight);


 for (const Muon & thismu : *event.muons){
   double pttmp = thismu.pt();
   double etatmp = thismu.eta();
   ((TH2F*)hist("met_ptmu"))->Fill(met, pttmp, weight);
   ((TH2F*)hist("met_etamu"))->Fill(met, etatmp, weight);
  }
 for (const Jet & thisjet : *event.jets){
   double pttmp = thisjet.pt();
   double etatmp = thisjet.eta();
   ((TH2F*)hist("met_ptjet"))->Fill(met, pttmp, weight);
   ((TH2F*)hist("met_etajet"))->Fill(met, etatmp, weight);
   }


} // end of "fill"

MET2dHists::~MET2dHists(){}
