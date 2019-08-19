#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

HT2dHists::HT2dHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  book<TH2F>("ht_ptmu", "H_{T} vs p_{T}^{#mu}", 100, 0, 5000, 200, 0, 1000);  
  book<TH2F>("ht_ptjet", "H_{T} vs p_{T}^{jet}", 100, 0, 5000, 40, 0, 1500);
  book<TH2F>("ht_ptmu1", "H_{T} vs p_{T}^{leading #mu}", 100, 0, 5000, 200, 0, 1000);
  book<TH2F>("ht_ptjet1", "H_{T} vs p_{T}^{leading jet}", 100, 0, 5000, 40, 0, 1500);
  book<TH2F>("ht_etamu", "H_{T} vs #eta^{#mu}", 100, 0, 5000, 40, -2.5, 2.5);
  book<TH2F>("ht_etajet", "H_{T} vs #eta^{jet}", 100, 0, 5000, 40, -2.5, 2.5);
  book<TH2F>("ht_met", "H_{T} vs missing E_{T}", 100, 0, 5000, 50, 0, 1000);
  book<TH2F>("ht_Nmu", "H_{T} vs missing N_{#mu}", 100, 0, 5000, 6, -0.5, 5.5);

} // end of constructor


void HT2dHists::fill(const Event & event){

  double weight = event.weight;
  double met = event.met->pt();
  double ptmu1 = 0;
  if(event.muons->size()>0){
    ptmu1 = event.muons->at(0).pt();
  }
  double ptjet1 = event.jets->at(0).pt();

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
  /*for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
  }*/
  ht = ht_lep + ht_jets + met;

 
   if(event.muons->size()>0){
     ((TH2F*)hist("ht_ptmu1"))->Fill(ht, ptmu1, weight);
   }
  ((TH2F*)hist("ht_ptjet1"))->Fill(ht, ptjet1, weight);
  ((TH2F*)hist("ht_met"))->Fill(ht, met, weight);
  ((TH2F*)hist("ht_Nmu"))->Fill(ht, event.muons->size(), weight);


 for (const Muon & thismu : *event.muons){
   double pttmp = thismu.pt();
   double etatmp = thismu.eta();
   ((TH2F*)hist("ht_ptmu"))->Fill(ht, pttmp, weight);
   ((TH2F*)hist("ht_etamu"))->Fill(ht, etatmp, weight);
  }
 for (const Jet & thisjet : *event.jets){
   double pttmp = thisjet.pt();
   double etatmp = thisjet.eta();
   ((TH2F*)hist("ht_ptjet"))->Fill(ht, pttmp, weight);
   ((TH2F*)hist("ht_etajet"))->Fill(ht, etatmp, weight);
   }


} // end of "fill"

HT2dHists::~HT2dHists(){}
