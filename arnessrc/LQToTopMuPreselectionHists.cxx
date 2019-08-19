#include "UHH2/LQToTopMu/include/LQToTopMuPreselectionHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include <math.h>

#include "TH1F.h"
#include "TH1D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMuPreselectionHists::LQToTopMuPreselectionHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1D>("N_jets", "N_{jets}", 20, 0, 20);  
  book<TH1D>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1D>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1D>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1D>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);
  book<TH1D>("pt_jet1", "#p_{T}^{jet 1}", 100, 0, 3000);
  book<TH1D>("pt_jet2", "#p_{T}^{jet 2}", 100, 0, 3000);
  book<TH1D>("pt_jet3", "#p_{T}^{jet 3}", 100, 0, 3000);

  // leptons
  book<TH1D>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1D>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  book<TH1D>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1D>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  // general
  book<TH1D>("N_pv", "N_{PV}", 50, 0, 50);
  book<TH1D>("H_T", "H_{T}", 100, 0, 5000);
  book<TH1D>("sum_event_weights", "BinContent = sum(eventweights)", 1, 0.5, 1.5);
  book<TH1D>("mindr_mujet", "min(#Delta R) between a muon and a jet", 50,0,5);
  book<TH1D>("mindr_elejet", "min(#Delta R) between an electron and a jet", 50,0,5);
}


void LQToTopMuPreselectionHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;

  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  hist("N_jets")->Fill(Njets, weight);
  
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



  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      hist("pt_mu")->Fill(thismu.pt(), weight);
      hist("eta_mu")->Fill(thismu.eta(), weight);
      hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }
  
  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

  //HT
  auto met = event.met->pt();
  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
  }
  //Bedeutung der for-Schleife
  /*const auto jets = event.jets;
    for(unsigned int i=0; i<jets.size();i++){
    Jet jet=jets[i];
    ht +=jet.pt();
    }*/
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }

  ht = ht_lep + ht_jets + met;
  hist("H_T")->Fill(ht, weight);

  hist("sum_event_weights")->Fill(1., weight);

  double mindr_mujet = 999;
  double mindr_elejet = 999;
  for(const Muon & thismu : *event.muons){
    for(const Jet & thisjet : *event.jets){
      double dr = deltaR(thismu,thisjet);
      if (dr < mindr_mujet) mindr_mujet = dr;
    }
    hist("mindr_mujet")->Fill(mindr_mujet,weight);
  }

  for(const Electron & thisele : *event.electrons){
    for(const Jet & thisjet : *event.jets){
      double dr = deltaR(thisele,thisjet);
      if (dr < mindr_elejet) mindr_elejet = dr;
    }
    hist("mindr_elejet")->Fill(mindr_elejet,weight);
  }





} //Methode


LQToTopMuPreselectionHists::~LQToTopMuPreselectionHists(){}
