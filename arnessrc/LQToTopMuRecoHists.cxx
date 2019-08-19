#include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include <math.h>

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

//!only fill histograms of this class after all nessecary recomodules have processed the event in the analysis cycle and if all requirements for reconstrucions (e.g. >= 1 electron etc.) are met, if not required explicitly before filling!

LQToTopMuRecoHists::LQToTopMuRecoHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  

book<TH1F>("MLQ_HT_Mix", "M_{LQ,mean} & H_{T}", 100, 0, 5000);


}

void LQToTopMuRecoHists::fill(const Event & event){

  double weight = event.weight;
  hist("MLQ_HT_Mix")->Fill(1,weight);
 
}




LQToTopMuRecoHists::~LQToTopMuRecoHists(){}



