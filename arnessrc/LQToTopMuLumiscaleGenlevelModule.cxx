#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuPreselectionHists.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuLumiscaleGenlevelModule: public AnalysisModule {
  public:

    explicit LQToTopMuLumiscaleGenlevelModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    unique_ptr<CommonModules> common;
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts,
    h_catA, h_jets_catA, h_ele_catA, h_mu_catA, h_event_catA, h_topjets_catA, h_lumi_catA,
    h_catB, h_jets_catB, h_ele_catB, h_mu_catB, h_event_catB, h_topjets_catB, h_lumi_catB,
    h_2muons, h_jets_2muons, h_ele_2muons, h_mu_2muons, h_event_2muons, h_topjets_2muons, h_lumi_2muons,
    h_finalselection, h_jets_finalselection, h_ele_finalselection, h_mu_finalselection, h_event_finalselection, h_topjets_finalselection, h_lumi_finalselection;

    double muonpt, muoneta, jetpt, jeteta, electronpt, electroneta;

  };


  LQToTopMuLumiscaleGenlevelModule::LQToTopMuLumiscaleGenlevelModule(Context & ctx){

    cout << "Hello from LQToTopMuLumiscaleGenlevelModule!" << endl;


    muonpt = 30;
    muoneta = 2.4;
    jetpt = 30;
    jeteta = 2.4;
    electronpt = 30;
    electroneta = 2.4;




    common.reset(new CommonModules());
    common->init(ctx);

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuPreselectionHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "Topjets_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));


    h_2muons.reset(new LQToTopMuPreselectionHists(ctx, "2Muons"));
    h_jets_2muons.reset(new JetHists(ctx, "Jets_2Muons"));
    h_ele_2muons.reset(new ElectronHists(ctx, "Ele_2Muons"));
    h_mu_2muons.reset(new MuonHists(ctx, "Mu_2Muons"));
    h_event_2muons.reset(new EventHists(ctx, "Event_2Muons"));
    h_topjets_2muons.reset(new TopJetHists(ctx, "Topjets_2Muons"));
    h_lumi_2muons.reset(new LuminosityHists(ctx, "Lumi_2Muons"));



    h_finalselection.reset(new LQToTopMuPreselectionHists(ctx, "FinalSelection"));
    h_jets_finalselection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalselection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalselection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_event_finalselection.reset(new EventHists(ctx, "Event_FinalSelection"));
    h_topjets_finalselection.reset(new TopJetHists(ctx, "Topjets_FinalSelection"));
    h_lumi_finalselection.reset(new LuminosityHists(ctx, "Lumi_FinalSelection"));

    h_catA.reset(new LQToTopMuPreselectionHists(ctx, "CatA"));
    h_jets_catA.reset(new JetHists(ctx, "Jets_CatA"));
    h_ele_catA.reset(new ElectronHists(ctx, "Ele_CatA"));
    h_mu_catA.reset(new MuonHists(ctx, "Mu_CatA"));
    h_event_catA.reset(new EventHists(ctx, "Event_CatA"));
    h_topjets_catA.reset(new TopJetHists(ctx, "Topjets_CatA"));
    h_lumi_catA.reset(new LuminosityHists(ctx, "Lumi_CatA"));


    h_catB.reset(new LQToTopMuPreselectionHists(ctx, "CatB"));
    h_jets_catB.reset(new JetHists(ctx, "Jets_CatB"));
    h_ele_catB.reset(new ElectronHists(ctx, "Ele_CatB"));
    h_mu_catB.reset(new MuonHists(ctx, "Mu_CatB"));
    h_event_catB.reset(new EventHists(ctx, "Event_CatB"));
    h_topjets_catB.reset(new TopJetHists(ctx, "Topjets_CatB"));
    h_lumi_catB.reset(new LuminosityHists(ctx, "Lumi_CatB"));

  }


  bool LQToTopMuLumiscaleGenlevelModule::process(Event & event) {

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    common->process(event);

    // 2 Muons
    int n_muons = 0, n_electrons = 0;
    for(unsigned int i=0; i<event.genparticles->size(); i++){
      auto gp = event.genparticles->at(i);
      if(abs(gp.pdgId()) == 13){

        // only gen-muons left
        if(gp.pt() < muonpt) continue;
        if(gp.eta() > muoneta) continue;
        n_muons++;
      }
      else if(abs(gp.pdgId()) == 11){

        // only gen-electrons left
        if(gp.pt() < electronpt) continue;
        if(gp.eta() > electroneta) continue;
        n_electrons++;

      }
    }
    if(n_muons < 2) return false;
    cout << "Number of genmuons: " << n_muons << endl;

    h_2muons->fill(event);
    h_jets_2muons->fill(event);
    h_ele_2muons->fill(event);
    h_mu_2muons->fill(event);
    h_event_2muons->fill(event);
    h_topjets_2muons->fill(event);
    h_lumi_2muons->fill(event);

    // Require m_ll > 111
    for(unsigned int i=0; i<event.genparticles->size(); i++){
      for(unsigned int j=0; j<event.genparticles->size(); j++){
        if(i >= j) continue;
        auto gp1 = event.genparticles->at(i);
        auto gp2 = event.genparticles->at(j);
        if(!( (abs(gp1.pdgId()) == 13) && (abs(gp2.pdgId()) == 13) )) continue;
        cout << "These two muons -- pt: " << gp1.pt() << " and " << gp2.pt() << endl;
        cout << "eta: " << gp1.eta() << " and " << gp2.eta() << endl;
        if(gp1.pt() < muonpt || gp2.pt() < muonpt) continue;
        if(gp1.eta() > muoneta || gp2.eta() > muoneta) continue;

        cout << "survived kinematic thresholds! Invariant mass: " << (gp1.v4() + gp2.v4()).M() << endl;

        //Here we have 2 muons. As soon as we find a pair with M_ll < 111: return false
        if((gp1.v4() + gp2.v4()).M() < 111) return false;
        cout << "survived" << endl;
      }
    }


    // 2 Jets
    int n_jets = 0;
    for(unsigned int i=0; i<event.genjets->size(); i++){
      auto gp = event.genjets->at(i);

      if(gp.pt() < jetpt) continue;
      if(gp.eta() > jeteta) continue;
      n_jets++;
    }
    if(n_jets < 2) return false;


    h_finalselection->fill(event);
    h_jets_finalselection->fill(event);
    h_ele_finalselection->fill(event);
    h_mu_finalselection->fill(event);
    h_event_finalselection->fill(event);
    h_topjets_finalselection->fill(event);
    h_lumi_finalselection->fill(event);


    // categorize for A, B

    if(n_muons >= 3 || n_electrons >= 1){
      //fill histos for A

      h_catA->fill(event);
      h_jets_catA->fill(event);
      h_ele_catA->fill(event);
      h_mu_catA->fill(event);
      h_event_catA->fill(event);
      h_topjets_catA->fill(event);
      h_lumi_catA->fill(event);

    }
    else{
      //fill histos for B

      h_catB->fill(event);
      h_jets_catB->fill(event);
      h_ele_catB->fill(event);
      h_mu_catB->fill(event);
      h_event_catB->fill(event);
      h_topjets_catB->fill(event);
      h_lumi_catB->fill(event);
    }


    return false;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuLumiscaleGenlevelModule)

}
