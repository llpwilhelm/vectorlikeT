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

  class LQToTopMuPreselectionCrossCheckModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuPreselectionCrossCheckModule(Context & ctx);
    virtual bool process(Event & event) override;
  
  private:
  
    unique_ptr<CommonModules> common;
    unique_ptr<JetCleaner> jetcleaner;
  
    // declare the Selections to use.
    unique_ptr<Selection> muon_2loose_nopt_sel, muon_2tight_nopt_sel, muon_2tight_sel, trigger_singlemu1, trigger_singlemu2, trigger_doublemu1, trigger_doublemu2, trigger_doublemu3, trigger_doublemu4;
  
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts,
      h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_lumi_cleaner, 
      h_2mu30, h_jets_2mu30, h_ele_2mu30, h_mu_2mu30, h_event_2mu30, h_topjets_2mu30, h_lumi_2mu30, 
      h_trigger_singlemu, h_jets_trigger_singlemu, h_ele_trigger_singlemu, h_mu_trigger_singlemu, h_event_trigger_singlemu, h_topjets_trigger_singlemu, h_lumi_trigger_singlemu,
      h_trigger_combination, h_jets_trigger_combination, h_ele_trigger_combination, h_mu_trigger_combination, h_event_trigger_combination, h_topjets_trigger_combination, h_lumi_trigger_combination,
      h_trigger_singlemu_2mu30, h_jets_trigger_singlemu_2mu30, h_ele_trigger_singlemu_2mu30, h_mu_trigger_singlemu_2mu30, h_event_trigger_singlemu_2mu30, h_topjets_trigger_singlemu_2mu30, h_lumi_trigger_singlemu_2mu30,
      h_trigger_combination_2mu30, h_jets_trigger_combination_2mu30, h_ele_trigger_combination_2mu30, h_mu_trigger_combination_2mu30, h_event_trigger_combination_2mu30, h_topjets_trigger_combination_2mu30, h_lumi_trigger_combination_2mu30,
      h_2mu_loose_nopt, h_jets_2mu_loose_nopt, h_ele_2mu_loose_nopt, h_mu_2mu_loose_nopt, h_event_2mu_loose_nopt, h_topjets_2mu_loose_nopt, h_lumi_2mu_loose_nopt,
      h_2mu_tight_nopt, h_jets_2mu_tight_nopt, h_ele_2mu_tight_nopt, h_mu_2mu_tight_nopt, h_event_2mu_tight_nopt, h_topjets_2mu_tight_nopt, h_lumi_2mu_tight_nopt;

    MuonId MuId_Loose_NoPt, MuId_Tight_NoPt, MuId_Tight;
    JetId Jet_ID;

  };


  LQToTopMuPreselectionCrossCheckModule::LQToTopMuPreselectionCrossCheckModule(Context & ctx){
    
    cout << "Hello from LQToTopMuPreselectionCrossCheckModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    MuId_Loose_NoPt = AndId<Muon>(MuonIDLoose(), PtEtaCut(0.,2.4), MuonIso(0.15));
    MuId_Tight_NoPt = AndId<Muon>(MuonIDTight(), PtEtaCut(0.,2.4), MuonIso(0.15));
    MuId_Tight = AndId<Muon>(MuonIDTight(), PtEtaCut(30.,2.4), MuonIso(0.15));
    Jet_ID = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));

    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_muon_id(MuId_Loose_NoPt);
    common->init(ctx);
    jetcleaner.reset(new JetCleaner(ctx,Jet_ID));

    // 2. set up selections

    //Selection
    muon_2loose_nopt_sel.reset(new NMuonSelection(2, -1, MuId_Loose_NoPt)); 
    muon_2tight_nopt_sel.reset(new NMuonSelection(2, -1, MuId_Tight_NoPt)); 
    muon_2tight_sel.reset(new NMuonSelection(2, -1, MuId_Tight)); 

    trigger_singlemu1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trigger_singlemu2.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    trigger_doublemu1.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*"));
    trigger_doublemu2.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"));
    trigger_doublemu3.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"));
    trigger_doublemu4.reset(new TriggerSelection("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"));

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuPreselectionHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "Topjets_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_cleaner.reset(new LQToTopMuPreselectionHists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_Cleaner"));
    h_topjets_cleaner.reset(new TopJetHists(ctx, "Topjets_Cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_Cleaner"));

    h_2mu30.reset(new LQToTopMuPreselectionHists(ctx, "2mu30"));
    h_jets_2mu30.reset(new JetHists(ctx, "Jets_2mu30"));
    h_ele_2mu30.reset(new ElectronHists(ctx, "Ele_2mu30"));
    h_mu_2mu30.reset(new MuonHists(ctx, "Mu_2mu30"));
    h_event_2mu30.reset(new EventHists(ctx, "Event_2mu30"));
    h_topjets_2mu30.reset(new TopJetHists(ctx, "Topjets_2mu30"));
    h_lumi_2mu30.reset(new LuminosityHists(ctx, "Lumi_2mu30"));

    h_trigger_singlemu.reset(new LQToTopMuPreselectionHists(ctx, "trigger_singlemu"));
    h_jets_trigger_singlemu.reset(new JetHists(ctx, "Jets_trigger_singlemu"));
    h_ele_trigger_singlemu.reset(new ElectronHists(ctx, "Ele_trigger_singlemu"));
    h_mu_trigger_singlemu.reset(new MuonHists(ctx, "Mu_trigger_singlemu"));
    h_event_trigger_singlemu.reset(new EventHists(ctx, "Event_trigger_singlemu"));
    h_topjets_trigger_singlemu.reset(new TopJetHists(ctx, "Topjets_trigger_singlemu"));
    h_lumi_trigger_singlemu.reset(new LuminosityHists(ctx, "Lumi_trigger_singlemu"));

    h_trigger_combination.reset(new LQToTopMuPreselectionHists(ctx, "trigger_combination"));
    h_jets_trigger_combination.reset(new JetHists(ctx, "Jets_trigger_combination"));
    h_ele_trigger_combination.reset(new ElectronHists(ctx, "Ele_trigger_combination"));
    h_mu_trigger_combination.reset(new MuonHists(ctx, "Mu_trigger_combination"));
    h_event_trigger_combination.reset(new EventHists(ctx, "Event_trigger_combination"));
    h_topjets_trigger_combination.reset(new TopJetHists(ctx, "Topjets_trigger_combination"));
    h_lumi_trigger_combination.reset(new LuminosityHists(ctx, "Lumi_trigger_combination"));

    h_trigger_singlemu_2mu30.reset(new LQToTopMuPreselectionHists(ctx, "trigger_singlemu_2mu30"));
    h_jets_trigger_singlemu_2mu30.reset(new JetHists(ctx, "Jets_trigger_singlemu_2mu30"));
    h_ele_trigger_singlemu_2mu30.reset(new ElectronHists(ctx, "Ele_trigger_singlemu_2mu30"));
    h_mu_trigger_singlemu_2mu30.reset(new MuonHists(ctx, "Mu_trigger_singlemu_2mu30"));
    h_event_trigger_singlemu_2mu30.reset(new EventHists(ctx, "Event_trigger_singlemu_2mu30"));
    h_topjets_trigger_singlemu_2mu30.reset(new TopJetHists(ctx, "Topjets_trigger_singlemu_2mu30"));
    h_lumi_trigger_singlemu_2mu30.reset(new LuminosityHists(ctx, "Lumi_trigger_singlemu_2mu30"));

    h_trigger_combination_2mu30.reset(new LQToTopMuPreselectionHists(ctx, "trigger_combination_2mu30"));
    h_jets_trigger_combination_2mu30.reset(new JetHists(ctx, "Jets_trigger_combination_2mu30"));
    h_ele_trigger_combination_2mu30.reset(new ElectronHists(ctx, "Ele_trigger_combination_2mu30"));
    h_mu_trigger_combination_2mu30.reset(new MuonHists(ctx, "Mu_trigger_combination_2mu30"));
    h_event_trigger_combination_2mu30.reset(new EventHists(ctx, "Event_trigger_combination_2mu30"));
    h_topjets_trigger_combination_2mu30.reset(new TopJetHists(ctx, "Topjets_trigger_combination_2mu30"));
    h_lumi_trigger_combination_2mu30.reset(new LuminosityHists(ctx, "Lumi_trigger_combination_2mu30"));

    h_2mu_loose_nopt.reset(new LQToTopMuPreselectionHists(ctx, "2mu_loose_nopt"));
    h_jets_2mu_loose_nopt.reset(new JetHists(ctx, "Jets_2mu_loose_nopt"));
    h_ele_2mu_loose_nopt.reset(new ElectronHists(ctx, "Ele_2mu_loose_nopt"));
    h_mu_2mu_loose_nopt.reset(new MuonHists(ctx, "Mu_2mu_loose_nopt"));
    h_event_2mu_loose_nopt.reset(new EventHists(ctx, "Event_2mu_loose_nopt"));
    h_topjets_2mu_loose_nopt.reset(new TopJetHists(ctx, "Topjets_2mu_loose_nopt"));
    h_lumi_2mu_loose_nopt.reset(new LuminosityHists(ctx, "Lumi_2mu_loose_nopt"));

    h_2mu_tight_nopt.reset(new LQToTopMuPreselectionHists(ctx, "2mu_tight_nopt"));
    h_jets_2mu_tight_nopt.reset(new JetHists(ctx, "Jets_2mu_tight_nopt"));
    h_ele_2mu_tight_nopt.reset(new ElectronHists(ctx, "Ele_2mu_tight_nopt"));
    h_mu_2mu_tight_nopt.reset(new MuonHists(ctx, "Mu_2mu_tight_nopt"));
    h_event_2mu_tight_nopt.reset(new EventHists(ctx, "Event_2mu_tight_nopt"));
    h_topjets_2mu_tight_nopt.reset(new TopJetHists(ctx, "Topjets_2mu_tight_nopt"));
    h_lumi_2mu_tight_nopt.reset(new LuminosityHists(ctx, "Lumi_2mu_tight_nopt"));
  }


  bool LQToTopMuPreselectionCrossCheckModule::process(Event & event) {

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_lumi_nocuts->fill(event);



    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);

    h_cleaner->fill(event);
    h_jets_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_event_cleaner->fill(event);
    h_topjets_cleaner->fill(event);
    h_lumi_cleaner->fill(event);

    if(muon_2tight_sel->passes(event)){
      h_2mu30->fill(event);
      h_jets_2mu30->fill(event);
      h_ele_2mu30->fill(event);
      h_mu_2mu30->fill(event);
      h_event_2mu30->fill(event);
      h_topjets_2mu30->fill(event);
      h_lumi_2mu30->fill(event);
    }

    bool singlemu = trigger_singlemu1->passes(event) || trigger_singlemu2->passes(event);
    bool combination = singlemu || trigger_doublemu1->passes(event) || trigger_doublemu2->passes(event) || trigger_doublemu3->passes(event) || trigger_doublemu4->passes(event);

    if(singlemu){
      h_trigger_singlemu->fill(event);
      h_jets_trigger_singlemu->fill(event);
      h_ele_trigger_singlemu->fill(event);
      h_mu_trigger_singlemu->fill(event);
      h_event_trigger_singlemu->fill(event);
      h_topjets_trigger_singlemu->fill(event);
      h_lumi_trigger_singlemu->fill(event);
    }

    if(combination){
      h_trigger_combination->fill(event);
      h_jets_trigger_combination->fill(event);
      h_ele_trigger_combination->fill(event);
      h_mu_trigger_combination->fill(event);
      h_event_trigger_combination->fill(event);
      h_topjets_trigger_combination->fill(event);
    }

    if(muon_2tight_sel->passes(event)){

      if(singlemu){
	h_trigger_singlemu_2mu30->fill(event);
	h_jets_trigger_singlemu_2mu30->fill(event);
	h_ele_trigger_singlemu_2mu30->fill(event);
	h_mu_trigger_singlemu_2mu30->fill(event);
	h_event_trigger_singlemu_2mu30->fill(event);
	h_topjets_trigger_singlemu_2mu30->fill(event);
	h_lumi_trigger_singlemu_2mu30->fill(event);
      }

      if(combination){
	h_trigger_combination_2mu30->fill(event);
	h_jets_trigger_combination_2mu30->fill(event);
	h_ele_trigger_combination_2mu30->fill(event);
	h_mu_trigger_combination_2mu30->fill(event);
	h_event_trigger_combination_2mu30->fill(event);
	h_topjets_trigger_combination_2mu30->fill(event);
      }

    }
    /*
    if(muon_2loose_nopt_sel->passes(event)){
      h_2mu_loose_nopt->fill(event);
      h_jets_2mu_loose_nopt->fill(event);
      h_ele_2mu_loose_nopt->fill(event);
      h_mu_2mu_loose_nopt->fill(event);
      h_event_2mu_loose_nopt->fill(event);
      h_topjets_2mu_loose_nopt->fill(event);
      h_lumi_2mu_loose_nopt->fill(event);
    }

    if(muon_2tight_nopt_sel->passes(event)){
      h_2mu_tight_nopt->fill(event);
      h_jets_2mu_tight_nopt->fill(event);
      h_ele_2mu_tight_nopt->fill(event);
      h_mu_2mu_tight_nopt->fill(event);
      h_event_2mu_tight_nopt->fill(event);
      h_topjets_2mu_tight_nopt->fill(event);
      h_lumi_2mu_tight_nopt->fill(event);
    }
    */
    return false;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuPreselectionCrossCheckModule)
  
}
