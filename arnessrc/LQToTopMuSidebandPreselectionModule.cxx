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

  class LQToTopMuSidebandPreselectionModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuSidebandPreselectionModule(Context & ctx);
    virtual bool process(Event & event) override;
  
  private:
  
    std::unique_ptr<CommonModules> common;
    //std::unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer;
  
    std::unique_ptr<JetCleaner> jetcleaner;
    std::unique_ptr<JetLeptonCleaner> jetleptoncleaner;
    std::unique_ptr<MuonCleaner> muoncleaner;
    std::unique_ptr<MuonCleaner> muoncleaner_tight;
    std::unique_ptr<ElectronCleaner> electroncleaner;

    std::unique_ptr<AnalysisModule> SF_muonID;

  
    // declare the Selections to use.
    std::unique_ptr<Selection> njet_sel, ht_sel, lumi_sel, mu_trigger_sel1, mu_trigger_sel2, ele_trigger_sel1, ele_trigger_sel2, lep1_sel, lep2_sel, mu_veto;
  
    // store the Hists collection as member variables. 
    std::unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts, 
      h_trigger, h_jets_trigger, h_ele_trigger, h_mu_trigger, h_event_trigger, h_topjets_trigger, h_lumi_trigger, 
      h_lumi, h_jets_lumi, h_ele_lumi, h_mu_lumi, h_event_lumi, h_topjets_lumi, h_lumi_lumi, 
      h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_lumi_cleaner, 
      h_1lep, h_jets_1lep, h_ele_1lep, h_mu_1lep, h_event_1lep, h_topjets_1lep, h_lumi_1lep, 
      h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets, h_event_2jets, h_topjets_2jets, h_lumi_2jets,
      h_ht350, h_jets_ht350, h_ele_ht350, h_mu_ht350, h_event_ht350, h_topjets_ht350, h_lumi_ht350, 
      h_muveto, h_jets_muveto, h_ele_muveto, h_mu_muveto, h_event_muveto, h_topjets_muveto, h_lumi_muveto,
      h_2lep, h_jets_2lep, h_ele_2lep, h_mu_2lep, h_event_2lep, h_topjets_2lep, h_lumi_2lep,

      h_ht2d, h_eff_cleaner, h_eff_1lep, h_eff_2jets, h_eff_ht350, h_eff_muveto, h_eff_2lep;

    MuonId MuId, MuLoose, MuTight;
    ElectronId EleId, EleLoose, EleTight;

    bool is_mc, is_mu_e, is_e_e;
    TString Sys_EleFakeRate, path_EleFakeRate, Sys_MuFakeRate, path_MuFakeRate;

    uhh2::Event::Handle<vector<Jet>> h_raw_jets_ele, h_raw_jets_mu;   
    uhh2::Event::Handle<vector<Particle>> h_raw_genjets_ele, h_raw_genjets_mu;
  };


  LQToTopMuSidebandPreselectionModule::LQToTopMuSidebandPreselectionModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuSidebandPreselectionModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    EleId = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(30.0, 2.4));               //only tight electrons are used
    MuId = AndId<Muon>(MuonIDLoose(),PtEtaCut(30.0, 2.4),MuonIso(0.15));                   //loose ID for cleaner
    MuLoose = MuonIDLoose();
    MuTight = MuonIDTight();
    EleLoose = ElectronID_Spring16_loose;
    EleTight = ElectronID_Spring16_tight;

    is_mc = ctx.get("dataset_type") == "MC";
    is_mu_e = (ctx.get("channel") == "mu_e" || ctx.get("channel") == "e_mu");
    is_e_e = ctx.get("channel") == "e_e";
    if((!is_mu_e && !is_e_e)) throw runtime_error("In SidebandPreselectionModule: Invalid definition of 'channel' in config file, must be 'mu_e', 'e_mu', or 'e_e'.");

    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->init(ctx);
    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));


    // Muon triggers for mu_e sideband
    mu_trigger_sel1.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    mu_trigger_sel2.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    
    
    //switch to appropriate electron triggers for e_e sideband
    ele_trigger_sel1.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    ele_trigger_sel2.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")); //105
    
    njet_sel.reset(new NJetSelection(2, -1));
    ht_sel.reset(new HtSelection(350)); 
    lumi_sel.reset(new LumiSelection(ctx));

    if(is_mu_e){
      lep1_sel.reset(new NMuonSelection(1,-1, MuTight));      // at least 1 tight muon (standard WP)
      mu_veto.reset(new NMuonSelection(1, 1, MuLoose));       // exactly 1 loose muon: this is the tight one from above + there are no additional muons (only against signal)
      lep2_sel.reset(new NElectronSelection(1, 1, EleTight)); // + ==1 loose electron (standard WP)
    }
    else{
      lep1_sel.reset(new NElectronSelection(1,-1, EleTight)); // at least 1 loose electron  (standard WP)
      mu_veto.reset(new NMuonSelection(0, 0, MuLoose));       // no loose muons
      lep2_sel.reset(new NElectronSelection(2,-1, EleTight)); // at least 2 loose electrons (standard WP)
    }

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuPreselectionHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "Topjets_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_trigger.reset(new LQToTopMuPreselectionHists(ctx, "Trigger"));
    h_jets_trigger.reset(new JetHists(ctx, "Jets_Trigger"));
    h_ele_trigger.reset(new ElectronHists(ctx, "Ele_Trigger"));
    h_mu_trigger.reset(new MuonHists(ctx, "Mu_Trigger"));
    h_event_trigger.reset(new EventHists(ctx, "Event_Trigger"));
    h_topjets_trigger.reset(new TopJetHists(ctx, "Topjets_Trigger"));
    h_lumi_trigger.reset(new LuminosityHists(ctx, "Lumi_Trigger"));

    h_cleaner.reset(new LQToTopMuPreselectionHists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_Cleaner"));
    h_topjets_cleaner.reset(new TopJetHists(ctx, "Topjets_Cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_Cleaner"));
    h_eff_cleaner.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_Cleaner"));

    h_1lep.reset(new LQToTopMuPreselectionHists(ctx, "1Lepton"));
    h_jets_1lep.reset(new JetHists(ctx, "Jets_1Lepton"));
    h_ele_1lep.reset(new ElectronHists(ctx, "Ele_1Lepton"));
    h_mu_1lep.reset(new MuonHists(ctx, "Mu_1Lepton"));
    h_event_1lep.reset(new EventHists(ctx, "Event_1Lepton"));
    h_topjets_1lep.reset(new TopJetHists(ctx, "Topjets_1Lepton"));
    h_lumi_1lep.reset(new LuminosityHists(ctx, "Lumi_1Lepton"));
    h_eff_1lep.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_1Lepton"));

    h_2jets.reset(new LQToTopMuPreselectionHists(ctx, "2Jets"));
    h_jets_2jets.reset(new JetHists(ctx, "Jets_2Jets"));
    h_ele_2jets.reset(new ElectronHists(ctx, "Ele_2Jets"));
    h_mu_2jets.reset(new MuonHists(ctx, "Mu_2Jets"));
    h_event_2jets.reset(new EventHists(ctx, "Event_2Jets"));
    h_topjets_2jets.reset(new TopJetHists(ctx, "Topjets_2Jets"));
    h_lumi_2jets.reset(new LuminosityHists(ctx, "Lumi_2Jets"));
    h_eff_2jets.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_2Jets"));

    h_ht350.reset(new LQToTopMuPreselectionHists(ctx, "HT350"));
    h_jets_ht350.reset(new JetHists(ctx, "Jets_HT350"));
    h_ele_ht350.reset(new ElectronHists(ctx, "Ele_HT350"));
    h_mu_ht350.reset(new MuonHists(ctx, "Mu_HT350"));
    h_event_ht350.reset(new EventHists(ctx, "Event_HT350"));
    h_topjets_ht350.reset(new TopJetHists(ctx, "Topjets_HT350"));
    h_lumi_ht350.reset(new LuminosityHists(ctx, "Lumi_HT350"));
    h_eff_ht350.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_HT350"));
    h_ht2d.reset(new HT2dHists(ctx, "HT2d_HT350"));

    h_2lep.reset(new LQToTopMuPreselectionHists(ctx, "2Leptons"));
    h_jets_2lep.reset(new JetHists(ctx, "Jets_2Leptons"));
    h_ele_2lep.reset(new ElectronHists(ctx, "Ele_2Leptons"));
    h_mu_2lep.reset(new MuonHists(ctx, "Mu_2Leptons"));
    h_event_2lep.reset(new EventHists(ctx, "Event_2Leptons"));
    h_topjets_2lep.reset(new TopJetHists(ctx, "Topjets_2Leptons"));
    h_lumi_2lep.reset(new LuminosityHists(ctx, "Lumi_2Leptons"));
    h_eff_2lep.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_2Leptons"));

    h_muveto.reset(new LQToTopMuPreselectionHists(ctx, "NoAdditionalMuons"));
    h_jets_muveto.reset(new JetHists(ctx, "Jets_NoAdditionalMuons"));
    h_ele_muveto.reset(new ElectronHists(ctx, "Ele_NoAdditionalMuons"));
    h_mu_muveto.reset(new MuonHists(ctx, "Mu_NoAdditionalMuons"));
    h_event_muveto.reset(new EventHists(ctx, "Event_NoAdditionalMuons"));
    h_topjets_muveto.reset(new TopJetHists(ctx, "Topjets_NoAdditionalMuons"));
    h_lumi_muveto.reset(new LuminosityHists(ctx, "Lumi_NoAdditionalMuons"));
    h_eff_muveto.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_NoAdditionalMuons"));


  }


  bool LQToTopMuSidebandPreselectionModule::process(Event & event) {

    
    // 1. run all modules other modules.

    //lumi selection
    if(!is_mc){
      if(!lumi_sel->passes(event)) return false;
    }

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_lumi_nocuts->fill(event);
    
    // trigger
    if(is_mu_e){
      if(!(mu_trigger_sel1->passes(event) || mu_trigger_sel2->passes(event))) return false;
    }
    else{
      if(!(ele_trigger_sel1->passes(event) || ele_trigger_sel2->passes(event))) return false;
    }
    
   
    h_trigger->fill(event);
    h_jets_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_mu_trigger->fill(event);
    h_event_trigger->fill(event);
    h_topjets_trigger->fill(event);
    h_lumi_trigger->fill(event);


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
    h_eff_cleaner->fill(event);

    // >=1 mu (e_mu)/ >=1 ele (e_e)
    if(!lep1_sel->passes(event)) return false;
    h_1lep->fill(event);
    h_jets_1lep->fill(event);
    h_ele_1lep->fill(event);
    h_mu_1lep->fill(event);
    h_event_1lep->fill(event);
    h_topjets_1lep->fill(event);
    h_lumi_1lep->fill(event);
    h_eff_1lep->fill(event);

    // >= 2 jets
    if (!njet_sel->passes(event)) return false;
    h_2jets->fill(event);
    h_jets_2jets->fill(event);
    h_ele_2jets->fill(event);
    h_mu_2jets->fill(event);
    h_event_2jets->fill(event);
    h_topjets_2jets->fill(event);
    h_lumi_2jets->fill(event);
    h_eff_2jets->fill(event);

    //HT > 350
    if (!ht_sel->passes(event)) return false;
    h_ht350->fill(event);
    h_jets_ht350->fill(event);
    h_ele_ht350->fill(event);
    h_mu_ht350->fill(event);
    h_event_ht350->fill(event);
    h_topjets_ht350->fill(event);
    h_lumi_ht350->fill(event);
    h_eff_ht350->fill(event);
    h_ht2d->fill(event);

    // ==1 ele (e_mu) / >2 ele (e_e)
    if(!lep2_sel->passes(event)) return false;
    h_2lep->fill(event);
    h_jets_2lep->fill(event);
    h_ele_2lep->fill(event);
    h_mu_2lep->fill(event);
    h_event_2lep->fill(event);
    h_topjets_2lep->fill(event);
    h_lumi_2lep->fill(event);
    h_eff_2lep->fill(event);

    // veto on additional muons
    if(!mu_veto->passes(event)) return false;
    h_muveto->fill(event);
    h_jets_muveto->fill(event);
    h_ele_muveto->fill(event);
    h_mu_muveto->fill(event);
    h_event_muveto->fill(event);
    h_lumi_muveto->fill(event);
    h_topjets_muveto->fill(event);
    h_eff_muveto->fill(event);



    return true;
  }
  
  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the LQToTopMuSidebandPreselectionModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuSidebandPreselectionModule)
  
}
