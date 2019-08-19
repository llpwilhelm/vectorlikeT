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

  class LQToTopMuPreselectionLumiscaleModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuPreselectionLumiscaleModule(Context & ctx);
    virtual bool process(Event & event) override;
  
  private:
  
    unique_ptr<CommonModules> common;
    //unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer;
  
    unique_ptr<JetCleaner> jetcleaner;
    unique_ptr<JetLeptonCleaner> jetleptoncleaner;
    unique_ptr<MuonCleaner> muoncleaner;
    unique_ptr<MuonCleaner> muoncleaner_iso;
    unique_ptr<ElectronCleaner> electroncleaner;

    unique_ptr<AnalysisModule> syst_module;
    unique_ptr<ElectronFakeRateWeights> SF_eleFakeRate;
    unique_ptr<MuonFakeRateWeights> SF_muFakeRate;
  
    // declare the Selections to use.
    unique_ptr<Selection> njet_sel, nmuon_sel, n_gen_muon_sel, nele_sel, n_gen_ele_sel, ht_sel, lumi_sel, mu1_sel, trigger_sel, trigger_sel1, trigger_sel2, mttbargen_sel;
  
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts, 
      h_trigger, h_jets_trigger, h_ele_trigger, h_mu_trigger, h_event_trigger, h_topjets_trigger, h_lumi_trigger, 
      h_lumi, h_jets_lumi, h_ele_lumi, h_mu_lumi, h_event_lumi, h_topjets_lumi, h_lumi_lumi, 
      h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_lumi_cleaner, 
      h_1mu, h_jets_1mu, h_ele_1mu, h_mu_1mu, h_event_1mu, h_topjets_1mu, h_lumi_1mu, 
      h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets, h_event_2jets, h_topjets_2jets, h_lumi_2jets,
      h_ht350, h_jets_ht350, h_ele_ht350, h_mu_ht350, h_event_ht350, h_topjets_ht350, h_lumi_ht350, 
      h_2mu, h_jets_2mu, h_ele_2mu, h_mu_2mu, h_event_2mu, h_topjets_2mu, h_lumi_2mu,

      h_ht2d, h_eff_2jets, h_eff_cleaner, h_eff_ht350, h_eff_1mu;

    MuonId MuId;
    ElectronId EleId;
    JetId Btag_loose, Jet_ID;

    bool is_mc, is_foreff;
    TString Sys_EleFakeRate, path_EleFakeRate, Sys_MuFakeRate, path_MuFakeRate;

    uhh2::Event::Handle<vector<Jet>> h_raw_jets_ele, h_raw_jets_mu;   
    uhh2::Event::Handle<vector<Particle>> h_raw_genjets_ele, h_raw_genjets_mu;


    uhh2::Event::Handle<double> FakeRateWeightEle;
    uhh2::Event::Handle<double> FakeRateWeightEleUp;
    uhh2::Event::Handle<double> FakeRateWeightEleDown;
    uhh2::Event::Handle<double> FakeRateWeightMu;
    uhh2::Event::Handle<double> FakeRateWeightMuUp;
    uhh2::Event::Handle<double> FakeRateWeightMuDown;
  };


  LQToTopMuPreselectionLumiscaleModule::LQToTopMuPreselectionLumiscaleModule(Context & ctx){
    
    cout << "Hello from LQToTopMuPreselectionLumiscaleModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    is_foreff = ctx.get("IsForEff") == "true";

    // Because this is Lumiscale: all eta-acceptances up to 4.0!

    if(!is_foreff)EleId = AndId<Electron>(ElectronID_Spring16_loose, PtEtaCut(30.0, 4.0));
    else          EleId = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(30.0, 4.0));
    MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 4.0), MuonIso(0.15));
    //MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4));
    Jet_ID = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 4.0));
    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);

    is_mc = ctx.get("dataset_type") == "MC";
    path_EleFakeRate = ctx.get("path_EleFakeRate");
    Sys_EleFakeRate = ctx.get("Systematic_EleFakeRate");
    path_MuFakeRate = ctx.get("path_MuFakeRate");
    Sys_MuFakeRate = ctx.get("Systematic_MuFakeRate");

    common.reset(new CommonModules());
   
    common->switch_jetlepcleaner(true);
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->init(ctx);
    jetcleaner.reset(new JetCleaner(ctx,Jet_ID));
    electroncleaner.reset(new ElectronCleaner(EleId));
    muoncleaner.reset(new MuonCleaner(MuId));
    syst_module.reset(new MCScaleVariation(ctx));


    h_raw_jets_ele = ctx.get_handle<vector<Jet>>("RawJetsEle");
    h_raw_genjets_ele = ctx.get_handle<vector<Particle>>("RawGenJetsEle");
    h_raw_jets_mu = ctx.get_handle<vector<Jet>>("RawJetsMu");
    h_raw_genjets_mu = ctx.get_handle<vector<Particle>>("RawGenJetsMu");


    SF_eleFakeRate.reset(new ElectronFakeRateWeights(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC, path_EleFakeRate, Sys_EleFakeRate, "RawJetsEle", "RawGenJetsEle"));    
    SF_muFakeRate.reset(new MuonFakeRateWeights(ctx, path_MuFakeRate, Sys_MuFakeRate));

    // 2. set up selections

    //Preselection
    trigger_sel1.reset(new TriggerSelection("HLT_IsoMu24_v*")); //original: IsoMu24
    trigger_sel2.reset(new TriggerSelection("HLT_IsoTkMu24_v*")); //original: IsoMu24
    njet_sel.reset(new NJetSelection(2, -1));
    mu1_sel.reset(new NMuonSelection(1, -1));
    nmuon_sel.reset(new NMuonSelection(2, -1)); 
    n_gen_muon_sel.reset(new NGenMuonSelection(2, -1)); 
    nele_sel.reset(new NElectronSelection(2, -1)); 
    n_gen_ele_sel.reset(new NGenElectronSelection(2, -1)); 
    ht_sel.reset(new HtSelection(350));  //350
    lumi_sel.reset(new LumiSelection(ctx));
    mttbargen_sel.reset(new MttbarGenSelection(0,700));

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
    
    h_1mu.reset(new LQToTopMuPreselectionHists(ctx, "1Mu"));
    h_jets_1mu.reset(new JetHists(ctx, "Jets_1Mu"));
    h_ele_1mu.reset(new ElectronHists(ctx, "Ele_1Mu"));
    h_mu_1mu.reset(new MuonHists(ctx, "Mu_1Mu"));
    h_event_1mu.reset(new EventHists(ctx, "Event_1Mu"));
    h_topjets_1mu.reset(new TopJetHists(ctx, "Topjets_1Mu"));
    h_lumi_1mu.reset(new LuminosityHists(ctx, "Lumi_1Mu"));
    h_eff_1mu.reset(new LQToTopMuEfficiencyHists(ctx,"Eff_1Mu"));

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

    h_2mu.reset(new LQToTopMuPreselectionHists(ctx, "2Mu"));
    h_jets_2mu.reset(new JetHists(ctx, "Jets_2Mu"));
    h_ele_2mu.reset(new ElectronHists(ctx, "Ele_2Mu"));
    h_mu_2mu.reset(new MuonHists(ctx, "Mu_2Mu"));
    h_event_2mu.reset(new EventHists(ctx, "Event_2Mu"));
    h_topjets_2mu.reset(new TopJetHists(ctx, "Topjets_2Mu"));
    h_lumi_2mu.reset(new LuminosityHists(ctx, "Lumi_2Mu"));



    FakeRateWeightEle = ctx.declare_event_output<double>("FakeRateWeightEle");
    FakeRateWeightEleUp = ctx.declare_event_output<double>("FakeRateWeightEleUp");
    FakeRateWeightEleDown = ctx.declare_event_output<double>("FakeRateWeightEleDown");
    FakeRateWeightMu = ctx.declare_event_output<double>("FakeRateWeightMu");
    FakeRateWeightMuUp = ctx.declare_event_output<double>("FakeRateWeightMuUp");
    FakeRateWeightMuDown = ctx.declare_event_output<double>("FakeRateWeightMuDown");

  }


  bool LQToTopMuPreselectionLumiscaleModule::process(Event & event) {


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
    if(!(trigger_sel1->passes(event) || trigger_sel2->passes(event) || is_foreff)) return false;
    

    h_trigger->fill(event);
    h_jets_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_mu_trigger->fill(event);
    h_event_trigger->fill(event);
    h_topjets_trigger->fill(event);
    h_lumi_trigger->fill(event);
    
    // apply fake rate weights for electrons and muons
    if(is_mc){

      //electrons
      vector<Jet> raw_jets;
      for(unsigned int i=0; i<event.jets->size(); i++){
	Jet jet = event.jets->at(i);
	raw_jets.push_back(jet);
      }
      event.set(h_raw_jets_ele,move(raw_jets));

      vector<Particle> gen_jets;
      for(unsigned int i=0; i<event.genjets->size(); i++){
	Particle jet = event.genjets->at(i);
	gen_jets.push_back(jet);
      }
      event.set(h_raw_genjets_ele,move(gen_jets));

      electroncleaner->process(event);
      SF_eleFakeRate->process(event);

      //muons
      muoncleaner->process(event);
      SF_muFakeRate->process(event);
    }
    else{
      event.set(FakeRateWeightEle,1.);
      event.set(FakeRateWeightEleUp,1.);
      event.set(FakeRateWeightEleDown,1.);
      event.set(FakeRateWeightMu,1.);
      event.set(FakeRateWeightMuUp,1.);
      event.set(FakeRateWeightMuDown,1.);
    }


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

    
    if(!(mu1_sel->passes(event) || is_foreff)) return false;
    h_1mu->fill(event);
    h_jets_1mu->fill(event);
    h_ele_1mu->fill(event);
    h_mu_1mu->fill(event);
    h_event_1mu->fill(event);
    h_topjets_1mu->fill(event);
    h_lumi_1mu->fill(event);
    h_eff_1mu->fill(event);
  
    if(!njet_sel->passes(event)) return false;
    h_2jets->fill(event);
    h_jets_2jets->fill(event);
    h_ele_2jets->fill(event);
    h_mu_2jets->fill(event);
    h_event_2jets->fill(event);
    h_topjets_2jets->fill(event);
    h_lumi_2jets->fill(event);
    h_eff_2jets->fill(event);


    if(!ht_sel->passes(event)) return false;
    h_ht350->fill(event);
    h_jets_ht350->fill(event);
    h_ele_ht350->fill(event);
    h_mu_ht350->fill(event);
    h_event_ht350->fill(event);
    h_topjets_ht350->fill(event);
    h_lumi_ht350->fill(event);
    h_eff_ht350->fill(event);
    h_ht2d->fill(event);
    

    if(!is_foreff){
      if(!nmuon_sel->passes(event)) return false;
    }
    else{
      if(!(nmuon_sel->passes(event) || nele_sel->passes(event) || n_gen_muon_sel->passes(event) || n_gen_ele_sel->passes(event))) return false;
    }
    h_2mu->fill(event);
    h_jets_2mu->fill(event);
    h_ele_2mu->fill(event);
    h_mu_2mu->fill(event);
    h_event_2mu->fill(event);
    h_topjets_2mu->fill(event);
    h_lumi_2mu->fill(event);
    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuPreselectionLumiscaleModule)
  
}
