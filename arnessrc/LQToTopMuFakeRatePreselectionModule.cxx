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

  class LQToTopMuFakeRatePreselectionModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuFakeRatePreselectionModule(Context & ctx);
    virtual bool process(Event & event) override;
  
  private:
  
    unique_ptr<CommonModules> common;  
    unique_ptr<JetCleaner> jetcleaner;
    unique_ptr<MuonCleaner> muoncleaner;
  
    // declare the Selections to use.
    unique_ptr<Selection> njet_sel, nmuon_sel, nele_sel, lumi_sel, m_ee_sel;
  
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts, 
      h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_lumi_cleaner, 
      h_nmu, h_jets_nmu, h_ele_nmu, h_mu_nmu, h_event_nmu, h_topjets_nmu, h_lumi_nmu, 
      h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets, h_event_2jets, h_topjets_2jets, h_lumi_2jets,
      h_2ele, h_jets_2ele, h_ele_2ele, h_mu_2ele, h_event_2ele, h_topjets_2ele, h_lumi_2ele,
      h_mee, h_jets_mee, h_ele_mee, h_mu_mee, h_event_mee, h_topjets_mee, h_lumi_mee;

    MuonId MuId;
    ElectronId EleId;

    bool is_mc, is_ele_channel;


  };


  LQToTopMuFakeRatePreselectionModule::LQToTopMuFakeRatePreselectionModule(Context & ctx){
    
    cout << "Hello from LQToTopMuFakeRatePreselectionModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    EleId = AndId<Electron>(ElectronID_Spring16_loose, PtEtaCut(30.0, 2.4));
    MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4),MuonIso(0.15));
    is_mc = ctx.get("dataset_type") == "MC";
    is_ele_channel = ctx.get("EleChannel") == "true";

    common.reset(new CommonModules());
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->switch_metcorrection(true);
    common->init(ctx);
    muoncleaner.reset(new MuonCleaner(MuId)); //to cut on ==0 muons before using met-corrections in commonmodules
    jetcleaner.reset(new JetCleaner(ctx, 10.0, 2.5));

    // 2. set up selections

    //Preselection
    if(is_ele_channel) nmuon_sel.reset(new NMuonSelection(0, 0)); 
    nele_sel.reset(new NElectronSelection(2, -1)); 
    lumi_sel.reset(new LumiSelection(ctx));
    m_ee_sel.reset(new InvMassEleEleSelection(71.,111.));

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuPreselectionHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "Topjets_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_nmu.reset(new LQToTopMuPreselectionHists(ctx, "NMu"));
    h_jets_nmu.reset(new JetHists(ctx, "Jets_NMu"));
    h_ele_nmu.reset(new ElectronHists(ctx, "Ele_NMu"));
    h_mu_nmu.reset(new MuonHists(ctx, "Mu_NMu"));
    h_event_nmu.reset(new EventHists(ctx, "Event_NMu"));
    h_topjets_nmu.reset(new TopJetHists(ctx, "Topjets_NMu"));
    h_lumi_nmu.reset(new LuminosityHists(ctx, "Lumi_NMu"));

    h_cleaner.reset(new LQToTopMuPreselectionHists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_Cleaner"));
    h_topjets_cleaner.reset(new TopJetHists(ctx, "Topjets_Cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_Cleaner"));

    h_2ele.reset(new LQToTopMuPreselectionHists(ctx, "2Ele"));
    h_jets_2ele.reset(new JetHists(ctx, "Jets_2Ele"));
    h_ele_2ele.reset(new ElectronHists(ctx, "Ele_2Ele"));
    h_mu_2ele.reset(new MuonHists(ctx, "Mu_2Ele"));
    h_event_2ele.reset(new EventHists(ctx, "Event_2Ele"));
    h_topjets_2ele.reset(new TopJetHists(ctx, "Topjets_2Ele"));
    h_lumi_2ele.reset(new LuminosityHists(ctx, "Lumi_2Ele"));

    h_mee.reset(new LQToTopMuPreselectionHists(ctx, "Mee"));
    h_jets_mee.reset(new JetHists(ctx, "Jets_Mee"));
    h_ele_mee.reset(new ElectronHists(ctx, "Ele_Mee"));
    h_mu_mee.reset(new MuonHists(ctx, "Mu_Mee"));
    h_event_mee.reset(new EventHists(ctx, "Event_Mee"));
    h_topjets_mee.reset(new TopJetHists(ctx, "Topjets_Mee"));
    h_lumi_mee.reset(new LuminosityHists(ctx, "Lumi_Mee"));


  }


  bool LQToTopMuFakeRatePreselectionModule::process(Event & event) {


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

    muoncleaner->process(event);
    if(is_ele_channel){
      if(!nmuon_sel->passes(event)) return false;
    }
    h_nmu->fill(event);
    h_jets_nmu->fill(event);
    h_ele_nmu->fill(event);
    h_mu_nmu->fill(event);
    h_event_nmu->fill(event);
    h_topjets_nmu->fill(event);
    h_lumi_nmu->fill(event);

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

    if(!nele_sel->passes(event)) return false;
    h_2ele->fill(event);
    h_jets_2ele->fill(event);
    h_ele_2ele->fill(event);
    h_mu_2ele->fill(event);
    h_event_2ele->fill(event);
    h_topjets_2ele->fill(event);
    h_lumi_2ele->fill(event);

    if(!m_ee_sel->passes(event)) return false;
    h_mee->fill(event);
    h_jets_mee->fill(event);
    h_ele_mee->fill(event);
    h_mu_mee->fill(event);
    h_event_mee->fill(event);
    h_topjets_mee->fill(event);
    h_lumi_mee->fill(event);


    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuFakeRatePreselectionModule)
  
}
