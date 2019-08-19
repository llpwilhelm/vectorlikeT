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

  class SingleVLQStephaniePreselectionModule: public AnalysisModule {
  public:

    explicit SingleVLQStephaniePreselectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    unique_ptr<CommonModules> common;

    // declare the Selections to use.
    unique_ptr<Selection> jet_sel, muon_sel, ele_sel;

    // store the Hists collection as member variables.
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_lumi_cleaner;
    unique_ptr<Hists> h_muon, h_jets_muon, h_ele_muon, h_mu_muon, h_event_muon, h_topjets_muon, h_lumi_muon;
    unique_ptr<Hists> h_electron, h_jets_electron, h_ele_electron, h_mu_electron, h_event_electron, h_topjets_electron, h_lumi_electron;
    unique_ptr<Hists> h_jets, h_jets_jets, h_ele_jets, h_mu_jets, h_event_jets, h_topjets_jets, h_lumi_jets;
    unique_ptr<Hists> h_met, h_jets_met, h_ele_met, h_mu_met, h_event_met, h_topjets_met, h_lumi_met;

    MuonId MuId;
    ElectronId EleId;
    JetId Jet_ID;

    bool is_mc;

  };


  SingleVLQStephaniePreselectionModule::SingleVLQStephaniePreselectionModule(Context & ctx){

    cout << "Hello from SingleVLQStephaniePreselectionModule!" << endl;

    for(auto & kv : ctx.get_all()) cout << " " << kv.first << " = " << kv.second << endl;


    EleId = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(30.0, 2.4));
    MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4), MuonIso(0.15));
    Jet_ID = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));

    is_mc = ctx.get("dataset_type") == "MC";

    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->set_jet_id(Jet_ID);
    common->switch_jetPtSorter();
    common->init(ctx);


    // 2. set up selections

    //Preselection
    muon_sel.reset(new NMuonSelection(1, 1));
    ele_sel.reset(new NElectronSelection(0, 0));
    jet_sel.reset(new NJetSelection(3, -1));

    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuPreselectionHists(ctx, "nocuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_nocuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_nocuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_nocuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_nocuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "Topjets_nocuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_nocuts"));

    h_cleaner.reset(new LQToTopMuPreselectionHists(ctx, "cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_cleaner"));
    h_topjets_cleaner.reset(new TopJetHists(ctx, "Topjets_cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_cleaner"));

    h_muon.reset(new LQToTopMuPreselectionHists(ctx, "muon"));
    h_jets_muon.reset(new JetHists(ctx, "Jets_muon"));
    h_ele_muon.reset(new ElectronHists(ctx, "Ele_muon"));
    h_mu_muon.reset(new MuonHists(ctx, "Mu_muon"));
    h_event_muon.reset(new EventHists(ctx, "Event_muon"));
    h_topjets_muon.reset(new TopJetHists(ctx, "Topjets_muon"));
    h_lumi_muon.reset(new LuminosityHists(ctx, "Lumi_muon"));

    h_electron.reset(new LQToTopMuPreselectionHists(ctx, "electron"));
    h_jets_electron.reset(new JetHists(ctx, "Jets_electron"));
    h_ele_electron.reset(new ElectronHists(ctx, "Ele_electron"));
    h_mu_electron.reset(new MuonHists(ctx, "Mu_electron"));
    h_event_electron.reset(new EventHists(ctx, "Event_electron"));
    h_topjets_electron.reset(new TopJetHists(ctx, "Topjets_electron"));
    h_lumi_electron.reset(new LuminosityHists(ctx, "Lumi_electron"));

    h_jets.reset(new LQToTopMuPreselectionHists(ctx, "jets"));
    h_jets_jets.reset(new JetHists(ctx, "Jets_jets"));
    h_ele_jets.reset(new ElectronHists(ctx, "Ele_jets"));
    h_mu_jets.reset(new MuonHists(ctx, "Mu_jets"));
    h_event_jets.reset(new EventHists(ctx, "Event_jets"));
    h_topjets_jets.reset(new TopJetHists(ctx, "Topjets_jets"));
    h_lumi_jets.reset(new LuminosityHists(ctx, "Lumi_jets"));

    h_met.reset(new LQToTopMuPreselectionHists(ctx, "met"));
    h_jets_met.reset(new JetHists(ctx, "Jets_met"));
    h_ele_met.reset(new ElectronHists(ctx, "Ele_met"));
    h_mu_met.reset(new MuonHists(ctx, "Mu_met"));
    h_event_met.reset(new EventHists(ctx, "Event_met"));
    h_topjets_met.reset(new TopJetHists(ctx, "Topjets_met"));
    h_lumi_met.reset(new LuminosityHists(ctx, "Lumi_met"));



  }


  bool SingleVLQStephaniePreselectionModule::process(Event & event) {


    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_lumi_nocuts->fill(event);


    bool pass_common = common->process(event);
    if(!pass_common) return false;

    h_cleaner->fill(event);
    h_jets_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_event_cleaner->fill(event);
    h_lumi_cleaner->fill(event);


    if(!(muon_sel->passes(event))) return false;
    h_muon->fill(event);
    h_jets_muon->fill(event);
    h_ele_muon->fill(event);
    h_mu_muon->fill(event);
    h_event_muon->fill(event);
    h_lumi_muon->fill(event);


    if(!(ele_sel->passes(event))) return false;
    h_electron->fill(event);
    h_jets_electron->fill(event);
    h_ele_electron->fill(event);
    h_mu_electron->fill(event);
    h_event_electron->fill(event);
    h_lumi_electron->fill(event);

    if(!jet_sel->passes(event)) return false;
    h_jets->fill(event);
    h_jets_jets->fill(event);
    h_ele_jets->fill(event);
    h_mu_jets->fill(event);
    h_event_jets->fill(event);
    h_lumi_jets->fill(event);


    if(event.met->pt() < 30) return false;
    h_met->fill(event);
    h_jets_met->fill(event);
    h_ele_met->fill(event);
    h_mu_met->fill(event);
    h_event_met->fill(event);
    h_lumi_met->fill(event);
    return true;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(SingleVLQStephaniePreselectionModule)

}
