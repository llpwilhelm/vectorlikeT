#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/SingleVLQHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuMLQRecoHists.h"
#include "UHH2/LQToTopMu/include/MET2dHists.h"
#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/LQToTopMu/include/SingleVLQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/SingleVLQReconstruction.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "LHAPDF/LHAPDF.h"
#include "UHH2/LQToTopMu/include/LQToTopMu14TeVCrossCheckHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class SingleVLQStephanieAnalysisModule: public AnalysisModule {
  public:

    explicit SingleVLQStephanieAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    unique_ptr<CommonModules> common;
    unique_ptr<AnalysisModule> SF_muonID, SF_muonIso, SF_btag, scale_variation_module;
    unique_ptr<MCMuonScaleFactor> SF_muonTrigger;

    // declare the Selections to use.
    unique_ptr<Selection>  btag_loose_sel, btag_2medium_sel, btag_3medium_sel, btag_2tight_sel, btag_3tight_sel, trigger1_sel, trigger2_sel;

    unique_ptr<HighMassSingleVLQReconstruction> tprime_reco;
    unique_ptr<SingleVLQChi2Discriminator> tprime_chi2;

    // store the Hists collection as member variables.
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_lumi_cleaner;
    unique_ptr<Hists> h_trigger, h_jets_trigger, h_ele_trigger, h_mu_trigger, h_event_trigger, h_lumi_trigger, h_btageff_trigger;
    unique_ptr<Hists> h_btag, h_jets_btag, h_ele_btag, h_mu_btag, h_event_btag, h_lumi_btag;
    unique_ptr<Hists> h_btag_2m, h_jets_btag_2m, h_ele_btag_2m, h_mu_btag_2m, h_event_btag_2m, h_lumi_btag_2m;
    unique_ptr<Hists> h_btag_3m, h_jets_btag_3m, h_ele_btag_3m, h_mu_btag_3m, h_event_btag_3m, h_lumi_btag_3m;
    unique_ptr<Hists> h_btag_2t, h_jets_btag_2t, h_ele_btag_2t, h_mu_btag_2t, h_event_btag_2t, h_lumi_btag_2t;
    unique_ptr<Hists> h_btag_3t, h_jets_btag_3t, h_ele_btag_3t, h_mu_btag_3t, h_event_btag_3t, h_lumi_btag_3t;
    unique_ptr<Hists> h_reco, h_jets_reco, h_ele_reco, h_mu_reco, h_event_reco, h_lumi_reco;
    unique_ptr<Hists> h_dRbb, h_jets_dRbb, h_ele_dRbb, h_mu_dRbb, h_event_dRbb, h_lumi_dRbb;
    unique_ptr<Hists> h_dRbw, h_jets_dRbw, h_ele_dRbw, h_mu_dRbw, h_event_dRbw, h_lumi_dRbw;
    unique_ptr<Hists> h_chi2_10, h_jets_chi2_10, h_ele_chi2_10, h_mu_chi2_10, h_event_chi2_10, h_lumi_chi2_10;
    unique_ptr<Hists> h_chi2h_5, h_jets_chi2h_5, h_ele_chi2h_5, h_mu_chi2h_5, h_event_chi2h_5, h_lumi_chi2h_5;
    unique_ptr<Hists> h_dRbb_11, h_jets_dRbb_11, h_ele_dRbb_11, h_mu_dRbb_11, h_event_dRbb_11, h_lumi_dRbb_11;
    unique_ptr<Hists> h_dRbw_15, h_jets_dRbw_15, h_ele_dRbw_15, h_mu_dRbw_15, h_event_dRbw_15, h_lumi_dRbw_15;
    unique_ptr<Hists> h_chi2h_2, h_jets_chi2h_2, h_ele_chi2h_2, h_mu_chi2h_2, h_event_chi2h_2, h_lumi_chi2h_2;
    unique_ptr<Hists> h_dRbb_10, h_jets_dRbb_10, h_ele_dRbb_10, h_mu_dRbb_10, h_event_dRbb_10, h_lumi_dRbb_10;
    unique_ptr<Hists> h_dRth, h_jets_dRth, h_ele_dRth, h_mu_dRth, h_event_dRth, h_lumi_dRth;
    unique_ptr<Hists> h_final, h_jets_final, h_ele_final, h_mu_final, h_event_final, h_lumi_final;

    JetId Btag_loose, Btag_medium, Btag_tight;
    CSVBTag::wp wp_btag_loose, wp_btag_medium, wp_btag_tight;

    bool is_mc, do_scale_variation;
    TString dataset_version, scalevariation_process, region;

    uhh2::Event::Handle<bool> h_is_tprime_reco;
    uhh2::Event::Handle<std::vector<SingleVLQReconstructionHypothesis>> h_hyps;
  };


  SingleVLQStephanieAnalysisModule::SingleVLQStephanieAnalysisModule(Context & ctx){

    cout << "Hello World from SingleVLQStephanieAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()) cout << " " << kv.first << " = " << kv.second << endl;

    is_mc = ctx.get("dataset_type") == "MC";
    dataset_version = ctx.get("dataset_version");
    region = ctx.get("Region");
    if(region != "SR" && region != "HighChi2") throw runtime_error("Invalid value of 'Region' in xml file.");
    scalevariation_process   = ctx.get("ScaleVariationProcess");
    do_scale_variation       = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    if(do_scale_variation && (scalevariation_process != "ttbar") && (scalevariation_process != "dy") && (scalevariation_process != "st") && (scalevariation_process != "wj") && (scalevariation_process != "db") && (scalevariation_process != "ttv")) throw runtime_error("In LQToTopMuAnalysisModule.cxx: Invalid process specified for 'ScaleVariationProcess'.");
    wp_btag_loose = CSVBTag::WP_LOOSE;
    wp_btag_medium = CSVBTag::WP_MEDIUM;
    wp_btag_tight = CSVBTag::WP_TIGHT;
    Btag_loose = CSVBTag(wp_btag_loose);
    Btag_medium = CSVBTag(wp_btag_medium);
    Btag_tight = CSVBTag(wp_btag_tight);

    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->init(ctx);

    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, "nominal"));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, "nominal"));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, "nominal"));
    SF_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_tight,"jets","central"));



    h_is_tprime_reco = ctx.get_handle<bool>("is_tprime_reco");
    h_hyps = ctx.get_handle<vector<SingleVLQReconstructionHypothesis>>("TprimeHypotheses");

    // 2. set up selections
    //Selection
    btag_loose_sel.reset(new NJetSelection(3, -1, Btag_loose));
    btag_2medium_sel.reset(new NJetSelection(2,-1,Btag_medium));
    btag_3medium_sel.reset(new NJetSelection(3,-1,Btag_medium));
    btag_2tight_sel.reset(new NJetSelection(2,-1,Btag_tight));
    btag_3tight_sel.reset(new NJetSelection(3,-1,Btag_tight));
    trigger1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trigger2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));

    tprime_reco.reset(new HighMassSingleVLQReconstruction(ctx, SingleVLQNeutrinoReconstruction, wp_btag_loose));
    tprime_chi2.reset(new SingleVLQChi2Discriminator(ctx));

    scale_variation_module.reset(new MCScaleVariation(ctx));

    // 3. Set up Hists classes:
    h_nocuts.reset(new SingleVLQHists(ctx, "nocuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_nocuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_nocuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_nocuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_nocuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_nocuts"));

    h_cleaner.reset(new SingleVLQHists(ctx, "cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_cleaner"));

    h_trigger.reset(new SingleVLQHists(ctx, "trigger"));
    h_jets_trigger.reset(new JetHists(ctx, "Jets_trigger"));
    h_ele_trigger.reset(new ElectronHists(ctx, "Ele_trigger"));
    h_mu_trigger.reset(new MuonHists(ctx, "Mu_trigger"));
    h_event_trigger.reset(new EventHists(ctx, "Event_trigger"));
    h_lumi_trigger.reset(new LuminosityHists(ctx, "Lumi_trigger"));

    h_btageff_trigger.reset(new BTagMCEfficiencyHists(ctx, "btageff_trigger",wp_btag_tight));

    h_btag.reset(new SingleVLQHists(ctx, "btag"));
    h_jets_btag.reset(new JetHists(ctx, "Jets_btag"));
    h_ele_btag.reset(new ElectronHists(ctx, "Ele_btag"));
    h_mu_btag.reset(new MuonHists(ctx, "Mu_btag"));
    h_event_btag.reset(new EventHists(ctx, "Event_btag"));
    h_lumi_btag.reset(new LuminosityHists(ctx, "Lumi_btag"));

    h_btag_2m.reset(new SingleVLQHists(ctx, "btag_2m"));
    h_jets_btag_2m.reset(new JetHists(ctx, "Jets_btag_2m"));
    h_ele_btag_2m.reset(new ElectronHists(ctx, "Ele_btag_2m"));
    h_mu_btag_2m.reset(new MuonHists(ctx, "Mu_btag_2m"));
    h_event_btag_2m.reset(new EventHists(ctx, "Event_btag_2m"));
    h_lumi_btag_2m.reset(new LuminosityHists(ctx, "Lumi_btag_2m"));

    h_btag_3m.reset(new SingleVLQHists(ctx, "btag_3m"));
    h_jets_btag_3m.reset(new JetHists(ctx, "Jets_btag_3m"));
    h_ele_btag_3m.reset(new ElectronHists(ctx, "Ele_btag_3m"));
    h_mu_btag_3m.reset(new MuonHists(ctx, "Mu_btag_3m"));
    h_event_btag_3m.reset(new EventHists(ctx, "Event_btag_3m"));
    h_lumi_btag_3m.reset(new LuminosityHists(ctx, "Lumi_btag_3m"));

    h_btag_2t.reset(new SingleVLQHists(ctx, "btag_2t"));
    h_jets_btag_2t.reset(new JetHists(ctx, "Jets_btag_2t"));
    h_ele_btag_2t.reset(new ElectronHists(ctx, "Ele_btag_2t"));
    h_mu_btag_2t.reset(new MuonHists(ctx, "Mu_btag_2t"));
    h_event_btag_2t.reset(new EventHists(ctx, "Event_btag_2t"));
    h_lumi_btag_2t.reset(new LuminosityHists(ctx, "Lumi_btag_2t"));

    h_btag_3t.reset(new SingleVLQHists(ctx, "btag_3t"));
    h_jets_btag_3t.reset(new JetHists(ctx, "Jets_btag_3t"));
    h_ele_btag_3t.reset(new ElectronHists(ctx, "Ele_btag_3t"));
    h_mu_btag_3t.reset(new MuonHists(ctx, "Mu_btag_3t"));
    h_event_btag_3t.reset(new EventHists(ctx, "Event_btag_3t"));
    h_lumi_btag_3t.reset(new LuminosityHists(ctx, "Lumi_btag_3t"));

    h_reco.reset(new SingleVLQHists(ctx, "reco"));
    h_jets_reco.reset(new JetHists(ctx, "Jets_reco"));
    h_ele_reco.reset(new ElectronHists(ctx, "Ele_reco"));
    h_mu_reco.reset(new MuonHists(ctx, "Mu_reco"));
    h_event_reco.reset(new EventHists(ctx, "Event_reco"));
    h_lumi_reco.reset(new LuminosityHists(ctx, "Lumi_reco"));

    h_dRbb.reset(new SingleVLQHists(ctx, "dRbb"));
    h_jets_dRbb.reset(new JetHists(ctx, "Jets_dRbb"));
    h_ele_dRbb.reset(new ElectronHists(ctx, "Ele_dRbb"));
    h_mu_dRbb.reset(new MuonHists(ctx, "Mu_dRbb"));
    h_event_dRbb.reset(new EventHists(ctx, "Event_dRbb"));
    h_lumi_dRbb.reset(new LuminosityHists(ctx, "Lumi_dRbb"));

    h_dRbw.reset(new SingleVLQHists(ctx, "dRbw"));
    h_jets_dRbw.reset(new JetHists(ctx, "Jets_dRbw"));
    h_ele_dRbw.reset(new ElectronHists(ctx, "Ele_dRbw"));
    h_mu_dRbw.reset(new MuonHists(ctx, "Mu_dRbw"));
    h_event_dRbw.reset(new EventHists(ctx, "Event_dRbw"));
    h_lumi_dRbw.reset(new LuminosityHists(ctx, "Lumi_dRbw"));

    h_chi2_10.reset(new SingleVLQHists(ctx, "chi2_10"));
    h_jets_chi2_10.reset(new JetHists(ctx, "Jets_chi2_10"));
    h_ele_chi2_10.reset(new ElectronHists(ctx, "Ele_chi2_10"));
    h_mu_chi2_10.reset(new MuonHists(ctx, "Mu_chi2_10"));
    h_event_chi2_10.reset(new EventHists(ctx, "Event_chi2_10"));
    h_lumi_chi2_10.reset(new LuminosityHists(ctx, "Lumi_chi2_10"));

    h_chi2h_5.reset(new SingleVLQHists(ctx, "chi2h_5"));
    h_jets_chi2h_5.reset(new JetHists(ctx, "Jets_chi2h_5"));
    h_ele_chi2h_5.reset(new ElectronHists(ctx, "Ele_chi2h_5"));
    h_mu_chi2h_5.reset(new MuonHists(ctx, "Mu_chi2h_5"));
    h_event_chi2h_5.reset(new EventHists(ctx, "Event_chi2h_5"));
    h_lumi_chi2h_5.reset(new LuminosityHists(ctx, "Lumi_chi2h_5"));

    h_dRbb_11.reset(new SingleVLQHists(ctx, "dRbb_11"));
    h_jets_dRbb_11.reset(new JetHists(ctx, "Jets_dRbb_11"));
    h_ele_dRbb_11.reset(new ElectronHists(ctx, "Ele_dRbb_11"));
    h_mu_dRbb_11.reset(new MuonHists(ctx, "Mu_dRbb_11"));
    h_event_dRbb_11.reset(new EventHists(ctx, "Event_dRbb_11"));
    h_lumi_dRbb_11.reset(new LuminosityHists(ctx, "Lumi_dRbb_11"));

    h_dRbw_15.reset(new SingleVLQHists(ctx, "dRbw_15"));
    h_jets_dRbw_15.reset(new JetHists(ctx, "Jets_dRbw_15"));
    h_ele_dRbw_15.reset(new ElectronHists(ctx, "Ele_dRbw_15"));
    h_mu_dRbw_15.reset(new MuonHists(ctx, "Mu_dRbw_15"));
    h_event_dRbw_15.reset(new EventHists(ctx, "Event_dRbw_15"));
    h_lumi_dRbw_15.reset(new LuminosityHists(ctx, "Lumi_dRbw_15"));

    h_chi2h_2.reset(new SingleVLQHists(ctx, "chi2h_2"));
    h_jets_chi2h_2.reset(new JetHists(ctx, "Jets_chi2h_2"));
    h_ele_chi2h_2.reset(new ElectronHists(ctx, "Ele_chi2h_2"));
    h_mu_chi2h_2.reset(new MuonHists(ctx, "Mu_chi2h_2"));
    h_event_chi2h_2.reset(new EventHists(ctx, "Event_chi2h_2"));
    h_lumi_chi2h_2.reset(new LuminosityHists(ctx, "Lumi_chi2h_2"));

    h_dRbb_10.reset(new SingleVLQHists(ctx, "dRbb_10"));
    h_jets_dRbb_10.reset(new JetHists(ctx, "Jets_dRbb_10"));
    h_ele_dRbb_10.reset(new ElectronHists(ctx, "Ele_dRbb_10"));
    h_mu_dRbb_10.reset(new MuonHists(ctx, "Mu_dRbb_10"));
    h_event_dRbb_10.reset(new EventHists(ctx, "Event_dRbb_10"));
    h_lumi_dRbb_10.reset(new LuminosityHists(ctx, "Lumi_dRbb_10"));

    h_dRth.reset(new SingleVLQHists(ctx, "dRth"));
    h_jets_dRth.reset(new JetHists(ctx, "Jets_dRth"));
    h_ele_dRth.reset(new ElectronHists(ctx, "Ele_dRth"));
    h_mu_dRth.reset(new MuonHists(ctx, "Mu_dRth"));
    h_event_dRth.reset(new EventHists(ctx, "Event_dRth"));
    h_lumi_dRth.reset(new LuminosityHists(ctx, "Lumi_dRth"));

    h_final.reset(new SingleVLQHists(ctx, "final"));
    h_jets_final.reset(new JetHists(ctx, "Jets_final"));
    h_ele_final.reset(new ElectronHists(ctx, "Ele_final"));
    h_mu_final.reset(new MuonHists(ctx, "Mu_final"));
    h_event_final.reset(new EventHists(ctx, "Event_final"));
    h_lumi_final.reset(new LuminosityHists(ctx, "Lumi_final"));

  }


  bool SingleVLQStephanieAnalysisModule::process(Event & event) {

    event.set(h_is_tprime_reco, false);

    if(do_scale_variation){
      if((dataset_version.Contains("TTbar") && scalevariation_process == "ttbar") || (dataset_version.Contains("DYJets") && scalevariation_process == "dy") || (dataset_version.Contains("SingleTop") && scalevariation_process == "st") || (dataset_version.Contains("WJets") && scalevariation_process == "wj") || (dataset_version.Contains("Diboson") && scalevariation_process == "db") || (dataset_version.Contains("TTV") && scalevariation_process == "ttv")){
        if(event.genInfo->systweights().size() < 10 && dataset_version.Contains("Diboson")) cout << "SystWeight size: " << event.genInfo->systweights().size() << endl;
        else{
          scale_variation_module->process(event);
        }
      }
    }

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_lumi_nocuts->fill(event);


    bool pass_common = common->process(event);
    if(!pass_common) return false;
    SF_muonID->process(event);
    SF_muonIso->process(event);

    h_cleaner->fill(event);
    h_jets_cleaner->fill(event);
    h_ele_cleaner->fill(event);
    h_mu_cleaner->fill(event);
    h_event_cleaner->fill(event);
    h_lumi_cleaner->fill(event);

    // Trigger
    if(!(trigger1_sel->passes(event) || trigger2_sel->passes(event)) ) return false;
    SF_muonTrigger->process_onemuon(event,0);
    SF_btag->process(event);

    h_trigger->fill(event);
    h_jets_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_mu_trigger->fill(event);
    h_event_trigger->fill(event);
    h_lumi_trigger->fill(event);

    h_btageff_trigger->fill(event);

    // 3 b-tags
    if(!btag_loose_sel->passes(event)) return false;
    h_btag->fill(event);
    h_jets_btag->fill(event);
    h_ele_btag->fill(event);
    h_mu_btag->fill(event);
    h_event_btag->fill(event);
    h_lumi_btag->fill(event);

    if(!(btag_2medium_sel->passes(event))) return false;
    h_btag_2m->fill(event);
    h_jets_btag_2m->fill(event);
    h_ele_btag_2m->fill(event);
    h_mu_btag_2m->fill(event);
    h_event_btag_2m->fill(event);
    h_lumi_btag_2m->fill(event);

    if(!(btag_3medium_sel->passes(event))) return false;
    h_btag_3m->fill(event);
    h_jets_btag_3m->fill(event);
    h_ele_btag_3m->fill(event);
    h_mu_btag_3m->fill(event);
    h_event_btag_3m->fill(event);
    h_lumi_btag_3m->fill(event);

    if(!(btag_2tight_sel->passes(event))) return false;
    h_btag_2t->fill(event);
    h_jets_btag_2t->fill(event);
    h_ele_btag_2t->fill(event);
    h_mu_btag_2t->fill(event);
    h_event_btag_2t->fill(event);
    h_lumi_btag_2t->fill(event);


    if(!(btag_3tight_sel->passes(event))) return false;
    h_btag_3t->fill(event);
    h_jets_btag_3t->fill(event);
    h_ele_btag_3t->fill(event);
    h_mu_btag_3t->fill(event);
    h_event_btag_3t->fill(event);
    h_lumi_btag_3t->fill(event);

    tprime_reco->process(event);
    tprime_chi2->process(event);
    bool is_tprime_reco = event.get(h_is_tprime_reco);
    if(!is_tprime_reco) throw runtime_error("After reconstruction, the T still isn't reconstructed. How?");

    h_reco->fill(event);
    h_jets_reco->fill(event);
    h_ele_reco->fill(event);
    h_mu_reco->fill(event);
    h_event_reco->fill(event);
    h_lumi_reco->fill(event);

    std::vector<SingleVLQReconstructionHypothesis> hyps = event.get(h_hyps);
    const SingleVLQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );

    float chi2 = hyp->discriminator("Chi2");
    float chi2h = hyp->discriminator("Chi2_higgs");
    float dR_bH_bH = deltaR(hyp->higgs_jets().at(0),hyp->higgs_jets().at(1));
    LorentzVector W = (hyp->lepton_v4() + hyp->neutrino_v4());
    float dR_bt_Wt = deltaR(W, hyp->toplep_jets().at(0));
    float dR_t_H = deltaR(hyp->toplep_v4(), hyp->higgs_v4());

    if(!(dR_bH_bH < 1.5)) return false;
    h_dRbb->fill(event);
    h_jets_dRbb->fill(event);
    h_ele_dRbb->fill(event);
    h_mu_dRbb->fill(event);
    h_event_dRbb->fill(event);
    h_lumi_dRbb->fill(event);

    if(!(dR_bt_Wt < 2.)) return false;
    h_dRbw->fill(event);
    h_jets_dRbw->fill(event);
    h_ele_dRbw->fill(event);
    h_mu_dRbw->fill(event);
    h_event_dRbw->fill(event);
    h_lumi_dRbw->fill(event);

    if(region == "SR"){
      if(!(chi2 < 10.)) return false;
    }
    else if(region == "HighChi2"){
      if(!(chi2 > 10.)) return false;
    }
    h_chi2_10->fill(event);
    h_jets_chi2_10->fill(event);
    h_ele_chi2_10->fill(event);
    h_mu_chi2_10->fill(event);
    h_event_chi2_10->fill(event);
    h_lumi_chi2_10->fill(event);

    if(region == "SR"){
      if(!(chi2h < 5.)) return false;
    }
    else if(region == "HighChi2"){
      //do nothing, effective cut comes later
    }
    h_chi2h_5->fill(event);
    h_jets_chi2h_5->fill(event);
    h_ele_chi2h_5->fill(event);
    h_mu_chi2h_5->fill(event);
    h_event_chi2h_5->fill(event);
    h_lumi_chi2h_5->fill(event);

    if(!(dR_bH_bH < 1.1)) return false;
    h_dRbb_11->fill(event);
    h_jets_dRbb_11->fill(event);
    h_ele_dRbb_11->fill(event);
    h_mu_dRbb_11->fill(event);
    h_event_dRbb_11->fill(event);
    h_lumi_dRbb_11->fill(event);

    if(!(dR_bt_Wt < 1.5)) return false;
    h_dRbw_15->fill(event);
    h_jets_dRbw_15->fill(event);
    h_ele_dRbw_15->fill(event);
    h_mu_dRbw_15->fill(event);
    h_event_dRbw_15->fill(event);
    h_lumi_dRbw_15->fill(event);

    if(region == "SR"){
      if(!(chi2h < 2.)) return false;
    }
    else if(region == "HighChi2"){
      if(!(chi2h > 2. || chi2 > 5.)) return false;      
    }
    h_chi2h_2->fill(event);
    h_jets_chi2h_2->fill(event);
    h_ele_chi2h_2->fill(event);
    h_mu_chi2h_2->fill(event);
    h_event_chi2h_2->fill(event);
    h_lumi_chi2h_2->fill(event);

    if(!(dR_bH_bH < 1.0)) return false;
    h_dRbb_10->fill(event);
    h_jets_dRbb_10->fill(event);
    h_ele_dRbb_10->fill(event);
    h_mu_dRbb_10->fill(event);
    h_event_dRbb_10->fill(event);
    h_lumi_dRbb_10->fill(event);

    if(!(dR_t_H > 1.5)) return false;
    h_dRth->fill(event);
    h_jets_dRth->fill(event);
    h_ele_dRth->fill(event);
    h_mu_dRth->fill(event);
    h_event_dRth->fill(event);
    h_lumi_dRth->fill(event);












    h_final->fill(event);
    h_jets_final->fill(event);
    h_ele_final->fill(event);
    h_mu_final->fill(event);
    h_event_final->fill(event);
    h_lumi_final->fill(event);

    return true;
  }


  UHH2_REGISTER_ANALYSIS_MODULE(SingleVLQStephanieAnalysisModule)
}
