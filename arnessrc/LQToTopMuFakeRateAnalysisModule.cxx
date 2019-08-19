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
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuFakeRateHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuRecoHists.h"
#include "UHH2/LQToTopMu/include/MET2dHists.h"
#include "UHH2/LQToTopMu/include/HT2dHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstruction.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/AdditionalSelections.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"


using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class LQToTopMuFakeRateAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuFakeRateAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    unique_ptr<CommonModules>  common;
    unique_ptr<AnalysisModule> SF_btag, SF_eleReco, SF_eleID, SF_muonID, SF_muonIso;
    unique_ptr<MuonTrkWeights> SF_muonTrk;
    unique_ptr<ElectronTriggerWeights> SF_eleTrigger;
    unique_ptr<DibosonScaleFactors> SF_Diboson;
    unique_ptr<ZEEFinder> ZEE_finder;
    unique_ptr<ElectronJetOverlapCleaner> elejetcleaner;
    
    unique_ptr<JetCleaner> jetcleaner;
    
    // declare the Selections to use.
    unique_ptr<Selection>  ele_trigger_sel1, ele_trigger_sel2, njet_sel, nbtag_loose_sel, m_ee_sel, met_sel, nele_sel, ht_sel, ele2_sel, Mee_sel, nmuon_sel;
    
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts, h_eff_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner, h_topjets_cleaner, h_tau_cleaner, h_eff_cleaner, h_lumi_cleaner;
    unique_ptr<Hists> h_2ele, h_jets_2ele, h_ele_2ele, h_mu_2ele, h_event_2ele, h_topjets_2ele, h_tau_2ele, h_eff_2ele, h_lumi_2ele;
    unique_ptr<Hists> h_mee, h_jets_mee, h_ele_mee, h_mu_mee, h_event_mee, h_topjets_mee, h_tau_mee, h_eff_mee, h_lumi_mee;
    unique_ptr<Hists> h_2jets, h_jets_2jets, h_ele_2jets, h_mu_2jets, h_event_2jets, h_topjets_2jets, h_tau_2jets, h_eff_2jets, h_lumi_2jets;
    unique_ptr<Hists> h_ht, h_jets_ht, h_ele_ht, h_mu_ht, h_event_ht, h_topjets_ht, h_tau_ht, h_eff_ht, h_lumi_ht;
    unique_ptr<Hists> h_btag, h_jets_btag, h_ele_btag, h_mu_btag, h_event_btag, h_topjets_btag, h_tau_btag, h_eff_btag, h_lumi_btag;
    unique_ptr<Hists> h_met, h_jets_met, h_ele_met, h_mu_met, h_event_met, h_topjets_met, h_tau_met, h_eff_met, h_lumi_met;
    unique_ptr<Hists> h_trigger, h_jets_trigger, h_ele_trigger, h_mu_trigger, h_event_trigger, h_topjets_trigger, h_tau_trigger, h_eff_trigger, h_lumi_trigger;
    unique_ptr<Hists> h_triggerSF, h_jets_triggerSF, h_ele_triggerSF, h_mu_triggerSF, h_event_triggerSF, h_topjets_triggerSF, h_tau_triggerSF, h_eff_triggerSF, h_lumi_triggerSF, h_btageff_triggerSF;
    unique_ptr<Hists> h_dibosonsf, h_jets_dibosonsf, h_ele_dibosonsf, h_mu_dibosonsf, h_event_dibosonsf, h_topjets_dibosonsf, h_tau_dibosonsf, h_eff_dibosonsf, h_lumi_dibosonsf, h_fakerate_dibosonsf;
    unique_ptr<Hists> h_3ele, h_jets_3ele, h_ele_3ele, h_mu_3ele, h_event_3ele, h_topjets_3ele, h_tau_3ele, h_eff_3ele, h_lumi_3ele, h_fakerate_3ele;
    unique_ptr<Hists> h_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection, h_tau_finalSelection, h_eff_finalSelection, h_lumi_finalSelection, h_fakerate_finalSelection;

    ElectronId EleId;
    JetId Btag_loose;
    CSVBTag::wp wp_btag_loose;


    bool is_mc, do_pdf_variations, is_dy, is_ele_channel;
    string Sys_BTag, Sys_PU, Sys_EleID, Sys_EleReco, Sys_EleTrigger, Sys_MuonID, Sys_MuonIso, Sys_MuonTrk;
    TString Sys_DibosonXSec, Sys_DibosonBTag;


    vector<unique_ptr<AnalysisModule>> recomodules;
    uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps;
    
  };
  
  
  LQToTopMuFakeRateAnalysisModule::LQToTopMuFakeRateAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuFakeRateAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    is_mc = ctx.get("dataset_type") == "MC";
    is_dy = ctx.get("Channel") == "DY";
    is_ele_channel = ctx.get("EleChannel") == "true";
    if(ctx.get("Channel") != "DY" && ctx.get("Channel") != "Diboson") throw runtime_error("In LQToTopMuFakeRateAnalysisModule.cxx: Channel is neither 'DY' nor 'Diboson'.");
    Sys_EleID = ctx.get("Systematic_EleID");
    Sys_EleReco = ctx.get("Systematic_EleReco");
    Sys_EleTrigger = ctx.get("Systematic_EleTrigger");
    Sys_MuonIso = ctx.get("Systematic_MuonIso");
    Sys_MuonTrk = ctx.get("Systematic_MuonTrk");
    Sys_MuonID = ctx.get("Systematic_MuonID");
    Sys_BTag = ctx.get("Systematic_BTag");
    Sys_PU = ctx.get("Systematic_PU");
    Sys_DibosonXSec = ctx.get("Systematic_DibosonXSec");
    Sys_DibosonBTag = ctx.get("Systematic_DibosonBTag");

    // 1. setup other modules. CommonModules and the JetCleaner:
    if(is_dy) EleId = AndId<Electron>(ElectronID_Spring16_loose,PtEtaCut(30.0, 2.4));
    else      EleId = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(30.0, 2.4));
    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
    wp_btag_loose = CSVBTag::WP_LOOSE;
    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4)); //10

    common.reset(new CommonModules());    
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(EleId);
    common->init(ctx,Sys_PU);

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", Sys_EleReco));
    if(is_dy) SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", Sys_EleID));
    else      SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "", Sys_EleID));
    SF_eleTrigger.reset(new ElectronTriggerWeights(ctx, "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TagProbe/Optimization/35867fb_Iso27_NonIso115/ElectronEfficiencies.root", Sys_EleTrigger));
    SF_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_loose,"jets",Sys_BTag));
    //if(is_dy) SF_Diboson.reset(new DibosonScaleFactors(ctx, "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/ElectronFakeRate/Optimization/35867fb_CorrMET_Diboson_Pt30/DibosonSF.root", Sys_DibosonXSec, Sys_DibosonBTag));
    //if(is_dy) SF_Diboson.reset(new DibosonScaleFactors(ctx, "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/ElectronFakeRate/Optimization/35867fb_FakeRateNoHF_Diboson/DibosonSF.root", Sys_DibosonXSec, Sys_DibosonBTag));
    if(is_dy) SF_Diboson.reset(new DibosonScaleFactors(ctx, "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/ElectronFakeRate/Optimization/35867fb_FakeRateNoHF_Diboson/WithPromptHistograms/DibosonSF_NLO.root", Sys_DibosonXSec, Sys_DibosonBTag));

    if(is_dy && !is_ele_channel){
      SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, Sys_MuonID));
      SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, Sys_MuonIso));
      SF_muonTrk.reset(new MuonTrkWeights(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", Sys_MuonTrk));
    }

    ZEE_finder.reset(new ZEEFinder());
    elejetcleaner.reset(new ElectronJetOverlapCleaner());


    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));

    
    // 2. set up selections
    if(!is_ele_channel && is_dy) nmuon_sel.reset(new NMuonSelection(0,1)); //DY mu, to be able to use ==1 muon for fake rate calculations
    ht_sel.reset(new HtSelection(350., -1));                         //DY
    nbtag_loose_sel.reset(new NJetSelection(0, 2, Btag_loose));      //both
    njet_sel.reset(new NJetSelection(2, -1));                        //DY
    m_ee_sel.reset(new InvMassEleEleSelection(71.,111.));            //both
    if(is_dy) met_sel.reset(new METSelection(0, 60));                //DY
    else      met_sel.reset(new METSelection(60,-1));                //Diboson
    if(is_ele_channel) nele_sel.reset(new NElectronSelection(3, 3)); //DY ele + Diboson ele
    else nele_sel.reset(new NMuonSelection(1, 1));                   //DY mu  + Diboson mu
    if(is_ele_channel){
      if(!is_dy) ele2_sel.reset(new NElectronSelection(2, 3));        //Diboson
      else       ele2_sel.reset(new NElectronSelection(2, 3));         //DY
    }
    else{
      if(!is_dy) ele2_sel.reset(new NElectronSelection(2, 2));        //Diboson
      else       ele2_sel.reset(new NElectronSelection(2, 2));         //DY
    }
    Mee_sel.reset(new InvMassEleEleSelection(71.,111.));


    ele_trigger_sel1.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    ele_trigger_sel2.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
    h_eff_nocuts.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_cleaner.reset(new LQToTopMuHists(ctx, "Cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_Cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_Cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_Cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_Cleaner"));
    h_topjets_cleaner.reset(new TopJetHists(ctx, "TopJets_Cleaner"));
    h_eff_cleaner.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Cleaner"));
    h_lumi_cleaner.reset(new LuminosityHists(ctx, "Lumi_Cleaner"));

    h_2ele.reset(new LQToTopMuHists(ctx, "2Ele"));
    h_jets_2ele.reset(new JetHists(ctx, "Jets_2Ele"));
    h_ele_2ele.reset(new ElectronHists(ctx, "Ele_2Ele"));
    h_mu_2ele.reset(new MuonHists(ctx, "Mu_2Ele"));
    h_event_2ele.reset(new EventHists(ctx, "Event_2Ele"));
    h_topjets_2ele.reset(new TopJetHists(ctx, "TopJets_2Ele"));
    h_eff_2ele.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_2Ele"));
    h_lumi_2ele.reset(new LuminosityHists(ctx, "Lumi_2Ele"));

    h_mee.reset(new LQToTopMuHists(ctx, "Mee"));
    h_jets_mee.reset(new JetHists(ctx, "Jets_Mee"));
    h_ele_mee.reset(new ElectronHists(ctx, "Ele_Mee"));
    h_mu_mee.reset(new MuonHists(ctx, "Mu_Mee"));
    h_event_mee.reset(new EventHists(ctx, "Event_Mee"));
    h_topjets_mee.reset(new TopJetHists(ctx, "TopJets_Mee"));
    h_eff_mee.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Mee"));
    h_lumi_mee.reset(new LuminosityHists(ctx, "Lumi_Mee"));

    h_2jets.reset(new LQToTopMuHists(ctx, "2Jets"));
    h_jets_2jets.reset(new JetHists(ctx, "Jets_2Jets"));
    h_ele_2jets.reset(new ElectronHists(ctx, "Ele_2Jets"));
    h_mu_2jets.reset(new MuonHists(ctx, "Mu_2Jets"));
    h_event_2jets.reset(new EventHists(ctx, "Event_2Jets"));
    h_topjets_2jets.reset(new TopJetHists(ctx, "TopJets_2Jets"));
    h_eff_2jets.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_2Jets"));
    h_lumi_2jets.reset(new LuminosityHists(ctx, "Lumi_2Jets"));

    h_ht.reset(new LQToTopMuHists(ctx, "Ht"));
    h_jets_ht.reset(new JetHists(ctx, "Jets_Ht"));
    h_ele_ht.reset(new ElectronHists(ctx, "Ele_Ht"));
    h_mu_ht.reset(new MuonHists(ctx, "Mu_Ht"));
    h_event_ht.reset(new EventHists(ctx, "Event_Ht"));
    h_topjets_ht.reset(new TopJetHists(ctx, "TopJets_Ht"));
    h_eff_ht.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Ht"));
    h_lumi_ht.reset(new LuminosityHists(ctx, "Lumi_Ht"));

    h_btag.reset(new LQToTopMuHists(ctx, "Btag"));
    h_jets_btag.reset(new JetHists(ctx, "Jets_Btag"));
    h_ele_btag.reset(new ElectronHists(ctx, "Ele_Btag"));
    h_mu_btag.reset(new MuonHists(ctx, "Mu_Btag"));
    h_event_btag.reset(new EventHists(ctx, "Event_Btag"));
    h_topjets_btag.reset(new TopJetHists(ctx, "TopJets_Btag"));
    h_eff_btag.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Btag"));
    h_lumi_btag.reset(new LuminosityHists(ctx, "Lumi_Btag"));

    h_met.reset(new LQToTopMuHists(ctx, "Met"));
    h_jets_met.reset(new JetHists(ctx, "Jets_Met"));
    h_ele_met.reset(new ElectronHists(ctx, "Ele_Met"));
    h_mu_met.reset(new MuonHists(ctx, "Mu_Met"));
    h_event_met.reset(new EventHists(ctx, "Event_Met"));
    h_topjets_met.reset(new TopJetHists(ctx, "TopJets_Met"));
    h_eff_met.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Met"));
    h_lumi_met.reset(new LuminosityHists(ctx, "Lumi_Met"));

    h_trigger.reset(new LQToTopMuHists(ctx, "Trigger"));
    h_jets_trigger.reset(new JetHists(ctx, "Jets_Trigger"));
    h_ele_trigger.reset(new ElectronHists(ctx, "Ele_Trigger"));
    h_mu_trigger.reset(new MuonHists(ctx, "Mu_Trigger"));
    h_event_trigger.reset(new EventHists(ctx, "Event_Trigger"));
    h_topjets_trigger.reset(new TopJetHists(ctx, "TopJets_Trigger"));
    h_eff_trigger.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Trigger"));
    h_lumi_trigger.reset(new LuminosityHists(ctx, "Lumi_Trigger"));

    h_triggerSF.reset(new LQToTopMuHists(ctx, "TriggerSF"));
    h_jets_triggerSF.reset(new JetHists(ctx, "Jets_TriggerSF"));
    h_ele_triggerSF.reset(new ElectronHists(ctx, "Ele_TriggerSF"));
    h_mu_triggerSF.reset(new MuonHists(ctx, "Mu_TriggerSF"));
    h_event_triggerSF.reset(new EventHists(ctx, "Event_TriggerSF"));
    h_topjets_triggerSF.reset(new TopJetHists(ctx, "TopJets_TriggerSF"));
    h_eff_triggerSF.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_TriggerSF"));
    h_lumi_triggerSF.reset(new LuminosityHists(ctx, "Lumi_TriggerSF"));
    h_btageff_triggerSF.reset(new BTagMCEfficiencyHists(ctx, "BTagEff_TriggerSF",wp_btag_loose));

    h_dibosonsf.reset(new LQToTopMuHists(ctx, "DibosonSF"));
    h_jets_dibosonsf.reset(new JetHists(ctx, "Jets_DibosonSF"));
    h_ele_dibosonsf.reset(new ElectronHists(ctx, "Ele_DibosonSF"));
    h_mu_dibosonsf.reset(new MuonHists(ctx, "Mu_DibosonSF"));
    h_event_dibosonsf.reset(new EventHists(ctx, "Event_DibosonSF"));
    h_topjets_dibosonsf.reset(new TopJetHists(ctx, "TopJets_DibosonSF"));
    h_eff_dibosonsf.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_DibosonSF"));
    h_lumi_dibosonsf.reset(new LuminosityHists(ctx, "Lumi_DibosonSF"));
    h_fakerate_dibosonsf.reset(new LQToTopMuFakeRateHists(ctx, "FakeRate_DibosonSF", is_ele_channel));

    h_3ele.reset(new LQToTopMuHists(ctx, "3Ele"));
    h_jets_3ele.reset(new JetHists(ctx, "Jets_3Ele"));
    h_ele_3ele.reset(new ElectronHists(ctx, "Ele_3Ele"));
    h_mu_3ele.reset(new MuonHists(ctx, "Mu_3Ele"));
    h_event_3ele.reset(new EventHists(ctx, "Event_3Ele"));
    h_topjets_3ele.reset(new TopJetHists(ctx, "TopJets_3Ele"));
    h_eff_3ele.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_3Ele"));
    h_lumi_3ele.reset(new LuminosityHists(ctx, "Lumi_3Ele"));
    h_fakerate_3ele.reset(new LQToTopMuFakeRateHists(ctx, "FakeRate_3Ele", is_ele_channel));

    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection")); 
    h_eff_finalSelection.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_FinalSelection"));
    h_lumi_finalSelection.reset(new LuminosityHists(ctx, "Lumi_FinalSelection"));
    h_fakerate_finalSelection.reset(new LQToTopMuFakeRateHists(ctx, "FakeRate_FinalSelection", is_ele_channel));
    
  }
  
  
  bool LQToTopMuFakeRateAnalysisModule::process(Event & event) {
    //cout << endl << endl <<"+++++ NEW EVENT +++++" << endl;

    SF_eleReco->process(event);
    SF_eleID->process(event);
    if(event.electrons->size() >= 1 && event.muons->size() >= 2 && event.jets->size() >= 2){
      for(auto & m : recomodules){
	m->process(event);
      }
    }

    h_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_eff_nocuts->fill(event);
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
    h_eff_cleaner->fill(event);
    h_lumi_cleaner->fill(event);

    if(!ele2_sel->passes(event)) return false;
    if(!is_ele_channel && is_dy && !nmuon_sel->passes(event)) return false;
    if(!is_ele_channel && is_dy && nele_sel->passes(event)){
      SF_muonID->process(event);
      SF_muonIso->process(event);
      SF_muonTrk->process(event);
    }
    h_2ele->fill(event);
    h_jets_2ele->fill(event);
    h_ele_2ele->fill(event);
    h_mu_2ele->fill(event);
    h_event_2ele->fill(event);
    h_topjets_2ele->fill(event);
    h_eff_2ele->fill(event);
    h_lumi_2ele->fill(event);

    if(!Mee_sel->passes(event) && !is_dy) return false;
    h_mee->fill(event);
    h_jets_mee->fill(event);
    h_ele_mee->fill(event);
    h_mu_mee->fill(event);
    h_event_mee->fill(event);
    h_topjets_mee->fill(event);
    h_eff_mee->fill(event);
    h_lumi_mee->fill(event);

    //find 2 electrons closest to Z mass and delete jets overlapping with them
    pair<int,int> best_ele = ZEE_finder->search(event);
    elejetcleaner->process(event, best_ele.first, best_ele.second);
   
    //2Jets
    if(!njet_sel->passes(event) && is_dy) return false;
    h_2jets->fill(event);
    h_jets_2jets->fill(event);
    h_ele_2jets->fill(event);
    h_mu_2jets->fill(event);
    h_event_2jets->fill(event);
    h_topjets_2jets->fill(event);
    h_eff_2jets->fill(event);
    h_lumi_2jets->fill(event);

    //HT
    if(!ht_sel->passes(event) && is_dy) return false;
    h_ht->fill(event);
    h_jets_ht->fill(event);
    h_ele_ht->fill(event);
    h_mu_ht->fill(event);
    h_event_ht->fill(event);
    h_topjets_ht->fill(event);
    h_eff_ht->fill(event);
    h_lumi_ht->fill(event);

    //BTag
    if(!nbtag_loose_sel->passes(event)) return false;
    SF_btag->process(event);

    h_btag->fill(event);
    h_jets_btag->fill(event);
    h_ele_btag->fill(event);
    h_mu_btag->fill(event);
    h_event_btag->fill(event);
    h_topjets_btag->fill(event);
    h_eff_btag->fill(event);
    h_lumi_btag->fill(event);

    //MET
    if(!met_sel->passes(event)) return false;
    h_met->fill(event);
    h_jets_met->fill(event);
    h_ele_met->fill(event);
    h_mu_met->fill(event);
    h_event_met->fill(event);
    h_topjets_met->fill(event);
    h_eff_met->fill(event);
    h_lumi_met->fill(event);

    //Trigger
    if(!(ele_trigger_sel1->passes(event) || ele_trigger_sel2->passes(event)) ) return false;
    h_trigger->fill(event);
    h_jets_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_mu_trigger->fill(event);
    h_event_trigger->fill(event);
    h_topjets_trigger->fill(event);
    h_eff_trigger->fill(event);
    h_lumi_trigger->fill(event);

    //TriggerSF
    SF_eleTrigger->process(event);
    h_triggerSF->fill(event);
    h_jets_triggerSF->fill(event);
    h_ele_triggerSF->fill(event);
    h_mu_triggerSF->fill(event);
    h_event_triggerSF->fill(event);
    h_topjets_triggerSF->fill(event);
    h_eff_triggerSF->fill(event);
    h_lumi_triggerSF->fill(event);
    h_btageff_triggerSF->fill(event);

    //DibosonSF
    if(is_dy) SF_Diboson->process(event);
    h_dibosonsf->fill(event);
    h_jets_dibosonsf->fill(event);
    h_ele_dibosonsf->fill(event);
    h_mu_dibosonsf->fill(event);
    h_event_dibosonsf->fill(event);
    h_topjets_dibosonsf->fill(event);
    h_eff_dibosonsf->fill(event);
    h_lumi_dibosonsf->fill(event);
    if(is_dy) h_fakerate_dibosonsf->fill(event); //only needed for DY region

    //3 Ele
    if(!nele_sel->passes(event)) return false;
    h_3ele->fill(event);
    h_jets_3ele->fill(event);
    h_ele_3ele->fill(event);
    h_mu_3ele->fill(event);
    h_event_3ele->fill(event);
    h_topjets_3ele->fill(event);
    h_eff_3ele->fill(event);
    h_lumi_3ele->fill(event);
    h_fakerate_3ele->fill(event);

    //Final Selection
    h_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);
    h_eff_finalSelection->fill(event);
    h_lumi_finalSelection->fill(event);
    h_fakerate_finalSelection->fill(event);

    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuFakeRateAnalysisModule)
} 

