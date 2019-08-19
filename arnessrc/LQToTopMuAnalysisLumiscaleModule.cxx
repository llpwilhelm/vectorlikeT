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
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/LQToTopMuEfficiencyHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuMLQRecoHists.h"
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

  class LQToTopMuAnalysisLumiscaleModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuAnalysisLumiscaleModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    unique_ptr<CommonModules> common;
    unique_ptr<AnalysisModule> Muon_printer, Electron_printer, Jet_printer, GenParticles_printer, syst_module, SF_muonID, SF_muonIso, SF_btag, SF_eleReco, SF_eleID;
    unique_ptr<AnalysisModule> ttgenprod, LQgenprod;    
    unique_ptr<AnalysisModule> sf_top_pt_reweight;
    unique_ptr<MCMuonScaleFactor> SF_muonTrigger;
    unique_ptr<ElectronFakeRateWeights> SF_eleFakeRate;
    unique_ptr<MuonTrkWeights> SF_muonTrk;

    unique_ptr<JetCleaner> jetcleaner;
    
    // declare the Selections to use.
    unique_ptr<Selection>  nele_sel, njet_sel, nbtag_loose_sel, m_mumu_veto, htlept_sel, mttbar_gen_sel, m_muele_veto, m_ee_veto, ht_sel, dr_lepjet_sel, genlvl_TopDilepton_sel, ele_trigger_sel1, ele_trigger_sel2, genlvl_ZMuMu_sel, genlvl_ZEE_sel, trig1_sel, trig2_sel, trig3_sel, trig4_sel, lq_matchable_sel, ttbar_matchable_sel;
    
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_mlqreco_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts, h_eff_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_hyphists, h_1ele, h_jets_1ele, h_ele_1ele, h_mu_1ele, h_event_1ele, h_topjets_1ele, h_tau_1ele, h_eff_1ele, h_lumi_1ele;
    unique_ptr<Hists> h_1bJetLoose, h_mlqreco_1bJetLoose, h_jets_1bJetLoose, h_mu_1bJetLoose, h_ele_1bJetLoose, h_event_1bJetLoose, h_topjets_1bJetLoose, h_tau_1bJetLoose, h_eff_1bJetLoose, h_lumi_1bJetLoose;
    unique_ptr<Hists> h_InvMassVeto, h_mlqreco_InvMassVeto, h_jets_InvMassVeto, h_mu_InvMassVeto, h_ele_InvMassVeto, h_event_InvMassVeto, h_topjets_InvMassVeto, h_tau_InvMassVeto, h_eff_InvMassVeto, h_lumi_InvMassVeto;
    unique_ptr<Hists> h_htlept200, h_mlqreco_htlept200, h_jets_htlept200, h_ele_htlept200, h_mu_htlept200, h_event_htlept200, h_topjets_htlept200, h_tau_htlept200, h_btageff_htlept200, h_eff_htlept200, h_lumi_htlept200;
    unique_ptr<Hists> h_ttbarmatchable, h_mlqreco_ttbarmatchable, h_jets_ttbarmatchable, h_ele_ttbarmatchable, h_mu_ttbarmatchable, h_event_ttbarmatchable, h_topjets_ttbarmatchable, h_tau_ttbarmatchable, h_eff_ttbarmatchable, h_lumi_ttbarmatchable;
    unique_ptr<Hists> h_lqmatchable, h_mlqreco_lqmatchable, h_jets_lqmatchable, h_ele_lqmatchable, h_mu_lqmatchable, h_event_lqmatchable, h_topjets_lqmatchable, h_tau_lqmatchable, h_eff_lqmatchable, h_lumi_lqmatchable;
    unique_ptr<Hists> h_correctlymatched, h_mlqreco_correctlymatched, h_jets_correctlymatched, h_ele_correctlymatched, h_mu_correctlymatched, h_event_correctlymatched, h_topjets_correctlymatched, h_tau_correctlymatched, h_eff_correctlymatched, h_lumi_correctlymatched;
    unique_ptr<Hists> h_finalSelection, h_mlqreco_finalSelection, h_jets_finalSelection, h_ele_finalSelection, h_mu_finalSelection, h_event_finalSelection, h_topjets_finalSelection, h_tau_finalSelection, h_eff_finalSelection, h_lumi_finalSelection;
    unique_ptr<Hists> h_ht_InvMassVeto, h_ht_finalSelection;
    unique_ptr<Hists> h_Sideband;
    unique_ptr<Hists> h_PDF_variations;

    
    MuonId MuId;
    ElectronId EleId;
    JetId Btag_loose, Btag_medium, Btag_tight;
    CSVBTag::wp wp_btag_loose;
    int n_misid;


    Event::Handle<TTbarGen> h_ttbargen;
    Event::Handle<LQGen> h_LQLQbargen;

    //additional variables to be written to the output
    Event::Handle<double> handle_evtwgt, handle_st, handle_stlep;
    Event::Handle<vector<double>> handle_mmumu;
    Event::Handle<int> handle_nbjets, handle_njets;


    bool do_scale_variation, is_mc, do_pdf_variations, is_foreff, forML;
    string Sys_MuonID, Sys_BTag, Sys_MuonTrigger, Sys_MuonIso, Sys_PU, Sys_EleID, Sys_EleReco, Sys_WJ, Sys_TTbar, Sys_DY, Sys_ST, Sys_DB, Sys_QCD, Sys_TTV, Sys_EleFake, Sys_MuFake, Sys_MuonTrk;
    TString dataset_version, scalevariation_process;
 
    
    vector<unique_ptr<AnalysisModule>> recomodules, muonic_recomodules, inclusive_recomodules;
    uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps, h_muonic_hyps, h_inclusive_hyps;

    uhh2::Event::Handle<double> h_FakeRateWeightEle, h_FakeRateWeightMu;
    
  };
  
  
  LQToTopMuAnalysisLumiscaleModule::LQToTopMuAnalysisLumiscaleModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuAnalysisLumiscaleModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    n_misid = 0;
    scalevariation_process   = ctx.get("ScaleVariationProcess");
    do_scale_variation       = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    if(do_scale_variation && (scalevariation_process != "ttbar") && (scalevariation_process != "dy") && (scalevariation_process != "st") && (scalevariation_process != "wj") && (scalevariation_process != "db") && (scalevariation_process != "ttv")) throw runtime_error("In LQToTopMuAnalysisLumiscaleModule.cxx: Invalid process specified for 'ScaleVariationProcess'.");
    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    dataset_version = ctx.get("dataset_version");
    is_mc = ctx.get("dataset_type") == "MC";
    is_foreff = ctx.get("IsForEff") == "true";
    forML = ctx.get("ForML") == "true";
    Sys_MuonID = ctx.get("Systematic_MuonID");
    Sys_MuonTrigger = ctx.get("Systematic_MuonTrigger");
    Sys_MuonIso = ctx.get("Systematic_MuonIso");
    Sys_MuonTrk = ctx.get("Systematic_MuonTrk");
    Sys_MuFake = ctx.get("Systematic_MuFake");
    Sys_EleID = ctx.get("Systematic_EleID");
    Sys_EleReco = ctx.get("Systematic_EleReco");
    Sys_EleFake = ctx.get("Systematic_EleFake");
    Sys_BTag = ctx.get("Systematic_BTag");
    Sys_PU = ctx.get("Systematic_PU");
    Sys_TTbar = ctx.get("Systematic_TTbar");
    Sys_DY = ctx.get("Systematic_DY");
    Sys_ST = ctx.get("Systematic_ST");
    Sys_DB = ctx.get("Systematic_DB");
    Sys_QCD = ctx.get("Systematic_QCD");
    Sys_TTV = ctx.get("Systematic_TTV");
    Sys_WJ = ctx.get("Systematic_WJ");


    // 1. setup other modules. CommonModules and the JetCleaner:
    Jet_printer.reset(new JetPrinter("Jet-Printer", 0));
    Electron_printer.reset(new ElectronPrinter("Electron-Printer"));
    Muon_printer.reset(new MuonPrinter("Muon-Printer"));
    GenParticles_printer.reset(new GenParticlesPrinter(ctx));
    MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 4.0),MuonIso(0.15));
    //MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4));
    EleId = AndId<Electron>(ElectronID_Spring16_loose,PtEtaCut(30.0, 4.0));

    Btag_loose = CSVBTag(CSVBTag::WP_LOOSE);
    wp_btag_loose = CSVBTag::WP_LOOSE;
    jetcleaner.reset(new JetCleaner(ctx,30.0, 4.0));


    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->init(ctx,Sys_PU);

    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, Sys_MuonID));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, Sys_MuonTrigger));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, Sys_MuonIso));
    SF_muonTrk.reset(new MuonTrkWeights(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", Sys_MuonTrk));

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", Sys_EleReco));
    SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", Sys_EleID));

    SF_btag.reset(new MCBTagScaleFactor(ctx,wp_btag_loose,"jets",Sys_BTag));

    //sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.0615, -0.0005, "", "", true, 1.0)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting

    
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
    h_muonic_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassMuonicLQReconstruction");
    h_inclusive_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassInclusiveLQReconstruction");
    
    if(is_mc){
      ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
      h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
      LQgenprod.reset(new LQGenProducer(ctx, "LQLQbargen", false));
      h_LQLQbargen = ctx.get_handle<LQGen>("LQLQbargen");
    }
    

    if(Sys_EleFake == "nominal")   h_FakeRateWeightEle =ctx.get_handle<double>("FakeRateWeightEle");
    else if(Sys_EleFake == "up")   h_FakeRateWeightEle =ctx.get_handle<double>("FakeRateWeightEleUp");
    else if(Sys_EleFake == "down") h_FakeRateWeightEle =ctx.get_handle<double>("FakeRateWeightEleDown");
    else throw runtime_error("Sys_EleFake is not one of the following: ['up', 'down', 'nominal']");
    if(Sys_MuFake == "nominal")   h_FakeRateWeightMu =ctx.get_handle<double>("FakeRateWeightMu");
    else if(Sys_MuFake == "up")   h_FakeRateWeightMu =ctx.get_handle<double>("FakeRateWeightMuUp");
    else if(Sys_MuFake == "down") h_FakeRateWeightMu =ctx.get_handle<double>("FakeRateWeightMuDown");
    else throw runtime_error("Sys_MuFake is not one of the following: ['up', 'down', 'nominal']");

    
    // 2. set up selections
    //Selection
    njet_sel.reset(new NJetSelection(2, -1));  
    nbtag_loose_sel.reset(new NJetSelection(1, -1, Btag_loose));  //1, -1 
    m_mumu_veto.reset(new InvMass2MuVeto(0, 111));
    m_ee_veto.reset(new InvMassEleEleVeto(0.,111.));
    nele_sel.reset(new NElectronSelection(1, -1));
    htlept_sel.reset(new HTLeptSelection(200., -1));
    ht_sel.reset(new HtSelection(350., -1));
    genlvl_ZMuMu_sel.reset(new GenLvlZMuMuSelection());
    genlvl_ZEE_sel.reset(new GenLvlZEESelection());
    lq_matchable_sel.reset(new LQSemiLepMatchable());
    ttbar_matchable_sel.reset(new TTbarSemiLepMatchable());


    //ele_trigger_sel1.reset(new TriggerSelection("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*"));
    ele_trigger_sel1.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    ele_trigger_sel2.reset(new TriggerSelection("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"));
    trig1_sel.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trig2_sel.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    //trig3_sel.reset(new TriggerSelection("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v*"));
    trig3_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
    trig4_sel.reset(new TriggerSelection("HLT_TkMu50_v*"));
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));
    if(is_mc) recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassLQReconstruction"));
    muonic_recomodules.emplace_back(new HighMassMuonicLQReconstruction(ctx,LQNeutrinoReconstruction));
    muonic_recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassMuonicLQReconstruction"));
    if(is_mc) muonic_recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassMuonicLQReconstruction"));
    inclusive_recomodules.emplace_back(new HighMassInclusiveLQReconstruction(ctx,LQNeutrinoReconstruction));
    inclusive_recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassInclusiveLQReconstruction"));
    if(is_mc) inclusive_recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassInclusiveLQReconstruction"));

    //systematics modules
    syst_module.reset(new MCScaleVariation(ctx));


    //additional variables to be written to the AnalysisTree
    handle_evtwgt = ctx.declare_event_output<double>("EventWeight");
    handle_nbjets = ctx.declare_event_output<int>("NBJets");
    handle_njets  = ctx.declare_event_output<int>("NJets");
    handle_st     = ctx.declare_event_output<double>("ST");
    handle_stlep  = ctx.declare_event_output<double>("STLep");
    handle_mmumu  = ctx.declare_event_output<vector<double>>("MMuMu");



    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_mlqreco_nocuts.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
    h_eff_nocuts.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));

    h_hyphists.reset(new HypothesisHistsOwn(ctx, "CorrectMatch_Hists", "HighMassLQReconstruction", "CorrectMatch"));
 
    h_1ele.reset(new LQToTopMuHists(ctx, "1Ele"));
    h_jets_1ele.reset(new JetHists(ctx, "Jets_1Ele"));
    h_ele_1ele.reset(new ElectronHists(ctx, "Ele_1Ele"));
    h_mu_1ele.reset(new MuonHists(ctx, "Mu_1Ele"));  
    h_event_1ele.reset(new EventHists(ctx, "Event_1Ele"));
    h_topjets_1ele.reset(new TopJetHists(ctx, "TopJets_1Ele"));
    h_eff_1ele.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_1Ele"));
    h_lumi_1ele.reset(new LuminosityHists(ctx, "Lumi_1Ele"));

    h_1bJetLoose.reset(new LQToTopMuHists(ctx, "1bJetLoose"));
    h_mlqreco_1bJetLoose.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_1bJetLoose"));
    h_jets_1bJetLoose.reset(new JetHists(ctx, "Jets_1bJetLoose"));
    h_ele_1bJetLoose.reset(new ElectronHists(ctx, "Ele_1bJetLoose"));
    h_mu_1bJetLoose.reset(new MuonHists(ctx, "Mu_1bJetLoose"));
    h_event_1bJetLoose.reset(new EventHists(ctx, "Event_1bJetLoose"));
    h_topjets_1bJetLoose.reset(new TopJetHists(ctx, "TopJets_1bJetLoose"));
    h_eff_1bJetLoose.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_1bJetLoose"));
    h_lumi_1bJetLoose.reset(new LuminosityHists(ctx, "Lumi_1bJetLoose"));


    h_InvMassVeto.reset(new LQToTopMuHists(ctx, "InvMassVeto"));
    h_mlqreco_InvMassVeto.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_InvMassVeto"));
    h_jets_InvMassVeto.reset(new JetHists(ctx, "Jets_InvMassVeto"));
    h_ele_InvMassVeto.reset(new ElectronHists(ctx, "Ele_InvMassVeto"));
    h_mu_InvMassVeto.reset(new MuonHists(ctx, "Mu_InvMassVeto"));
    h_event_InvMassVeto.reset(new EventHists(ctx, "Event_InvMassVeto"));
    h_topjets_InvMassVeto.reset(new TopJetHists(ctx, "TopJets_InvMassVeto"));
    h_eff_InvMassVeto.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_InvMassVeto"));
    h_lumi_InvMassVeto.reset(new LuminosityHists(ctx, "Lumi_InvMassVeto"));


    h_htlept200.reset(new LQToTopMuHists(ctx, "HTLept200"));
    h_mlqreco_htlept200.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_HTLept200"));
    h_jets_htlept200.reset(new JetHists(ctx, "Jets_HTLept200"));
    h_ele_htlept200.reset(new ElectronHists(ctx, "Ele_HTLept200"));
    h_mu_htlept200.reset(new MuonHists(ctx, "Mu_HTLept200"));
    h_topjets_htlept200.reset(new TopJetHists(ctx, "TopJets_HTLept200"));
    h_event_htlept200.reset(new EventHists(ctx, "Event_HTLept200"));
    h_eff_htlept200.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_HTLept200"));
    h_lumi_htlept200.reset(new LuminosityHists(ctx, "Lumi_HTLept200"));
    h_btageff_htlept200.reset(new BTagMCEfficiencyHists(ctx, "BTagEff_HTLept200",wp_btag_loose));


    h_ttbarmatchable.reset(new LQToTopMuHists(ctx, "TTbar_Matchable"));
    h_mlqreco_ttbarmatchable.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_TTbar_Matchable"));
    h_jets_ttbarmatchable.reset(new JetHists(ctx, "Jets_TTbar_Matchable"));
    h_ele_ttbarmatchable.reset(new ElectronHists(ctx, "Ele_TTbar_Matchable"));
    h_mu_ttbarmatchable.reset(new MuonHists(ctx, "Mu_TTbar_Matchable"));
    h_topjets_ttbarmatchable.reset(new TopJetHists(ctx, "TopJets_TTbar_Matchable"));
    h_event_ttbarmatchable.reset(new EventHists(ctx, "Event_TTbar_Matchable")); 
    h_eff_ttbarmatchable.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_TTbar_Matchable"));
    h_lumi_ttbarmatchable.reset(new LuminosityHists(ctx, "Lumi_TTbar_Matchable"));


    h_lqmatchable.reset(new LQToTopMuHists(ctx, "LQ_Matchable"));
    h_mlqreco_lqmatchable.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_LQ_Matchable"));
    h_jets_lqmatchable.reset(new JetHists(ctx, "Jets_LQ_Matchable"));
    h_ele_lqmatchable.reset(new ElectronHists(ctx, "Ele_LQ_Matchable"));
    h_mu_lqmatchable.reset(new MuonHists(ctx, "Mu_LQ_Matchable"));
    h_topjets_lqmatchable.reset(new TopJetHists(ctx, "TopJets_LQ_Matchable"));
    h_event_lqmatchable.reset(new EventHists(ctx, "Event_LQ_Matchable")); 
    h_eff_lqmatchable.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_LQ_Matchable"));
    h_lumi_lqmatchable.reset(new LuminosityHists(ctx, "Lumi_LQ_Matchable"));


    h_correctlymatched.reset(new LQToTopMuHists(ctx, "CorrectlyMatched"));
    h_mlqreco_correctlymatched.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_CorrectlyMatched"));
    h_jets_correctlymatched.reset(new JetHists(ctx, "Jets_CorrectlyMatched"));
    h_ele_correctlymatched.reset(new ElectronHists(ctx, "Ele_CorrectlyMatched"));
    h_mu_correctlymatched.reset(new MuonHists(ctx, "Mu_CorrectlyMatched"));
    h_topjets_correctlymatched.reset(new TopJetHists(ctx, "TopJets_CorrectlyMatched"));
    h_event_correctlymatched.reset(new EventHists(ctx, "Event_CorrectlyMatched")); 
    h_eff_correctlymatched.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_CorrectlyMatched"));
    h_lumi_correctlymatched.reset(new LuminosityHists(ctx, "Lumi_CorrectlyMatched"));


    h_finalSelection.reset(new LQToTopMuHists(ctx, "FinalSelection"));
    h_mlqreco_finalSelection.reset(new LQToTopMuMLQRecoHists(ctx, "MLQReco_FinalSelection"));
    h_jets_finalSelection.reset(new JetHists(ctx, "Jets_FinalSelection"));
    h_ele_finalSelection.reset(new ElectronHists(ctx, "Ele_FinalSelection"));
    h_mu_finalSelection.reset(new MuonHists(ctx, "Mu_FinalSelection"));
    h_topjets_finalSelection.reset(new TopJetHists(ctx, "TopJets_FinalSelection"));
    h_event_finalSelection.reset(new EventHists(ctx, "Event_FinalSelection")); 
    h_eff_finalSelection.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_FinalSelection"));
    h_lumi_finalSelection.reset(new LuminosityHists(ctx, "Lumi_FinalSelection"));
    h_ht_finalSelection.reset(new HT2dHists(ctx, "HT2d_FinalSelection"));
    h_PDF_variations.reset(new LQToTopMuPDFHists(ctx, "PDF_variations", true, do_pdf_variations));
    h_Sideband.reset(new LQToTopMuHists(ctx, "Sideband_weights_applied"));
    
    
  }
  
  
  bool LQToTopMuAnalysisLumiscaleModule::process(Event & event) {
    //cout << endl << endl << "+++NEW EVENT+++" << endl;

    // int ntags = 0;
    // int upper = event.jets->size();
    // if(upper > 7) upper = 7;
    // for(unsigned int i=0; i<upper; i++) {
    //   if(Btag_loose(event.jets->at(i),event)){
    // 	ntags++;
    //   }
    // }
    // if(ntags < 2) return false;

    if(is_mc){
      double factor_xsec = -1;
      int control = (dataset_version.Contains("TTbar") && Sys_TTbar != "nominal") + (dataset_version.Contains("DYJets") && Sys_DY != "nominal") + (dataset_version.Contains("SingleTop") && Sys_ST != "nominal") + (dataset_version.Contains("WJets") && Sys_WJ != "nominal") + (dataset_version.Contains("Diboson") && Sys_DB != "nominal") + (dataset_version.Contains("QCD") && Sys_QCD != "nominal")+ (dataset_version.Contains("TTV") && Sys_TTV != "nominal");
      if(!(control == 0 || control == 1)) throw runtime_error("In LQToTopMuSidebandAnalysisModule.cxx: More than one rate systematic is set to something different than 'nominal'");

      if(control == 0) factor_xsec = 0;
      else if(control == 1){
	if(dataset_version.Contains("TTbar") && Sys_TTbar != "nominal")       factor_xsec = 0.056;
	else if(dataset_version.Contains("DYJets") && Sys_DY != "nominal")    factor_xsec = 0.1;
	else if(dataset_version.Contains("SingleTop") && Sys_ST != "nominal") factor_xsec = 0.1;
	else if(dataset_version.Contains("WJets") && Sys_WJ != "nominal")     factor_xsec = 0.1;
	else if(dataset_version.Contains("Diboson") && Sys_DB != "nominal")   factor_xsec = 0.2;
	else if(dataset_version.Contains("QCD") && Sys_QCD != "nominal")      factor_xsec = 1;
	else if(dataset_version.Contains("TTV") && Sys_TTV != "nominal")      factor_xsec = 0.25;
	else if(dataset_version.Contains("LQ"))                               factor_xsec = 0;
      }
      double sf_xsec = 1;
      if(Sys_TTbar == "up" || Sys_DY == "up"|| Sys_ST == "up"|| Sys_DB == "up"|| Sys_WJ == "up"|| Sys_QCD == "up"|| Sys_TTV == "up") sf_xsec += factor_xsec;
      else if(Sys_TTbar == "down" || Sys_DY == "down"|| Sys_ST == "down"|| Sys_DB == "down"|| Sys_WJ == "down"|| Sys_QCD == "down"|| Sys_TTV == "down") sf_xsec -= factor_xsec;
      else if(control != 0) throw runtime_error("In LQToTopMuAnalysisLumiscaleModule.cxx: Invalid direction for 'Sys_Rate_YYY' specified.");

      event.weight *= sf_xsec;
    }


    //Jet_printer->process(event);
    //Muon_printer->process(event);
    //Electron_printer->process(event);
    //GenParticles_printer->process(event);

    
    //Apply fake-rate scale factors
    double SF_ele = 1, SF_mu = 1;
    double n_matched_to_taus = 0, n_matched_to_muons = 0;
    if(is_mc){
      //cout <<endl << "NEW EVENT" << endl << "Before applying SF: weight = " << event.weight << endl;
      SF_ele = event.get(h_FakeRateWeightEle);
      if(event.electrons->size() == 0 && SF_ele != 1) throw runtime_error("There are no electrons in the event, still the fake-rate SF is != 1...");
      event.weight *= SF_ele;    
      SF_mu = event.get(h_FakeRateWeightMu);
      if(event.muons->size() == 0 && SF_mu != 1) throw runtime_error("There are no muons in the event, still the fake-rate SF is != 1...");     
      event.weight *= SF_mu;
    }
    
    if(!is_foreff){
      //apply muon SFs as in preselection
      SF_muonTrigger->process_onemuon(event,0);
      SF_muonID->process(event);
      SF_muonIso->process(event);
      SF_muonTrk->process(event);

      if(event.electrons->size() >= 1){
	SF_eleReco->process(event);
	SF_eleID->process(event);
      }
    }

    //Only for NoLep selection
    if(is_foreff){
      if(!(ele_trigger_sel1->passes(event) || ele_trigger_sel2->passes(event)) ) return false;
      //if(!ht_sel->passes(event)) return false;
    }
     
    if(do_scale_variation){
      if((dataset_version.Contains("TTbar") && scalevariation_process == "ttbar") || (dataset_version.Contains("DYJets") && scalevariation_process == "dy") || (dataset_version.Contains("SingleTop") && scalevariation_process == "st") || (dataset_version.Contains("WJets") && scalevariation_process == "wj") || (dataset_version.Contains("Diboson") && scalevariation_process == "db") || (dataset_version.Contains("TTV") && scalevariation_process == "ttv")){
	if(event.genInfo->systweights().size() < 10 && dataset_version.Contains("Diboson")) cout << "SystWeight size: " << event.genInfo->systweights().size() << endl;
	else syst_module->process(event);    
      }
      //Signals not needed for scale variation!
      else if(dataset_version.Contains("LQtoT")) return false;
    }

    

    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);
    if(!njet_sel->passes(event)) return false;
    if(!ht_sel->passes(event)) return false;



    //calculate additional vairables
    int nbjets = 0;
    int njets = 0;
    double st = 0;
    double stlep = 0;
    vector<double> mmumu;
    for(unsigned int i=0; i<event.jets->size(); i++){
      if(Btag_loose(event.jets->at(i), event)) nbjets++;
      njets++;
    }
    
    st    = event.met->pt();
    for(const auto & jet : *event.jets){
      st += jet.pt();
    }
    for(const auto & electron : *event.electrons){
      st += electron.pt();
      stlep += electron.pt();
    }
    for(const auto & muon : *event.muons){
      st += muon.pt();
      stlep += muon.pt();
    }

    for(unsigned int i=0; i<event.muons->size(); i++){
      for(unsigned int j=0; j<event.muons->size(); j++){
	if(j > i){
	  mmumu.emplace_back((event.muons->at(i).v4() + event.muons->at(j).v4()).M());
	}
      }
    }


    //fill additional variables for ML
    if(forML){
      event.set(handle_evtwgt, event.weight);
      event.set(handle_nbjets, nbjets);
      event.set(handle_njets, njets);
      event.set(handle_st, st);
      event.set(handle_stlep, stlep);
      event.set(handle_mmumu, mmumu);
      return true;
    }






    //check for at least 1 muon pair with opposite charge
    bool charge_opposite = false;
    for(unsigned int i=0; i<event.muons->size(); i++){
      for(unsigned int j=0; j<event.muons->size(); j++){
	if(j>i){
	  if(event.muons->at(i).charge() != event.muons->at(j).charge()) {
	    charge_opposite = true;
	  }
	}
      }
    }
    double sum_mu_charge = 0, sum_mu_charge_first3=0;
    int idx=0;
    for(const auto & mu : *event.muons){
      sum_mu_charge += mu.charge();
      if(idx < 3) sum_mu_charge_first3 += mu.charge();
      idx++;
    }
    
    bool reconstruct_mlq_mu = (event.electrons->size() == 0 && event.muons->size() == 3 && fabs(sum_mu_charge) == 1);
    bool reconstruct_mlq_ele = (event.electrons->size() >= 1 && event.muons->size() >= 2 && charge_opposite);
    bool reconstruct_mlq = reconstruct_mlq_ele || reconstruct_mlq_mu;

    
    //if(!reconstruct_mlq) return false;
    //if(event.muons->size() + event.electrons->size() < 3) return false;
    

    bool reconstruct_mlq_inclusive = (event.muons->size() >= 2 && event.electrons->size() >= 1 && charge_opposite) || (event.electrons->size() == 0 && event.muons->size() >= 3 && charge_opposite && fabs(sum_mu_charge_first3) == 1);
    
    if(is_mc){
      ttgenprod->process(event);
      LQgenprod->process(event);
    }
    
    double CorrectMatchDiscriminator = 99;
    // MLQ reco
    if(reconstruct_mlq_ele){
      for(auto & m : recomodules){
	m->process(event);
      }
      
      if(is_mc){
	std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps); 
	const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
	CorrectMatchDiscriminator = hyp->discriminator("CorrectMatch");
	//if(!hyp->is_real_neutrino()) return false;
      }
      
    }

    
    if(reconstruct_mlq_mu){
      for(auto & m : muonic_recomodules){
	m->process(event);
      }
      if(is_mc){
	std::vector<LQReconstructionHypothesis> hyps = event.get(h_muonic_hyps); 
	const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
	CorrectMatchDiscriminator = hyp->discriminator("CorrectMatch");
	//if(!hyp->is_real_neutrino()) return false;
      }
    }

    if(reconstruct_mlq_inclusive){
      
      for(auto & m : inclusive_recomodules){
	m->process(event);
      }
      
      if(is_mc){
	std::vector<LQReconstructionHypothesis> hyps = event.get(h_inclusive_hyps); 
	const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
	CorrectMatchDiscriminator = hyp->discriminator("CorrectMatch");
      }


    }

    /*
    if(!(trig1_sel->passes(event) || trig2_sel->passes(event) || trig3_sel->passes(event))) return false;   
    if(trig3_sel->passes(event) && !trig2_sel->passes(event) && !trig1_sel->passes(event)){
      if(event.met->pt() < 200) return false;
    }
    */
    /*
    if(!(trig1_sel->passes(event) || trig2_sel->passes(event) || trig3_sel->passes(event) || trig4_sel->passes(event))) return false;   
    if((trig3_sel->passes(event) || trig4_sel->passes(event))  && !trig2_sel->passes(event) && !trig1_sel->passes(event)){
      if(event.muons->at(0).pt() < 55) return false;
    }
    */

    //sf_top_pt_reweight->process(event);  

    

  
    h_nocuts->fill(event);
    h_mlqreco_nocuts->fill(event);
    h_jets_nocuts->fill(event);
    h_ele_nocuts->fill(event);
    h_mu_nocuts->fill(event);
    h_event_nocuts->fill(event);
    h_topjets_nocuts->fill(event);
    h_eff_nocuts->fill(event);
    h_lumi_nocuts->fill(event);

    if(reconstruct_mlq_ele && is_mc) h_hyphists->fill(event);


    //1 bTag1 loose
    if(!nbtag_loose_sel->passes(event)) return false; 
    SF_btag->process(event);
    //LOESCHEMICH
    h_jets_1bJetLoose->fill(event);
    h_1bJetLoose->fill(event);
    h_mlqreco_1bJetLoose->fill(event);
    h_ele_1bJetLoose->fill(event);
    h_mu_1bJetLoose->fill(event);
    h_event_1bJetLoose->fill(event);
    h_topjets_1bJetLoose->fill(event);
    h_eff_1bJetLoose->fill(event);
    h_lumi_1bJetLoose->fill(event);

    //InvMassVeto
    if(!(m_mumu_veto->passes(event) || is_foreff)) return false; 
    h_jets_InvMassVeto->fill(event);
    h_InvMassVeto->fill(event);
    h_mlqreco_InvMassVeto->fill(event);
    h_ele_InvMassVeto->fill(event);
    h_mu_InvMassVeto->fill(event);
    h_event_InvMassVeto->fill(event);
    h_topjets_InvMassVeto->fill(event);
    h_eff_InvMassVeto->fill(event);
    h_lumi_InvMassVeto->fill(event);

    //HT Lept
    if(!htlept_sel->passes(event)) return false;
    h_jets_htlept200->fill(event);
    h_htlept200->fill(event);
    h_mlqreco_htlept200->fill(event);
    h_ele_htlept200->fill(event);
    h_mu_htlept200->fill(event);
    h_event_htlept200->fill(event);
    h_topjets_htlept200->fill(event);    
    h_eff_htlept200->fill(event);
    h_lumi_htlept200->fill(event);
    h_btageff_htlept200->fill(event);

    if(is_mc && reconstruct_mlq && !dataset_version.Contains("LQtoTTau")&& !dataset_version.Contains("LQtoTLep")){
      if(ttbar_matchable_sel->passes(event)){
	h_ttbarmatchable->fill(event);
	h_mlqreco_ttbarmatchable->fill(event);
	h_jets_ttbarmatchable->fill(event);
	h_ele_ttbarmatchable->fill(event);
	h_mu_ttbarmatchable->fill(event);
	h_event_ttbarmatchable->fill(event);
	h_topjets_ttbarmatchable->fill(event);
	h_eff_ttbarmatchable->fill(event);
	h_lumi_ttbarmatchable->fill(event);	
      }

      if(lq_matchable_sel->passes(event)){
	h_lqmatchable->fill(event);
	h_mlqreco_lqmatchable->fill(event);
	h_jets_lqmatchable->fill(event);
	h_ele_lqmatchable->fill(event);
	h_mu_lqmatchable->fill(event);
	h_event_lqmatchable->fill(event);
	h_topjets_lqmatchable->fill(event);
	h_eff_lqmatchable->fill(event);
	h_lumi_lqmatchable->fill(event);
	
	if(CorrectMatchDiscriminator < 99){
	  h_correctlymatched->fill(event);
	  h_mlqreco_correctlymatched->fill(event);
	  h_jets_correctlymatched->fill(event);
	  h_ele_correctlymatched->fill(event);
	  h_mu_correctlymatched->fill(event);
	  h_event_correctlymatched->fill(event);
	  h_topjets_correctlymatched->fill(event);
	  h_eff_correctlymatched->fill(event);
	  h_lumi_correctlymatched->fill(event);
	}
	
      }
    }

    //Final Selection
    h_finalSelection->fill(event);
    h_mlqreco_finalSelection->fill(event);
    h_jets_finalSelection->fill(event);
    h_ele_finalSelection->fill(event);
    h_mu_finalSelection->fill(event);
    h_event_finalSelection->fill(event);
    h_topjets_finalSelection->fill(event);
    h_eff_finalSelection->fill(event);
    h_lumi_finalSelection->fill(event);
    h_ht_finalSelection->fill(event);

    h_PDF_variations->fill(event);
    h_Sideband->fill(event);



    event.set(handle_evtwgt, event.weight);
    event.set(handle_nbjets, nbjets);
    event.set(handle_njets, njets);
    event.set(handle_st, st);
    event.set(handle_stlep, stlep);
    event.set(handle_mmumu, mmumu);

    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuAnalysisLumiscaleModule)
} 

