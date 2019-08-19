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
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/LQToTopMu/include/LQToTopMuSelections.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
#include "UHH2/LQToTopMu/include/LQToTopMuTagProbeHists.h"
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

  class LQToTopMuTagProbeAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuTagProbeAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    unique_ptr<CommonModules> common;
    unique_ptr<AnalysisModule> SF_muonID, SF_muonTrigger, SF_muonIso, SF_btag, SF_eleReco, SF_eleID;
    unique_ptr<ElectronTriggerWeights> SF_eleTrigger;
    unique_ptr<JetCleaner> jetcleaner;
    unique_ptr<ElectronCleaner> elecleaner;
    
    // declare the Selections to use.
    unique_ptr<Selection>  trigger_sel1, trigger_sel2, nele_sel_trig1, nele_sel_trig2, nele_sel_trig50;
    
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_eff_nocuts, h_tagprobe_nocuts,
      h_notfromiso, h_jets_notfromiso, h_ele_notfromiso, h_mu_notfromiso, h_event_notfromiso, h_topjets_notfromiso, h_eff_notfromiso, h_tagprobe_notfromiso,
      h_notfromiso_ele110, h_jets_notfromiso_ele110, h_ele_notfromiso_ele110, h_mu_notfromiso_ele110, h_event_notfromiso_ele110, h_topjets_notfromiso_ele110, h_eff_notfromiso_ele110, h_tagprobe_notfromiso_ele110,
      h_ele30, h_jets_ele30, h_ele_ele30, h_mu_ele30, h_event_ele30, h_topjets_ele30, h_eff_ele30, h_tagprobe_ele30,
      h_ele30to110, h_jets_ele30to110, h_ele_ele30to110, h_mu_ele30to110, h_event_ele30to110, h_topjets_ele30to110, h_eff_ele30to110, h_tagprobe_ele30to110,
      h_ele30to50, h_jets_ele30to50, h_ele_ele30to50, h_mu_ele30to50, h_event_ele30to50, h_topjets_ele30to50, h_eff_ele30to50, h_tagprobe_ele30to50,
      h_ele50to110, h_jets_ele50to110, h_ele_ele50to110, h_mu_ele50to110, h_event_ele50to110, h_topjets_ele50to110, h_eff_ele50to110, h_tagprobe_ele50to110,
      h_ele110, h_jets_ele110, h_ele_ele110, h_mu_ele110, h_event_ele110, h_topjets_ele110, h_eff_ele110, h_tagprobe_ele110,
      h_trigger, h_jets_trigger, h_ele_trigger, h_mu_trigger, h_event_trigger, h_topjets_trigger, h_eff_trigger, h_tagprobe_trigger,
      h_isotrigger, h_jets_isotrigger, h_ele_isotrigger, h_mu_isotrigger, h_event_isotrigger,h_topjets_isotrigger, h_eff_isotrigger, h_tagprobe_isotrigger,
      h_nonisotrigger, h_jets_nonisotrigger, h_ele_nonisotrigger, h_mu_nonisotrigger, h_event_nonisotrigger, h_topjets_nonisotrigger, h_eff_nonisotrigger, h_tagprobe_nonisotrigger,
      h_isotrigger_plateau, h_jets_isotrigger_plateau, h_ele_isotrigger_plateau, h_mu_isotrigger_plateau, h_event_isotrigger_plateau,h_topjets_isotrigger_plateau, h_eff_isotrigger_plateau, h_tagprobe_isotrigger_plateau,
      h_nonisotrigger_plateau, h_jets_nonisotrigger_plateau, h_ele_nonisotrigger_plateau, h_mu_nonisotrigger_plateau, h_event_nonisotrigger_plateau, h_topjets_nonisotrigger_plateau, h_eff_nonisotrigger_plateau, h_tagprobe_nonisotrigger_plateau,     
      h_30plateau, h_jets_30plateau, h_ele_30plateau, h_mu_30plateau, h_event_30plateau, h_topjets_30plateau, h_eff_30plateau, h_tagprobe_30plateau,
      h_30to110plateau, h_jets_30to110plateau, h_ele_30to110plateau, h_mu_30to110plateau, h_event_30to110plateau, h_topjets_30to110plateau, h_eff_30to110plateau, h_tagprobe_30to110plateau,
      h_30to50plateau, h_jets_30to50plateau, h_ele_30to50plateau, h_mu_30to50plateau, h_event_30to50plateau, h_topjets_30to50plateau, h_eff_30to50plateau, h_tagprobe_30to50plateau,
      h_50to110plateau, h_jets_50to110plateau, h_ele_50to110plateau, h_mu_50to110plateau, h_event_50to110plateau, h_topjets_50to110plateau, h_eff_50to110plateau, h_tagprobe_50to110plateau,
      h_110plateau, h_jets_110plateau, h_ele_110plateau, h_mu_110plateau, h_event_110plateau, h_topjets_110plateau, h_eff_110plateau, h_tagprobe_110plateau,
      h_plateau, h_jets_plateau, h_ele_plateau, h_mu_plateau, h_event_plateau, h_topjets_plateau, h_eff_plateau, h_tagprobe_plateau;


    
    MuonId MuId;
    ElectronId EleId, EleId_trig1, EleId_trig2, EleId_trig50;


    bool is_mc, closuretest;
    string Sys_MuonID, Sys_MuonTrigger, Sys_MuonIso, Sys_PU, Sys_EleID, Sys_EleReco, Sys_EleTrigger;
 
    
    vector<unique_ptr<AnalysisModule>> recomodules;
    uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps;
    
  };
  
  
  LQToTopMuTagProbeAnalysisModule::LQToTopMuTagProbeAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuTagProbeAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
    is_mc = ctx.get("dataset_type") == "MC";
    closuretest = ctx.get("Closuretest") == "true" || ctx.get("Closuretest") == "True";
    Sys_MuonID = ctx.get("Systematic_MuonID");
    Sys_MuonTrigger = ctx.get("Systematic_MuonTrigger");
    Sys_MuonIso = ctx.get("Systematic_MuonIso");
    Sys_EleID = ctx.get("Systematic_EleID");
    Sys_EleReco = ctx.get("Systematic_EleReco");
    Sys_EleTrigger = ctx.get("Systematic_EleTrigger");
    Sys_PU = ctx.get("Systematic_PU");


    // 1. setup other modules. CommonModules and the JetCleaner:
    MuId = AndId<Muon>(MuonIDTight(),PtEtaCut(30.0, 2.4),MuonIso(0.15));
    EleId = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(10.0, 2.4));
    EleId_trig1 = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(30.0, 2.4)); //30
    EleId_trig2 = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(120.0, 2.4)); 
    EleId_trig50 = AndId<Electron>(ElectronID_Spring16_tight,PtEtaCut(50.0, 2.4)); 
    jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4));


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

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", Sys_EleReco));
    SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", Sys_EleID));
    SF_eleTrigger.reset(new ElectronTriggerWeights(ctx, "/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/TagProbe/Optimization/35867fb_Iso27_NonIso115/ElectronEfficiencies.root", Sys_EleTrigger));
    
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");

    
    // 2. set up selections
    //Selection
    trigger_sel1.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    trigger_sel2.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")); //105
    nele_sel_trig1.reset(new NElectronSelection(1, 1, EleId_trig1));
    nele_sel_trig50.reset(new NElectronSelection(1, 1, EleId_trig50));
    nele_sel_trig2.reset(new NElectronSelection(1, 1, EleId_trig2));
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
    h_eff_nocuts.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoCuts"));    
    h_tagprobe_nocuts.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_NoCuts"));  

    h_notfromiso.reset(new LQToTopMuHists(ctx, "NotfromIso"));
    h_jets_notfromiso.reset(new JetHists(ctx, "Jets_NotfromIso"));
    h_ele_notfromiso.reset(new ElectronHists(ctx, "Ele_NotfromIso"));
    h_mu_notfromiso.reset(new MuonHists(ctx, "Mu_NotfromIso"));
    h_event_notfromiso.reset(new EventHists(ctx, "Event_NotfromIso"));
    h_topjets_notfromiso.reset(new TopJetHists(ctx, "TopJets_NotfromIso"));
    h_eff_notfromiso.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NotfromIso"));    
    h_tagprobe_notfromiso.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_NotfromIso"));  

    h_notfromiso_ele110.reset(new LQToTopMuHists(ctx, "NotfromIso_Ele110"));
    h_jets_notfromiso_ele110.reset(new JetHists(ctx, "Jets_NotfromIso_Ele110"));
    h_ele_notfromiso_ele110.reset(new ElectronHists(ctx, "Ele_NotfromIso_Ele110"));
    h_mu_notfromiso_ele110.reset(new MuonHists(ctx, "Mu_NotfromIso_Ele110"));
    h_event_notfromiso_ele110.reset(new EventHists(ctx, "Event_NotfromIso_Ele110"));
    h_topjets_notfromiso_ele110.reset(new TopJetHists(ctx, "TopJets_NotfromIso_Ele110"));
    h_eff_notfromiso_ele110.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NotfromIso_Ele110"));    
    h_tagprobe_notfromiso_ele110.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_NotfromIso_Ele110"));  

    h_ele30.reset(new LQToTopMuHists(ctx, "Ele30"));
    h_jets_ele30.reset(new JetHists(ctx, "Jets_Ele30"));
    h_ele_ele30.reset(new ElectronHists(ctx, "Ele_Ele30"));
    h_mu_ele30.reset(new MuonHists(ctx, "Mu_Ele30"));
    h_event_ele30.reset(new EventHists(ctx, "Event_Ele30"));
    h_topjets_ele30.reset(new TopJetHists(ctx, "TopJets_Ele30"));
    h_eff_ele30.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Ele30"));    
    h_tagprobe_ele30.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Ele30"));  

    h_ele30to110.reset(new LQToTopMuHists(ctx, "Ele30to110"));
    h_jets_ele30to110.reset(new JetHists(ctx, "Jets_Ele30to110"));
    h_ele_ele30to110.reset(new ElectronHists(ctx, "Ele_Ele30to110"));
    h_mu_ele30to110.reset(new MuonHists(ctx, "Mu_Ele30to110"));
    h_event_ele30to110.reset(new EventHists(ctx, "Event_Ele30to110"));
    h_topjets_ele30to110.reset(new TopJetHists(ctx, "TopJets_Ele30to110"));
    h_eff_ele30to110.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Ele30to110"));    
    h_tagprobe_ele30to110.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Ele30to110"));  

    h_ele30to50.reset(new LQToTopMuHists(ctx, "Ele30to50"));
    h_jets_ele30to50.reset(new JetHists(ctx, "Jets_Ele30to50"));
    h_ele_ele30to50.reset(new ElectronHists(ctx, "Ele_Ele30to50"));
    h_mu_ele30to50.reset(new MuonHists(ctx, "Mu_Ele30to50"));
    h_event_ele30to50.reset(new EventHists(ctx, "Event_Ele30to50"));
    h_topjets_ele30to50.reset(new TopJetHists(ctx, "TopJets_Ele30to50"));
    h_eff_ele30to50.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Ele30to50"));    
    h_tagprobe_ele30to50.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Ele30to50")); 

    h_ele50to110.reset(new LQToTopMuHists(ctx, "Ele50to110"));
    h_jets_ele50to110.reset(new JetHists(ctx, "Jets_Ele50to110"));
    h_ele_ele50to110.reset(new ElectronHists(ctx, "Ele_Ele50to110"));
    h_mu_ele50to110.reset(new MuonHists(ctx, "Mu_Ele50to110"));
    h_event_ele50to110.reset(new EventHists(ctx, "Event_Ele50to110"));
    h_topjets_ele50to110.reset(new TopJetHists(ctx, "TopJets_Ele50to110"));
    h_eff_ele50to110.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Ele50to110"));    
    h_tagprobe_ele50to110.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Ele50to110"));  

    h_ele110.reset(new LQToTopMuHists(ctx, "Ele110"));
    h_jets_ele110.reset(new JetHists(ctx, "Jets_Ele110"));
    h_ele_ele110.reset(new ElectronHists(ctx, "Ele_Ele110"));
    h_mu_ele110.reset(new MuonHists(ctx, "Mu_Ele110"));
    h_event_ele110.reset(new EventHists(ctx, "Event_Ele110"));
    h_topjets_ele110.reset(new TopJetHists(ctx, "TopJets_Ele110"));
    h_eff_ele110.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Ele110"));    
    h_tagprobe_ele110.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Ele110"));  

    h_trigger.reset(new LQToTopMuHists(ctx, "Trigger"));
    h_jets_trigger.reset(new JetHists(ctx, "Jets_Trigger"));
    h_ele_trigger.reset(new ElectronHists(ctx, "Ele_Trigger"));
    h_mu_trigger.reset(new MuonHists(ctx, "Mu_Trigger"));
    h_event_trigger.reset(new EventHists(ctx, "Event_Trigger"));
    h_topjets_trigger.reset(new TopJetHists(ctx, "TopJets_Trigger"));
    h_eff_trigger.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Trigger")); 
    h_tagprobe_trigger.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Trigger")); 

    h_isotrigger.reset(new LQToTopMuHists(ctx, "IsoTrigger"));
    h_jets_isotrigger.reset(new JetHists(ctx, "Jets_IsoTrigger"));
    h_ele_isotrigger.reset(new ElectronHists(ctx, "Ele_IsoTrigger"));
    h_mu_isotrigger.reset(new MuonHists(ctx, "Mu_IsoTrigger"));
    h_event_isotrigger.reset(new EventHists(ctx, "Event_IsoTrigger"));
    h_topjets_isotrigger.reset(new TopJetHists(ctx, "TopJets_IsoTrigger"));
    h_eff_isotrigger.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_IsoTrigger"));    
    h_tagprobe_isotrigger.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_IsoTrigger"));  

    h_nonisotrigger.reset(new LQToTopMuHists(ctx, "NonIsoTrigger"));
    h_jets_nonisotrigger.reset(new JetHists(ctx, "Jets_NonIsoTrigger"));
    h_ele_nonisotrigger.reset(new ElectronHists(ctx, "Ele_NonIsoTrigger"));
    h_mu_nonisotrigger.reset(new MuonHists(ctx, "Mu_NonIsoTrigger"));
    h_event_nonisotrigger.reset(new EventHists(ctx, "Event_NonIsoTrigger"));
    h_topjets_nonisotrigger.reset(new TopJetHists(ctx, "TopJets_NonIsoTrigger"));
    h_eff_nonisotrigger.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NonIsoTrigger"));    
    h_tagprobe_nonisotrigger.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_NonIsoTrigger"));  

    h_isotrigger_plateau.reset(new LQToTopMuHists(ctx, "IsoTrigger_Plateau"));
    h_jets_isotrigger_plateau.reset(new JetHists(ctx, "Jets_IsoTrigger_Plateau"));
    h_ele_isotrigger_plateau.reset(new ElectronHists(ctx, "Ele_IsoTrigger_Plateau"));
    h_mu_isotrigger_plateau.reset(new MuonHists(ctx, "Mu_IsoTrigger_Plateau"));
    h_event_isotrigger_plateau.reset(new EventHists(ctx, "Event_IsoTrigger_Plateau"));
    h_topjets_isotrigger_plateau.reset(new TopJetHists(ctx, "TopJets_IsoTrigger_Plateau"));
    h_eff_isotrigger_plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_IsoTrigger_Plateau"));    
    h_tagprobe_isotrigger_plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_IsoTrigger_Plateau"));  

    h_nonisotrigger_plateau.reset(new LQToTopMuHists(ctx, "NonIsoTrigger_Plateau"));
    h_jets_nonisotrigger_plateau.reset(new JetHists(ctx, "Jets_NonIsoTrigger_Plateau"));
    h_ele_nonisotrigger_plateau.reset(new ElectronHists(ctx, "Ele_NonIsoTrigger_Plateau"));
    h_mu_nonisotrigger_plateau.reset(new MuonHists(ctx, "Mu_NonIsoTrigger_Plateau"));
    h_event_nonisotrigger_plateau.reset(new EventHists(ctx, "Event_NonIsoTrigger_Plateau"));
    h_topjets_nonisotrigger_plateau.reset(new TopJetHists(ctx, "TopJets_NonIsoTrigger_Plateau"));
    h_eff_nonisotrigger_plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NonIsoTrigger_Plateau"));    
    h_tagprobe_nonisotrigger_plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_NonIsoTrigger_Plateau")); 
 
    h_30plateau.reset(new LQToTopMuHists(ctx, "30Plateau"));
    h_jets_30plateau.reset(new JetHists(ctx, "Jets_30Plateau"));
    h_ele_30plateau.reset(new ElectronHists(ctx, "Ele_30Plateau"));
    h_mu_30plateau.reset(new MuonHists(ctx, "Mu_30Plateau"));
    h_event_30plateau.reset(new EventHists(ctx, "Event_30Plateau"));
    h_topjets_30plateau.reset(new TopJetHists(ctx, "TopJets_30Plateau"));
    h_eff_30plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_30Plateau"));    
    h_tagprobe_30plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_30Plateau"));  

    h_30to110plateau.reset(new LQToTopMuHists(ctx, "30to110Plateau"));
    h_jets_30to110plateau.reset(new JetHists(ctx, "Jets_30to110Plateau"));
    h_ele_30to110plateau.reset(new ElectronHists(ctx, "Ele_30to110Plateau"));
    h_mu_30to110plateau.reset(new MuonHists(ctx, "Mu_30to110Plateau"));
    h_event_30to110plateau.reset(new EventHists(ctx, "Event_30to110Plateau"));
    h_topjets_30to110plateau.reset(new TopJetHists(ctx, "TopJets_30to110Plateau"));
    h_eff_30to110plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_30to110Plateau")); 
    h_tagprobe_30to110plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_30to110Plateau")); 

    h_30to50plateau.reset(new LQToTopMuHists(ctx, "30to50Plateau"));
    h_jets_30to50plateau.reset(new JetHists(ctx, "Jets_30to50Plateau"));
    h_ele_30to50plateau.reset(new ElectronHists(ctx, "Ele_30to50Plateau"));
    h_mu_30to50plateau.reset(new MuonHists(ctx, "Mu_30to50Plateau"));
    h_event_30to50plateau.reset(new EventHists(ctx, "Event_30to50Plateau"));
    h_topjets_30to50plateau.reset(new TopJetHists(ctx, "TopJets_30to50Plateau"));
    h_eff_30to50plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_30to50Plateau")); 
    h_tagprobe_30to50plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_30to50Plateau")); 

    h_50to110plateau.reset(new LQToTopMuHists(ctx, "50to110Plateau"));
    h_jets_50to110plateau.reset(new JetHists(ctx, "Jets_50to110Plateau"));
    h_ele_50to110plateau.reset(new ElectronHists(ctx, "Ele_50to110Plateau"));
    h_mu_50to110plateau.reset(new MuonHists(ctx, "Mu_50to110Plateau"));
    h_event_50to110plateau.reset(new EventHists(ctx, "Event_50to110Plateau"));
    h_topjets_50to110plateau.reset(new TopJetHists(ctx, "TopJets_50to110Plateau"));
    h_eff_50to110plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_50to110Plateau")); 
    h_tagprobe_50to110plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_50to110Plateau")); 

    h_110plateau.reset(new LQToTopMuHists(ctx, "110Plateau"));
    h_jets_110plateau.reset(new JetHists(ctx, "Jets_110Plateau"));
    h_ele_110plateau.reset(new ElectronHists(ctx, "Ele_110Plateau"));
    h_mu_110plateau.reset(new MuonHists(ctx, "Mu_110Plateau"));
    h_event_110plateau.reset(new EventHists(ctx, "Event_110Plateau"));
    h_topjets_110plateau.reset(new TopJetHists(ctx, "TopJets_110Plateau"));
    h_eff_110plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_110Plateau")); 
    h_tagprobe_110plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_110Plateau")); 

    h_plateau.reset(new LQToTopMuHists(ctx, "Plateau"));
    h_jets_plateau.reset(new JetHists(ctx, "Jets_Plateau"));
    h_ele_plateau.reset(new ElectronHists(ctx, "Ele_Plateau"));
    h_mu_plateau.reset(new MuonHists(ctx, "Mu_Plateau"));
    h_event_plateau.reset(new EventHists(ctx, "Event_Plateau"));
    h_topjets_plateau.reset(new TopJetHists(ctx, "TopJets_Plateau"));
    h_eff_plateau.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_Plateau")); 
    h_tagprobe_plateau.reset(new LQToTopMuTagProbeHists(ctx, "TagProbe_Plateau")); 
    
  }
  
  
  bool LQToTopMuTagProbeAnalysisModule::process(Event & event) {
    //cout << endl << endl << "+++NEW EVENT+++" << endl;
    
    //apply muon SFs as in preselection
    if(event.muons->size() >= 1){
      SF_muonTrigger->process(event);
      SF_muonID->process(event);
      SF_muonIso->process(event);
    }

    if(event.electrons->size() >= 1){
      SF_eleReco->process(event);
      SF_eleID->process(event);
    }





    bool pass_common = common->process(event);
    if(!pass_common) return false;
    jetcleaner->process(event);

    if(event.electrons->size() >= 1 && event.muons->size() >= 2){
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
    h_tagprobe_nocuts->fill(event);

    if(!trigger_sel1->passes(event)){
      h_notfromiso->fill(event);
      h_jets_notfromiso->fill(event);
      h_ele_notfromiso->fill(event);
      h_mu_notfromiso->fill(event);
      h_event_notfromiso->fill(event);
      h_topjets_notfromiso->fill(event);
      h_eff_notfromiso->fill(event);
      h_tagprobe_notfromiso->fill(event);
    }

    if(!trigger_sel1->passes(event) && nele_sel_trig2->passes(event)){
      h_notfromiso_ele110->fill(event);
      h_jets_notfromiso_ele110->fill(event);
      h_ele_notfromiso_ele110->fill(event);
      h_mu_notfromiso_ele110->fill(event);
      h_event_notfromiso_ele110->fill(event);
      h_topjets_notfromiso_ele110->fill(event);
      h_eff_notfromiso_ele110->fill(event);
      h_tagprobe_notfromiso_ele110->fill(event);
    }

    if(nele_sel_trig1->passes(event)){
      h_ele30->fill(event);
      h_jets_ele30->fill(event);
      h_ele_ele30->fill(event);
      h_mu_ele30->fill(event);
      h_event_ele30->fill(event);
      h_topjets_ele30->fill(event);
      h_eff_ele30->fill(event);
      h_tagprobe_ele30->fill(event);
    }

    if(nele_sel_trig1->passes(event) && !nele_sel_trig2->passes(event)){
      h_ele30to110->fill(event);
      h_jets_ele30to110->fill(event);
      h_ele_ele30to110->fill(event);
      h_mu_ele30to110->fill(event);
      h_event_ele30to110->fill(event);
      h_topjets_ele30to110->fill(event);
      h_eff_ele30to110->fill(event);
      h_tagprobe_ele30to110->fill(event);
    }

    if(nele_sel_trig1->passes(event) && !nele_sel_trig50->passes(event)){
      h_ele30to50->fill(event);
      h_jets_ele30to50->fill(event);
      h_ele_ele30to50->fill(event);
      h_mu_ele30to50->fill(event);
      h_event_ele30to50->fill(event);
      h_topjets_ele30to50->fill(event);
      h_eff_ele30to50->fill(event);
      h_tagprobe_ele30to50->fill(event);
    }

    if(nele_sel_trig50->passes(event) && !nele_sel_trig2->passes(event)){
      h_ele50to110->fill(event);
      h_jets_ele50to110->fill(event);
      h_ele_ele50to110->fill(event);
      h_mu_ele50to110->fill(event);
      h_event_ele50to110->fill(event);
      h_topjets_ele50to110->fill(event);
      h_eff_ele50to110->fill(event);
      h_tagprobe_ele50to110->fill(event);
    }

    if(nele_sel_trig2->passes(event)){
      h_ele110->fill(event);
      h_jets_ele110->fill(event);
      h_ele_ele110->fill(event);
      h_mu_ele110->fill(event);
      h_event_ele110->fill(event);
      h_topjets_ele110->fill(event);
      h_eff_ele110->fill(event);
      h_tagprobe_ele110->fill(event);
    }

    if(!(trigger_sel1->passes(event) || trigger_sel2->passes(event))) return false;

    if(closuretest){
      if(event.electrons->at(0).pt() < 30) return false;      
      SF_eleTrigger->process(event);
    }

    if(event.electrons->at(0).pt() < 30) cout << "Hello, this ele has pt < 30" << endl;



    h_trigger->fill(event);
    h_jets_trigger->fill(event);
    h_ele_trigger->fill(event);
    h_mu_trigger->fill(event);
    h_event_trigger->fill(event);
    h_topjets_trigger->fill(event);
    h_eff_trigger->fill(event);
    h_tagprobe_trigger->fill(event);

    if(trigger_sel1->passes(event)){
      h_isotrigger->fill(event);
      h_jets_isotrigger->fill(event);
      h_ele_isotrigger->fill(event);
      h_mu_isotrigger->fill(event);
      h_event_isotrigger->fill(event);
      h_topjets_isotrigger->fill(event);
      h_eff_isotrigger->fill(event);
      h_tagprobe_isotrigger->fill(event);
    }

    if(!trigger_sel1->passes(event) && trigger_sel2->passes(event)){
      h_nonisotrigger->fill(event);
      h_jets_nonisotrigger->fill(event);
      h_ele_nonisotrigger->fill(event);
      h_mu_nonisotrigger->fill(event);
      h_event_nonisotrigger->fill(event);
      h_topjets_nonisotrigger->fill(event);
      h_eff_nonisotrigger->fill(event);
      h_tagprobe_nonisotrigger->fill(event);
    }

    if(trigger_sel1->passes(event) && nele_sel_trig1->passes(event)){
      h_isotrigger_plateau->fill(event);
      h_jets_isotrigger_plateau->fill(event);
      h_ele_isotrigger_plateau->fill(event);
      h_mu_isotrigger_plateau->fill(event);
      h_event_isotrigger_plateau->fill(event);
      h_topjets_isotrigger_plateau->fill(event);
      h_eff_isotrigger_plateau->fill(event);
      h_tagprobe_isotrigger_plateau->fill(event);
    }

    if(!trigger_sel1->passes(event) && trigger_sel2->passes(event) && nele_sel_trig2->passes(event)){
      h_nonisotrigger_plateau->fill(event);
      h_jets_nonisotrigger_plateau->fill(event);
      h_ele_nonisotrigger_plateau->fill(event);
      h_mu_nonisotrigger_plateau->fill(event);
      h_event_nonisotrigger_plateau->fill(event);
      h_topjets_nonisotrigger_plateau->fill(event);
      h_eff_nonisotrigger_plateau->fill(event);
      h_tagprobe_nonisotrigger_plateau->fill(event);
    }

    if(nele_sel_trig1->passes(event)){
      h_30plateau->fill(event);
      h_jets_30plateau->fill(event);
      h_ele_30plateau->fill(event);
      h_mu_30plateau->fill(event);
      h_event_30plateau->fill(event);
      h_topjets_30plateau->fill(event);
      h_eff_30plateau->fill(event);
      h_tagprobe_30plateau->fill(event);
    }

    if(nele_sel_trig1->passes(event) && !nele_sel_trig2->passes(event)){
      h_30to110plateau->fill(event);
      h_jets_30to110plateau->fill(event);
      h_ele_30to110plateau->fill(event);
      h_mu_30to110plateau->fill(event);
      h_event_30to110plateau->fill(event);
      h_topjets_30to110plateau->fill(event);
      h_eff_30to110plateau->fill(event);
      h_tagprobe_30to110plateau->fill(event);
    }

    if(nele_sel_trig1->passes(event) && !nele_sel_trig50->passes(event)){
      h_30to50plateau->fill(event);
      h_jets_30to50plateau->fill(event);
      h_ele_30to50plateau->fill(event);
      h_mu_30to50plateau->fill(event);
      h_event_30to50plateau->fill(event);
      h_topjets_30to50plateau->fill(event);
      h_eff_30to50plateau->fill(event);
      h_tagprobe_30to50plateau->fill(event);
    }

    if(nele_sel_trig50->passes(event) && !nele_sel_trig2->passes(event)){
      h_50to110plateau->fill(event);
      h_jets_50to110plateau->fill(event);
      h_ele_50to110plateau->fill(event);
      h_mu_50to110plateau->fill(event);
      h_event_50to110plateau->fill(event);
      h_topjets_50to110plateau->fill(event);
      h_eff_50to110plateau->fill(event);
      h_tagprobe_50to110plateau->fill(event);
    }

    if(nele_sel_trig2->passes(event)){
      h_110plateau->fill(event);
      h_jets_110plateau->fill(event);
      h_ele_110plateau->fill(event);
      h_mu_110plateau->fill(event);
      h_event_110plateau->fill(event);
      h_topjets_110plateau->fill(event);
      h_eff_110plateau->fill(event);
      h_tagprobe_110plateau->fill(event);
    }

    if(!((trigger_sel1->passes(event) && nele_sel_trig1->passes(event)) || (trigger_sel2->passes(event) && nele_sel_trig2->passes(event))) ) return false;
    h_plateau->fill(event);
    h_jets_plateau->fill(event);
    h_ele_plateau->fill(event);
    h_mu_plateau->fill(event);
    h_event_plateau->fill(event);
    h_topjets_plateau->fill(event);
    h_eff_plateau->fill(event);
    h_tagprobe_plateau->fill(event);


    return false;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuTagProbeAnalysisModule)
} 

