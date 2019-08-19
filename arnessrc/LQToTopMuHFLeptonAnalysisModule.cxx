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
#include "UHH2/LQToTopMu/include/LQToTopMuHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuHFLeptonHists.h"
#include "UHH2/LQToTopMu/include/LQToTopMuModules.h"
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

  class LQToTopMuHFLeptonAnalysisModule: public AnalysisModule {
  public:
    
    explicit LQToTopMuHFLeptonAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;
    
  private:
    
    unique_ptr<CommonModules> common;
    unique_ptr<AnalysisModule> SF_muonID, SF_muonIso, SF_eleReco, SF_eleID;
    unique_ptr<AnalysisModule> ttgenprod, LQgenprod;
    unique_ptr<MCMuonScaleFactor> SF_muonTrigger;
    unique_ptr<MuonTrkWeights> SF_muonTrk;
    
    // declare the Selections to use.
    unique_ptr<Selection>  nele_sel, nmu_sel;
    
    // store the Hists collection as member variables. 
    unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts, h_topjets_nocuts, h_tau_nocuts, h_eff_nocuts, h_lumi_nocuts;
    unique_ptr<Hists> h_hyphists;
    unique_ptr<Hists> h_noaddlepton, h_jets_noaddlepton, h_ele_noaddlepton, h_mu_noaddlepton, h_event_noaddlepton, h_topjets_noaddlepton, h_tau_noaddlepton, h_eff_noaddlepton, h_lumi_noaddlepton;
    unique_ptr<Hists> h_addlepton, h_jets_addlepton, h_ele_addlepton, h_mu_addlepton, h_event_addlepton, h_topjets_addlepton, h_tau_addlepton, h_eff_addlepton, h_lumi_addlepton;
    unique_ptr<Hists> h_hflepton_addlepton;

    
    MuonId MuId, MuId_Sel;
    ElectronId EleId, EleId_Sel;
    JetId Btag_loose, Btag_medium, Btag_tight, JetID;


    Event::Handle<TTbarGen> h_ttbargen;
    Event::Handle<LQGen> h_LQLQbargen;


    bool is_mc, do_pdf_variations, is_ele;
    string Sys_MuonID, Sys_MuonTrigger, Sys_MuonIso, Sys_PU, Sys_EleID, Sys_EleReco, Sys_MuonTrk;
    TString dataset_version;
 
    
    vector<unique_ptr<AnalysisModule>> recomodules, muonic_recomodules;
    uhh2::Event::Handle<vector<LQReconstructionHypothesis>> h_hyps, h_muonic_hyps;
    
  };
  
  
  LQToTopMuHFLeptonAnalysisModule::LQToTopMuHFLeptonAnalysisModule(Context & ctx){
    
    cout << "Hello World from LQToTopMuHFLeptonAnalysisModule!" << endl;

    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }

    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    dataset_version = ctx.get("dataset_version");
    is_mc = ctx.get("dataset_type") == "MC";
    is_ele = ctx.get("is_ele_channel") == "true";
    Sys_MuonID = ctx.get("Systematic_MuonID");
    Sys_MuonTrigger = ctx.get("Systematic_MuonTrigger");
    Sys_MuonIso = ctx.get("Systematic_MuonIso");
    Sys_MuonTrk = ctx.get("Systematic_MuonTrk");
    Sys_EleID = ctx.get("Systematic_EleID");
    Sys_EleReco = ctx.get("Systematic_EleReco");
    Sys_PU = ctx.get("Systematic_PU");


    // 1. setup other modules. CommonModules and the JetCleaner:
    EleId = AndId<Electron>(ElectronID_Spring16_tight_noIso, PtEtaCut(30.0, 2.4)); //IDs without any Iso statement for cleaners. 
    MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4));
    EleId_Sel = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(30.0, 2.4)); 
    MuId_Sel = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4),MuonIso(0.15));
    JetID = AndId<Jet>(PtEtaCut(30.0, 2.4), JetPFID(JetPFID::WP_LOOSE));


    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->set_jet_id(JetID);
    common->init(ctx,Sys_PU);

    SF_muonID.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "tightID", true, Sys_MuonID));
    SF_muonTrigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, Sys_MuonTrigger));
    SF_muonIso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1., "iso", true, Sys_MuonIso));
    SF_muonTrk.reset(new MuonTrkWeights(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", Sys_MuonTrk));

    SF_eleReco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", Sys_EleReco));
    SF_eleID.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Loose_ID.root", 1, "", Sys_EleID));

    
    h_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
    h_muonic_hyps = ctx.get_handle<vector<LQReconstructionHypothesis>>("HighMassMuonicLQReconstruction");
    
    if(is_mc){
      ttgenprod.reset(new TTbarGenProducer(ctx, "ttbargen", false));
      h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
      LQgenprod.reset(new LQGenProducer(ctx, "LQLQbargen", false));
      h_LQLQbargen = ctx.get_handle<LQGen>("LQLQbargen");
    }

    
    // 2. set up selections
    //In PreSel we required == 1 iso ele + mu
    //now we want to require one additional non-isolated lepton, either ele or mu
    
    //We cleaned without any iso statement. Therefore, if there is an additional ele/mu in the event, it must be non-isolated (we required ==1 iso in PreSel). The isolated one also fulfills the selection, so we have to require ==2
    nele_sel.reset(new NElectronSelection(2, 2));
    nmu_sel.reset(new NMuonSelection(2, 2));
    
    //make reconstruction hypotheses
    recomodules.emplace_back(new LQPrimaryLepton(ctx));
    recomodules.emplace_back(new HighMassLQReconstruction(ctx,LQNeutrinoReconstruction));
    recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassLQReconstruction"));
    if(is_mc) recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassLQReconstruction"));
    muonic_recomodules.emplace_back(new HighMassMuonicLQReconstruction(ctx,LQNeutrinoReconstruction));
    muonic_recomodules.emplace_back(new LQChi2Discriminator(ctx,"HighMassMuonicLQReconstruction"));
    if(is_mc) muonic_recomodules.emplace_back(new LQCorrectMatchDiscriminator(ctx,"HighMassMuonicLQReconstruction"));
    
    // 3. Set up Hists classes:
    h_nocuts.reset(new LQToTopMuHists(ctx, "NoCuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_NoCuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_NoCuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_NoCuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_NoCuts"));
    h_topjets_nocuts.reset(new TopJetHists(ctx, "TopJets_NoCuts"));
    h_eff_nocuts.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoCuts"));
    h_lumi_nocuts.reset(new LuminosityHists(ctx, "Lumi_NoCuts"));
    h_hyphists.reset(new HypothesisHistsOwn(ctx, "CorrectMatch_Hists", "HighMassLQReconstruction", "CorrectMatch"));

    h_noaddlepton.reset(new LQToTopMuHists(ctx, "NoAdditionalLepton"));
    h_jets_noaddlepton.reset(new JetHists(ctx, "Jets_NoAdditionalLepton"));
    h_ele_noaddlepton.reset(new ElectronHists(ctx, "Ele_NoAdditionalLepton"));
    h_mu_noaddlepton.reset(new MuonHists(ctx, "Mu_NoAdditionalLepton"));
    h_event_noaddlepton.reset(new EventHists(ctx, "Event_NoAdditionalLepton"));
    h_topjets_noaddlepton.reset(new TopJetHists(ctx, "TopJets_NoAdditionalLepton"));
    h_eff_noaddlepton.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_NoAdditionalLepton"));
    h_lumi_noaddlepton.reset(new LuminosityHists(ctx, "Lumi_NoAdditionalLepton"));

    h_addlepton.reset(new LQToTopMuHists(ctx, "AdditionalLepton"));
    h_hflepton_addlepton.reset(new LQToTopMuHFLeptonHists(ctx, "HFLepton_AdditionalLepton", EleId_Sel, MuId_Sel));
    h_jets_addlepton.reset(new JetHists(ctx, "Jets_AdditionalLepton"));
    h_ele_addlepton.reset(new ElectronHists(ctx, "Ele_AdditionalLepton"));
    h_mu_addlepton.reset(new MuonHists(ctx, "Mu_AdditionalLepton"));
    h_event_addlepton.reset(new EventHists(ctx, "Event_AdditionalLepton"));
    h_topjets_addlepton.reset(new TopJetHists(ctx, "TopJets_AdditionalLepton"));
    h_eff_addlepton.reset(new LQToTopMuEfficiencyHists(ctx, "Eff_AdditionalLepton"));
    h_lumi_addlepton.reset(new LuminosityHists(ctx, "Lumi_AdditionalLepton"));
    
    
  }
  
  
  bool LQToTopMuHFLeptonAnalysisModule::process(Event & event) {
    //apply muon SFs as in preselection
    SF_muonTrigger->process_onemuon(event,0);
    SF_muonID->process(event);
    SF_muonIso->process(event);
    SF_muonTrk->process(event);

    if(event.electrons->size() >= 1){
      SF_eleReco->process(event);
      SF_eleID->process(event);
    }

    if(is_mc) event.weight *= 0.9626 /* * 0.9392*/;

    bool pass_common = common->process(event);
    if(!pass_common) return false;


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
    double sum_mu_charge = 0;
    for(const auto & mu : *event.muons) sum_mu_charge += mu.charge();
    bool reconstruct_mlq_mu = (event.electrons->size() == 0 && event.muons->size() == 3 && fabs(sum_mu_charge) == 1);
    bool reconstruct_mlq_ele = (event.electrons->size() >= 1 && event.muons->size() >= 2 && charge_opposite);
    bool reconstruct_mlq = reconstruct_mlq_ele || reconstruct_mlq_mu;

    
    if(is_mc){
      ttgenprod->process(event);
      LQgenprod->process(event);
    }
    
    // MLQ reco
    if(reconstruct_mlq_ele){
      for(auto & m : recomodules){
	m->process(event);
      }
      
    }

    
    if(reconstruct_mlq_mu){
      for(auto & m : muonic_recomodules){
	m->process(event);
      }
      if(is_mc){
	std::vector<LQReconstructionHypothesis> hyps = event.get(h_muonic_hyps); 
	const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
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
    if(reconstruct_mlq_ele && is_mc) h_hyphists->fill(event);

    //if there are no additional leptons, fill hists to be able to scale normalization, which otherwise could interfere with HF-investigation
    if(event.electrons->size() == 1 && event.muons->size() == 1){
      h_noaddlepton->fill(event);
      h_jets_noaddlepton->fill(event);
      h_ele_noaddlepton->fill(event);
      h_mu_noaddlepton->fill(event);
      h_event_noaddlepton->fill(event);
      h_topjets_noaddlepton->fill(event);
      h_eff_noaddlepton->fill(event);
      h_lumi_noaddlepton->fill(event);
    }


    //require one additional lepton, but not 2 ele AND 2 muons
    if(event.electrons->size() > 2 || event.muons->size() > 2) return false;
    if(event.electrons->size() < 1 || event.muons->size() < 1) return false;
    if(!(event.electrons->size() == 2 && event.muons->size() == 1) && !(event.electrons->size() == 1 && event.muons->size() == 2)) return false;
    if(is_ele){
      //if(!(event.electrons->size() == 2 && event.muons->size() == 1)) return false;
    }
    else{
      //if(!(event.electrons->size() == 1 && event.muons->size() == 2)) return false;
    }
	
    h_addlepton->fill(event);
    h_hflepton_addlepton->fill(event);
    h_jets_addlepton->fill(event);
    h_ele_addlepton->fill(event);
    h_mu_addlepton->fill(event);
    h_event_addlepton->fill(event);
    h_topjets_addlepton->fill(event);
    h_eff_addlepton->fill(event);
    h_lumi_addlepton->fill(event);

    return true;
  }
  

  UHH2_REGISTER_ANALYSIS_MODULE(LQToTopMuHFLeptonAnalysisModule)
} 

