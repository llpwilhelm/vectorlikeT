#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/NSelections.h"
//#include "UHH2/vectorlikeT/include/vectorlikeTSelections.h"
//#include "UHH2/vectorlikeT/include/vectorlikeTHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/EventHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class vectorlikeTtestpreselectionModule: public AnalysisModule {
public:
    
  explicit vectorlikeTtestpreselectionModule(Context & ctx);
  virtual bool process(Event & event) override;

private:
    
  std::unique_ptr<CommonModules> common;
    
  //std::unique_ptr<JetCleaner> jetcleaner;
   
  // declare the Selections to use. 
  //Use unique_ptr to ensure automatic call of delete in the destructor,
  // to avoid memory leaks.
    unique_ptr<Selection> jet_sel, muon_sel, ele_sel;

  
  // store the Hists collection as member variables.
  //Again, use unique_ptr to avoid memory leaks.
  // std::unique_ptr<Hists> h_nocuts, h_njet, h_ele;   
  unique_ptr<Hists> h_nocuts, h_jets_nocuts, h_ele_nocuts, h_mu_nocuts, h_event_nocuts;
  unique_ptr<Hists> h_cleaner, h_jets_cleaner, h_ele_cleaner, h_mu_cleaner, h_event_cleaner;
    unique_ptr<Hists> h_muon, h_jets_muon, h_ele_muon, h_mu_muon, h_event_muon;
    unique_ptr<Hists> h_electron, h_jets_electron, h_ele_electron, h_mu_electron, h_event_electron;
    unique_ptr<Hists> h_jets, h_jets_jets, h_ele_jets, h_mu_jets, h_event_jets;
    unique_ptr<Hists> h_met, h_jets_met, h_ele_met, h_mu_met, h_event_met;

  MuonId MuId;
  ElectronId EleId;
  JetId Jet_ID;
  bool lep_is_mu;
  bool is_mc;

};


vectorlikeTtestpreselectionModule::vectorlikeTtestpreselectionModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello folks, from vectorlikeTtestpreselectionModule!" << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    EleId = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(30.0, 2.4));
    MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4), MuonIso(0.15));
    Jet_ID = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0, 2.4));
    
    is_mc = ctx.get("dataset_type") == "MC";
    
    // 1. setup other modules. CommonModules and the JetCleaner:
    //common.reset(new CommonModules());
    // TODO: configure common here, e.g. by 
    // calling common->set_*_id or common->disable_*   
    
    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->set_jet_id(Jet_ID);
    common->switch_jetPtSorter();
    common->init(ctx);
  
    // the cleaning can also be achieved with less code via CommonModules with:
    // common->set_jet_id(PtEtaCut(30.0, 2.4));
    // before the 'common->init(ctx)' line.
    
    // 2. set up selections
    // njet_sel.reset(new NJetSelection(2)); // see common/include/NSelections.h
    // dijet_sel.reset(new DijetSelection()); // see vectorlikeTSelections

    
    lep_is_mu = true;
    
    muon_sel.reset(new NMuonSelection(1, 1));
    ele_sel.reset(new NElectronSelection(1, 1));
    
    
    jet_sel.reset(new NJetSelection(3, -1));


    // 3. Set up Hists classes:
    //h_nocuts.reset(new vectorlikeTHists(ctx, "NoCuts"));
    //h_njet.reset(new vectorlikeTHists(ctx, "Njet"));
    //h_dijet.reset(new vectorlikeTHists(ctx, "Dijet"));
    // h_ele.reset(new ElectronHists(ctx, "ele_nocuts"));   

    //h_nocuts.reset(new vectorlikeTHists(ctx, "nocuts"));
    h_jets_nocuts.reset(new JetHists(ctx, "Jets_nocuts"));
    h_ele_nocuts.reset(new ElectronHists(ctx, "Ele_nocuts"));
    h_mu_nocuts.reset(new MuonHists(ctx, "Mu_nocuts"));
    h_event_nocuts.reset(new EventHists(ctx, "Event_nocuts"));

    // h_cleaner.reset(new vectorlikeTHists(ctx, "cleaner"));
    h_jets_cleaner.reset(new JetHists(ctx, "Jets_cleaner"));
    h_ele_cleaner.reset(new ElectronHists(ctx, "Ele_cleaner"));
    h_mu_cleaner.reset(new MuonHists(ctx, "Mu_cleaner"));
    h_event_cleaner.reset(new EventHists(ctx, "Event_cleaner"));

    // h_muon.reset(new vectorlikeTHists(ctx, "muon")); 
    h_jets_muon.reset(new JetHists(ctx, "Jets_muon"));
    h_ele_muon.reset(new ElectronHists(ctx, "Ele_muon"));
    h_mu_muon.reset(new MuonHists(ctx, "Mu_muon"));
    h_event_muon.reset(new EventHists(ctx, "Event_muon"));

    // h_electron.reset(new vectorlikeTHists(ctx, "electron"));
    h_jets_electron.reset(new JetHists(ctx, "Jets_electron"));
    h_ele_electron.reset(new ElectronHists(ctx, "Ele_electron"));
    h_mu_electron.reset(new MuonHists(ctx, "Mu_electron"));
    h_event_electron.reset(new EventHists(ctx, "Event_electron"));

    // h_jets.reset(new vectorlikeTHists(ctx, "jets"));
    h_jets_jets.reset(new JetHists(ctx, "Jets_jets"));
    h_ele_jets.reset(new ElectronHists(ctx, "Ele_jets"));
    h_mu_jets.reset(new MuonHists(ctx, "Mu_jets"));
    h_event_jets.reset(new EventHists(ctx, "Event_jets"));


}


bool vectorlikeTtestpreselectionModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
  // cout << "vectorlikeTtestpreselectionModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
    // 1. run all modules other modules.
    //common->process(event);
    //jetcleaner->process(event);
    
    // 2. test selections and fill histograms
    // h_ele->fill(event);
    //h_nocuts->fill(event);
    
    //bool njet_selection = njet_sel->passes(event);
    //if(njet_selection){
  //  h_njet->fill(event);
	//}
    //bool dijet_selection = dijet_sel->passes(event);
    //if(dijet_selection){
    //  h_dijet->fill(event);
	// }
    // 3. decide whether or not to keep the current event in the output:
    // return njet_selection && dijet_selection; h_nocuts->fill(event);

  // here hist fill as stolen from SingleVLQStephaniePreselectionModule
  // h_nocuts->fill(event);
  h_jets_nocuts->fill(event);
  h_ele_nocuts->fill(event);
  h_mu_nocuts->fill(event);
  h_event_nocuts->fill(event);

  bool pass_common = common->process(event);
  if(!pass_common) return false;

  // h_cleaner->fill(event);
  h_jets_cleaner->fill(event);
  h_ele_cleaner->fill(event);
  h_mu_cleaner->fill(event);
  h_event_cleaner->fill(event);
  
  
    if(!((muon_sel->passes(event)) or (ele_sel->passes(event)))) return false;
    if(muon_sel->passes(event)){
    // h_muon->fill(event);
    h_jets_muon->fill(event);
    h_ele_muon->fill(event);
    h_mu_muon->fill(event);
    h_event_muon->fill(event);
    }
    else{
    //lep_is_mu = false;
    // h_electron->fill(event);
    h_jets_electron->fill(event);
    h_ele_electron->fill(event);
    h_mu_electron->fill(event);
    h_event_electron->fill(event);
    }
  
/*  if(!lep_is_mu){
    if(!(muon_sel->passes(event))) return false;
    // h_muon->fill(event);
    h_jets_muon->fill(event);
    h_ele_muon->fill(event);
    h_mu_muon->fill(event);
    h_event_muon->fill(event);
  
    if(!(ele_sel->passes(event))) return false;
    // h_electron->fill(event);
    h_jets_electron->fill(event);
    h_ele_electron->fill(event);
    h_mu_electron->fill(event);
    h_event_electron->fill(event);
  } */  
      
  if(!jet_sel->passes(event)) return false;
  // h_jets->fill(event);
  h_jets_jets->fill(event);
  h_ele_jets->fill(event);
  h_mu_jets->fill(event);
  h_event_jets->fill(event);


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the vectorlikeTvectorlikeTtestpreselectionModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(vectorlikeTtestpreselectionModule)

}
