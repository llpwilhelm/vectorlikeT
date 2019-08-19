#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/JetIds.h"
#include <TFile.h>
#include <TGraphAsymmErrors.h>
//#include "LHAPDF/LHAPDF.h"
#include "TSystem.h"


class WeightsTo14TeV : public uhh2::AnalysisModule {
public:

  explicit WeightsTo14TeV(uhh2::Context & ctx, TString pdfname = "NNPDF30_lo_as_0130");
  virtual bool process(uhh2::Event & event) override {double dummy=event.weight;dummy+=1;throw std::runtime_error("In LQToTopMuModules.cxx: WeightsTo14TeV::process does not do anything. Use ::calculateWeight(Event) instead.");return false;};
  double calculateWeight(uhh2::Event & event);
  
 private:

  bool m_libvalid;
  double xmin, xmax, qmin, qmax;

  //std::vector<LHAPDF::PDF*> m_pdfs;
  LHAPDF::PDF* pdf;
  //Handles related to 13->14TeV scaling
  uhh2::Event::Handle<double> h_x13_1, h_x13_2, h_x14_1, h_x14_2, h_Q, h_xf13_1, h_xf13_2, h_xf14_1, h_xf14_2, h_weight13, h_weight14, h_sf;
  uhh2::Event::Handle<int> h_f1, h_f2;
  
};

class JetLeptonOverlapCleaner: public uhh2::AnalysisModule {

 public:
  explicit JetLeptonOverlapCleaner(double RJet = 0.4);
  virtual bool process(uhh2::Event & event) override;

 private:
  double RJet;

};

class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> Eff_lowpt_MC, Eff_lowpt_DATA, Eff_highpt_MC, Eff_highpt_DATA;

};

class JetCorrectorVariable: public JetCorrector{

 public: 
  explicit JetCorrectorVariable(uhh2::Context & ctx, const std::vector<std::string> & JEC_files);
  bool correct_collection(uhh2::Event & event, std::vector<Jet> & jets);

  
};
/*
class JetSmearerVariable: public GenericJetResolutionSmearer{

 public: 
  explicit JetSmearerVariable(uhh2::Context & ctx, const std::string label_jets, const std::string label_genjets, const JERSmearing::SFtype1& JER_sf=JERSmearing::SF_13TeV_2016);
  bool process(uhh2::Event & event);

  
};*/

/*
class JetSmearerVariable: public JetResolutionSmearer{

 public: 
  explicit JetSmearerVariable(uhh2::Context & ctx, const JERSmearing::SFtype1& JER_sf=JERSmearing::SF_13TeV_2016);
  bool smear_collection(uhh2::Event & event, std::vector<Jet> & jets, std::vector<Particle> & genjets);

  
};*/

class ElectronFakeRateWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronFakeRateWeights(uhh2::Context & ctx, const std::vector<std::string> & JEC_files, TString path_, TString SysDirection_, const std::string label_jets, const std::string label_genjets);
  virtual bool process(uhh2::Event & event) override;

 protected:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> SF;
  std::vector<double> x_low, x_high;
  int n_points;
  std::unique_ptr<JetCorrectorVariable> jet_corrector;
  std::unique_ptr<GenericJetResolutionSmearer> jet_smearer;
  JetId jet_id;
  uhh2::Event::Handle<double> FakeRateWeightEle;
  uhh2::Event::Handle<double> FakeRateWeightEleUp;
  uhh2::Event::Handle<double> FakeRateWeightEleDown;
  uhh2::Event::Handle<std::vector<Jet>> h_jets;

};

class MuonFakeRateWeights: public uhh2::AnalysisModule{

 public:
  explicit MuonFakeRateWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 protected:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> SF;
  uhh2::Event::Handle<double> FakeRateWeightMu;
  uhh2::Event::Handle<double> FakeRateWeightMuUp;
  uhh2::Event::Handle<double> FakeRateWeightMuDown;

};

class ZEEFinder: public uhh2::AnalysisModule{

 public:
  explicit ZEEFinder();
  virtual bool process(uhh2::Event & event) override {double dummy=event.weight;dummy+=1;throw std::runtime_error("In LQToTopMuModules.cxx: ZEEFinder::process does not do anything. Use ::search(Event) instead.");return false;};
  std::pair<int,int> search(uhh2::Event & event);
  std::pair<int,int> search(const uhh2::Event & event);

};

class ElectronJetOverlapCleaner: public uhh2::AnalysisModule{

 public:
  explicit ElectronJetOverlapCleaner();
  bool process(uhh2::Event & event) override{double dummy=event.weight;dummy+=1;throw std::runtime_error("In LQToTopMuModules.cxx: ElectronJetOverlapCleaner::process does not do anything. Use ::process(Event,int,int) instead.");return false;};
  bool process(uhh2::Event & event, int idx1, int idx2);

};

class DibosonScaleFactors: public uhh2::AnalysisModule{

 public:
  explicit DibosonScaleFactors(uhh2::Context & ctx, TString path_, TString SysDirectionXSec_, TString SysDirectionBTag_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirectionXSec, SysDirectionBTag;
  bool is_diboson;

};


class MuonTrkWeights: public uhh2::AnalysisModule{

 public:
  explicit MuonTrkWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> Trk_SF;

  // protected:


};
