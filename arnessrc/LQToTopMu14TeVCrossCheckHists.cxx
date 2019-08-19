#include "../include/LQToTopMu14TeVCrossCheckHists.h"
#include "../include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

LQToTopMu14TeVCrossCheckHists::LQToTopMu14TeVCrossCheckHists(Context & ctx, const string & dirname): Hists(ctx, dirname){

  //Histograms
  book<TH1F>("x13_1", "x_{1}^{13TeV}", 100, 0, 1);  
  book<TH1F>("x13_2", "x_{2}^{13TeV}", 100, 0, 1);  
  book<TH1F>("x14_1", "x_{1}^{14TeV}", 100, 0, 1);  
  book<TH1F>("x14_2", "x_{2}^{14TeV}", 100, 0, 1);   
  book<TH1F>("q", "Momentum transfer Q [GeV]", 50, 0, 1000); 
  book<TH1F>("f_1", "flavor (pdgId) of parton 1", 51, -25.5, 25.5);  
  book<TH1F>("f_2", "flavor (pdgId) of parton 2", 51, -25.5, 25.5);     
  book<TH1F>("xf13_1", "x_{1}^{13TeV} #times f(x_{1}^{13TeV}, Q, f_{1})", 200, 0, 50); 
  book<TH1F>("xf13_2", "x_{2}^{13TeV} #times f(x_{2}^{13TeV}, Q, f_{2})", 200, 0, 50); 
  book<TH1F>("xf14_1", "x_{1}^{14TeV} #times f(x_{1}^{14TeV}, Q, f_{1})", 200, 0, 50); 
  book<TH1F>("xf14_2", "x_{2}^{14TeV} #times f(x_{2}^{14TeV}, Q, f_{2})", 200, 0, 50); 
  book<TH1F>("weight13", "PDF weight (13TeV)", 100, 0, 10); 
  book<TH1F>("weight14", "PDF weight (14TeV)", 100, 0, 10); 
  book<TH1F>("sf", "PDF scale factor (13TeV->14TeV)", 100, 0, 10); 

  book<TH1F>("x13_1_noweight", "x_{1}^{13TeV} , unweighted events", 20, 0, 1);  
  book<TH1F>("x13_2_noweight", "x_{2}^{13TeV} , unweighted events", 20, 0, 1);  
  book<TH1F>("x14_1_noweight", "x_{1}^{14TeV} , unweighted events", 20, 0, 1);  
  book<TH1F>("x14_2_noweight", "x_{2}^{14TeV} , unweighted events", 20, 0, 1);   
  book<TH1F>("q_noweight", "Momentum transfer Q [GeV] , unweighted events", 50, 0, 1000); 
  book<TH1F>("f_1_noweight", "flavor (pdgId) of parton 1 , unweighted events", 51, -25.5, 25.5);  
  book<TH1F>("f_2_noweight", "flavor (pdgId) of parton 2 , unweighted events", 51, -25.5, 25.5);     
  book<TH1F>("xf13_1_noweight", "x_{1}^{13TeV} #times f(x_{1}^{13TeV}, Q, f_{1}) , unweighted events", 200, 0, 50); 
  book<TH1F>("xf13_2_noweight", "x_{2}^{13TeV} #times f(x_{2}^{13TeV}, Q, f_{2}) , unweighted events", 200, 0, 50); 
  book<TH1F>("xf14_1_noweight", "x_{1}^{14TeV} #times f(x_{1}^{14TeV}, Q, f_{1}) , unweighted events", 200, 0, 50); 
  book<TH1F>("xf14_2_noweight", "x_{2}^{14TeV} #times f(x_{2}^{14TeV}, Q, f_{2}) , unweighted events", 200, 0, 50); 
  book<TH1F>("weight13_noweight", "PDF weight (13TeV) , unweighted events", 100, 0, 10); 
  book<TH1F>("weight14_noweight", "PDF weight (14TeV) , unweighted events", 100, 0, 10); 
  book<TH1F>("sf_noweight", "PDF scale factor (13TeV->14TeV) , unweighted events", 100, 0, 10); 


  //Handles containing the values
  h_x13_1 = ctx.get_handle<double>("x13_1");
  h_x13_2 = ctx.get_handle<double>("x13_2");
  h_x14_1 = ctx.get_handle<double>("x14_1");
  h_x14_2 = ctx.get_handle<double>("x14_2");
  h_Q = ctx.get_handle<double>("Q");
  h_xf13_1 = ctx.get_handle<double>("xf13_1");
  h_xf13_2 = ctx.get_handle<double>("xf13_1");
  h_xf14_1 = ctx.get_handle<double>("xf14_1");
  h_xf14_2 = ctx.get_handle<double>("xf14_2");
  h_f1 = ctx.get_handle<int>("f1");
  h_f2 = ctx.get_handle<int>("f2");
  h_weight13 = ctx.get_handle<double>("weight13");
  h_weight14 = ctx.get_handle<double>("weight14");
  h_sf = ctx.get_handle<double>("sf");

  is_mc = ctx.get("dataset_type") == "MC";

}


void LQToTopMu14TeVCrossCheckHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  

  double weight = event.weight;

  double x13_1    = event.get(h_x13_1);
  double x13_2    = event.get(h_x13_2);
  double x14_1    = event.get(h_x14_1);
  double x14_2    = event.get(h_x14_2);
  double q        = event.get(h_Q);
  double xf13_1   = event.get(h_xf13_1);
  double xf13_2   = event.get(h_xf13_2);
  double xf14_1   = event.get(h_xf14_1);
  double xf14_2   = event.get(h_xf14_2);
  double f_1       = event.get(h_f1);
  double f_2       = event.get(h_f2);
  double weight13 = event.get(h_weight13);
  double weight14 = event.get(h_weight14);
  double sf       = event.get(h_sf);

  // Where these histograms are filled, the 13->14TeV weight has already been applied
  double fillweight = weight / sf;

  hist("x13_1")->Fill(x13_1, fillweight);
  hist("x13_2")->Fill(x13_2, fillweight);
  hist("x14_1")->Fill(x14_1, fillweight);
  hist("x14_2")->Fill(x14_2, fillweight);
  hist("q")->Fill(q, fillweight);
  hist("xf13_1")->Fill(xf13_1, fillweight);
  hist("xf13_2")->Fill(xf13_2, fillweight);
  hist("xf14_1")->Fill(xf14_1, fillweight);
  hist("xf14_2")->Fill(xf14_2, fillweight);
  hist("f_1")->Fill(f_1, fillweight);
  hist("f_2")->Fill(f_2, fillweight);
  hist("weight13")->Fill(weight13, fillweight);
  hist("weight14")->Fill(weight14, fillweight);
  hist("sf")->Fill(sf, fillweight);

  hist("x13_1_noweight")->Fill(x13_1, 1.);
  hist("x13_2_noweight")->Fill(x13_2, 1.);
  hist("x14_1_noweight")->Fill(x14_1, 1.);
  hist("x14_2_noweight")->Fill(x14_2, 1.);
  hist("q_noweight")->Fill(q, 1.);
  hist("xf13_1_noweight")->Fill(xf13_1, 1.);
  hist("xf13_2_noweight")->Fill(xf13_2, 1.);
  hist("xf14_1_noweight")->Fill(xf14_1, 1.);
  hist("xf14_2_noweight")->Fill(xf14_2, 1.);
  hist("f_1_noweight")->Fill(f_1, 1.);
  hist("f_2_noweight")->Fill(f_2, 1.);
  hist("weight13_noweight")->Fill(weight13, 1.);
  hist("weight14_noweight")->Fill(weight14, 1.);
  hist("sf_noweight")->Fill(sf, 1.);



} //Methode



LQToTopMu14TeVCrossCheckHists::~LQToTopMu14TeVCrossCheckHists(){}
