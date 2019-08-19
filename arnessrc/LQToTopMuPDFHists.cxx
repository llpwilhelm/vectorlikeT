#include "UHH2/LQToTopMu/include/LQToTopMuPDFHists.h"
#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <sstream>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;


LQToTopMuPDFHists::LQToTopMuPDFHists(Context & ctx, const string & dirname, bool use_ntupleweights_, bool use_pdf_weights_): Hists(ctx, dirname), use_ntupleweights(use_ntupleweights_), use_pdf_weights(use_pdf_weights_){  

  is_mc = ctx.get("dataset_type") == "MC";
  //For MLQ reconstruction
  h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassLQReconstruction");
  h_muonic_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>("HighMassMuonicLQReconstruction");
  m_discriminator_name ="Chi2"; 
  m_oname = ctx.get("dataset_version");
  is_LO = m_oname.Contains("LQtoT") || m_oname.Contains("Diboson") || m_oname.Contains("DYJets") || m_oname.Contains("QCD");
  TString m_pdfname = "NNPDF30_lo_as_0130";
  //if(!is_LO && !m_oname.Contains("SingleTop")) m_pdfname = "PDF4LHC15_nlo_mc";
  TString weightpath = ctx.get("PDFWeightPath");  cout << "File: " << weightpath+m_oname << endl; 

  //take ntupleweights if
  //1) use_ntupleweights = true and the sample has ntupleweights stored

  //take weights from txt files if
  //1) the sample is LO (and doesn't have ntupleweights for that reason) (this assumption is protected by a runtime_error later)
  //2) the sample is NLO and yet doesn't have ntupleweights 

  take_ntupleweights =  use_ntupleweights && (!is_LO || m_oname.Contains("DYJets")) && (m_oname != "SingleTop_T_tWch" && m_oname != "SingleTop_Tbar_tWch");

  cout << "For this sample '" << m_oname << "' is_LO is set to " << is_LO << endl;
  cout << "Are ntupleweights taken for this sample?: " << take_ntupleweights << endl;

  if(is_mc && !take_ntupleweights){
    if(m_oname.Contains("LQtoT")) m_pdfweights.reset(new PDFWeights(m_pdfname,weightpath+m_oname)); 
    else m_pdfweights.reset(new PDFWeights(m_pdfname)); 
  }

  for(int i=0; i<100; i++){
    stringstream ss_name;
    ss_name << "H_T_from350_rebin_PDF_"  << i+1 ;
    stringstream ss_name2;
    ss_name2 << "H_T_from350_all_filled_rebin_PDF_"  << i+1 ;
    stringstream ss_name3;
    ss_name3 << "M_LQ_MuEle_comb_all_filled_PDF_"  << i+1 ;
    stringstream ss_name4;
    ss_name4 << "H_T_comb_NoMLQ_all_filled_rebin_PDF_"  << i+1 ;

    stringstream ss_title;
    ss_title << "H_{T} [GeV] in alpha-binning for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title2;
    ss_title2 << "H_{T} [GeV] (from 350) for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title3;
    ss_title3 << "M_{LQ} [GeV] for PDF No. "  << i+1 << " out of 100" ;
    stringstream ss_title4;
    ss_title4 << "H_{T} [GeV] No Ele for PDF No. "  << i+1 << " out of 100" ;

    string s_name = ss_name.str();
    string s_title = ss_title.str();
    string s_name2 = ss_name2.str();
    string s_title2 = ss_title2.str();
    string s_name3 = ss_name3.str();
    string s_title3 = ss_title3.str();
    string s_name4 = ss_name4.str();
    string s_title4 = ss_title4.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();
    const char* char_name2 = s_name2.c_str();
    const char* char_title2 = s_title2.c_str();
    const char* char_name3 = s_name3.c_str();
    const char* char_title3 = s_title3.c_str();
    const char* char_name4 = s_name4.c_str();
    const char* char_title4 = s_title4.c_str();
    histo_names[i] = s_name;
    histo_names2[i] = s_name2;
    histo_names3[i] = s_name3;
    histo_names4[i] = s_name4;


    double bins_low_NoEle2[11] = {0,350,500,650,800,950,1100,1250,1450,1750,2050};
    double bins_mlq_low3[7] = {0,250,350,450,600,750,1000};
   double bins_from350_allfilled[16] = {0,175,350,525,700,875,1050,1225,1400,1575,1750,1925,2100,2450,2800,3000}; //same binning as _from350_all_filled up to 2100, then \
larger bins 
    book<TH1F>(char_name, char_title, 48,0,4200);
    book<TH1F>(char_name2, char_title2, 15,bins_from350_allfilled);  
    book<TH1F>(char_name3, char_title3,6, bins_mlq_low3);
    book<TH1F>(char_name4, char_title4, 15,bins_from350_allfilled);
  }

}

void LQToTopMuPDFHists::fill(const Event & event){
  double weight = event.weight;

  if(is_mc){
    //HT
    auto met = event.met->pt();
    double ht = 0.0;
    double ht_jets = 0.0;
    double ht_lep = 0.0;
    for(const auto & jet : *event.jets){
      ht_jets += jet.pt();
    }
    for(const auto & electron : *event.electrons){
      ht_lep += electron.pt();
    }
    for(const auto & muon : *event.muons){
      ht_lep += muon.pt();
    }
    
    ht = ht_lep + ht_jets + met;


    //collect information on wether MLQ is reconstructable or it is not
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

    bool reconstruct_mlq_ele = (event.electrons->size() >= 1 && event.muons->size() >= 2 && charge_opposite);
    bool reconstruct_mlq_mu = (event.electrons->size() == 0 && event.muons->size() == 3 && fabs(sum_mu_charge) == 1);

    if(event.genInfo->systweights().size() == 0 && take_ntupleweights) throw runtime_error("In LQToTopMuPDFHists.cxx: Systweights in event.genInfo() is empty but ntupleweights shall be taken. Is this correct? In this case add exception to take_ntupleweights.");    
    if(event.genInfo->systweights().size() != 0 && (is_LO && !m_oname.Contains("DYJets") && !m_oname.Contains("Diboson"))) throw runtime_error("In LQToTopMuPDFHists.cxx: Systweights in event.genInfo() is NOT empty but this IS a LO sample. Is this correct? In this case Thomas says the genInfo weight should be used. Add this sample to take_ntupleweights");

    //Fill HT (2 cases)
    if(take_ntupleweights){
      for(int i=0; i<100; i++){
	if(use_pdf_weights){
 
	  double pdf_weight = event.genInfo->systweights().at(i+9);
	  double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	  const char* name = histo_names[i].c_str();
	  const char* name2 = histo_names2[i].c_str();
	  const char* name4 = histo_names4[i].c_str();

	  hist(name)->Fill(ht,fillweight);

	  if(ht <= 2900) hist(name2)->Fill(ht,fillweight);
	  else hist(name2)->Fill(2900,fillweight);

	  if(!reconstruct_mlq_ele && !reconstruct_mlq_mu){
	    if(ht <= 2900) hist(name4)->Fill(ht,fillweight);
	    else hist(name4)->Fill(2900,fillweight);
	  }
	}
      }
    }
    else{ 
      std::vector<double> weights = m_pdfweights->GetWeightList(event);
      for(int i=0; i<100; i++){
	if(use_pdf_weights){
	  double fillweight = weight*weights[i];
	  const char* name = histo_names[i].c_str();
	  const char* name2 = histo_names2[i].c_str();
	  const char* name4 = histo_names4[i].c_str();

	  hist(name)->Fill(ht,fillweight);

	  if(ht <= 2900) hist(name2)->Fill(ht,fillweight);
	  else hist(name2)->Fill(2900,fillweight);

	  if(!reconstruct_mlq_ele && !reconstruct_mlq_mu){
	    if(ht <= 2900) hist(name4)->Fill(ht,fillweight);
	    else hist(name4)->Fill(2900,fillweight);
	  }
	}
      }
    }



    if(reconstruct_mlq_ele){   
      std::vector<LQReconstructionHypothesis> hyps = event.get(h_hyps); 
      const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );

      double mLQlep_rec = 0;
      double mLQhad_rec = 0;
      double mLQmed_rec = 0;


      if( (hyp->LQlep_v4()).isTimelike() ) mLQlep_rec = (hyp->LQlep_v4()).M();
      else mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());
      if( (hyp->LQhad_v4()).isTimelike() ) mLQhad_rec = (hyp->LQhad_v4()).M();
      else mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());

      mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;

      if(take_ntupleweights){
	for(int i=0; i<100; i++){
	  if(use_pdf_weights){
	    double pdf_weight = event.genInfo->systweights().at(i+9);
	    double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	    const char* name3 = histo_names3[i].c_str();
	    if(mLQmed_rec < 900) hist(name3)->Fill(mLQmed_rec,fillweight);
	    else hist(name3)->Fill(900,fillweight);
	  }
	}
      } //systweights filled
      else{
	std::vector<double> weights = m_pdfweights->GetWeightList(event);
	for(int i=0; i<100; i++){
          if(use_pdf_weights){
	    double fillweight = weight*weights[i]; // LQ PDF
	    const char* name3 = histo_names3[i].c_str();
	    if(mLQmed_rec < 900) hist(name3)->Fill(mLQmed_rec,fillweight);
	    else hist(name3)->Fill(900,fillweight);
          }
	}
      }
    } //LQ reconstructable
    else if(reconstruct_mlq_mu){
      std::vector<LQReconstructionHypothesis> muonic_hyps = event.get(h_muonic_hyps); 
      const LQReconstructionHypothesis* hyp = get_best_hypothesis( muonic_hyps, m_discriminator_name );
    
      double mLQlep_rec = 0;
      double mLQhad_rec = 0;
      double mLQmed_rec = 0;


      if( (hyp->LQlep_v4()).isTimelike() ) mLQlep_rec = (hyp->LQlep_v4()).M();
      else mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());
      if( (hyp->LQhad_v4()).isTimelike() ) mLQhad_rec = (hyp->LQhad_v4()).M();
      else mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());

      mLQmed_rec = (mLQhad_rec + mLQlep_rec) / 2;

      if(take_ntupleweights){
	for(int i=0; i<100; i++){
	  if(use_pdf_weights){
	    double pdf_weight = event.genInfo->systweights().at(i+9);
	    double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
	    const char* name3 = histo_names3[i].c_str();
	    if(mLQmed_rec < 900) hist(name3)->Fill(mLQmed_rec,fillweight);
	    else hist(name3)->Fill(900,fillweight);
	  }
	}
      } //systweights filled
      else{
	std::vector<double> weights = m_pdfweights->GetWeightList(event);
	for(int i=0; i<100; i++){
          if(use_pdf_weights){
	    double fillweight = weight*weights[i]; // LQ PDF
	    const char* name3 = histo_names3[i].c_str();
	    if(mLQmed_rec < 900) hist(name3)->Fill(mLQmed_rec,fillweight);
	    else hist(name3)->Fill(900,fillweight);
          }
	}
      }


    }

  } //is_mc
}

LQToTopMuPDFHists::~LQToTopMuPDFHists(){}














