#include "UHH2/LQToTopMu/include/HypothesisHistsOwn.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;
using namespace std;

HypothesisHistsOwn::HypothesisHistsOwn(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  TString name = discriminator_name;
  double min=0;
  double max=500;
  
  if(discriminator_name=="Chi2"){
    name = "#Chi^{2}";
  }
  else{
    name += " discriminator";
  }

  if( discriminator_name=="CorrectMatch"){
    min=0;
    max=2;
  }


    Discriminator = book<TH1F>("Discriminator",name,100,min,max);
    Discriminator_2 = book<TH1F>("Discriminator_2",name,50,0,10);
    Discriminator_3 = book<TH1F>("Discriminator_3",name,300,0,30); 
 
    M_LQlep_rec  = book<TH1F>("M_LQlep_rec", "M_{LQ,lep}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQhad_rec  = book<TH1F>("M_LQhad_rec", "M_{LQ,had}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmax_rec  = book<TH1F>("M_LQmax_rec", "M_{LQmax}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmean_rec = book<TH1F>("M_LQmean_rec", "M_{LQmean}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQ_rec_diff = book<TH1F>("M_LQ_rec_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV/c^2]", 50, -500, 500);
    M_LQ_rec_diff_rel = book<TH1F>("M_LQ_rec_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep}) / M_{LQ}^{mean} [GeV/c^]", 50, -0.5, 0.5);
    M_ttbar_rec = book<TH1F>( "M_ttbar_rec", "M_{t#bar{t}}^{rec} [GeV/c^{2}]", 100, 0, 5000 ) ;
    M_toplep_rec = book<TH1F>( "M_toplep_rec", "M^{top,lep} [GeV/c^{2}]", 70, 0, 700 ) ;
    M_tophad_rec = book<TH1F>( "M_tophad_rec", "M^{top,had} [GeV/c^{2}]", 70, 0, 700 ) ;
    Pt_ttbar_rec = book<TH1F>( "Pt_ttbar_rec", "P_{T,t#bar{t}}^{rec} [GeV/c]", 60, 0, 600 ) ;
    Pt_diff_LQ = book<TH1F>( "Pt_diff_LQ", "#Delta p_{T}^{LQ} [GeV]", 50, -500, 500 ) ;
    Pt_diff_LQ_rel = book<TH1F>( "Pt_diff_LQ_rel", "#Delta p_{T}^{LQ} [GeV]", 50, -0.5, 0.5 ) ;
    Pt_ratio = book<TH1F>( "Pt_ratio", "p_{T}^{LQ had} / p_{T}^{LQ lep} [GeV]", 100, -5, 5 ) ;
    DeltaPhi_LQ_LQ = book<TH1F>( "DeltaPhi_LQ_LQ", "#Delta #phi (LQ_{lep}, LQ_{had})", 126, 0, 6.3 ) ;


    M_LQlep_rec_discriminator_cut  = book<TH1F>("M_LQlep_rec_discriminator_cut", "M_{LQ,lep}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQhad_rec_discriminator_cut  = book<TH1F>("M_LQhad_rec_discriminator_cut", "M_{LQ,had}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmax_rec_discriminator_cut  = book<TH1F>("M_LQmax_rec_discriminator_cut", "M_{LQmax}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQmean_rec_discriminator_cut = book<TH1F>("M_LQmean_rec_discriminator_cut", "M_{LQmean}^{rec} [GeV/c^{2}]", 40, 0, 2000 );
    M_LQ_rec_discriminator_cut_diff = book<TH1F>("M_LQ_rec_discriminator_cut_diff", "M_{LQ}^{had} - M_{LQ}^{lep} [GeV/c^2]", 50, -500, 500);
    M_LQ_rec_discriminator_cut_diff_rel = book<TH1F>("M_LQ_rec_discriminator_cut_diff_rel", "(M_{LQ}^{had} - M_{LQ}^{lep}) / M_{LQ}^{mean} [GeV/c^]", 50, -0.5, 0.5);
    M_ttbar_rec_discriminator_cut = book<TH1F>( "M_ttbar_rec_discriminator_cut", "M_{t#bar{t}}^{rec} [GeV/c^{2}]", 100, 0, 5000 ) ;
    M_toplep_rec_discriminator_cut = book<TH1F>( "M_toplep_rec_discriminator_cut", "M^{top,lep} [GeV/c^{2}]", 70, 0, 700 ) ;
    M_tophad_rec_discriminator_cut = book<TH1F>( "M_tophad_rec_discriminator_cut", "M^{top,had} [GeV/c^{2}]", 70, 0, 700 ) ;
    Pt_ttbar_rec_discriminator_cut = book<TH1F>( "Pt_ttbar_rec_discriminator_cut", "P_{T,t#bar{t}}^{rec} [GeV/c]", 60, 0, 600 ) ;
    Pt_diff_LQ_discriminator_cut = book<TH1F>( "Pt_diff_LQ_discriminator_cut", "#Delta p_{T}^{LQ} [GeV]", 50, -500, 500 ) ;
    Pt_diff_LQ_rel_discriminator_cut = book<TH1F>( "Pt_diff_LQ_rel_discriminator_cut", "#Delta p_{T}^{LQ} [GeV]", 50, -0.5, 0.5 ) ;
    Pt_ratio_discriminator_cut = book<TH1F>( "Pt_ratio_discriminator_cut", "p_{T}^{LQ had} / p_{T}^{LQ lep} [GeV]", 100, -5, 5 ) ;
    DeltaPhi_LQ_LQ_discriminator_cut = book<TH1F>( "DeltaPhi_LQ_LQ_discriminator_cut", "#Delta #phi (LQ_{lep}, LQ_{had})", 126, 0, 6.3 ) ;

    
    
    h_hyps = ctx.get_handle<std::vector<LQReconstructionHypothesis>>(hyps_name);
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    m_discriminator_name = discriminator_name;
}


void HypothesisHistsOwn::fill(const uhh2::Event & e){
  bool charge_opposite = false;
  for(unsigned int i=0; i<e.muons->size(); i++){
    for(unsigned int j=0; j<e.muons->size(); j++){
      if(j>i){
	if(e.muons->at(i).charge() != e.muons->at(j).charge()) {
	  charge_opposite = true;
	}
      }
    }
  }
  if(charge_opposite){
    std::vector<LQReconstructionHypothesis> hyps = e.get(h_hyps);
    const LQReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
    double weight = e.weight;

    double mttbar_rec = 0;


    if( (hyp->top_v4()+hyp->antitop_v4()).isTimelike() )
      mttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).M();
    else{
      mttbar_rec = sqrt( -(hyp->top_v4()+hyp->antitop_v4()).mass2());
    }
    double ptttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).Pt();

    M_ttbar_rec->Fill(mttbar_rec, weight);
    Pt_ttbar_rec->Fill ( ptttbar_rec, weight);
  
    double mtoplep=0;
    double mtophad=0;

    if(hyp->toplep_v4().isTimelike()) mtoplep = hyp->toplep_v4().M();
    else mtoplep = sqrt( -(hyp->toplep_v4()).mass2());
    if(hyp->tophad_v4().isTimelike()) mtophad = hyp->tophad_v4().M();
    else mtophad = sqrt( -(hyp->tophad_v4()).mass2());
    M_toplep_rec->Fill(mtoplep,weight);
    M_tophad_rec->Fill(mtophad,weight);

    Discriminator->Fill(hyp->discriminator(m_discriminator_name) ,weight);
    Discriminator_2->Fill(hyp->discriminator(m_discriminator_name) ,weight); 
    Discriminator_3->Fill(hyp->discriminator(m_discriminator_name) ,weight);

  
    //Combine Top and Muon (Electron and Muon Charge have to be opposite)
    double mLQlep_rec = 0;
    double mLQhad_rec = 0;
    double mLQmax_rec = 0;
    double mLQmed_rec = 0;
    double mLQ_rec_diff = 0;
    double mLQ_rec_diff_rel = 0;

    if( (hyp->LQlep_v4()).isTimelike() ) {mLQlep_rec = (hyp->LQlep_v4()).M();}
    else {mLQlep_rec = sqrt( -(hyp->LQlep_v4()).mass2());}
    if( (hyp->LQhad_v4()).isTimelike() ) {mLQhad_rec = (hyp->LQhad_v4()).M();}
    else {mLQhad_rec = sqrt( -(hyp->LQhad_v4()).mass2());}
  
    if(mLQhad_rec>mLQlep_rec){
      mLQmax_rec = mLQhad_rec;
    }
    else{
      mLQmax_rec = mLQlep_rec;
    }
  
    mLQmed_rec = (mLQhad_rec + mLQlep_rec)/2;
    mLQ_rec_diff = mLQhad_rec - mLQlep_rec;
    mLQ_rec_diff_rel = mLQ_rec_diff/mLQmed_rec;

    double deltaPt_LQ = hyp->LQhad_v4().Pt() - hyp->LQlep_v4().Pt();
    double deltaPt_LQ_rel = deltaPt_LQ/((hyp->LQhad_v4().Pt() + hyp->LQlep_v4().Pt())/2);
    double pt_ratio = hyp->LQhad_v4().Pt() / hyp->LQlep_v4().Pt();

    //calculate dPhi in [0, 2pi]
    double dphi = (hyp->LQhad_v4().Phi()) - (hyp->LQlep_v4().Phi());
    if(dphi < 0) dphi += 2*M_PI;
    
    M_LQlep_rec->Fill(mLQlep_rec, weight);
    M_LQhad_rec->Fill(mLQhad_rec, weight);
    M_LQmax_rec->Fill(mLQmax_rec, weight);
    M_LQmean_rec->Fill(mLQmed_rec, weight);
    M_LQ_rec_diff->Fill(mLQ_rec_diff, weight);
    M_LQ_rec_diff_rel->Fill(mLQ_rec_diff_rel,weight);
    Pt_diff_LQ->Fill(deltaPt_LQ,weight);
    Pt_diff_LQ_rel->Fill(deltaPt_LQ_rel,weight);
    Pt_ratio->Fill(pt_ratio, weight);
    DeltaPhi_LQ_LQ->Fill(dphi, weight);
    if(hyp->discriminator(m_discriminator_name) < 20){ // 999999 is set as discr. value on CorrectMatchDiscr, if one of the required particles is not matched
      M_LQ_rec_discriminator_cut_diff->Fill(mLQ_rec_diff, weight);
      M_LQ_rec_discriminator_cut_diff_rel->Fill(mLQ_rec_diff_rel,weight);
      M_LQlep_rec_discriminator_cut->Fill(mLQlep_rec, weight);
      M_LQhad_rec_discriminator_cut->Fill(mLQhad_rec, weight);
      M_LQmax_rec_discriminator_cut->Fill(mLQmax_rec, weight);
      M_LQmean_rec_discriminator_cut->Fill(mLQmed_rec, weight);
      M_toplep_rec_discriminator_cut->Fill(mtoplep,weight);
      M_tophad_rec_discriminator_cut->Fill(mtophad,weight);
      M_ttbar_rec_discriminator_cut->Fill(mttbar_rec, weight);
      Pt_ttbar_rec_discriminator_cut->Fill ( ptttbar_rec, weight);
      Pt_diff_LQ_discriminator_cut->Fill(deltaPt_LQ,weight);
      Pt_diff_LQ_rel_discriminator_cut->Fill(deltaPt_LQ_rel,weight);
      Pt_ratio_discriminator_cut->Fill(pt_ratio, weight);
      DeltaPhi_LQ_LQ_discriminator_cut->Fill(dphi, weight);
    }
  }





  
}
