#include "UHH2/LQToTopMu/include/LQGenHists.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;

LQGenHists::LQGenHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

    MLQLQbar_gen =  book<TH1F>( "M_LQLQbar_gen", "M_{LQ#bar{LQ}} [GeV/c^{2}]", 1000, 0, 5000 ) ;
    Pt_LQLQbar_gen = book< TH1F>( "Pt_LQLQbar_gen", "P_{T,LQ#bar{LQ}} [GeV/c]", 600, 0, 600 ) ;
    shat = book< TH1F>( "shat", "#hat{s} [GeV]", 1000, 0, 5000 ) ;
    //DecayChannel = book< TH1F>( "DecayChannel", "decay channel", 11, 0, 11 ) ;

    M_LQ = book< TH1F>("M_LQ", "M_{LQ} [GeV/c^{2}]", 30000, 0, 1650) ;
    M_antiLQ = book< TH1F>("M_antiLQ", "M_{#bar{LQ}} [GeV/c^{2}]", 30000, 0, 1650) ;  
    M_LQs = book< TH1F>("M_LQs", "M_{LQs} [GeV/c^{2}]", 30000, 0, 1650) ;


    Pt_LQ = book< TH1F>( "Pt_LQ", "P_{T,LQ} [GeV/c]", 1000, 0, 2000 ) ;
    Pt_antiLQ = book< TH1F>( "Pt_antiLQ", "P_{T,#bar{LQ}} [GeV/c]", 1000, 0, 2000 ) ;
    Pt_LQs = book< TH1F>( "Pt_LQs", "P_{T,LQs} [GeV/c]", 1000, 0, 2000 ) ;

    phi_LQ = book< TH1F>( "phi_LQ", "#phi_{LQ} [GeV/c]", 25, -M_PI, M_PI ) ;
    phi_antiLQ = book< TH1F>( "phi_antiLQ", "#phi_{#bar{LQ}} [GeV/c]", 25, -M_PI, M_PI ) ;
    phi_LQs = book< TH1F>( "phi_LQs", "#phi_{LQs} [GeV/c]", 25, -M_PI, M_PI ) ;
    
    Pt_LQ_over_shat = book< TH1F>( "Pt_LQ_over_shat", "P_{T,LQ}/#hat{s}", 1000, 0, 1 ) ;
    Pt_antiLQ_over_shat = book< TH1F>( "Pt_antiLQ_over_shat", "P_{T,#bar{LQ}}/#hat{s}", 1000, 0, 1 ) ;
    Pt_LQ_over_M_LQLQbar = book< TH1F>( "Pt_LQ_over_M_LQLQbar", "P_{T,LQ}/M_{LQ#bar{LQ}}", 1000, 0, 1 ) ;
    Pt_antiLQ_over_M_LQLQbar = book< TH1F>( "Pt_antiLQ_over_M_LQLQbar", "P_{T,#bar{LQ}}/M_{LQ#bar{LQ}}", 1000, 0, 1 ) ;
    
    M_LQLQbar_vs_shat = book< TH2F>( "M_LQLQbar_vs_shat","M_{LQ#bar{LQ}} [GeV/c^{2}] vs #hat{s} [GeV]", 500,0,5000,500,0,5000);
    
    eta_LQ = book< TH1F>( "eta_LQ", "#eta_{LQ}", 1000, -5, 5 ) ;
    eta_antiLQ = book< TH1F>( "eta_antiLQ", "#eta_{#bar{LQ}}", 1000, -5, 5 ) ;
    eta_LQs = book< TH1F>( "eta_LQs", "#eta_{LQs}", 1000, -5, 5 ) ;

    y_LQ = book< TH1F>( "y_LQ", "y_{LQ}", 1000, -5, 5 ) ;
    y_antiLQ = book< TH1F>( "y_antiLQ", "y_{#bar{LQ}}", 1000, -5, 5 ) ;
    y_LQs = book< TH1F>( "y_LQs", "y_{LQs}", 1000, -5, 5 ) ;    

    diffabseta = book<TH1F>( "diffabseta", "|#eta_{LQ}|-|#eta_{#bar{LQ}}|",1000,-5,5);
    diffabsy = book<TH1F>( "diffabsy", "|y_{LQ}|-|y_{#bar{LQ}}|",1000,-5,5);
    
    deltaR_LQ_decays = book<TH1F>( "deltaR_LQ_decays", "#DeltaR(LQ decay prod.)",1000,0,5);
    deltaR_antiLQ_decays = book<TH1F>( "deltaR_antiLQ_decays", "#DeltaR(#bar{LQ} decay prod.)",1000,0,5);  
    M_LQLQbar_vs_deltaR_LQ = book<TH2F>( "M_LQLQbar_vs_deltaR_LQ", "M_{LQ#bar{LQ}} [GeV/c^{2}] vs #DeltaR(LQ decay prod.)",500,0,5000,500,0,5);
    M_LQLQbar_vs_deltaR_antiLQ = book<TH2F>( "M_LQLQbar_vs_deltaR_antiLQ", "M_{LQ#bar{LQ}} [GeV/c^{2}] vs #DeltaR(#bar{LQ} decay prod.)",500,0,5000,500,0,5); 
    
    shat_vs_deltaR_LQ = book<TH2F>( "shat_vs_deltaR_LQ", "#hat{s} [GeV] vs #DeltaR(LQ decay prod.)",500,0,5000,500,0,5);
    shat_vs_deltaR_antiLQ = book<TH2F>( "shat_vs_deltaR_antiLQ", "#hat{s} [GeV] vs #DeltaR(#bar{LQ} decay prod.)",500,0,5000,500,0,5); 

    Pt_LQ_vs_deltaR_LQ = book<TH2F>( "Pt_LQ_vs_deltaR_LQ", "P_{T,LQ} [GeV/c] vs #DeltaR(LQ decay prod.)",500,0,2000,500,0,5);
    Pt_antiLQ_vs_deltaR_antiLQ = book<TH2F>( "Pt_antiLQ_vs_deltaR_antiLQ", "P_{T,#bar{LQ}} [GeV/c] vs #DeltaR(#bar{LQ} decay prod.)",500,0,2000,500,0,5);

    /*deltaR_Wplus_decays = book<TH1F>( "deltaR_Wplus_decays", "#DeltaR(W^{+} decay prod.)",1000,0,5);
    deltaR_Wminus_decays = book<TH1F>( "deltaR_Wminus_decays", "#DeltaR(W^{-} decay prod.)",1000,0,5);  
    M_ttbar_vs_deltaR_Wplus = book<TH2F>( "M_ttbar_vs_deltaR_Wplus", "M_{t#bar{t}} [GeV/c^{2}] vs #DeltaR(W^{+} decay prod.)",500,0,5000,500,0,5);
    M_ttbar_vs_deltaR_Wminus = book<TH2F>( "M_ttbar_vs_deltaR_Wminus", "M_{t#bar{t}} [GeV/c^{2}] vs #DeltaR(W^{-} decay prod.)",500,0,5000,500,0,5); */

    M_LQLQbar_vs_Pt_LQ = book<TH2F>( "M_LQLQbar_vs_Pt_LQ", "M_{LQ#bar{LQ}} [GeV/c^{2}] vs P_{T,LQ} [GeV/c]",500,0,5000,500,0,2000); 
    M_LQLQbar_vs_Pt_antiLQ = book<TH2F>( "M_LQLQbar_vs_Pt_antiLQ", "M_{LQ#bar{LQ}} [GeV/c^{2}] vs P_{T,#bar{LQ}} [GeV/c]",500,0,5000,500,0,2000);
    shat_vs_Pt_LQ = book<TH2F>( "shat_vs_Pt_LQ", "#hat{s} [GeV] vs P_{T,LQ} [GeV/c]",500,0,5000,500,0,2000); 
    shat_vs_Pt_antiLQ = book<TH2F>( "shat_vs_Pt_antiLQ", "#hat{s} [GeV] vs P_{T,#bar{LQ}} [GeV/c]",500,0,5000,500,0,2000);
    Pt_LQ_vs_Pt_antiLQ = book<TH2F>( "Pt_LQ_vs_Pt_antiLQ", "P_{T,LQ} [GeV/c] vs P_{T,#bar{LQ}} [GeV/c]",500,0,2000,500,0,2000);
    
    Pt_mu = book<TH1F>( "Pt_mu","P_{T,mu} [GeV/c]",250,0,1000);
    Pt_antimu = book<TH1F>( "Pt_antimu","P_{T,#bar{mu}} [GeV/c]",250,0,1000);
    Pt_mus = book<TH1F>( "Pt_mus","P_{T,mus} [GeV/c]",250,0,1000);
    eta_mu = book<TH1F>( "eta_mu","#eta_{mu}",100,-5,5);
    eta_antimu = book<TH1F>( "eta_antimu","#eta_{#bar{mu}}",100,-5,5);
    eta_mus = book<TH1F>( "eta_mus","#eta_{mus}",100,-5,5);

    y_mu = book<TH1F>( "y_mu","y_{mu}",1000,-5,5);
    y_antimu = book<TH1F>( "y_antimu","y_{#bar{mu}}",1000,-5,5);
    y_mus = book<TH1F>( "y_mus","y_{mus}",1000,-5,5);

    phi_mu = book< TH1F>( "phi_mu", "#phi_{#mu}", 25, -M_PI, M_PI ) ;
    phi_antimu = book< TH1F>( "phi_antimu", "#phi_{#bar{#mu}}", 25, -M_PI, M_PI ) ;
    phi_mus = book< TH1F>( "phi_mus", "#phi_{#mus}", 25, -M_PI, M_PI ) ;

    M_mu = book<TH1F>( "M_mu","M_{mu} [GeV/c^{2}]",100,0,1);
    M_antimu = book<TH1F>( "M_antimu","M_{#bar{mu}} [GeV/c^{2}]",100,0,1);
    M_mus = book<TH1F>( "M_mus","M_{mus} [GeV/c^{2}]",100,0,1);


    Pt_top = book<TH1F>( "Pt_top","P_{T,t} [GeV/c]",500,0,2000);
    Pt_antitop = book<TH1F>( "Pt_antitop","P_{T,#bar{t}} [GeV/c]",500,0,2000);
    Pt_tops = book<TH1F>( "Pt_tops","P_{T,ts} [GeV/c]",500,0,2000);

    eta_top = book<TH1F>( "eta_top","#eta_{t}",100,-5,5);
    eta_antitop = book<TH1F>( "eta_antitop","#eta_{#bar{t}}",100,-5,5);
    eta_tops = book<TH1F>( "eta_tops","#eta_{ts}",100,-5,5);

    y_top = book<TH1F>( "y_top","y_{t}",1000,-5,5);
    y_antitop = book<TH1F>( "y_antitop","y_{#bar{t}}",1000,-5,5);
    y_tops = book<TH1F>( "y_tops","y_{ts}",1000,-5,5);

    phi_top = book< TH1F>( "phi_top", "#phi_{t}", 25, -M_PI, M_PI ) ;
    phi_antitop = book< TH1F>( "phi_antitop", "#phi_{#bar{t}}", 25, -M_PI, M_PI ) ;
    phi_tops = book< TH1F>( "phi_tops", "#phi_{ts}", 25, -M_PI, M_PI ) ;

    M_top = book<TH1F>( "M_top","M_{t} [GeV/c^{2}]",1000,150,200);
    M_antitop = book<TH1F>( "M_antitop","M_{#bar{t}} [GeV/c^{2}]",1000,150,200);
    M_tops = book<TH1F>( "M_tops","M_{ts} [GeV/c^{2}]",1000,150,200);

  
    /*cosThetastar_top_ttframe = book< TH1F>( "cosThetastar_top_ttframe", "cos(#Theta*)_{t}",1000,-1,1);
    cosThetastar_antitop_ttframe = book< TH1F>( "cosThetastar_antitop_ttframe", "cos(#Theta*)_{#bar{t}}",1000,-1,1);
    Pt_top_ttframe = book< TH1F>( "Pt_top_ttframe", "P_{T,t}* [GeV/c]",1000,0,2000);
    Pt_antitop_ttframe = book< TH1F>( "Pt_antitop_ttframe", "P_{T,#bar{t}}* [GeV/c]",1000,0,2000);
    M_ttbar_vs_Pt_top_ttframe = book<TH2F>( "M_ttbar_vs_Pt_top_ttframe", "M_{t#bar{t}} [GeV/c^{2}] vs P_{T,t}* [GeV/c]",500,0,5000,500,0,2000); 
    M_ttbar_vs_Pt_antitop_ttframe = book<TH2F>( "M_ttbar_vs_Pt_antitop_ttframe", "M_{t#bar{t}} [GeV/c^{2}] vs P_{T,#bar{t}}* [GeV/c]",500,0,5000,500,0,2000);*/

    M_LQLQbar_vs_eta_LQ = book<TH2F>( "M_LQLQbar_vs_eta_LQ", "M_{LQ#bar{LQ}} [GeV/c^{2}] vs #eta_{LQ}",500,0,5000,500,-5,5); 
    M_LQLQbar_vs_eta_antiLQ = book<TH2F>( "M_LQLQbar_vs_eta_antiLQ", "M_{LQ#bar{LQ}} [GeV/c^{2}] vs #eta_{#bar{LQ}}",500,0,5000,500,-5,5);

    h_LQLQbargen = ctx.get_handle<LQGen>("LQLQbargen");
}


void LQGenHists::fill(const uhh2::Event & e){
    //do not fill histograms if LQLQbargen information has not been filled
    if(!e.is_valid(h_LQLQbargen)){
      return;
    }
    const auto & LQLQbargen = e.get(h_LQLQbargen);
    
    LorentzVector LQ = LQLQbargen.LQ().v4();
    LorentzVector antiLQ = LQLQbargen.AntiLQ().v4();

 
    double mLQLQbar_gen = (LQ+antiLQ).M();
    double  ptLQLQbar_gen = ( LQ + antiLQ ).Pt();
    double sh = (e.genparticles->at(0).v4()+ e.genparticles->at(1).v4()).M();

    //DecayChannel->Fill(LQLQbargen.DecayChannel(), e.weight);
    //if(LQ.M() == 0 || antiLQ.M() == 0){return;}

    //std::cout << "MLQLQ = " << mLQLQbar_gen << " & eventweight = " << e.weight << std::endl;
    //std::cout << "MLQ = " << LQ.M() << " & MAntiLQ = " << antiLQ.M() << std::endl << std::endl;
    //    if (LQ.M() > 0 && antiLQ.M() > 0){
    MLQLQbar_gen->Fill( mLQLQbar_gen,e.weight);   //}
    Pt_LQLQbar_gen->Fill ( ptLQLQbar_gen, e.weight);
    shat->Fill(sh, e.weight);
    
    //   if (LQ.M() > 0 && antiLQ.M() > 0){
    M_LQ->Fill( LQ.M(), e.weight);
    M_LQs->Fill( LQ.M(), e.weight);
    M_antiLQ->Fill(antiLQ.M(), e.weight);
    M_LQs->Fill(antiLQ.M(), e.weight);
    Pt_LQ->Fill( LQ.Pt(), e.weight);
    Pt_LQs->Fill( LQ.Pt(), e.weight);
    Pt_antiLQ->Fill(antiLQ.Pt(), e.weight);
    Pt_LQs->Fill( antiLQ.Pt(), e.weight); //}

    Pt_LQ_over_shat->Fill( LQ.Pt()/sh, e.weight);
    Pt_antiLQ_over_shat->Fill( antiLQ.Pt()/sh, e.weight);
    Pt_LQ_over_M_LQLQbar->Fill( LQ.Pt()/mLQLQbar_gen, e.weight);
    Pt_antiLQ_over_M_LQLQbar->Fill( antiLQ.Pt()/mLQLQbar_gen, e.weight);

    eta_LQ->Fill( LQ.eta(), e.weight);
    eta_LQs->Fill( LQ.eta(), e.weight);
    eta_antiLQ->Fill(antiLQ.eta(), e.weight);
    eta_LQs->Fill( antiLQ.eta(), e.weight);
    y_LQ->Fill( LQ.Rapidity(), e.weight);
    y_LQs->Fill( LQ.Rapidity(), e.weight);
    y_antiLQ->Fill(antiLQ.Rapidity(), e.weight);
    y_LQs->Fill( antiLQ.Rapidity(), e.weight);
    phi_LQ->Fill( LQ.phi(), e.weight);
    phi_LQs->Fill( LQ.phi(), e.weight);
    phi_antiLQ->Fill(antiLQ.phi(), e.weight);
    phi_LQs->Fill( antiLQ.phi(), e.weight);

    double difabseta = fabs( LQ.eta()) - fabs( antiLQ.eta());
    double difabsy = fabs( LQ.Rapidity()) - fabs( antiLQ.Rapidity());

    diffabseta->Fill(difabseta, e.weight);
    diffabsy->Fill(difabsy, e.weight);

    double deltaR_LQ = std::max (std::max( uhh2::deltaR(LQLQbargen.muLQ(), LQLQbargen.Topdecay1() ), 
				  uhh2::deltaR(LQLQbargen.muLQ(), LQLQbargen.Topdecay2() ) )
			     , uhh2::deltaR(LQLQbargen.Topdecay1(), LQLQbargen.Topdecay2() ) );

    double deltaR_antiLQ = std::max (std::max( uhh2::deltaR(LQLQbargen.muAntiLQ(), LQLQbargen.Antitopdecay1() ), 
				      uhh2::deltaR(LQLQbargen.muAntiLQ(), LQLQbargen.Antitopdecay2() ) )
				 , uhh2::deltaR(LQLQbargen.Antitopdecay1(), LQLQbargen.Antitopdecay2() ) );

    deltaR_LQ_decays->Fill(deltaR_LQ,e.weight);
    deltaR_antiLQ_decays->Fill(deltaR_antiLQ,e.weight);
    
    M_LQLQbar_vs_deltaR_LQ->Fill(mLQLQbar_gen, deltaR_LQ, e.weight);
    M_LQLQbar_vs_deltaR_antiLQ->Fill(mLQLQbar_gen, deltaR_antiLQ, e.weight);

    /* double deltaR_Top = uhh2::deltaR(LQLQbargen.Topdecay1(), LQLQbargen.Topdecay2());
    double deltaR_Antitop = uhh2::deltaR(LQLQbargen.Antitopdecay1(), LQLQbargen.Antitopdecay2());
    deltaR_Top_decays->Fill(deltaR_Top,e.weight);
    deltaR_Antitop_decays->Fill(deltaR_Antitop,e.weight);
    M_LQLQbar_vs_deltaR_Top->Fill(mLQLQbar_gen, deltaR_Top, e.weight);
    M_LQLQbar_vs_deltaR_Antitop->Fill(mLQLQbar_gen, deltaR_Antitop, e.weight);*/

    Pt_LQ_vs_deltaR_LQ->Fill(LQ.Pt(), deltaR_LQ, e.weight);
    Pt_antiLQ_vs_deltaR_antiLQ->Fill(antiLQ.Pt(), deltaR_antiLQ, e.weight);
    M_LQLQbar_vs_Pt_LQ->Fill(mLQLQbar_gen, LQ.Pt(), e.weight);
    M_LQLQbar_vs_Pt_antiLQ->Fill(mLQLQbar_gen, antiLQ.Pt(), e.weight);
    Pt_LQ_vs_Pt_antiLQ->Fill(LQ.Pt(), antiLQ.Pt(), e.weight);
    M_LQLQbar_vs_shat->Fill(mLQLQbar_gen,sh,e.weight);

    shat_vs_deltaR_LQ->Fill(sh, deltaR_LQ, e.weight);
    shat_vs_deltaR_antiLQ->Fill(sh, deltaR_antiLQ, e.weight);
    shat_vs_Pt_LQ->Fill(sh, LQ.Pt(),  e.weight);
    shat_vs_Pt_antiLQ->Fill(sh, antiLQ.Pt(),  e.weight);

    Pt_mu->Fill( LQLQbargen.muLQ().pt(), e.weight);
    Pt_mus->Fill( LQLQbargen.muLQ().pt(), e.weight);
    Pt_antimu->Fill( LQLQbargen.muAntiLQ().pt(), e.weight); 
    Pt_mus->Fill( LQLQbargen.muAntiLQ().pt(), e.weight); 
    eta_mu->Fill( LQLQbargen.muLQ().eta(), e.weight);
    eta_mus->Fill( LQLQbargen.muLQ().eta(), e.weight);
    eta_antimu->Fill( LQLQbargen.muAntiLQ().eta(), e.weight); 
    eta_mus->Fill( LQLQbargen.muAntiLQ().eta(), e.weight); 
    y_mu->Fill( LQLQbargen.muLQ().v4().Rapidity(), e.weight);
    y_mus->Fill( LQLQbargen.muLQ().v4().Rapidity(), e.weight);
    y_antimu->Fill( LQLQbargen.muAntiLQ().v4().Rapidity(), e.weight); 
    y_mus->Fill( LQLQbargen.muAntiLQ().v4().Rapidity(), e.weight); 
    phi_mu->Fill( LQLQbargen.muLQ().phi(), e.weight);
    phi_mus->Fill( LQLQbargen.muLQ().phi(), e.weight);
    phi_antimu->Fill( LQLQbargen.muAntiLQ().phi(), e.weight); 
    phi_mus->Fill( LQLQbargen.muAntiLQ().phi(), e.weight); 
    M_mu->Fill( LQLQbargen.muLQ().v4().M(), e.weight);
    M_mus->Fill( LQLQbargen.muLQ().v4().M(), e.weight);
    M_antimu->Fill( LQLQbargen.muAntiLQ().v4().M(), e.weight); 
    M_mus->Fill( LQLQbargen.muAntiLQ().v4().M(), e.weight); 

    Pt_top->Fill( LQLQbargen.TopLQ().pt(), e.weight);
    Pt_tops->Fill( LQLQbargen.TopLQ().pt(), e.weight);
    Pt_antitop->Fill( LQLQbargen.TopAntiLQ().pt(), e.weight); 
    Pt_tops->Fill( LQLQbargen.TopAntiLQ().pt(), e.weight); 
    eta_top->Fill( LQLQbargen.TopLQ().eta(), e.weight);
    eta_tops->Fill( LQLQbargen.TopLQ().eta(), e.weight);
    eta_antitop->Fill( LQLQbargen.TopAntiLQ().eta(), e.weight);
    eta_tops->Fill( LQLQbargen.TopAntiLQ().eta(), e.weight);
    y_top->Fill( LQLQbargen.TopLQ().v4().Rapidity(), e.weight);
    y_tops->Fill( LQLQbargen.TopLQ().v4().Rapidity(), e.weight);
    y_antitop->Fill( LQLQbargen.TopAntiLQ().v4().Rapidity(), e.weight); 
    y_tops->Fill( LQLQbargen.TopAntiLQ().v4().Rapidity(), e.weight); 
    phi_top->Fill( LQLQbargen.TopLQ().phi(), e.weight);
    phi_tops->Fill( LQLQbargen.TopLQ().phi(), e.weight);
    phi_antitop->Fill( LQLQbargen.TopAntiLQ().phi(), e.weight);
    phi_tops->Fill( LQLQbargen.TopAntiLQ().phi(), e.weight);
    if(LQLQbargen.TopLQ().v4().isTimelike()){
      M_top->Fill( LQLQbargen.TopLQ().v4().M(), e.weight);
      M_tops->Fill( LQLQbargen.TopLQ().v4().M(), e.weight);
    }
    if(LQLQbargen.TopAntiLQ().v4().isTimelike()){
      M_antitop->Fill( LQLQbargen.TopAntiLQ().v4().M(), e.weight); 
      M_tops->Fill( LQLQbargen.TopAntiLQ().v4().M(), e.weight); 
    }
    /*  TLorentzVector LQ_LQLQframe(0,0,0,0);
    LQ_LQLQframe.SetPtEtaPhiE(LQ.pt(), LQ.eta(), LQ.phi(), LQ.E());
    TLorentzVector antiLQ_LQLQframe(0,0,0,0);
    antiLQ_LQLQframe.SetPtEtaPhiE(antiLQ.pt(), antiLQ.eta(), antiLQ.phi(), antiLQ.E());
    TLorentzVector LQLQbar(0,0,0,0);
    LQLQbar.SetPtEtaPhiE((LQ+antiLQ).pt(), (LQ+antiLQ).eta(), (LQ+antiLQ).phi(), (LQ+antiLQ).E()); 
    
    LQ_LQLQframe.Boost(-LQLQbar.BoostVector());
    antiLQ_LQLQframe.Boost(-LQLQbar.BoostVector());

    cosThetastar_LQ_LQLQframe->Fill(cos(LQ_LQLQframe.Theta()) ,e.weight);
    cosThetastar_antiLQ_LQLQframe->Fill(cos(antiLQ_LQLQframe.Theta()) ,e.weight);
    Pt_LQ_LQLQframe->Fill(LQ_LQLQframe.Pt(),e.weight);
    Pt_antiLQ_LQLQframe->Fill(antiLQ_LQLQframe.Pt(),e.weight); 

    M_LQLQbar_vs_Pt_LQ_LQLQframe->Fill(mLQLQbar_gen, LQ_LQLQframe.Pt(),e.weight );
    M_LQLQbar_vs_Pt_antiLQ_LQLQframe->Fill(mLQLQbar_gen, antiLQ_LQLQframe.Pt(),e.weight ); */
    M_LQLQbar_vs_eta_LQ->Fill(mLQLQbar_gen, LQ.eta(),  e.weight);
    M_LQLQbar_vs_eta_antiLQ->Fill(mLQLQbar_gen, antiLQ.eta(),  e.weight);
}
