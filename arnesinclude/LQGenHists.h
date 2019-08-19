#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/LQToTopMu/include/LQGen.h"

/** \brief Histograms for LQLQbar quantities on generator (parton) level
 * 
 * LQGen container has to be filled before calling this histogram class
 */
class LQGenHists: public uhh2::Hists {
public:
    LQGenHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F* MLQLQbar_gen, *Pt_LQLQbar_gen, *shat, *M_LQ, *M_antiLQ, *M_LQs, *Pt_LQ, *Pt_antiLQ, *Pt_LQs, *Pt_LQ_over_shat, *Pt_antiLQ_over_shat, *Pt_LQ_over_M_LQLQbar, *Pt_antiLQ_over_M_LQLQbar, *eta_LQ, *eta_antiLQ, *eta_LQs, *y_LQ, *y_antiLQ, *y_LQs, *phi_LQ, *phi_antiLQ, *phi_LQs, *diffabseta, *diffabsy, *deltaR_LQ_decays, *deltaR_antiLQ_decays, *deltaR_top_decays, *deltaR_antitop_decays, *Pt_mu, *Pt_antimu, *Pt_mus, *eta_mu, *eta_antimu, *eta_mus, *y_mu, *y_antimu, *y_mus, *phi_mu, *phi_antimu, *phi_mus, *M_mu, *M_antimu, *M_mus, *Pt_top, *Pt_antitop, *Pt_tops, *eta_top, *eta_antitop, *eta_tops, *y_top, *y_antitop, *y_tops, *M_top, *M_antitop, *M_tops, *phi_top, *phi_antitop, *phi_tops;

    TH2F* M_LQLQbar_vs_shat, *M_LQLQbar_vs_deltaR_LQ, *M_LQLQbar_vs_deltaR_antiLQ, *shat_vs_deltaR_LQ, *shat_vs_deltaR_antiLQ, *Pt_LQ_vs_deltaR_LQ, *Pt_antiLQ_vs_deltaR_antiLQ, *M_LQLQbar_vs_deltaR_top, *M_LQLQbar_vs_deltaR_antitop, *M_LQLQbar_vs_Pt_LQ, *M_LQLQbar_vs_Pt_antiLQ, *shat_vs_Pt_LQ, *shat_vs_Pt_antiLQ, *Pt_LQ_vs_Pt_antiLQ,/* *M_LQLQbar_vs_Pt_LQ_LQLQframe, *M_LQLQbar_vs_Pt_antiLQ_LQLQframe,*/ *M_LQLQbar_vs_eta_LQ, *M_LQLQbar_vs_eta_antiLQ;
    uhh2::Event::Handle<LQGen> h_LQLQbargen;
};
