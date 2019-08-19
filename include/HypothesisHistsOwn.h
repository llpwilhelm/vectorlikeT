#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
/*
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesis.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
*/

/** \brief Common histograms for reconstruction hypotheses
 *
 * hyps_name is the name of the reconstruction hypothesis collection, for instance "HighMassReconstruction"
 * discriminator_name is the name of the discriminator used to choose the best reconstruction hypothesis, for instance "Chi2"
 */
class HypothesisHistsOwn: public uhh2::Hists {
public:
    HypothesisHistsOwn(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F *Discriminator, *Discriminator_2, *Discriminator_3;
    TH1F *M_LQlep_rec, *M_LQhad_rec, *M_LQmax_rec, *M_LQmean_rec, *M_LQmean_rec_rebin, *M_LQ_rec_diff, *M_LQ_rec_diff_rel;
    TH1F *Pt_diff_LQ, *Pt_diff_LQ_rel, *Pt_diff_LQ_discriminator_cut, *Pt_diff_LQ_rel_discriminator_cut;
    TH1F *Pt_ratio, *Pt_ratio_discriminator_cut, *DeltaPhi_LQ_LQ, *DeltaPhi_LQ_LQ_discriminator_cut;
    TH1F *M_ttbar_rec, *M_toplep_rec, *M_tophad_rec;
    TH1F *Pt_toplep_rec, *Pt_tophad_rec, *Pt_ttbar_rec;
    TH1F *M_LQlep_rec_discriminator_cut, *M_LQhad_rec_discriminator_cut, *M_LQmax_rec_discriminator_cut, *M_LQmean_rec_discriminator_cut, *M_LQmean_rec_discriminator_cut_rebin, *M_LQ_rec_discriminator_cut_diff, *M_LQ_rec_discriminator_cut_diff_rel;
    TH1F *M_ttbar_rec_discriminator_cut, *M_toplep_rec_discriminator_cut, *M_tophad_rec_discriminator_cut;
    TH1F *Pt_toplep_rec_discriminator_cut, *Pt_tophad_rec_discriminator_cut, *Pt_ttbar_rec_discriminator_cut;

 //   uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    std::string m_discriminator_name;
    
};
