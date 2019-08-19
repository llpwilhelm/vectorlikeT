#include "UHH2/vectorlikeT/include/vectorlikeTReconstructionHypothesisDiscriminators.h"
#include "UHH2/core/include/Utils.h"

#include <set>

using namespace uhh2;
using namespace std;

namespace {

  // invariant mass of a lorentzVector, but save for timelike / spacelike vectors
  float inv_mass(const LorentzVector & p4){
    if(p4.isTimelike()) return p4.mass();
    else return -sqrt(-p4.mass2());
  }
}


const vectorlikeTReconstructionHypothesis * get_best_hypothesis(const std::vector<vectorlikeTReconstructionHypothesis> & hyps, const std::string & label){
  const vectorlikeTReconstructionHypothesis * best = nullptr;
  float current_best_disc = numeric_limits<float>::infinity();
  for(const auto & hyp : hyps){
    if(!hyp.has_discriminator(label)) continue;
    auto disc = hyp.discriminator(label);
    if(disc < current_best_disc){
      best = &hyp;
      current_best_disc = disc;
    }
  }
  if(std::isfinite(current_best_disc)){
    return best;
  }
  else{
    return nullptr;
  }
}

vectorlikeTChi2Discriminator::vectorlikeTChi2Discriminator(Context & ctx, const cfg & config_): config(config_){
  h_hyps = ctx.get_handle<vector<vectorlikeTReconstructionHypothesis>>("TprimeHypotheses");
}


bool vectorlikeTChi2Discriminator::process(uhh2::Event & event){
  auto & hyps = event.get(h_hyps);

  const double mass_tlep = 173.;
  const double mass_tlep_sigma = 22.;
  const double mass_higgs = 121.;
  const double mass_higgs_sigma = 13.;


  for(auto & hyp: hyps){
    double mass_tlep_rec = inv_mass(hyp.toplep_v4());
    double mass_higgs_rec = inv_mass(hyp.higgs_v4());
    double chi2_top = pow((mass_tlep_rec - mass_tlep)/(mass_tlep_sigma),2);
    double chi2_higgs = pow((mass_higgs_rec - mass_higgs)/(mass_higgs_sigma),2);
    hyp.set_discriminator(config.discriminator_label + "_top", chi2_top);
    hyp.set_discriminator(config.discriminator_label + "_higgs", chi2_higgs);
    hyp.set_discriminator(config.discriminator_label, chi2_higgs + chi2_top); // modified
  }
  return true;
}
