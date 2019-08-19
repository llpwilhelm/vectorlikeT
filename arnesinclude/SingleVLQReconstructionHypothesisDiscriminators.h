#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/LQToTopMu/include/SingleVLQReconstructionHypothesis.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQGen.h"

const SingleVLQReconstructionHypothesis * get_best_hypothesis(const std::vector<SingleVLQReconstructionHypothesis> & hyps, const std::string & label);


class SingleVLQChi2Discriminator: public uhh2::AnalysisModule {
public:
    struct cfg {
        std::string discriminator_label;
        cfg(): discriminator_label("Chi2"){}
    };

    SingleVLQChi2Discriminator(uhh2::Context & ctx, const cfg & config = cfg());
    virtual bool process(uhh2::Event & event) override;

private:
    uhh2::Event::Handle<std::vector<SingleVLQReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<bool> h_is_tprime_reco;
    cfg config;
};
