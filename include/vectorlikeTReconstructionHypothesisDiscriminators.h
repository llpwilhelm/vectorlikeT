#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/vectorlikeT/include/vectorlikeTReconstructionHypothesis.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/vectorlikeT/arnesinclude/LQGen.h"

const vectorlikeTReconstructionHypothesis * get_best_hypothesis(const std::vector<vectorlikeTReconstructionHypothesis> & hyps, const std::string & label);


class vectorlikeTChi2Discriminator: public uhh2::AnalysisModule {
public:
    struct cfg {
        std::string discriminator_label;
        cfg(): discriminator_label("Chi2"){}
    };

    vectorlikeTChi2Discriminator(uhh2::Context & ctx, const cfg & config = cfg());
    virtual bool process(uhh2::Event & event) override;

private:
    uhh2::Event::Handle<std::vector<vectorlikeTReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<bool> h_is_tprime_reco;
    cfg config;
};
