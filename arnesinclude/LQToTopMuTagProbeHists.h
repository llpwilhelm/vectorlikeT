#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesis.h"
#include "UHH2/LQToTopMu/include/LQGen.h"


namespace uhh2examples {
 
class LQToTopMuTagProbeHists: public uhh2::Hists {

  public:
    LQToTopMuTagProbeHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & ev) override;

  protected:
    virtual ~LQToTopMuTagProbeHists();

};

}
