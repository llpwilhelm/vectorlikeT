#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/vectorlikeT/arnesinclude/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/vectorlikeT/arnesinclude/LQReconstructionHypothesis.h"
#include "UHH2/vectorlikeT/arnesinclude/LQGen.h"


namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class LQToTopMuMLQRecoHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    LQToTopMuMLQRecoHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

  protected:
    uhh2::Event::Handle<std::vector<LQReconstructionHypothesis>> h_inclusive_hyps;
    std::string m_discriminator_name;
    bool is_mc;


    virtual ~LQToTopMuMLQRecoHists();
};

}
