#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesis.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ObjectIdUtils.h"



namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class LQToTopMuHFLeptonHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    LQToTopMuHFLeptonHists(uhh2::Context & ctx, const std::string & dirname, ElectronId & EleIdIso_, MuonId & MuIdIso_);

    virtual void fill(const uhh2::Event & ev) override;

  protected:
    bool is_mc;
    ElectronId EleIdIso;
    MuonId MuIdIso;


    virtual ~LQToTopMuHFLeptonHists();
};

}
