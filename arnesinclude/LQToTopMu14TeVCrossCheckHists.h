#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesisDiscriminators.h"
#include "UHH2/LQToTopMu/include/LQReconstructionHypothesis.h"
#include "UHH2/LQToTopMu/include/LQGen.h"


namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class LQToTopMu14TeVCrossCheckHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    LQToTopMu14TeVCrossCheckHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

  protected:
    //Handles related to 13->14TeV scaling
    uhh2::Event::Handle<double> h_x13_1, h_x13_2, h_x14_1, h_x14_2, h_Q, h_xf13_1, h_xf13_2, h_xf14_1, h_xf14_2, h_weight13, h_weight14, h_sf;
    uhh2::Event::Handle<int> h_f1, h_f2;

    bool is_mc;


    virtual ~LQToTopMu14TeVCrossCheckHists();
};

}
