#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class LQToTopMuEfficiencyHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    LQToTopMuEfficiencyHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;



    virtual ~LQToTopMuEfficiencyHists();

  private:
    bool is_mc = false;
    double bins_from350[21] = {0,175,350,525,700,875,1050,1225,1400,1575,1750,1925,2100,2275,2450,2625,2800,2975,3325,3675,4200}; //same binning as _from350 up to 2975, then two double-size and one triple-size bin

};

}
