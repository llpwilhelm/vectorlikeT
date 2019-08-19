#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/LQToTopMu/include/LQGen.h"
//#include "UHH2/common/include/LQGen.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/LQToTopMu/include/LQGenHists.h"
//#include "UHH2/common/include/LQGenHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Example for calculating and accessing the TTbarGen interpretation
 * 
 */
class LQGenModule: public AnalysisModule {
public:
    
    explicit LQGenModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
  //std::unique_ptr<AnalysisModule> printer;
  std::unique_ptr<AnalysisModule> LQgenprod;
  std::unique_ptr<Hists> h_LQgenhists;
  Event::Handle<LQGen> h_LQLQbargen;
};


LQGenModule::LQGenModule(Context & ctx){


    //printer.reset(new GenParticlesPrinter(ctx));
  LQgenprod.reset(new LQGenProducer(ctx, "LQLQbargen", false));
  h_LQLQbargen = ctx.get_handle<LQGen>("LQLQbargen");
  h_LQgenhists.reset(new LQGenHists(ctx, "LQgenhists"));
}


bool LQGenModule::process(Event & event) {
  //printer->process(event);
  LQgenprod->process(event);

  const auto & LQLQbargen = event.get(h_LQLQbargen);
    
  h_LQgenhists->fill(event);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(LQGenModule)

}
