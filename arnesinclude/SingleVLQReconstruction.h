#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/vectorlikeT/arnesinclude/SingleVLQReconstructionHypothesis.h"
#include "TMinuit.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TopJetIds.h"

typedef std::function< std::vector<LorentzVector>  (const LorentzVector & lepton, const LorentzVector & met)> NeutrinoReconstructionMethod;


class HighMassSingleVLQReconstruction: public uhh2::AnalysisModule {
public:

  explicit HighMassSingleVLQReconstruction(uhh2::Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, CSVBTag::wp wp_btag);
  virtual bool process(uhh2::Event & event) override;
  virtual ~HighMassSingleVLQReconstruction();

private:
  NeutrinoReconstructionMethod m_neutrinofunction;
  uhh2::Event::Handle<std::vector<SingleVLQReconstructionHypothesis>> h_recohyps;
  uhh2::Event::Handle<bool> h_is_tprime_reco;
  CSVBTag::wp wp_btag;
};

std::vector<LorentzVector> SingleVLQNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met);
