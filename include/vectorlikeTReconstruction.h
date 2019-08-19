#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/vectorlikeT/include/vectorlikeTReconstructionHypothesis.h"
#include "TMinuit.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TopJetIds.h"

typedef std::function< std::vector<LorentzVector>  (const LorentzVector & lepton, const LorentzVector & met)> NeutrinoReconstructionMethod;


class HighMassvectorlikeTReconstruction: public uhh2::AnalysisModule {
public:

  explicit HighMassvectorlikeTReconstruction(uhh2::Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, CSVBTag::wp wp_btag);
  virtual bool process(uhh2::Event & event) override;
  virtual ~HighMassvectorlikeTReconstruction();

private:
  NeutrinoReconstructionMethod m_neutrinofunction;
  uhh2::Event::Handle<std::vector<vectorlikeTReconstructionHypothesis>> h_recohyps;
  uhh2::Event::Handle<bool> h_is_tprime_reco;
  CSVBTag::wp wp_btag;
};

std::vector<LorentzVector> vectorlikeTNeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met);
