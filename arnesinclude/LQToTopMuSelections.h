#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/GenParticle.h"

#include <vector>

namespace uhh2examples {
    
/* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
 * below 20% of the average of the leading two jets, where the minimum deltaphi and
 * maximum third jet pt fraction can be changed in the constructor.
 * The jets are assumed to be sorted in pt.
 */

  class DijetSelection: public uhh2::Selection {
  public:
    DijetSelection(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float dphi_min, third_frac_max;
  };
  
  class HtSelection: public uhh2::Selection {
  public:
    explicit HtSelection(double ht_min=0., double ht_max=-1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  class InvMass2MuVeto: public uhh2::Selection {
  public:
    explicit InvMass2MuVeto(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class InvMass2MuVetoInverted: public uhh2::Selection {
  public:
    explicit InvMass2MuVetoInverted(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class PtLeadingMuonSelection: public uhh2::Selection {
  public:
    explicit PtLeadingMuonSelection(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class Pt2ndMuonSelection: public uhh2::Selection {
  public:
    explicit Pt2ndMuonSelection(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class PtLeadingJetSelection: public uhh2::Selection {
  public:
    explicit PtLeadingJetSelection(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class Pt2ndJetSelection: public uhh2::Selection {
  public:
    explicit Pt2ndJetSelection(double pt_min = 30., double m_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double pt_min, pt_max;
  };

  class PtRelMuJetSelection : public uhh2::Selection{
  public:
    explicit PtRelMuJetSelection(double ptrel_min = 0, double ptrel_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ptrel_min, ptrel_max;
  };

  class METSelection : public uhh2::Selection{
  public:
    explicit METSelection(double met_min = 0, double met_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double met_min, met_max;
  };

  class GenLvlZMuMuSelection : public uhh2::Selection{
  public: 
    explicit GenLvlZMuMuSelection();
    virtual bool passes(const uhh2::Event & event);
  };

  class GenLvlZEESelection : public uhh2::Selection{
  public: 
    explicit GenLvlZEESelection();
    virtual bool passes(const uhh2::Event & event);
  };

  class GenLvlTopDileptonSelection : public uhh2::Selection{
  public: 
    explicit GenLvlTopDileptonSelection();
    virtual bool passes(const uhh2::Event & event);
  };

  class PtRelMu1JetSelection : public uhh2::Selection{
  public:
    explicit PtRelMu1JetSelection(double ptrel_min = 0., double ptrel_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ptrel_min, ptrel_max;
  };

  class HTJetsSelection : public uhh2::Selection{
  public:
    explicit HTJetsSelection(double ht_min = 0., double ht_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  class HTLeptSelection : public uhh2::Selection{
  public:
    explicit HTLeptSelection(double ht_min = 0., double ht_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double ht_min, ht_max;
  };

  class EtaLeadingJetSelection : public uhh2::Selection{
  public:
    explicit EtaLeadingJetSelection(double eta_max = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double eta_max;
  };

  class InvMassMuEleVeto: public uhh2::Selection {
  public:
    explicit InvMassMuEleVeto(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class InvMassEleEleVeto: public uhh2::Selection {
  public:
    explicit InvMassEleEleVeto(double m_min = 81., double m_max = 101.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class InvMassEleEleSelection: public uhh2::Selection {
  public:
    explicit InvMassEleEleSelection(double m_min = 71., double m_max = 111.);
    virtual bool passes(const uhh2::Event & event);
  private:
    double m_min, m_max;
  };

  class dRLeptonJetSelection : public uhh2::Selection{
  public:
    explicit dRLeptonJetSelection(double dRmin = 0., double dRmax = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double dRmin, dRmax;
  };

  class NGenElectronSelection : public uhh2::Selection{
  public:
    explicit NGenElectronSelection(int n_min_ = 0., int n_max_ = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double n_min, n_max;
  };

 class NGenMuonSelection : public uhh2::Selection{
  public:
    explicit NGenMuonSelection(int n_min_ = 0., int n_max_ = -1);
    virtual bool passes(const uhh2::Event & event);
  private:
    double n_min, n_max;
  };

  class LQSemiLepMatchable: public uhh2::Selection{
  public:
    LQSemiLepMatchable();
    ~LQSemiLepMatchable(){};
    virtual bool passes(const uhh2::Event & event);
  private:
  };

  class TTbarSemiLepMatchable: public uhh2::Selection{
  public:
    TTbarSemiLepMatchable();
    ~TTbarSemiLepMatchable(){};
    virtual bool passes(const uhh2::Event & event);
  private:
  };

}
