#ifndef RIVET_PseudoTop_HH
#define RIVET_PseudoTop_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Event.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  // @brief Pseudo top finder
  // 
  // Find top quark in the particle level.
  // The definition is based on the agreement at the LHC working group.
  class PseudoTop : public FinalState {
  public:
    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    PseudoTop(double minpt = 0.0*GeV,
              double lepR = 0.1, double lepMinPt = 20, double lepMaxEta = 2.4,
              double jetR = 0.4, double jetMinPt = 30, double jetMaxEta = 2.4)
      : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
        _lepR(lepR), _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta),
        _jetR(jetR), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta)
    {
      setName("PseudoTop");
    }

    enum TTbarMode { CH_NONE=-1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
    enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

    TTbarMode mode() const {
      if ( !_isValid ) return CH_NONE;
      if ( _mode1 == CH_HADRON && _mode2 == CH_HADRON ) return CH_FULLHADRON;
      else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON ) return CH_FULLLEPTON;
      else return CH_SEMILEPTON;
    }
    DecayMode mode1() const { return _mode1; }
    DecayMode mode2() const { return _mode2; }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new PseudoTop(*this);
    }

    //@}

  public:
    Particle t1() const { return _t1; }
    Particle t2() const { return _t2; }
    Particle b1() const { return _b1; }
    Particle b2() const { return _b2; }
    ParticleVector wDecays1() const { return _wDecays1; }
    ParticleVector wDecays2() const { return _wDecays2; }

  protected:
    // Apply the projection to the event
    void project(const Event& e) override;
    void cleanup(std::map<double, std::pair<size_t, size_t> >& v) const;

  private:
    const double _lepR, _lepMinPt, _lepMaxEta;
    const double _jetR, _jetMinPt, _jetMaxEta;

    constexpr static double _tMass = 172.5;
    constexpr static double _wMass = 80.4;

  private:
    bool _isValid;
    DecayMode _mode1, _mode2;

    Particle _t1, _t2;
    Particle _b1, _b2;
    ParticleVector _wDecays1, _wDecays2;

  };

}

#endif
