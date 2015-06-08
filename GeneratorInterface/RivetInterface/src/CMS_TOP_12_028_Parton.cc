#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

namespace Rivet {

class CMS_TOP_12_028_Parton : public Analysis {
public:
  CMS_TOP_12_028_Parton() : Analysis("CMS_TOP_12_028_Parton") {
  }

  void init() {
    // Parton level top quarks
    PartonTTbarState ttbarState;
    addProjection(ttbarState, "ttbar");

    _h14_diffXSecTopSemiLepPartontopPt         = bookHisto1D("h14_diffXSecTopSemiLepPartontopPt"        );
    _h15_diffXSecTopSemiLepPartontopPtTtbarSys = bookHisto1D("h15_diffXSecTopSemiLepPartontopPtTtbarSys");
    _h16_diffXSecTopSemiLepPartontopY          = bookHisto1D("h16_diffXSecTopSemiLepPartontopY"         );
    _h17_diffXSecTopSemiLepPartonttbarDelPhi   = bookHisto1D("h17_diffXSecTopSemiLepPartonttbarDelPhi"  );
    _h18_diffXSecTopSemiLepPartontopPtLead     = bookHisto1D("h18_diffXSecTopSemiLepPartontopPtLead"    );
    _h19_diffXSecTopSemiLepPartontopPtSubLead  = bookHisto1D("h19_diffXSecTopSemiLepPartontopPtSubLead" );
    _h20_diffXSecTopSemiLepPartonttbarPt       = bookHisto1D("h20_diffXSecTopSemiLepPartonttbarPt"      );
    _h21_diffXSecTopSemiLepPartonttbarY        = bookHisto1D("h21_diffXSecTopSemiLepPartonttbarY"       );
    _h22_diffXSecTopSemiLepPartonttbarMass     = bookHisto1D("h22_diffXSecTopSemiLepPartonttbarMass"    );

    _h23_diffXSecTopDiLepPartontopPt           = bookHisto1D("h23_diffXSecTopDiLepPartontopPt"          );
    _h24_diffXSecTopDiLepPartontopPtTtbarSys   = bookHisto1D("h24_diffXSecTopDiLepPartontopPtTtbarSys"  );
    _h25_diffXSecTopDiLepPartontopY            = bookHisto1D("h25_diffXSecTopDiLepPartontopY"           );
    _h26_diffXSecTopDiLepPartonttbarDelPhi     = bookHisto1D("h26_diffXSecTopDiLepPartonttbarDelPhi"    );
    _h27_diffXSecTopDiLepPartontopPtLead       = bookHisto1D("h27_diffXSecTopDiLepPartontopPtLead"      );
    _h28_diffXSecTopDiLepPartontopPtSubLead    = bookHisto1D("h28_diffXSecTopDiLepPartontopPtSubLead"   );
    _h29_diffXSecTopDiLepPartonttbarPt         = bookHisto1D("h29_diffXSecTopDiLepPartonttbarPt"        );
    _h30_diffXSecTopDiLepPartonttbarY          = bookHisto1D("h30_diffXSecTopDiLepPartonttbarY"         );
    _h31_diffXSecTopDiLepPartonttbarMass       = bookHisto1D("h31_diffXSecTopDiLepPartonttbarMass"      );
  };

  void analyze(const Event& event) {
    const double weight = event.weight();

    // Get the parton level ttbar candidate
    const PartonTTbarState& ttbarState = applyProjection<PartonTTbarState>(event, "ttbar");
    const Particle& tCand1 = ttbarState.t1();
    const Particle& tCand2 = ttbarState.t2();

    // Do the anlaysis only for semileptonic and full leptonic channels.
    // Veto tau decays
    if ( ttbarState.mode() != PartonTTbarState::CH_SEMILEPTON and
         ttbarState.mode() != PartonTTbarState::CH_FULLLEPTON ) vetoEvent;
    if ( ttbarState.mode1() >= PartonTTbarState::CH_TAU_HADRON ||
         ttbarState.mode2() >= PartonTTbarState::CH_TAU_HADRON ) vetoEvent;

    // Fill top quarks they are defined in the parton level, full phase space
    const FourMomentum& t1P4 = tCand1.momentum();
    const FourMomentum& t2P4 = tCand2.momentum();
    const double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
    const FourMomentum ttbarP4 = t1P4+t2P4;
    const FourMomentum t1P4AtCM = LorentzTransform(-ttbarP4.boostVector()).transform(t1P4);
    const double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

    if ( ttbarState.mode() == PartonTTbarState::CH_SEMILEPTON ) {
      _h14_diffXSecTopSemiLepPartontopPt->fill(t1Pt, weight); 
      _h14_diffXSecTopSemiLepPartontopPt->fill(t2Pt, weight); 
      _h15_diffXSecTopSemiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight); 
      _h16_diffXSecTopSemiLepPartontopY->fill(t1P4.rapidity(), weight); 
      _h16_diffXSecTopSemiLepPartontopY->fill(t2P4.rapidity(), weight); 
      _h17_diffXSecTopSemiLepPartonttbarDelPhi->fill(dPhi, weight); 
      _h18_diffXSecTopSemiLepPartontopPtLead->fill(std::max(t1Pt, t2Pt), weight); 
      _h19_diffXSecTopSemiLepPartontopPtSubLead->fill(std::min(t1Pt, t2Pt), weight); 
      _h20_diffXSecTopSemiLepPartonttbarPt->fill(ttbarP4.pT(), weight); 
      _h21_diffXSecTopSemiLepPartonttbarY->fill(ttbarP4.rapidity(), weight); 
      _h22_diffXSecTopSemiLepPartonttbarMass->fill(ttbarP4.mass(), weight); 
    }
    else if ( ttbarState.mode() == PartonTTbarState::CH_FULLLEPTON ) {
      _h23_diffXSecTopDiLepPartontopPt->fill(t1Pt, weight);
      _h23_diffXSecTopDiLepPartontopPt->fill(t2Pt, weight);
      _h24_diffXSecTopDiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t1P4.rapidity(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t2P4.rapidity(), weight);
      _h26_diffXSecTopDiLepPartonttbarDelPhi->fill(dPhi, weight);
      _h27_diffXSecTopDiLepPartontopPtLead->fill(std::max(t1Pt, t2Pt), weight);
      _h28_diffXSecTopDiLepPartontopPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
      _h29_diffXSecTopDiLepPartonttbarPt->fill(ttbarP4.pT(), weight);
      _h30_diffXSecTopDiLepPartonttbarY->fill(ttbarP4.rapidity(), weight);
      _h31_diffXSecTopDiLepPartonttbarMass->fill(ttbarP4.mass(), weight);
    }
  };

  void finalize() {
    normalize(_h14_diffXSecTopSemiLepPartontopPt        );
    normalize(_h15_diffXSecTopSemiLepPartontopPtTtbarSys);
    normalize(_h16_diffXSecTopSemiLepPartontopY         );
    normalize(_h17_diffXSecTopSemiLepPartonttbarDelPhi  );
    normalize(_h18_diffXSecTopSemiLepPartontopPtLead    );
    normalize(_h19_diffXSecTopSemiLepPartontopPtSubLead );
    normalize(_h20_diffXSecTopSemiLepPartonttbarPt      );
    normalize(_h21_diffXSecTopSemiLepPartonttbarY       );
    normalize(_h22_diffXSecTopSemiLepPartonttbarMass    );

    normalize(_h23_diffXSecTopDiLepPartontopPt        );
    normalize(_h24_diffXSecTopDiLepPartontopPtTtbarSys);
    normalize(_h25_diffXSecTopDiLepPartontopY         );
    normalize(_h26_diffXSecTopDiLepPartonttbarDelPhi  );
    normalize(_h27_diffXSecTopDiLepPartontopPtLead    );
    normalize(_h28_diffXSecTopDiLepPartontopPtSubLead );
    normalize(_h29_diffXSecTopDiLepPartonttbarPt      );
    normalize(_h30_diffXSecTopDiLepPartonttbarY       );
    normalize(_h31_diffXSecTopDiLepPartonttbarMass    );
  };

private:
  Histo1DPtr _h14_diffXSecTopSemiLepPartontopPt        ;
  Histo1DPtr _h15_diffXSecTopSemiLepPartontopPtTtbarSys;
  Histo1DPtr _h16_diffXSecTopSemiLepPartontopY         ;
  Histo1DPtr _h17_diffXSecTopSemiLepPartonttbarDelPhi  ;
  Histo1DPtr _h18_diffXSecTopSemiLepPartontopPtLead    ;
  Histo1DPtr _h19_diffXSecTopSemiLepPartontopPtSubLead ;
  Histo1DPtr _h20_diffXSecTopSemiLepPartonttbarPt      ;
  Histo1DPtr _h21_diffXSecTopSemiLepPartonttbarY       ;
  Histo1DPtr _h22_diffXSecTopSemiLepPartonttbarMass    ;

  Histo1DPtr _h23_diffXSecTopDiLepPartontopPt        ;
  Histo1DPtr _h24_diffXSecTopDiLepPartontopPtTtbarSys;
  Histo1DPtr _h25_diffXSecTopDiLepPartontopY         ;
  Histo1DPtr _h26_diffXSecTopDiLepPartonttbarDelPhi  ;
  Histo1DPtr _h27_diffXSecTopDiLepPartontopPtLead    ;
  Histo1DPtr _h28_diffXSecTopDiLepPartontopPtSubLead ;
  Histo1DPtr _h29_diffXSecTopDiLepPartonttbarPt      ;
  Histo1DPtr _h30_diffXSecTopDiLepPartonttbarY       ;
  Histo1DPtr _h31_diffXSecTopDiLepPartonttbarMass    ;

private:
  class PartonTTbarState : public FinalState {
  public:
    enum TTbarMode { CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON };
    enum DecayMode { CH_HADRON = 0, CH_MUON, CH_ELECTRON, CH_TAU_HADRON, CH_TAU_MUON, CH_TAU_ELECTRON };

    /// @name Standard constructors and destructors.
    //@{

    /// The default constructor. May specify the minimum and maximum
    /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    PartonTTbarState(double mineta = -MAXDOUBLE,
                    double maxeta =  MAXDOUBLE,
                    double minpt = 0.0*GeV)
      : FinalState(mineta, maxeta, minpt)
    {
      setName("PartonTop");
    }

    /// Clone on the heap.
    virtual const Projection* clone() const {
      return new PartonTTbarState(*this);
    }

    //@}

  public:
    TTbarMode mode() const { 
      if ( _mode1 == CH_HADRON && _mode2 == CH_HADRON ) return CH_FULLHADRON;
      else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON ) return CH_FULLLEPTON;
      else return CH_SEMILEPTON;
    }
    DecayMode mode1() const { return _mode1; }
    DecayMode mode2() const { return _mode2; }

    Particle t1() const { return _t1; }
    Particle t2() const { return _t2; }
    Particle b1() const { return _b1; }
    Particle b2() const { return _b2; }
    ParticleVector wDecays1() const { return _wDecays1; }
    ParticleVector wDecays2() const { return _wDecays2; }

  protected:
    // Apply the projection to the event
    void project(const Event& e);

  private:
    DecayMode _mode1, _mode2;
    Particle _t1, _t2;
    Particle _b1, _b2;
    ParticleVector _wDecays1, _wDecays2;
  };

};

void CMS_TOP_12_028_Parton::PartonTTbarState::project(const Event& e) {
  _theParticles.clear();
  _wDecays1.clear();
  _wDecays2.clear();
  _mode1 = _mode2 = CH_HADRON; // Set default decay mode to full-hadronic

  const double ptmin = 0;
  const double etamin = -MAXDOUBLE, etamax = MAXDOUBLE;

  int nTop = 0;
  bool isTau1 = false, isTau2 = false;
  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    if ( status == 1 ) continue;

    const int pdgId = p->pdg_id();
    if ( PID::isHadron(pdgId) ) continue; // skip hadrons
    if ( pdgId == 22 ) continue; // skip photons

    if ( isZero(p->momentum().perp()) || p->momentum().perp() < ptmin ) continue;
    if ( !inRange(p->momentum().eta(), etamin, etamax) ) continue;

    // Avoid double counting by skipping if particle ID == parent ID
    std::vector<GenParticle*> pps;
    if ( abs(pdgId) == 6 and p->end_vertex() != 0 ) {
      pps = Rivet::particles(p->end_vertex(), HepMC::children);
    }
    else if ( abs(pdgId) != 6 and p->production_vertex() != 0 )
    {
      pps = Rivet::particles(p->production_vertex(), HepMC::parents);
    }
    else continue;

    bool isDuplicated = false;
    foreach (GenParticle* pp, pps) {
      if ( p != pp && p->pdg_id() == pp->pdg_id() ) {
        isDuplicated = true;
        break;
      }
    }
    if ( isDuplicated ) continue;

    // Build Rivet::Particle
    Particle rp(*p);

    // Skip particles from hadronization (and keep tau decay)
    if ( rp.fromDecay() and !rp.hasAncestor(15) and !rp.hasAncestor(-15) ) continue;

    if      ( pdgId ==  6 ) { nTop++; _t1 = Particle(rp); }
    else if ( pdgId == -6 ) { nTop++; _t2 = Particle(rp); }
    else if ( pdgId ==  5 ) _b1 = Particle(rp);
    else if ( pdgId == -5 ) _b2 = Particle(rp);
    else if ( pdgId != 24 && rp.hasAncestor( 24) ) {
      if ( pdgId == -15 ) isTau1 = true;
      else if ( pdgId == -11 ) _mode1 = CH_ELECTRON;
      else if ( pdgId == -13 ) _mode1 = CH_MUON;
      _wDecays1.push_back(rp);
    }
    else if ( pdgId != -24 && rp.hasAncestor(-24) ) {
      if ( pdgId == 15 ) isTau2 = true;
      else if ( pdgId == 11 ) _mode2 = CH_ELECTRON;
      else if ( pdgId == 13 ) _mode2 = CH_MUON;
      _wDecays2.push_back(rp);
    }
  }
  if ( isTau1 ) _mode1 = static_cast<DecayMode>(_mode1+3);
  if ( isTau2 ) _mode2 = static_cast<DecayMode>(_mode2+3);
}

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_TOP_12_028_Parton> plugin_CMS_TOP_12_028_Parton;

}
