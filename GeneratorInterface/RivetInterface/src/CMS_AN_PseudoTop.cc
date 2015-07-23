#include "Rivet/Analysis.hh"
//#include "Rivet/AnalysisLoader.hh"
//#include "Rivet/Particle.fhh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

#include "GeneratorInterface/RivetInterface/interface/PseudoTop.hh"

namespace Rivet {

class CMS_AN_PseudoTop : public Analysis {
public:
  CMS_AN_PseudoTop() : Analysis("CMS_AN_PseudoTop") {
  }

  void init() {
    addProjection(PseudoTop(), "ttbar");

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
    const PseudoTop& ttbar = applyProjection<PseudoTop>(event, "ttbar");

    const FourMomentum& t1P4 = ttbar.t1().momentum();
    const FourMomentum& t2P4 = ttbar.t2().momentum();
    const double pt1 = std::max(t1P4.pT(), t2P4.pT());
    const double pt2 = std::min(t1P4.pT(), t2P4.pT());
    const double dPhi = deltaPhi(t1P4, t2P4);
    const FourMomentum ttP4 = t1P4+t2P4;
    const FourMomentum t1P4AtCM = LorentzTransform(-ttP4.boostVector()).transform(t1P4);

    if ( ttbar.mode() == PseudoTop::CH_SEMILEPTON ) {
      //const Particle lCand1 = ttbar.wDecays1()[0];
      //const Particle lCand2 = ttbar.wDecays2()[0];
      _h14_diffXSecTopSemiLepPartontopPt->fill(t1P4.pT(), weight);
      _h14_diffXSecTopSemiLepPartontopPt->fill(t2P4.pT(), weight);
      _h15_diffXSecTopSemiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _h16_diffXSecTopSemiLepPartontopY->fill(t1P4.rapidity(), weight);
      _h16_diffXSecTopSemiLepPartontopY->fill(t2P4.rapidity(), weight);
      _h17_diffXSecTopSemiLepPartonttbarDelPhi->fill(dPhi, weight);
      _h18_diffXSecTopSemiLepPartontopPtLead->fill(pt1, weight);
      _h19_diffXSecTopSemiLepPartontopPtSubLead->fill(pt2, weight);
      _h20_diffXSecTopSemiLepPartonttbarPt->fill(ttP4.pT(), weight);
      _h21_diffXSecTopSemiLepPartonttbarY->fill(ttP4.rapidity(), weight);
      _h22_diffXSecTopSemiLepPartonttbarMass->fill(ttP4.mass(), weight);
    }
    else if ( ttbar.mode() == PseudoTop::CH_FULLLEPTON ) {
      _h23_diffXSecTopDiLepPartontopPt->fill(t1P4.pT(), weight);
      _h23_diffXSecTopDiLepPartontopPt->fill(t2P4.pT(), weight);
      _h24_diffXSecTopDiLepPartontopPtTtbarSys->fill(t1P4AtCM.pT(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t1P4.rapidity(), weight);
      _h25_diffXSecTopDiLepPartontopY->fill(t2P4.rapidity(), weight);
      _h26_diffXSecTopDiLepPartonttbarDelPhi->fill(dPhi, weight);
      _h27_diffXSecTopDiLepPartontopPtLead->fill(pt1, weight);
      _h28_diffXSecTopDiLepPartontopPtSubLead->fill(pt2, weight);
      _h29_diffXSecTopDiLepPartonttbarPt->fill(ttP4.pT(), weight);
      _h30_diffXSecTopDiLepPartonttbarY->fill(ttP4.rapidity(), weight);
      _h31_diffXSecTopDiLepPartonttbarMass->fill(ttP4.mass(), weight);
    }
/*
    const FourMomentum& l1P4 = lCand1.momentum();
    const FourMomentum& l2P4 = lCand2.momentum();
    const FourMomentum& b1P4 = ttbar.b1().momentum();
    const FourMomentum& b2P4 = ttbar.b2().momentum();

    const double l1Pt = l1P4.pT(), l1Eta = l1P4.eta();
    const double l2Pt = l2P4.pT(), l2Eta = l2P4.eta();
    const double b1Pt = b1P4.pT(), b1Eta = b1P4.eta();
    const double b2Pt = b2P4.pT(), b2Eta = b2P4.eta();
    // Leptons and b jet observables are considered only within the particle level phase space
    if ( l1Pt > 20 and l2Pt > 20 and std::abs(l1Eta) < 2.4 and std::abs(l2Eta) < 2.4 and
         b1Pt > 30 and b2Pt > 30 and std::abs(b1Eta) < 2.4 and std::abs(b2Eta) < 2.4 ) {
      _h_lepton_pt->fill(l1Pt, weight);
      _h_lepton_pt->fill(l2Pt, weight);
      _h_lepton_eta->fill(l1Eta, weight);
      _h_lepton_eta->fill(l2Eta, weight);

      _h_bjet_pt->fill(b1Pt, weight);
      _h_bjet_pt->fill(b2Pt, weight);
      _h_bjet_eta->fill(b1Eta, weight);
      _h_bjet_eta->fill(b2Eta, weight);

      const FourMomentum dileptonP4 = l1P4+l2P4;
      _h_dilepton_mass->fill(dileptonP4.mass(), weight);
      _h_dilepton_pt->fill(dileptonP4.pT(), weight);

      const FourMomentum lb11P4 = l1P4 + b1P4;
      const FourMomentum lb22P4 = l2P4 + b2P4;
      const FourMomentum lb12P4 = l1P4 + b2P4;
      const FourMomentum lb21P4 = l2P4 + b1P4;
      _h_lb_mass->fill(lb11P4.mass(), weight);
      _h_lb_mass->fill(lb22P4.mass(), weight);
      _h_lb_mass->fill(lb12P4.mass(), weight);
      _h_lb_mass->fill(lb21P4.mass(), weight);

      const FourMomentum dijetP4 = b1P4 + b2P4;
      _h_dijet_mass->fill(dijetP4.mass(), weight);
      _h_dijet_pt->fill(dijetP4.pT(), weight);
    }
*/

  };

  void finalize() {
    // Correction functions for TOP-12-028 paper, (parton bin height)/(pseudotop bin height)
    const double ch14[] = {1.161660, 1.136480, 1.020996, 0.895649, 0.772136, 0.685911, 0.559711, 0.566430};
    const double ch16[] = {2.101211, 1.099831, 0.937698, 0.883005, 0.868135, 0.882153, 0.878180, 0.941096, 1.095958, 2.056497};
    const double ch20[] = {1.602612, 0.913407, 0.816876, 0.849766, 0.889415, 0.857082};
    const double ch21[] = {2.461665, 1.147150, 0.908031, 0.848166, 0.814687, 0.803214, 0.824948, 0.947269, 1.122359, 2.428979};
    const double ch22[] = {1.498358, 1.362128, 1.024490, 0.819021, 0.646227, 0.475925, 0.372441};

    const double ch23[] = {0.933825, 1.069645, 1.051336, 0.919932, 0.774565};
    const double ch25[] = {1.682022, 1.002849, 0.925246, 0.924734, 0.880097, 0.901330, 1.042041, 1.733911};
    const double ch29[] = {1.129278, 0.908123, 0.933110, 0.963850};
    const double ch30[] = {2.401265, 1.140515, 0.937143, 0.889803, 0.833903, 0.946386, 1.179555, 2.445021};
    const double ch31[] = {0.803342, 1.136017, 1.206834, 1.037619, 1.081579, 0.741247};

    applyCorrection(_h14_diffXSecTopSemiLepPartontopPt        , ch14);
    applyCorrection(_h16_diffXSecTopSemiLepPartontopY         , ch16);
    applyCorrection(_h20_diffXSecTopSemiLepPartonttbarPt      , ch20);
    applyCorrection(_h21_diffXSecTopSemiLepPartonttbarY       , ch21);
    applyCorrection(_h22_diffXSecTopSemiLepPartonttbarMass    , ch22);

    applyCorrection(_h23_diffXSecTopDiLepPartontopPt        , ch23);
    applyCorrection(_h25_diffXSecTopDiLepPartontopY         , ch25);
    applyCorrection(_h29_diffXSecTopDiLepPartonttbarPt      , ch29);
    applyCorrection(_h30_diffXSecTopDiLepPartonttbarY       , ch30);
    applyCorrection(_h31_diffXSecTopDiLepPartonttbarMass    , ch31);

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

  void applyCorrection(Histo1DPtr h, const double* cf) {
    std::vector<YODA::HistoBin1D>& bins = h->bins();
    for ( int i=0, n=bins.size(); i<n; ++i ) {
      const double s = cf[i];
      YODA::HistoBin1D& bin = bins[i];
      bin.scaleW(s);
    }
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

};

// This global object acts as a hook for the plugin system
AnalysisBuilder<CMS_AN_PseudoTop> plugin_CMS_AN_PseudoTop;

}

