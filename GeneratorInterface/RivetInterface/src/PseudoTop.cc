#include "GeneratorInterface/RivetInterface/interface/PseudoTop.hh"
//#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

using namespace Rivet;

struct GreaterByPt
{
  bool operator()(const Particle& a, const Particle& b) {
    return a.pt() > b.pt();
  }
};

void PseudoTop::cleanup(std::map<double, std::pair<size_t, size_t> >& v) const
{
  std::vector<std::map<double, std::pair<size_t, size_t> >::const_iterator> toErase;
  std::set<size_t> usedLeg1, usedLeg2;
  for ( auto key = v.begin(); key != v.end(); ++key )
  {
    const size_t leg1 = key->second.first;
    const size_t leg2 = key->second.second;
    if ( usedLeg1.find(leg1) == usedLeg1.end() and
         usedLeg2.find(leg2) == usedLeg2.end() )
    {
      usedLeg1.insert(leg1);
      usedLeg2.insert(leg2);
    }
    else
    {
      toErase.push_back(key);
    }
  }
  for ( auto& key : toErase ) v.erase(key);
}

void PseudoTop::project(const Event& e) {
  // Leptons : do the lepton clustering anti-kt R=0.1 using stable photons and leptons not from hadron decay
  // Neutrinos : neutrinos not from hadron decay
  // MET : vector sum of all invisible particles in x-y plane 
  // Jets : anti-kt R=0.4 using all particles excluding neutrinos and particles used in lepton clustering
  //        add ghost B hadrons during the jet clustering to identify B jets.

  // W->lv : dressed lepton and neutrino pairs
  // W->jj : light flavored dijet
  // W candidate : select lv or jj pairs which minimise |mW1-80.4|+|mW2-80.4|
  //               lepton-neutrino pair will be selected with higher priority

  // t->Wb : W candidate + b jet
  // t candidate : select Wb pairs which minimise |mtop1-172.5|+|mtop2-172.5|

  _isValid = false;
  _theParticles.clear();
  _wDecays1.clear();
  _wDecays2.clear();
  _mode1 = _mode2 = CH_HADRON;

  // Collect final state particles
  Particles pForLep, pForJet;
  Particles neutrinos; // Prompt neutrinos
  foreach (GenParticle* p, Rivet::particles(e.genEvent())) {
    const int status = p->status();
    const int pdgId = p->pdg_id();
    if ( status == 1 ) {
      Particle rp(*p);
      if ( !PID::isHadron(pdgId) && !rp.fromHadron() ) {
        // Collect particles not from hadron decay
        if ( rp.isNeutrino() ) {
          // Prompt neutrinos are kept in separate collection
          neutrinos.push_back(rp);
        }
        else if ( pdgId == 22 || rp.isLepton() ) {
          // Leptons and photons for the dressing
          pForLep.push_back(rp);
        }
      }
      else if ( !rp.isNeutrino() ) {
        // Use all particles from hadron decay
        pForJet.push_back(rp);
      }
    }
    else if ( PID::isHadron(pdgId) && PID::hasBottom(pdgId) ) {
      // NOTE: Consider B hadrons with pT > 5GeV - not in CMS proposal
      //if ( p->momentum().perp() < 5 ) continue; 

      // Do unstable particles, to be used in the ghost B clustering
      // Use last B hadrons only
      bool isLast = true;
      foreach ( GenParticle* pp, Rivet::particles(p->end_vertex(), HepMC::children) ) {
        if ( PID::hasBottom(pp->pdg_id()) ) {
          isLast = false;
          break;
        }
      }
      if ( !isLast ) continue;
      
      // Rescale momentum by 10^-20
      Particle ghost(pdgId, FourMomentum(p->momentum())*1e-20/p->momentum().rho());
      pForJet.push_back(ghost);
    }
  }

  // Start object building from trivial thing - prompt neutrinos
  std::sort(neutrinos.begin(), neutrinos.end(), GreaterByPt());

  // Proceed to lepton dressing
  FastJets fjLep(FastJets::ANTIKT, _lepR);
  fjLep.calc(pForLep);
  Jets leptons;
  std::vector<int> leptonsId;
  std::set<int> dressedIdxs;
  foreach ( const Jet& lep, fjLep.jetsByPt(_lepMinPt) ) {
    if ( std::abs(lep.eta()) > _lepMaxEta ) continue;

    double leadingPt = -1;
    int leptonId = 0;
    foreach ( const Particle& p, lep.particles() ) {
      dressedIdxs.insert(p.genParticle()->barcode());
      if ( p.isLepton() && p.momentum().pt() > leadingPt ) {
        leadingPt = p.momentum().pt();
        leptonId = p.pdgId();
      }
    }
    if ( leptonId == 0 ) continue;
    leptons.push_back(lep);
    leptonsId.push_back(leptonId);
  }

  // Re-use particles not used in lepton dressing
  foreach ( const Particle& rp, pForLep ) {
    const int barcode = rp.genParticle()->barcode();
    // Skip if the particle is used in dressing
    if ( dressedIdxs.find(barcode) != dressedIdxs.end() ) continue;

    // Put back to be used in jet clustering
    pForJet.push_back(rp);
  }

  // Then do the jet clustering
  FastJets fjJet(FastJets::ANTIKT, _jetR);
  //fjJet.useInvisibles(); // NOTE: CMS proposal to remove neutrinos
  fjJet.calc(pForJet);
  Jets bjets, ljets;
  foreach ( const Jet& jet, fjJet.jetsByPt(_jetMinPt) ) {
    if ( std::abs(jet.eta()) > _jetMaxEta ) continue;

    bool isBJet = false;
    foreach ( const Particle& rp, jet.particles() ) {
      if ( PID::hasBottom(rp.pdgId()) ) {
        isBJet = true;
        break;
      }
    }

    if ( isBJet ) bjets.push_back(jet);
    else ljets.push_back(jet);
  }

  // Every building blocks are ready. Continue to pseudo-W and pseudo-top combination

  if ( bjets.size() < 2 ) return; // Ignore single top for now
  std::map<double, std::pair<size_t, size_t> > wLepCandIdxs;
  std::map<double, std::pair<size_t, size_t> > wHadCandIdxs;

  // Collect leptonic-decaying W's
  for ( size_t iLep=0, nLep=leptons.size(); iLep<nLep; ++iLep )
  {
    const Jet& lep = leptons.at(iLep);
    for ( size_t iNu=0, nNu=neutrinos.size(); iNu<nNu; ++iNu ) {
      const Particle& nu = neutrinos.at(iNu);
      const double m = (lep.momentum()+nu.momentum()).mass();
      const double dm = std::abs(m-_wMass);
      wLepCandIdxs[dm] = make_pair(iLep, iNu);
    }
  }

  // Continue to hadronic decaying W's
  for ( size_t i=0, nLjet=ljets.size(); i<nLjet; ++i )
  {
    const Jet& ljet1 = ljets[i];
    for ( size_t j=i+1; j<nLjet; ++j )
    {
      const Jet& ljet2 = ljets[j];
      const double m = (ljet1.momentum()+ljet2.momentum()).mass();
      const double dm = std::abs(m-_wMass);
      wHadCandIdxs[dm] = make_pair(i, j);
    }
  }

  // Cleanup W candidate, choose pairs with minimum dm if they share decay products
  cleanup(wLepCandIdxs);
  cleanup(wHadCandIdxs);
  const size_t nWLepCand = wLepCandIdxs.size();
  const size_t nWHadCand = wHadCandIdxs.size();

  if ( nWLepCand + nWHadCand < 2 ) return; // We skip single top

  int w1Q = 1, w2Q = -1;
  int w1dau1Id = 1, w2dau1Id = 1;
  FourMomentum w1dau1LVec, w1dau2LVec;
  FourMomentum w2dau1LVec, w2dau2LVec;
  if ( nWLepCand == 0 ) // Full hadronic case
  {
    const auto& idPair1 = wHadCandIdxs.begin()->second;
    const auto& idPair2 = std::next(wHadCandIdxs.begin())->second;
    const auto& w1dau1 = ljets[idPair1.first];
    const auto& w1dau2 = ljets[idPair1.second];
    const auto& w2dau1 = ljets[idPair2.first];
    const auto& w2dau2 = ljets[idPair2.second];

    w1dau1LVec = w1dau1.momentum();
    w1dau2LVec = w1dau2.momentum();
    w2dau1LVec = w2dau1.momentum();
    w2dau2LVec = w2dau2.momentum();
  }
  else if ( nWLepCand == 1 ) // Semi-leptonic case
  {
    const auto& idPair1 = wLepCandIdxs.begin()->second;
    const auto& idPair2 = wHadCandIdxs.begin()->second;
    const auto& w1dau1 = leptons[idPair1.first];
    const auto& w1dau2 = neutrinos[idPair1.second];
    const auto& w2dau1 = ljets[idPair2.first];
    const auto& w2dau2 = ljets[idPair2.second];

    w1dau1LVec = w1dau1.momentum();
    w1dau2LVec = w1dau2.momentum();
    w2dau1LVec = w2dau1.momentum();
    w2dau2LVec = w2dau2.momentum();
    w1dau1Id = leptonsId[idPair1.first];
    w1Q = w1dau1Id > 0 ? -1 : 1;
    w2Q = -w1Q;

    switch ( w1dau1Id ) {
      case 13: case -13: _mode1 = CH_MUON; break;
      case 11: case -11: _mode1 = CH_ELECTRON; break;
    }
  }
  else // Full leptonic case
  {
    const auto& idPair1 = wLepCandIdxs.begin()->second;
    const auto& idPair2 = std::next(wLepCandIdxs.begin())->second;
    const auto& w1dau1 = leptons[idPair1.first];
    const auto& w1dau2 = neutrinos[idPair1.second];
    const auto& w2dau1 = leptons[idPair2.first];
    const auto& w2dau2 = neutrinos[idPair2.second];

    w1dau1LVec = w1dau1.momentum();
    w1dau2LVec = w1dau2.momentum();
    w2dau1LVec = w2dau1.momentum();
    w2dau2LVec = w2dau2.momentum();
    w1dau1Id = leptonsId[idPair1.first];
    w2dau1Id = leptonsId[idPair2.first];
    w1Q = w1dau1Id > 0 ? -1 : 1;
    w2Q = w2dau1Id > 0 ? -1 : 1;

    switch ( w1dau1Id ) {
      case 13: case -13: _mode1 = CH_MUON; break;
      case 11: case -11: _mode1 = CH_ELECTRON; break;
    }
    switch ( w2dau1Id ) {
      case 13: case -13: _mode2 = CH_MUON; break;
      case 11: case -11: _mode2 = CH_ELECTRON; break;
    }
  }
  const auto w1LVec = w1dau1LVec+w1dau2LVec;
  const auto w2LVec = w2dau1LVec+w2dau2LVec;

  // Combine b jets
  double sumDm = 1e9;
  FourMomentum b1LVec, b2LVec;
  for ( size_t i=0, n=bjets.size(); i<n; ++i ) {
    const Jet& bjet1 = bjets[i];
    const double mtop1 = (w1LVec+bjet1.momentum()).mass();
    const double dmtop1 = std::abs(mtop1-_tMass);
    for ( size_t j=0; j<n; ++j ) {
      if ( i == j ) continue;
      const Jet& bjet2 = bjets[j];
      const double mtop2 = (w2LVec+bjet2.momentum()).mass();
      const double dmtop2 = std::abs(mtop2-_tMass);

      if ( sumDm <= dmtop1+dmtop2 ) continue;

      sumDm = dmtop1+dmtop2;
      b1LVec = bjet1.momentum();
      b2LVec = bjet2.momentum();
    }
  }
  if ( sumDm >= 1e9 ) return; // Failed to make top, but this should not happen.

  const auto t1LVec = w1LVec + b1LVec;
  const auto t2LVec = w2LVec + b2LVec;

  // Put all of them into candidate collection
  _t1 = Particle(w1Q*6, t1LVec);
  _b1 = Particle(w1Q*5, b1LVec);
  _wDecays1.push_back(Particle(w1dau1Id, w1dau1LVec));
  _wDecays1.push_back(Particle(-w1dau1Id+w1Q, w1dau2LVec));

  _t2 = Particle(w2Q*6, t2LVec);
  _b2 = Particle(w2Q*5, b2LVec);
  _wDecays2.push_back(Particle(w2dau1Id, w2dau1LVec));
  _wDecays2.push_back(Particle(-w2dau1Id+w2Q, w2dau2LVec));

  _isValid = true;
}

