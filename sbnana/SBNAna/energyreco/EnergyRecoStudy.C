// Make a plot with cuts
#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

// ---- Constants ----
double mass_muon=0.1056583745; //GeV
double mass_proton=0.938272081; //GeV
double mass_pion=0.13957039; //GeV
double mass_pion_neut=0.1349768; //GeV

// ---- VARS -----
// Muons: energy measured by calorimeters
const Var kErecMuon([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  int bestplane=0;
  float kinetic=0;
  float energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen)
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      energy = kinetic + mass_muon;}

  return energy;});

const Var kEresMuon([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  int bestplane=0;
  float kinetic=0;
  float residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen)
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      residual = (kinetic+mass_muon) - trk.truth.p.genE;}

  return residual;});

// True muons: energy measured by calorimeters
const Var kErecTrueMuon([](const caf::SRSliceProxy* slc) -> float {
  int bestplane=0;
  float kinetic=0;
  float energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13)
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      energy = kinetic + mass_muon;}

  return energy;});

const Var kEresTrueMuon([](const caf::SRSliceProxy* slc) -> float {
  int bestplane=0;
  float kinetic=0;
  float residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13)
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      residual = (kinetic+mass_muon) - trk.truth.p.genE;}

  return residual;});

// Muons: energy calculated from momentum
const Var kErecMuon2([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  float energy=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen)
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));}

  return energy;});

const Var kEresMuon2([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  float energy=0;
  float residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen)
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));
      residual = energy - trk.truth.p.genE;}

  return residual;});

// True muons: energy calculated from momentum
const Var kErecTrueMuon2([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13)
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));}

  return energy;});

const Var kEresTrueMuon2([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  float residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13)
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));
      residual = energy - trk.truth.p.genE;}

  return residual;});

// Muons: true energy 
const Var kTrueErecMuon([](const caf::SRSliceProxy* slc) -> float {
  float energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13) {
      energy = trk.truth.p.genE;}}

  return energy;});

// Protons: KE measured by calorimeters
const MultiVar kErecProton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  int bestplane=0;
  double kinetic=0;
  std::vector<double> energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton < chi2_pion) {
          bestplane = trk.bestplane;
          if (bestplane == 0)
            kinetic = trk.calo0.ke/1000;
          if (bestplane == 1)
            kinetic = trk.calo1.ke/1000;
          if (bestplane == 2)
            kinetic = trk.calo2.ke/1000;
          energy.push_back(kinetic);}}}

  return energy;});

const MultiVar kEresProton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  int bestplane=0;
  float kinetic=0;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton < chi2_pion) {
          bestplane = trk.bestplane;
          if (bestplane == 0)
            kinetic = trk.calo0.ke/1000;
          if (bestplane == 1)
            kinetic = trk.calo1.ke/1000;
          if (bestplane == 2)
            kinetic = trk.calo2.ke/1000;
          residual.push_back(kinetic - trk.truth.p.genT);}}}

  return residual;});

// True protons: KE measured by calorimeters
const MultiVar kErecTrueProton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  int bestplane=0;
  double kinetic=0;
  std::vector<double> energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      energy.push_back(kinetic);}}

  return energy;});

const MultiVar kEresTrueProton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  int bestplane=0;
  float kinetic=0;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      residual.push_back(kinetic - trk.truth.p.genT);}}

  return residual;});

// Protons: KE calculated from momentum
const MultiVar kErecProton2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  double energy=0;
  std::vector<double> kinetic;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton < chi2_pion) {
          energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
          kinetic.push_back( energy - mass_proton );}}}

  return kinetic;});

const MultiVar kEresProton2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  double energy=0;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton < chi2_pion) {
          energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
          residual.push_back( energy - mass_proton - trk.truth.p.genT );}}}

  return residual;});

// True protons: KE calculated from momentum
const MultiVar kErecTrueProton2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double energy=0;
  std::vector<double> kinetic;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
      kinetic.push_back( energy - mass_proton );}}

  return kinetic;});

const MultiVar kEresTrueProton2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double energy=0;
  float kinetic;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
      kinetic = energy - mass_proton;
      residual.push_back(kinetic - trk.truth.p.genT);}}

  return residual;});

// Protons: true KE
const MultiVar kTrueErecProton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      energy.push_back(trk.truth.p.genT);}}

  return energy;});

// Proton momentum
//const MultiVar kPproton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  //double maxlen=0;
  //for (auto const& trk : slc->reco.trk) {
    //if (trk.len > maxlen) {
      //maxlen = trk.len;}}

  //double chi2_proton;
  //double chi2_pion;
  //std::vector<double> momentum;
  //for (auto const& trk : slc->reco.trk) {
    //if (trk.len < maxlen) {
      //chi2_proton = trk.chi2pid2.chi2_proton;
      //chi2_pion = trk.chi2pid2.chi2_pion;
        //if (chi2_proton < chi2_pion) {
          //momentum.push_back(trk.rangeP.p_proton);}}}

  //return momentum;});

// Pions: KE measured by calorimeters
const MultiVar kErecPion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  int bestplane=0;
  double kinetic=0;
  std::vector<double> energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton > chi2_pion) {
          bestplane = trk.bestplane;
          if (bestplane == 0)
            kinetic = trk.calo0.ke/1000;
          if (bestplane == 1)
            kinetic = trk.calo1.ke/1000;
          if (bestplane == 2)
            kinetic = trk.calo2.ke/1000;
          energy.push_back(kinetic);}}}

  return energy;});

const MultiVar kEresPion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  int bestplane=0;
  double kinetic=0;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton > chi2_pion) {
          bestplane = trk.bestplane;
          if (bestplane == 0)
            kinetic = trk.calo0.ke/1000;
          if (bestplane == 1)
            kinetic = trk.calo1.ke/1000;
          if (bestplane == 2)
            kinetic = trk.calo2.ke/1000;
          residual.push_back(kinetic - trk.truth.p.genT);}}}

  return residual;});

// True pions: KE measured by calorimeters
const MultiVar kErecTruePion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  int bestplane=0;
  double kinetic=0;
  std::vector<double> energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      energy.push_back(kinetic);}}

  return energy;});

const MultiVar kEresTruePion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  int bestplane=0;
  float kinetic=0;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      bestplane = trk.bestplane;
      if (bestplane == 0)
        kinetic = trk.calo0.ke/1000;
      if (bestplane == 1)
        kinetic = trk.calo1.ke/1000;
      if (bestplane == 2)
        kinetic = trk.calo2.ke/1000;
      residual.push_back(kinetic - trk.truth.p.genT);}}

  return residual;});

// Pions: KE calcuated from momentum
const MultiVar kErecPion2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  double energy=0;
  std::vector<double> kinetic;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton > chi2_pion) {
          energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
          kinetic.push_back( energy - mass_pion );}}}

  return kinetic;});

const MultiVar kEresPion2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen) {
      maxlen = trk.len;}}

  double chi2_proton;
  double chi2_pion;
  double energy=0;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton > chi2_pion) {
          energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
          residual.push_back( energy - mass_pion - trk.truth.p.genT );}}}

  return residual;});

// True pions: KE calcuated from momentum
const MultiVar kErecTruePion2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double energy=0;
  std::vector<double> kinetic;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
      kinetic.push_back( energy - mass_pion );}}

  return kinetic;});

const MultiVar kEresTruePion2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  double energy=0;
  float kinetic;
  std::vector<double> residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
      kinetic = energy - mass_pion;
      residual.push_back(kinetic - trk.truth.p.genT);}}

  return residual;});

// Pions: true KE 
const MultiVar kTrueErecPion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> energy;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      energy.push_back(trk.truth.p.genT);}}

  return energy;});

// Pion  momentum
//const MultiVar kPpion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  //double maxlen=0;
  //for (auto const& trk : slc->reco.trk) {
    //if (trk.len > maxlen) {
      //maxlen = trk.len;}}

  //double chi2_proton;
  //double chi2_pion;
  //std::vector<double> momentum;
  //for (auto const& trk : slc->reco.trk) {
    //if (trk.len < maxlen) {
      //chi2_proton = trk.chi2pid2.chi2_proton;
      //chi2_pion = trk.chi2pid2.chi2_pion;
        //if (chi2_proton > chi2_pion) {
          //momentum.push_back(trk.rangeP.p_pion);}}}

  //return momentum;});

// Neutral pion energy
const Var kErecPionNeut([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  for (auto const& shw : slc->reco.shw) {
    energy += shw.bestplane_energy/1000;}

  return energy;});

const Var kEresPionNeut([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  float true_energy=0;
  float residual=0;
  for (auto const& shw : slc->reco.shw) {
    energy += shw.bestplane_energy/1000;
    true_energy += shw.truth.p.genE;}
 
  residual = energy - true_energy;
  return residual;});

// True neutral pion energy
const Var kErecTruePionNeut([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++)
      if (shw.truth.p.daughters[i]==111)
        energy += shw.bestplane_energy/1000;
        break;}

  return energy;});

const Var kEresTruePionNeut([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  float true_energy=0;
  float residual=0;
  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++) {
      if (shw.truth.p.daughters[i]==111) {
        energy += shw.bestplane_energy/1000;
        true_energy += shw.truth.p.genE;
        break;}}}

  residual = energy - true_energy;
  return residual;});

// Neutral pion true energy
const MultiVar kTrueErecPionNeut([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> energy;
  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++) {
      if (shw.truth.p.daughters[i]==111)
        energy.push_back(shw.truth.p.genE);
        break;}}

  return energy;});

// Total energy (calorimeter)
const Var kErec([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  int bestplane=0;
  float kinetic=0;
  float energy=0;
  for (auto const& trk : slc->reco.trk) {
    bestplane = trk.bestplane;
    if (bestplane == 0)
      kinetic = trk.calo0.ke/1000;
    if (bestplane == 1)
      kinetic = trk.calo1.ke/1000;
    if (bestplane == 2)
      kinetic = trk.calo2.ke/1000;
    if (trk.len == maxlen) {
      energy += (kinetic + mass_muon);}
    if (trk.len < maxlen) {
      energy += kinetic;}}

  for (auto const& shw : slc->reco.shw) {
    energy += shw.bestplane_energy/1000;}

  return energy;});

const Var kEres([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  int bestplane=0;
  float true_energy=0;
  float kinetic=0;
  float energy=0;
  float residual=0;
  double chi2_proton;
  double chi2_pion;
  for (auto const& trk : slc->reco.trk) {
    bestplane = trk.bestplane;
    if (bestplane == 0)
      kinetic = trk.calo0.ke/1000;
    if (bestplane == 1)
      kinetic = trk.calo1.ke/1000;
    if (bestplane == 2)
      kinetic = trk.calo2.ke/1000;
    if (trk.len == maxlen) {
      true_energy += trk.truth.p.genE;
      energy += (kinetic + mass_muon);}
    if (trk.len < maxlen) {
      true_energy += trk.truth.p.genT;
      energy += kinetic;}}

  for (auto const& shw : slc->reco.shw) {
    true_energy += shw.truth.p.genE;
    energy += shw.bestplane_energy/1000;} 

  residual = energy - true_energy;
  return residual;});

// Total energy for true particles (calorimeter)
const Var kErecTrue([](const caf::SRSliceProxy* slc) -> float {
  int bestplane=0;
  float kinetic=0;
  float energy=0;
  for (auto const& trk : slc->reco.trk) {
    bestplane = trk.bestplane;
    if (bestplane == 0)
      kinetic = trk.calo0.ke/1000;
    if (bestplane == 1)
      kinetic = trk.calo1.ke/1000;
    if (bestplane == 2)
      kinetic = trk.calo2.ke/1000;
    if (trk.truth.p.pdg == 13) {
      energy += (kinetic + mass_muon);}
    else {
      energy += kinetic;}}

  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++) {
      if (shw.truth.p.daughters[i]==111) {
        energy += shw.bestplane_energy/1000;
        break;}}}

  return energy;});

const Var kEresTrue([](const caf::SRSliceProxy* slc) -> float {
  int bestplane=0;
  float true_energy=0;
  float kinetic=0;
  float energy=0;
  float residual=0;
  for (auto const& trk : slc->reco.trk) {
    bestplane = trk.bestplane;
    if (bestplane == 0)
      kinetic = trk.calo0.ke/1000;
    if (bestplane == 1)
      kinetic = trk.calo1.ke/1000;
    if (bestplane == 2)
      kinetic = trk.calo2.ke/1000;
    if (trk.truth.p.pdg == 13) {
      true_energy += trk.truth.p.genE;
      energy += (kinetic + mass_muon);}
    else {
      true_energy += trk.truth.p.genT;
      energy += kinetic;}}

  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++) {
      if (shw.truth.p.daughters[i]==111) {
        true_energy += shw.truth.p.genE;
        energy += shw.bestplane_energy/1000;
        break;}}}

  residual = energy - true_energy;
  return residual;});

// Total energy (momentum)
const Var kErec2([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  float energy;
  float total=0;
  double chi2_proton;
  double chi2_pion;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen) {
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));
      total += energy;} 
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
        if (chi2_proton < chi2_pion) {
          energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
          total += (energy - mass_proton);}
        if (chi2_proton > chi2_pion) {
          energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
          total += (energy - mass_pion);}}}     

  for (auto const& shw : slc->reco.shw) {
    total += shw.bestplane_energy/1000;} 

  return total;});

const Var kEres2([](const caf::SRSliceProxy* slc) -> float {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  float true_energy=0;
  float kinetic;
  float energy=0;
  float residual=0;
  double chi2_proton;
  double chi2_pion;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen) {
      true_energy += trk.truth.p.genE;
      energy += sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2)); 
    if (trk.len < maxlen) {
      chi2_proton = trk.chi2pid2.chi2_proton;
      chi2_pion = trk.chi2pid2.chi2_pion;
      if (chi2_proton < chi2_pion) {
        true_energy += trk.truth.p.genT;
        kinetic = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2)); 
        energy += kinetic;}
      if (chi2_proton > chi2_pion) {
        true_energy += trk.truth.p.genT;
        kinetic = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
        energy += kinetic;}}}}

  for (auto const& shw : slc->reco.shw) {
    true_energy += shw.truth.p.genE;
    energy += shw.bestplane_energy/1000;} 

  residual = energy - true_energy;
  return residual;});

// Total energy for true particles (momentum)
const Var kErecTrue2([](const caf::SRSliceProxy* slc) -> float {
  float energy;
  float total=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13) {
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));
      total += energy;}     
    if (trk.truth.p.pdg == 2212) {
      energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
      total += (energy - mass_proton);}
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
      total += (energy - mass_pion);}}

  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++) {
      if (shw.truth.p.daughters[i]==111) {
        total += shw.bestplane_energy/1000;
        break;}}}

  return total;});

const Var kEresTrue2([](const caf::SRSliceProxy* slc) -> float {
  float true_energy=0;
  float energy;
  float total=0;
  float residual;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13) {
      true_energy += trk.truth.p.genE;
      energy = sqrt(pow(trk.rangeP.p_muon,2) + pow(mass_muon,2));
      total += energy;}     
    if (trk.truth.p.pdg == 2212) {
      true_energy += trk.truth.p.genT;
      energy = sqrt(pow(trk.rangeP.p_proton,2) + pow(mass_proton,2));
      total += (energy - mass_proton);}
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      true_energy += trk.truth.p.genT;
      energy = sqrt(pow(trk.rangeP.p_pion,2) + pow(mass_pion,2));
      total += (energy - mass_pion);}}

  for (auto const& shw : slc->reco.shw) {
    int daughters=shw.truth.p.daughters.size();
    for (int i=0; i<daughters; i++) {
      if (shw.truth.p.daughters[i]==111) {
        true_energy += shw.truth.p.genE;
        total += shw.bestplane_energy/1000;
        break;}}}

  residual = total - true_energy;
  return residual;});

// True total energy
const Var kTrueErec([](const caf::SRSliceProxy* slc) -> float {
  float energy=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 13) {
      energy += trk.truth.p.genE;}
    else {
      energy += trk.truth.p.genT;}}

  for (auto const& shw : slc->reco.shw) {
    energy += shw.truth.p.genE;}

  return energy;});

//const MultiVar kBestplane([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  //std::vector<double> bestplane;
  //for (auto const& trk : slc->reco.trk) {
    //bestplane.push_back(trk.bestplane);}
//
  //return bestplane;});

const Var kPTrackInd([](const caf::SRSliceProxy* slc) -> int {
  // The (dis)qualification of a slice is based upon the track level information.
  float Atslc, Chi2Proton, Chi2Muon, Longest(0);
  bool AtSlice, Contained, MaybeMuonExiting, MaybeMuonContained;
  int PTrackInd(-1);
  for (std::size_t i(0); i < slc->reco.trk.size(); ++i) {
    auto const& trk = slc->reco.trk.at(i);
    // First we calculate the distance of each track to the slice vertex.
    Atslc = sqrt( pow( slc->vertex.x - trk.start.x, 2.0 ) +
		  pow( slc->vertex.y - trk.start.y, 2.0 ) +
		  pow( slc->vertex.z - trk.start.z, 2.0 ) );
    // We require that the distance of the track from the slice is less than
    // 10 cm and that the parent of the track has been marked as the primary.
    AtSlice = ( Atslc < 10.0 && trk.parent_is_primary);
	
    if (trk.bestplane == 0) {
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Muon = trk.chi2pid0.chi2_muon;}
    else if (trk.bestplane == 1) {
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Muon = trk.chi2pid1.chi2_muon;}
    else {
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Muon = trk.chi2pid2.chi2_muon;}

    Contained = ( !isnan(trk.end.x) &&
		( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		  !isnan(trk.end.y) &&
		( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		  !isnan(trk.end.z) &&
		      ( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );
    MaybeMuonExiting = ( !Contained && trk.len > 100);
    MaybeMuonContained = ( Contained && Chi2Proton > 60 && Chi2Muon < 30 && trk.len > 50 );
    if ( AtSlice && ( MaybeMuonExiting || MaybeMuonContained ) && trk.len > Longest ) {
      Longest = trk.len;
      PTrackInd = i;}}
  return PTrackInd;});

// chi2_proton for true protons
const MultiVar kProtonChi2_proton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> chi2;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      chi2.push_back(trk.chi2pid2.chi2_proton/trk.chi2pid2.pid_ndof);}}

  return chi2;});

// chi2_pion for true protons
const MultiVar kProtonChi2_pion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> chi2;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 2212) {
      chi2.push_back(trk.chi2pid2.chi2_pion/trk.chi2pid2.pid_ndof);}}

  return chi2;});

// chi2_proton for true pions
const MultiVar kPionChi2_proton([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> chi2;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      chi2.push_back(trk.chi2pid2.chi2_proton/trk.chi2pid2.pid_ndof);}}

  return chi2;});

// chi2_pion for true pions
const MultiVar kPionChi2_pion([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> chi2;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      chi2.push_back(trk.chi2pid2.chi2_pion/trk.chi2pid2.pid_ndof);}}

  return chi2;});

// ---- CUTS -----

//const Cut kIsNu([](const caf::SRSliceProxy* slc) {
  //return (!isnan(slc->truth.E) && slc->truth.iscc);});

const Cut kNotClearCosmic([](const caf::SRSliceProxy* slc) {
  return !slc->is_clear_cosmic;});

const Cut kFMScore([](const caf::SRSliceProxy* slc) {
  return ( !isnan(slc->fmatch.score) && slc->fmatch.score < 10 );});

const Cut kMaxLen50([](const caf::SRSliceProxy* slc) {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}
  return ( maxlen > 50 );});

const Cut kMaxLen100([](const caf::SRSliceProxy* slc) {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}
  return ( maxlen > 100 );});

const Cut kChi2Cut([](const caf::SRSliceProxy* slc) {
  float maxlen=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len > maxlen)
      maxlen = trk.len;}

  float chi2_muon=0;
  float chi2_proton=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.len == maxlen)
      chi2_muon=trk.chi2pid1.chi2_muon;
      chi2_proton=trk.chi2pid1.chi2_proton;}

  return (chi2_muon < 30 && chi2_proton > 60);});

const Cut kRFiducial([](const caf::SRSliceProxy* slc) {
  return ( !isnan(slc->vertex.x) &&
	   ( ( slc->vertex.x < -71.1 - 25 && slc->vertex.x > -369.33 + 25 ) ||
           ( slc->vertex.x > 71.1 + 25 && slc->vertex.x < 369.33 - 25 ) ) &&
	   ( !isnan(slc->vertex.y) ) &&
	   ( slc->vertex.y > -181.7 + 25 && slc->vertex.y < 134.8 - 25 ) &&
	   ( !isnan(slc->vertex.z) ) &&
	   ( slc->vertex.z > -895.95 + 30 && slc->vertex.z < 895.95 - 50 ) );});

const Cut kPTrackContained([](const caf::SRSliceProxy* slc) {
  int Ind = kPTrackInd(slc);
  bool Contained(false);
  if ( Ind >= 0 ) {
    auto const& trk = slc->reco.trk.at(Ind);
    Contained = ( !isnan(trk.end.x) &&
		( trk.end.x < -71.1 - 25 && trk.end.x > -369.33 + 25 ) &&
		  !isnan(trk.end.y) &&
		( trk.end.y > -181.7 + 25 && trk.end.y < 134.8 - 25 ) &&
		  !isnan(trk.end.z) &&
		( trk.end.z > -895.95 + 30 && trk.end.z < 895.95 - 50 ) );}
  return Contained;});

const Cut kNoShowers([](const caf::SRSliceProxy* slc) {
  int count=0;
  for (auto const& shw : slc->reco.shw) {
    count += 1;}

  return ( count == 0 );});

const Cut kNoPions([](const caf::SRSliceProxy* slc) {
  int count=0;
  for (auto const& trk : slc->reco.trk) {
    if (trk.truth.p.pdg == 211 or trk.truth.p.pdg == -211) {
      count += 1;}}

  return ( count == 0 );});

const Cut kContainedCut = (kNotClearCosmic && kFMScore && kMaxLen50 && kChi2Cut && kRFiducial && kPTrackContained);
const Cut kUncontainedCut = (kNotClearCosmic && kFMScore && kMaxLen100 && kRFiducial);
 
// one file
//void Muon(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/caf-01608068-6435-4a36-93b5-29ead574d963.root")

// twenty-one files
void Muon(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/*.root")

// 354 files
//void Muon(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/*/*.root")

// all files
//void Muon(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/*/*/*.root")
//void Muon(const std::string inputName = "workshop_SBNWorkshop0421_prodoverlay_corsika_cosmics_cosmics_proton_genie_nu_spill_gsimple-config_caf_icarus")

{
  SpectrumLoader loader(inputName);

  // ---- SPECTRA -----
  // A spectrum is a histogram with associated POT information
  const Binning bins1 = Binning::Simple(30, 0, 3);
  const Binning bins2 = Binning::Simple(100, 0, 10);
  const Binning bins3 = Binning::Simple(100, -5, 5);
  const Binning bins4 = Binning::Simple(10, -5, 5);
  const Binning bins_chi2 = Binning::Simple(1000, 0, 1000);

  // Spectrum(Spectrumloader, HistAxis, Cut)
  // calorimeter
  Spectrum sErecMuon("sErecMuon", bins1, loader, kErecMuon, kNoSpillCut, kUncontainedCut);
  Spectrum sEresMuon("sEresMuon", bins3, loader, kEresMuon, kNoSpillCut, kUncontainedCut);
  Spectrum sErecProton("sErecProton", bins1, loader, kErecProton, kNoSpillCut, kUncontainedCut);
  Spectrum sEresProton("sEresProton", bins3, loader, kEresProton, kNoSpillCut, kUncontainedCut);
  Spectrum sErecPion("sErecPion", bins1, loader, kErecPion, kNoSpillCut, kUncontainedCut);
  Spectrum sEresPion("sEresPion", bins3, loader, kEresPion, kNoSpillCut, kUncontainedCut);
  Spectrum sErecPionNeut("sErecPionNeut", bins1, loader, kErecPionNeut, kNoSpillCut, kUncontainedCut);
  Spectrum sEresPionNeut("sEresPionNeut", bins3, loader, kEresPionNeut, kNoSpillCut, kUncontainedCut);

  Spectrum sErec("sErec", bins2, loader, kErec, kNoSpillCut, kUncontainedCut);
  Spectrum sEres("sEres", bins3, loader, kEres, kNoSpillCut, kUncontainedCut);
  Spectrum sErecb("sErecb", bins2, loader, kErec, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);
  Spectrum sEresb("sEresb", bins3, loader, kEres, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);

  Spectrum sErecContain("sErecContain", bins2, loader, kErec, kNoSpillCut, kContainedCut);
  Spectrum sEresContain("sEresContain", bins3, loader, kEres, kNoSpillCut, kContainedCut);
  Spectrum sErecContainb("sErecContainb", bins2, loader, kErec, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);
  Spectrum sEresContainb("sEresContainb", bins3, loader, kEres, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);

  // momentum
  Spectrum sErecMuon2("sErecMuon2", bins1, loader, kErecMuon2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresMuon2("sEresMuon2", bins3, loader, kEresMuon2, kNoSpillCut, kUncontainedCut);
  Spectrum sErecProton2("sErecProton2", bins1, loader, kErecProton2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresProton2("sEresProton2", bins3, loader, kEresProton2, kNoSpillCut, kUncontainedCut);
  Spectrum sErecPion2("sErecPion2", bins1, loader, kErecPion2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresPion2("sEresPion2", bins3, loader, kEresPion2, kNoSpillCut, kUncontainedCut);

  Spectrum sErec2("sErec2", bins2, loader, kErec2, kNoSpillCut, kUncontainedCut);
  Spectrum sEres2("sEres2", bins3, loader, kEres2, kNoSpillCut, kUncontainedCut);
  Spectrum sErec2b("sErec2b", bins2, loader, kErec2, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);
  Spectrum sEres2b("sEres2b", bins3, loader, kEres2, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);

  Spectrum sErecContain2("sErecContain2", bins2, loader, kErec2, kNoSpillCut, kContainedCut);
  Spectrum sEresContain2("sEresContain2", bins3, loader, kEres2, kNoSpillCut, kContainedCut);
  Spectrum sErecContain2b("sErecContain2b", bins2, loader, kErec2, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);
  Spectrum sEresContain2b("sEresContain2b", bins3, loader, kEres2, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);

  // calorimeter for true particles
  Spectrum sErecTrueMuon("sErecTrueMuon", bins1, loader, kErecTrueMuon, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTrueMuon("sEresTrueMuon", bins3, loader, kEresTrueMuon, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTrueProton("sErecTrueProton", bins1, loader, kErecTrueProton, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTrueProton("sEresTrueProton", bins3, loader, kEresTrueProton, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTruePion("sErecTruePion", bins1, loader, kErecTruePion, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTruePion("sEresTruePion", bins3, loader, kEresTruePion, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTruePionNeut("sErecTruePionNeut", bins1, loader, kErecTruePionNeut, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTruePionNeut("sEresTruePionNeut", bins3, loader, kEresTruePionNeut, kNoSpillCut, kUncontainedCut);

  Spectrum sErecTrue("sErecTrue", bins2, loader, kErecTrue, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTrue("sEresTrue", bins3, loader, kEresTrue, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTrueb("sErecTrueb", bins2, loader, kErecTrue, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);
  Spectrum sEresTrueb("sEresTrueb", bins3, loader, kEresTrue, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);

  Spectrum sErecTrueContain("sErecTrueContain", bins2, loader, kErecTrue, kNoSpillCut, kContainedCut);
  Spectrum sEresTrueContain("sEresTrueContain", bins3, loader, kEresTrue, kNoSpillCut, kContainedCut);
  Spectrum sErecTrueContainb("sErecTrueContainb", bins2, loader, kErecTrue, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);
  Spectrum sEresTrueContainb("sEresTrueContainb", bins3, loader, kEresTrue, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);

  // momentum for true particles
  Spectrum sErecTrueMuon2("sErecTrueMuon2", bins1, loader, kErecTrueMuon2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTrueMuon2("sEresTrueMuon2", bins3, loader, kEresTrueMuon2, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTrueProton2("sErecTrueProton2", bins1, loader, kErecTrueProton2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTrueProton2("sEresTrueProton2", bins3, loader, kEresTrueProton2, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTruePion2("sErecTruePion2", bins1, loader, kErecTruePion2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTruePion2("sEresTruePion2", bins3, loader, kEresTruePion2, kNoSpillCut, kUncontainedCut);

  Spectrum sErecTrue2("sErecTrue2", bins2, loader, kErecTrue2, kNoSpillCut, kUncontainedCut);
  Spectrum sEresTrue2("sEresTrue2", bins3, loader, kEresTrue2, kNoSpillCut, kUncontainedCut);
  Spectrum sErecTrue2b("sErecTrue2b", bins2, loader, kErecTrue2, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);
  Spectrum sEresTrue2b("sEresTrue2b", bins3, loader, kEresTrue2, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);

  Spectrum sErecTrueContain2("sErecTrueContain2", bins2, loader, kErecTrue2, kNoSpillCut, kContainedCut);
  Spectrum sEresTrueContain2("sEresTrueContain2", bins3, loader, kEresTrue2, kNoSpillCut, kContainedCut);
  Spectrum sErecTrueContain2b("sErecTrueContain2b", bins2, loader, kErecTrue2, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);
  Spectrum sEresTrueContain2b("sEresTrueContain2b", bins3, loader, kEresTrue2, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);

  // true energy
  Spectrum sTrueErecMuon("sTrueErecMuon", bins1, loader, kTrueErecMuon, kNoSpillCut, kUncontainedCut);
  Spectrum sTrueErecProton("sTrueErecProton", bins1, loader, kTrueErecProton, kNoSpillCut, kUncontainedCut);
  Spectrum sTrueErecPion("sTrueErecPion", bins1, loader, kTrueErecPion, kNoSpillCut, kUncontainedCut);
  Spectrum sTrueErecPionNeut("sTrueErecPionNeut", bins1, loader, kTrueErecPionNeut, kNoSpillCut, kUncontainedCut);
  Spectrum sTrueErec("sTrueErec", bins2, loader, kTrueErec, kNoSpillCut, kUncontainedCut);
  Spectrum sTrueErecb("sTrueErecb", bins2, loader, kTrueErec, kNoSpillCut, kUncontainedCut && kNoShowers && kNoPions);
  Spectrum sTrueErecContain("sTrueErecContain", bins2, loader, kTrueErec, kNoSpillCut, kContainedCut);
  Spectrum sTrueErecContainb("sTrueErecContainb", bins2, loader, kTrueErec, kNoSpillCut, kContainedCut && kNoShowers && kNoPions);

  //Spectrum sPproton("sPproton", bins1, loader, kPproton, kNoSpillCut, kUncontainedCut);
  //Spectrum sPpion("sPpion", bins1, loader, kPpion, kNoSpillCut, kUncontainedCut);
  //Spectrum sBestplane("sBestplane", bins4, loader, kBestplane, kNoSpillCut, kUncontainedCut);

  //const HistAxis axProtonChi2_proton("chi2_proton", bins_chi2, kProtonChi2_proton);
  //const HistAxis axProtonChi2_pion("chi2_pion", bins_chi2, kProtonChi2_pion);
  //const HistAxis axPionChi2_proton("chi2_proton", bins_chi2, kPionChi2_proton);
  //const HistAxis axPionChi2_pion("chi2_pion", bins_chi2, kPionChi2_pion);

  //Spectrum sProtonChi2("Proton Chi2", loader, bins_chi2, kProtonChi2_proton, bins_chi2, kProtonChi2_pion, kNoSpillCut, kNotClearCosmic && kFMScore && kMaxLen && kRFiducial);
  //Spectrum sPionChi2("Pion Chi2", loader, bins_chi2, kPionChi2_proton, bins_chi2, kPionChi2_pion, kNoSpillCut, kNotClearCosmic && kFMScore && kMaxLen && kRFiducial);

  // This is the call that actually fills in the spectrum
  loader.Go();

  // ---- DRAW -----
  // For plotting purposes we can convert spectra to a TH1
  TH1* hErecMuon = sErecMuon.ToTH1(6.6e20);
  TH1* hErecProton = sErecProton.ToTH1(6.6e20);
  TH1* hErecPion = sErecPion.ToTH1(6.6e20);
  TH1* hErecPionNeut = sErecPionNeut.ToTH1(6.6e20);

  TH1* hErecMuon2 = sErecMuon2.ToTH1(6.6e20);
  TH1* hErecProton2 = sErecProton2.ToTH1(6.6e20);
  TH1* hErecPion2 = sErecPion2.ToTH1(6.6e20);

  TH1* hTrueErecMuon = sTrueErecMuon.ToTH1(6.6e20);
  TH1* hTrueErecProton = sTrueErecProton.ToTH1(6.6e20);
  TH1* hTrueErecPion = sTrueErecPion.ToTH1(6.6e20);
  TH1* hTrueErecPionNeut = sTrueErecPionNeut.ToTH1(6.6e20);

  TH1* hErecTrueMuon = sErecTrueMuon.ToTH1(6.6e20);
  TH1* hErecTrueProton = sErecTrueProton.ToTH1(6.6e20);
  TH1* hErecTruePion = sErecTruePion.ToTH1(6.6e20);
  TH1* hErecTruePionNeut = sErecTruePionNeut.ToTH1(6.6e20);

  TH1* hErecTrueMuon2 = sErecTrueMuon2.ToTH1(6.6e20);
  TH1* hErecTrueProton2 = sErecTrueProton2.ToTH1(6.6e20);
  TH1* hErecTruePion2 = sErecTruePion2.ToTH1(6.6e20);

  TH1* hErec = sErec.ToTH1(6.6e20);
  TH1* hErec2 = sErec2.ToTH1(6.6e20);
  TH1* hErecTrue = sErecTrue.ToTH1(6.6e20);
  TH1* hErecTrue2 = sErecTrue2.ToTH1(6.6e20);
  TH1* hTrueErec = sTrueErec.ToTH1(6.6e20);

  TH1* hErecb = sErecb.ToTH1(6.6e20);
  TH1* hErec2b = sErec2b.ToTH1(6.6e20);
  TH1* hErecTrueb = sErecTrueb.ToTH1(6.6e20);
  TH1* hErecTrue2b = sErecTrue2b.ToTH1(6.6e20);
  TH1* hTrueErecb = sTrueErecb.ToTH1(6.6e20);

  TH1* hEres = sEres.ToTH1(6.6e20);
  TH1* hEres2 = sEres2.ToTH1(6.6e20);
  TH1* hEresTrue = sEresTrue.ToTH1(6.6e20);
  TH1* hEresTrue2 = sEresTrue2.ToTH1(6.6e20);

  TH1* hEresb = sEresb.ToTH1(6.6e20);
  TH1* hEres2b = sEres2b.ToTH1(6.6e20);
  TH1* hEresTrueb = sEresTrueb.ToTH1(6.6e20);
  TH1* hEresTrue2b = sEresTrue2b.ToTH1(6.6e20);

  TH1* hErecContain = sErecContain.ToTH1(6.6e20);
  TH1* hErecContain2 = sErecContain2.ToTH1(6.6e20);
  TH1* hErecTrueContain = sErecTrueContain.ToTH1(6.6e20);
  TH1* hErecTrueContain2 = sErecTrueContain2.ToTH1(6.6e20);
  TH1* hTrueErecContain = sTrueErecContain.ToTH1(6.6e20);

  TH1* hErecContainb = sErecContainb.ToTH1(6.6e20);
  TH1* hErecContain2b = sErecContain2b.ToTH1(6.6e20);
  TH1* hErecTrueContainb = sErecTrueContainb.ToTH1(6.6e20);
  TH1* hErecTrueContain2b = sErecTrueContain2b.ToTH1(6.6e20);
  TH1* hTrueErecContainb = sTrueErecContainb.ToTH1(6.6e20);

  TH1* hEresContain = sEresContain.ToTH1(6.6e20);
  TH1* hEresContain2 = sEresContain2.ToTH1(6.6e20);
  TH1* hEresTrueContain = sEresTrueContain.ToTH1(6.6e20);
  TH1* hEresTrueContain2 = sEresTrueContain2.ToTH1(6.6e20);

  TH1* hEresContainb = sEresContainb.ToTH1(6.6e20);
  TH1* hEresContain2b = sEresContain2b.ToTH1(6.6e20);
  TH1* hEresTrueContainb = sEresTrueContainb.ToTH1(6.6e20);
  TH1* hEresTrueContain2b = sEresTrueContain2b.ToTH1(6.6e20);

  // Set the curve colors
  hErecMuon->SetLineColor(kRed);
  hErecProton->SetLineColor(kRed);
  hErecPion->SetLineColor(kRed);
  hErecPionNeut->SetLineColor(kRed);

  hErecMuon2->SetLineColor(kPink+9);
  hErecProton2->SetLineColor(kPink+9);
  hErecPion2->SetLineColor(kPink+9);

  hTrueErecMuon->SetLineColor(kBlue);
  hTrueErecProton->SetLineColor(kBlue);
  hTrueErecPion->SetLineColor(kBlue);
  hTrueErecPionNeut->SetLineColor(kBlue);

  hErecTrueMuon->SetLineColor(kViolet-6);
  hErecTrueProton->SetLineColor(kViolet-6);
  hErecTruePion->SetLineColor(kViolet-6);
  hErecTruePionNeut->SetLineColor(kViolet-6);

  hErecTrueMuon2->SetLineColor(kViolet-9);
  hErecTrueProton2->SetLineColor(kViolet-9);
  hErecTruePion2->SetLineColor(kViolet-9);

  hErec->SetLineColor(kRed);
  hErec2->SetLineColor(kPink+9);
  hErecTrue->SetLineColor(kViolet-6);
  hErecTrue2->SetLineColor(kViolet-9);
  hTrueErec->SetLineColor(kBlue);

  hErecb->SetLineColor(kRed);
  hErec2b->SetLineColor(kPink+9);
  hErecTrueb->SetLineColor(kViolet-6);
  hErecTrue2b->SetLineColor(kViolet-9);
  hTrueErecb->SetLineColor(kBlue);

  hEres->SetLineColor(kRed);
  hEres2->SetLineColor(kPink+9);
  hEresTrue->SetLineColor(kViolet-6);
  hEresTrue2->SetLineColor(kViolet-9);

  hEresb->SetLineColor(kRed);
  hEres2b->SetLineColor(kPink+9);
  hEresTrueb->SetLineColor(kViolet-6);
  hEresTrue2b->SetLineColor(kViolet-9);

  hErecContain->SetLineColor(kRed);
  hErecContain2->SetLineColor(kPink+9);
  hErecTrueContain->SetLineColor(kViolet-6);
  hErecTrueContain2->SetLineColor(kViolet-9);
  hTrueErecContain->SetLineColor(kBlue);

  hErecContainb->SetLineColor(kRed);
  hErecContain2b->SetLineColor(kPink+9);
  hErecTrueContainb->SetLineColor(kViolet-6);
  hErecTrueContain2b->SetLineColor(kViolet-9);
  hTrueErecContainb->SetLineColor(kBlue);

  hEresContain->SetLineColor(kRed);
  hEresContain2->SetLineColor(kPink+9);
  hEresTrueContain->SetLineColor(kViolet-6);
  hEresTrueContain2->SetLineColor(kViolet-9);

  hEresContainb->SetLineColor(kRed);
  hEresContain2b->SetLineColor(kPink+9);
  hEresTrueContainb->SetLineColor(kViolet-6);
  hEresTrueContain2b->SetLineColor(kViolet-9);

  // Plotted trk.bestplane and found that entries are 0, 1, or 2
  //TH1* hBestplane = sBestplane.ToTH1(6.6e20);
  //TCanvas* c0 = new TCanvas("c0", "c0");
  //hBestplane->Draw("hist");
  //c0->Print("Bestplane.png");

  // Plot muon energy
  TCanvas* c1 = new TCanvas("c1", "c1");
  hErecMuon->Draw("hist");
  hErecMuon2->Draw("same hist");
  hErecTrueMuon->Draw("same hist");
  hErecTrueMuon2->Draw("same hist");
  hTrueErecMuon->Draw("same hist");

  TLegend* leg = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg->AddEntry(hErecMuon, "KE for longest track particles (calorimeter)", "l");
  leg->AddEntry(hErecMuon2, "KE for longest track particles (momentum)", "l");
  leg->AddEntry(hErecTrueMuon, "KE for true muons (calorimeter)", "l");
  leg->AddEntry(hErecTrueMuon2, "KE for true muons (momentum)", "l"); 
  leg->AddEntry(hTrueErecMuon, "True KE for muons","l");
  leg->Draw();
  c1->Print("ErecMuon.png");

  // Plot proton KE 
  TCanvas* c2 = new TCanvas("c2", "c2");
  hErecProton->Draw("hist");
  hErecProton2->Draw("same hist");
  hErecTrueProton->Draw("same hist");
  hErecTrueProton2->Draw("same hist");
  hTrueErecProton->Draw("same hist");
  
  TLegend* leg2 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg2->AddEntry(hErecProton, "KE for smaller chi^2 particles (calorimeter)", "l");
  leg2->AddEntry(hErecProton2, "KE for smaller chi^2 particles (momentum)", "l");
  leg2->AddEntry(hErecTrueProton, "KE for true protons (calorimeter)", "l");
  leg2->AddEntry(hErecTrueProton2, "KE for true protons (momentum)", "l");
  leg2->AddEntry(hTrueErecProton, "True KE for protons", "l");
  leg2->Draw();
  c2->Print("ErecProton.png");
 
  // Plot pion KE 
  TCanvas* c3 = new TCanvas("c3", "c3");
  hErecPion->Draw("hist");
  hErecPion2->Draw("same hist");
  hErecTruePion->Draw("same hist");
  hErecTruePion2->Draw("same hist");
  hTrueErecPion->Draw("same hist");
  
  TLegend* leg3 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg3->AddEntry(hErecPion, "KE for smaller chi^2 particles (calorimeter)", "l");
  leg3->AddEntry(hErecPion2, "KE for smaller chi^2 particles (momentum)", "l");
  leg3->AddEntry(hErecTruePion, "KE for true pions (calorimeter)", "l");
  leg3->AddEntry(hErecTruePion2, "KE for true pions (momentum)", "l");
  leg3->AddEntry(hTrueErecPion, "True KE for pions", "l"); 
  leg3->Draw();
  c3->Print("ErecPion.png");

  // Plot neutral pion energy
  TCanvas* c9 = new TCanvas("c9", "c9");
  hErecPionNeut->Draw("hist");
  hErecTruePionNeut->Draw("same hist");
  hTrueErecPionNeut->Draw("same hist");
  
  TLegend* leg9 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg9->AddEntry(hErecPionNeut, "Energy of shower particles", "l");
  leg9->AddEntry(hErecTruePionNeut, "Energy of true neutral pions", "l");
  leg9->AddEntry(hTrueErecPionNeut, "True energy of neutral pions", "l"); 
  leg9->Draw();
  c9->Print("ErecPionNeut.png");

  // Plot Proton & Pion Momentum
  //TCanvas* c4 = new TCanvas("c4", "c4");
  //TH1* hPproton = sPproton.ToTH1(6.6e20);
  //hPproton->SetLineColor(kGreen);
  //hPproton->Draw("hist");
  //TH1* hPpion = sPpion.ToTH1(6.6e20);
  //hPpion->SetLineColor(kAzure+10);
  //hPpion->Draw("same hist");

  //TLegend* leg4 = new TLegend(0.6, 0.7, 0.85, 0.8);
  //leg4->AddEntry(hPproton, "Proton Momentum (GeV/c^2)", "l");
  //leg4->AddEntry(hPpion, "Pion Momentum (GeV/c^2)", "l");
  //leg4->Draw();
  //c4->Print("Momentum.png");

  // Plot residual muon energy
  TH1* hEresMuon = sEresMuon.ToTH1(6.6e20);
  TH1* hEresMuon2 = sEresMuon2.ToTH1(6.6e20);
  TH1* hEresTrueMuon = sEresTrueMuon.ToTH1(6.6e20);
  TH1* hEresTrueMuon2 = sEresTrueMuon2.ToTH1(6.6e20);

  hEresMuon->SetLineColor(kRed);
  hEresMuon2->SetLineColor(kPink+9);
  hEresTrueMuon->SetLineColor(kViolet-6);
  hEresTrueMuon2->SetLineColor(kViolet-9);

  TCanvas* c5 = new TCanvas("c5", "c5");
  hEresMuon->Draw("hist");
  hEresMuon2->Draw("same hist");
  hEresTrueMuon->Draw("same hist");
  hEresTrueMuon2->Draw("same hist");

  TLegend* leg5 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg5->AddEntry(hEresMuon, "Residual KE for longest track particles (calorimeter)", "l");
  leg5->AddEntry(hEresMuon2, "Residual KE for longest track particles (momentum)", "l");
  leg5->AddEntry(hEresTrueMuon, "Residual KE for true muons (calorimeter)", "l");
  leg5->AddEntry(hEresTrueMuon2, "Residual KE for true muons (momentum)", "l"); 
  leg5->Draw();
  c5->Print("EresMuon.png");

  // Plot residual proton KE
  TH1* hEresProton = sEresProton.ToTH1(6.6e20);
  TH1* hEresProton2 = sEresProton2.ToTH1(6.6e20);
  TH1* hEresTrueProton = sEresTrueProton.ToTH1(6.6e20);
  TH1* hEresTrueProton2 = sEresTrueProton2.ToTH1(6.6e20);

  hEresProton->SetLineColor(kRed);
  hEresProton2->SetLineColor(kPink+9);
  hEresTrueProton->SetLineColor(kViolet-6);
  hEresTrueProton2->SetLineColor(kViolet-9);

  TCanvas* c6 = new TCanvas("c6", "c6");
  hEresProton->Draw("hist");
  hEresProton2->Draw("same hist");
  hEresTrueProton->Draw("same hist");
  hEresTrueProton2->Draw("same hist");

  TLegend* leg6 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg6->AddEntry(hEresProton, "Residual KE for smaller chi^2 particles (calorimeter)", "l");
  leg6->AddEntry(hEresProton2, "Residual KE for smaller chi^2 particles (momentum)", "l");
  leg6->AddEntry(hEresTrueProton, "Residual KE for true protons (calorimeter)", "l");
  leg6->AddEntry(hEresTrueProton2, "Residual KE for true protons (momentum)", "l");
  leg6->Draw();
  c6->Print("EresProton.png");

  // Plot residual pion KE
  TH1* hEresPion = sEresPion.ToTH1(6.6e20);
  TH1* hEresPion2 = sEresPion2.ToTH1(6.6e20);
  TH1* hEresTruePion = sEresTruePion.ToTH1(6.6e20);
  TH1* hEresTruePion2 = sEresTruePion2.ToTH1(6.6e20);

  hEresPion->SetLineColor(kRed);
  hEresPion2->SetLineColor(kPink+9);
  hEresTruePion->SetLineColor(kViolet-6);
  hEresTruePion2->SetLineColor(kViolet-9);

  TCanvas* c7 = new TCanvas("c7", "c7");
  hEresPion->Draw("hist");
  hEresPion2->Draw("same hist");
  hEresTruePion->Draw("same hist");
  hEresTruePion2->Draw("same hist");

  TLegend* leg7 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg7->AddEntry(hEresPion, "Residual KE for smaller chi^2 particles (calorimeter)", "l");
  leg7->AddEntry(hEresPion2, "Residual KE for smaller chi^2 particles (momentum)", "l");
  leg7->AddEntry(hEresTruePion, "Residual KE for true pions (calorimeter)", "l");
  leg7->AddEntry(hEresTruePion2, "Residual KE for true pions (momentum)", "l");
  leg7->Draw();
  c7->Print("EresPion.png");

  // Plot residual neutral pion energy
  TH1* hEresPionNeut = sEresPionNeut.ToTH1(6.6e20);
  TH1* hEresTruePionNeut = sEresTruePionNeut.ToTH1(6.6e20);

  hEresPionNeut->SetLineColor(kRed);
  hEresTruePionNeut->SetLineColor(kViolet-6);

  TCanvas* c10 = new TCanvas("c10", "c10");
  hEresPionNeut->Draw("hist");
  hEresTruePionNeut->Draw("same hist");

  TLegend* leg10 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg10->AddEntry(hErecPionNeut, "Residual energy of shower particles", "l");
  leg10->AddEntry(hErecTruePionNeut, "Residual energy of true neutral pions", "l");
  leg10->Draw();
  c10->Print("EresPionNeut.png");

// Plot total reconstructed energy
  TCanvas* c11 = new TCanvas("c11", "c11");
  hErec->Draw("hist");
  hErec2->Draw("same hist");
  hErecTrue->Draw("same hist");
  hErecTrue2->Draw("same hist");
  hTrueErec->Draw("same hist");

  TLegend* leg11 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg11->AddEntry(hErec, "Reconstructed energy (calorimeter)", "l");
  leg11->AddEntry(hErec2, "Reconstructed energy (momentum)", "l");
  leg11->AddEntry(hErecTrue, "'True' reconstructed energy (calorimeter)", "l");
  leg11->AddEntry(hErecTrue2, "'True' reconstructed energy (momentum)", "l");
  leg11->AddEntry(hTrueErec, "True energy", "l");
  leg11->Draw();
  c11->Print("Erec.png");

// Plot total reconstructed energy in events with no pions
  TCanvas* c13 = new TCanvas("c13", "c13");
  hErecb->Draw("hist");
  hErec2b->Draw("same hist");
  hErecTrueb->Draw("same hist");
  hErecTrue2b->Draw("same hist");
  hTrueErecb->Draw("same hist");

  TLegend* leg13 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg13->AddEntry(hErecb, "Reconstructed energy (calorimeter)", "l");
  leg13->AddEntry(hErec2b, "Reconstructed energy (momentum)", "l");
  leg13->AddEntry(hErecTrueb, "'True' reconstructed energy (calorimeter)", "l");
  leg13->AddEntry(hErecTrue2b, "'True' reconstructed energy (momentum)", "l");
  leg13->AddEntry(hTrueErecb, "True energy", "l");
  leg13->Draw();
  c13->Print("Erecb.png");

// Plot total reconstructed energy (contained)
  TCanvas* c19 = new TCanvas("c19", "c19");
  hErecContain->Draw("hist");
  hErecContain2->Draw("same hist");
  hErecTrueContain->Draw("same hist");
  hErecTrueContain2->Draw("same hist");
  hTrueErecContain->Draw("same hist");

  TLegend* leg19 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg19->AddEntry(hErecContain, "Reconstructed energy (calorimeter)", "l");
  leg19->AddEntry(hErecContain2, "Reconstructed energy (momentum)", "l");
  leg19->AddEntry(hErecTrueContain, "'True' reconstructed energy (calorimeter)", "l");
  leg19->AddEntry(hErecTrueContain2, "'True' reconstructed energy (momentum)", "l");
  leg19->AddEntry(hTrueErecContain, "True energy", "l");
  leg19->Draw();
  c19->Print("ErecContain.png");

// Plot total reconstructed energy in events with no pions (contained)
  TCanvas* c20 = new TCanvas("c20", "c20");
  hErecContainb->Draw("hist");
  hErecContain2b->Draw("same hist");
  hErecTrueContainb->Draw("same hist");
  hErecTrueContain2b->Draw("same hist");
  hTrueErecContainb->Draw("same hist");

  TLegend* leg20 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg20->AddEntry(hErecContainb, "Reconstructed energy (calorimeter)", "l");
  leg20->AddEntry(hErecContain2b, "Reconstructed energy (momentum)", "l");
  leg20->AddEntry(hErecTrueContainb, "'True' reconstructed energy (calorimeter)", "l");
  leg20->AddEntry(hErecTrueContain2b, "'True' reconstructed energy (momentum)", "l");
  leg20->AddEntry(hTrueErecContainb, "True energy", "l");
  leg20->Draw();
  c20->Print("ErecContainb.png");

// Plot total residual energy
  TCanvas* c12 = new TCanvas("c12", "c12");
  hEres->Draw("hist");
  hEres2->Draw("same hist");
  hEresTrue->Draw("same hist");
  hEresTrue2->Draw("same hist");

  TLegend* leg12 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg12->AddEntry(hEres, "Residual energy (calorimeter)", "l");
  leg12->AddEntry(hEres2, "Residual energy (momentum)", "l");
  leg12->AddEntry(hEresTrue, "'True' residual energy (calorimeter)", "l");
  leg12->AddEntry(hEresTrue2, "'True' residual energy (momentum)", "l");
  leg12->Draw();
  c12->Print("Eres.png");

// Plot total residual energy in events with no pions
  TCanvas* c14 = new TCanvas("c14", "c14");
  hEresb->Draw("hist");
  hEres2b->Draw("same hist");
  hEresTrueb->Draw("same hist");
  hEresTrue2b->Draw("same hist");

  TLegend* leg14 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg14->AddEntry(hEresb, "Residual energy (calorimeter)", "l");
  leg14->AddEntry(hEres2b, "Residual energy (momentum)", "l");
  leg14->AddEntry(hEresTrueb, "'True' residual energy (calorimeter)", "l");
  leg14->AddEntry(hEresTrue2b, "'True' residual energy (momentum)", "l");
  leg14->Draw();
  c14->Print("Eresb.png");

// Plot total residual energy (contained)
  TCanvas* c17 = new TCanvas("c17", "c17");
  hEresContain->Draw("hist");
  hEresContain2->Draw("same hist");
  hEresTrueContain->Draw("same hist");
  hEresTrueContain2->Draw("same hist");

  TLegend* leg17 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg17->AddEntry(hEresContain, "Residual energy (calorimeter)", "l");
  leg17->AddEntry(hEresContain2, "Residual energy (momentum)", "l");
  leg17->AddEntry(hEresTrueContain, "'True' residual energy (calorimeter)", "l");
  leg17->AddEntry(hEresTrueContain2, "'True' residual energy (momentum)", "l");
  leg17->Draw();
  c17->Print("EresContain.png");

// Plot total residual energy in events with no pions (contained)
  TCanvas* c18 = new TCanvas("c18", "c18");
  hEresContainb->Draw("hist");
  hEresContain2b->Draw("same hist");
  hEresTrueContainb->Draw("same hist");
  hEresTrueContain2b->Draw("same hist");

  TLegend* leg18 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg18->AddEntry(hEresContainb, "Residual energy (calorimeter)", "l");
  leg18->AddEntry(hEresContain2b, "Residual energy (momentum)", "l");
  leg18->AddEntry(hEresTrueContainb, "'True' residual energy (calorimeter)", "l");
  leg18->AddEntry(hEresTrueContain2b, "'True' residual energy (momentum)", "l");
  leg18->Draw();
  c18->Print("EresContainb.png");

// Plot chi2_proton vs chi2_pion for true protons
  //TCanvas* c15 = new TCanvas("c15", "c15");
  //TH2* hProtonChi2 = sProtonChi2.ToTH2(6.6e20);
  //hProtonChi2->Draw("colz");
  //c15->Print("ProtonChi1.png");

// Plot chi2_proton vs chi2_pion for true pions
  //TCanvas* c16 = new TCanvas("c16", "c16");
  //TH2* hPionChi2 = sPionChi2.ToTH2(6.6e20);
  //hPionChi2->Draw("colz");
  //c16->Print("PionChi1.png");

}
