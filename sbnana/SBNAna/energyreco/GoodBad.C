#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

using namespace ana;

#include "myConstants.h"
#include "myEstimator.h"
#include "myEventSelection.h"
#include "myMuonSelection.h"
#include "myProtonSelection.h"
#include "myTruth.h"
#include "Chi2Plotter.h"
#include "varsTransverseP.h"
#include "varsErec.h"

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include <numeric>
#include <fstream>

// ---- VARS -----
const Var varEnergyReco([](const caf::SRSliceProxy* slc) -> double {
  // Muon energy
  double P_Muon = varMuonTrackCombinedP(slc);
  double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );
  // Proton KE
  std::vector<double> KE_Proton = varProtonRangeKE(slc);
  double Total_KE_Proton = std::accumulate(KE_Proton.begin(), KE_Proton.end(), 0.0);
  // Pion KE
  std::vector<double> E_Pion = varPionRangeE(slc);
  double Total_E_Pion = std::accumulate(E_Pion.begin(), E_Pion.end(), 0.0);
  // Shower energy
  double E_shws = varShwE(slc);
  // Neutrino energy
  double E_Nu = E_Muon + Total_KE_Proton + Total_E_Pion + E_shws;
  return E_Nu;
});

const Var varEnergyRes([](const caf::SRSliceProxy* slc) -> double {
  double E_Nu_Reco = varEnergyReco(slc);
  double E_Nu_True = varNeutrinoTruthE(slc);
  double E_Nu_Res = E_Nu_Reco - E_Nu_True;
  return E_Nu_Res;
});

// look at poorly reconstructed events with Eres below EresMax
double EresMax = -0.5; //GeV

int count_good = 0;
int countMCSP_good = 0;
int countRangep_good = 0;
int count_bad = 0;
int countMCSP_bad = 0;
int countRangep_bad = 0;

const Var varMuonPCounter([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex >= 0){
    int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
    double Eres = varEnergyRes(slc);
    if (Eres >= EresMax) {
      count_good += 1;
      if (mtc == 0) {
        countMCSP_good += 1;
      }
      if (mtc == 1) {
        countRangep_good += 1;
      }}
    if (Eres < EresMax) {
      count_bad += 1;
      if (mtc == 0) {
        countMCSP_bad += 1;
      }
      if (mtc == 1) {
        countRangep_bad += 1;
      }}}
  return 1.89; //the number of liters in Ocean Spray CranMango
});

const Var varEnergyTrue([](const caf::SRSliceProxy* slc) -> double {
  double E_Nu_True = varNeutrinoTruthE(slc);
  return E_Nu_True;
});

const Var varEmuon([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon = varMuonTrackCombinedP(slc);
  double  E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );
  return E_Muon;
});

const Var varTrueEmuon([](const caf::SRSliceProxy* slc) -> double {
  double E_Muon = varMuonTrackMatchedTruthE(slc); 
  return E_Muon;
});

const Var varEresMuon([](const caf::SRSliceProxy* slc) -> double {
  double Erec = varEmuon(slc);
  double Etrue = varTrueEmuon(slc);
  double Eres = Erec - Etrue;
  return Eres;
});

const Var varCosTheta([](const caf::SRSliceProxy* slc) -> double {
  double CosTheta=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    CosTheta = slc->reco.trk[muon_index].costh; 
  }
  return CosTheta;
});

const Var varTrueCosTheta([](const caf::SRSliceProxy* slc) -> double {
  double CosTheta=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    double px = slc->reco.trk[muon_index].truth.p.genp.x;
    double py = slc->reco.trk[muon_index].truth.p.genp.y;
    double pz = slc->reco.trk[muon_index].truth.p.genp.z;
    double p = sqrt( px*px + py*py + pz*pz );
    CosTheta = pz/p;
  }
  return CosTheta;
});

const Var varDeltaCosTheta([](const caf::SRSliceProxy* slc) -> double {
  double cos_theta = varCosTheta(slc);
  double true_cos_theta = varTrueCosTheta(slc);
  double delta_cos_theta = cos_theta - true_cos_theta;
  return delta_cos_theta;
});

const Var varTrueStartx([](const caf::SRSliceProxy* slc) -> double {
  double start_x=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_x = slc->reco.trk[muon_index].truth.p.start.x;
  }
  return start_x;
});

const Var varTrueStarty([](const caf::SRSliceProxy* slc) -> double {
  double start_y=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_y = slc->reco.trk[muon_index].truth.p.start.y;
  }
  return start_y;
});

const Var varTrueStartz([](const caf::SRSliceProxy* slc) -> double {
  double start_z=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_z = slc->reco.trk[muon_index].truth.p.start.z;
  }
  return start_z;
});

const Var varStartx([](const caf::SRSliceProxy* slc) -> double {
  double start_x=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_x = slc->reco.trk[muon_index].start.x;
  }
  return start_x;
});

const Var varDeltaStartx([](const caf::SRSliceProxy* slc) -> double {
  double start_x = varStartx(slc);
  double true_start_x = varTrueStartx(slc);
  double delta_start_x = start_x - true_start_x;
  return delta_start_x;
});

const Var varStarty([](const caf::SRSliceProxy* slc) -> double {
  double start_y=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_y = slc->reco.trk[muon_index].start.y;
  }
  return start_y;
});

const Var varDeltaStarty([](const caf::SRSliceProxy* slc) -> double {
  double start_y = varStarty(slc);
  double true_start_y = varTrueStarty(slc);
  double delta_start_y = start_y - true_start_y;
  return delta_start_y;
});

const Var varStartz([](const caf::SRSliceProxy* slc) -> double {
  double start_z=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_z = slc->reco.trk[muon_index].start.z;
  }
  return start_z;
});

const Var varDeltaStartz([](const caf::SRSliceProxy* slc) -> double {
  double start_z = varStartz(slc);
  double true_start_z = varTrueStartz(slc);
  double delta_start_z = start_z - true_start_z;
  return delta_start_z;
});

const Var varTrueEndx([](const caf::SRSliceProxy* slc) -> double {
  double end_x=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_x = slc->reco.trk[muon_index].truth.p.end.x;
  }
  return end_x;
});

const Var varTrueEndy([](const caf::SRSliceProxy* slc) -> double {
  double end_y=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_y = slc->reco.trk[muon_index].truth.p.end.y;
  }
  return end_y;
});

const Var varTrueEndz([](const caf::SRSliceProxy* slc) -> double {
  double end_z=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_z = slc->reco.trk[muon_index].truth.p.end.z;
  }
  return end_z;
});

const Var varEndx([](const caf::SRSliceProxy* slc) -> double {
  double end_x=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_x = slc->reco.trk[muon_index].end.x;
  }
  return end_x;
});

const Var varDeltaEndx([](const caf::SRSliceProxy* slc) -> double {
  double end_x = varEndx(slc);
  double true_end_x = varTrueEndx(slc);
  double delta_end_x = end_x - true_end_x;
  return delta_end_x;
});

const Var varEndy([](const caf::SRSliceProxy* slc) -> double {
  double end_y=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_y = slc->reco.trk[muon_index].end.y;
  }
  return end_y;
});

const Var varDeltaEndy([](const caf::SRSliceProxy* slc) -> double {
  double end_y = varEndy(slc);
  double true_end_y = varTrueEndy(slc);
  double delta_end_y = end_y - true_end_y;
  return delta_end_y;
});

const Var varEndz([](const caf::SRSliceProxy* slc) -> double {
  double end_z=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_z = slc->reco.trk[muon_index].end.z;
  }
  return end_z;
});

const Var varDeltaEndz([](const caf::SRSliceProxy* slc) -> double {
  double end_z = varEndz(slc);
  double true_end_z = varTrueEndz(slc);
  double delta_end_z = end_z - true_end_z;
  return delta_end_z;
});

const Var varDeltaEndStartx([](const caf::SRSliceProxy* slc) -> double {
  double end_x = varEndx(slc);
  double true_start_x = varTrueStartx(slc);
  double delta_end_start_x = end_x - true_start_x;
  return delta_end_start_x;
});

const Var varDeltaEndStarty([](const caf::SRSliceProxy* slc) -> double {
  double end_y = varEndy(slc);
  double true_start_y = varTrueStarty(slc);
  double delta_end_start_y = end_y - true_start_y;
  return delta_end_start_y;
});

const Var varDeltaEndStartz([](const caf::SRSliceProxy* slc) -> double {
  double end_z = varEndz(slc);
  double true_start_z = varTrueStartz(slc);
  double delta_end_start_z = end_z - true_start_z;
  return delta_end_start_z;
});

const Var varTrueMuonTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
  // -1: no trk found 
  // 0: exiting
  // 1: contained
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    double XMargin = 25.;
    double YMargin = 25.;
    double ZMarginUp = 30.;
    double ZMarginDown = 50.;
    bool isContained;
    // cryo 0
    if( trk_Muon.truth.p.end.x < 0 ){
      isContained = ( !isnan(trk_Muon.truth.p.end.x) &&
                    ( trk_Muon.truth.p.end.x < -71.1 - XMargin && trk_Muon.truth.p.end.x > -369.33 + XMargin ) &&
                    !isnan(trk_Muon.truth.p.end.y) &&
                    ( trk_Muon.truth.p.end.y > -181.7 + YMargin && trk_Muon.truth.p.end.y < 134.8 - YMargin ) &&
                    !isnan(trk_Muon.truth.p.end.z) &&
                    ( trk_Muon.truth.p.end.z > -895.95 + ZMarginUp && trk_Muon.truth.p.end.z < 895.95 - ZMarginDown ) );
    }
    // cryo 1
    else{
      isContained = ( !isnan(trk_Muon.truth.p.end.x) &&
                    ( trk_Muon.truth.p.end.x > 71.1 + XMargin && trk_Muon.truth.p.end.x < 369.33 - XMargin ) &&
                    !isnan(trk_Muon.truth.p.end.y) &&
                    ( trk_Muon.truth.p.end.y > -181.7 + YMargin && trk_Muon.truth.p.end.y < 134.8 - YMargin ) &&
                    !isnan(trk_Muon.truth.p.end.z) &&
                    ( trk_Muon.truth.p.end.z > -895.95 + ZMarginUp && trk_Muon.truth.p.end.z < 895.95 - ZMarginDown ) );
    }
    if(isContained) return 1;
    else return 0;
  }
  else{
    return -1;
  }
});

const Var varMCSP_good([](const caf::SRSliceProxy* slc) -> double {
  double P=-999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres >= EresMax && mtc == 0) {
    P = varMuonTrackMCSP(slc);
  }
  return P;
});

const Var varMCSP_bad([](const caf::SRSliceProxy* slc) -> double {
  double P=-999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres < EresMax && mtc == 0) {
    P = varMuonTrackMCSP(slc);
  }
  return P;
});

const Var varRangeP_good([](const caf::SRSliceProxy* slc) -> double {
  double P=-999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres >= EresMax && mtc == 1) {
    P = varMuonTrackRangeP(slc);
  }
  return P;
});

const Var varRangeP_bad([](const caf::SRSliceProxy* slc) -> double {
  double P=-999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres < EresMax && mtc == 1) {
    P = varMuonTrackRangeP(slc);
  }
  return P;
});

const Var varTrueP_good0([](const caf::SRSliceProxy* slc) -> double {
  double P=999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres >= EresMax && mtc == 0) {
    P = varMuonTrackMatchedTruthP(slc);
  }
  return P;
});

const Var varTrueP_good1([](const caf::SRSliceProxy* slc) -> double {
  double P=999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres >= EresMax && mtc == 1) {
    P = varMuonTrackMatchedTruthP(slc);
  }
  return P;
});

const Var varTrueP_bad0([](const caf::SRSliceProxy* slc) -> double {
  double P=999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres < EresMax && mtc == 0) {
    P = varMuonTrackMatchedTruthP(slc);
  }
  return P;
});

const Var varTrueP_bad1([](const caf::SRSliceProxy* slc) -> double {
  double P=999;
  double Eres = varEnergyRes(slc);
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  if (Eres < EresMax && mtc == 1) {
    P = varMuonTrackMatchedTruthP(slc);
  }
  return P;
});

const Var varPres_good0([](const caf::SRSliceProxy* slc) -> double {
  double Preco = varMCSP_good(slc);
  double Ptrue = varTrueP_good0(slc);
  double Pres = Preco - Ptrue;
  return Pres;
});

const Var varPres_good1([](const caf::SRSliceProxy* slc) -> double {
  double Preco = varRangeP_good(slc);
  double Ptrue = varTrueP_good1(slc);
  double Pres = Preco - Ptrue;
  return Pres;
});

const Var varPres_bad0([](const caf::SRSliceProxy* slc) -> double {
  double Preco = varMCSP_bad(slc);
  double Ptrue = varTrueP_bad0(slc);
  double Pres = Preco - Ptrue;
  return Pres;
});

const Var varPres_bad1([](const caf::SRSliceProxy* slc) -> double {
  double Preco = varRangeP_bad(slc);
  double Ptrue = varTrueP_bad1(slc);
  double Pres = Preco - Ptrue;
  return Pres;
});

// CUTS
const Cut mycutIsNuMuCC([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.iscc && slc->truth.pdg == 14 );
});

const Cut mycutIsAntiNuMuCC([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.iscc && slc->truth.pdg == -14 );
});

Cut cutNuMu = mycutIsNuMuCC && kRFiducial && cutHasMuonTrack && cutIsMuonTrackLong && cutExitingHadrons;
Cut cutAntiNuMu = mycutIsAntiNuMuCC && kRFiducial && cutHasMuonTrack && cutIsMuonTrackLong && cutExitingHadrons;

const Cut cutGoodReco([](const caf::SRSliceProxy* slc) {
  double Eres = varEnergyRes(slc);
  return (Eres>=EresMax);
});

const Cut cutBadReco([](const caf::SRSliceProxy* slc) {
  double Eres = varEnergyRes(slc);
  return (Eres<EresMax);
});

const Cut cutRecoContained([](const caf::SRSliceProxy* slc) {
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  return mtc==1;
});

const Cut cutRecoNotContained([](const caf::SRSliceProxy* slc) {
  int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
  return mtc==0;
});

double MuonEresMax = -0.25;

const Cut cutGoodMuonReco([](const caf::SRSliceProxy* slc) {
  double Eres = varEresMuon(slc);
  return (Eres>=MuonEresMax);
});

const Cut cutBadMuonReco([](const caf::SRSliceProxy* slc) {
  double Eres = varEresMuon(slc);
  return (Eres<MuonEresMax);
});

// DATA FILES
void GoodBad(const std::string inputName = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf")

{
  SpectrumLoader loader(inputName);

// SPECTRA
  const Binning bins1 = Binning::Simple(300, 0, 3);
  const Binning bins1b = Binning::Simple(300, -0.5, 0.5);
  const Binning bins1c = Binning::Simple(300, -3, 0);
  const Binning bins2 = Binning::Simple(300, 0, 2);
  const Binning bins3 = Binning::Simple(300, -2, 1);
  const Binning bins4 = Binning::Simple(300, -2, 1);
  const Binning bins4b = Binning::Simple(50, -2, 1);
  const Binning bins5 = Binning::Simple(50, -400, 400);
  const Binning bins5b = Binning::Simple(300, -5, 5);
  const Binning bins5c = Binning::Simple(300, -250, 250);
  const Binning bins6 = Binning::Simple(50, -200, 200);
  const Binning bins6b = Binning::Simple(300, -5, 5);
  const Binning bins6c = Binning::Simple(300, -250, 250);
  const Binning bins7 = Binning::Simple(50, -1000, 1000);
  const Binning bins7b = Binning::Simple(300, -20, 20);
  const Binning bins7c = Binning::Simple(300, -200, 1300);
  const Binning bins9 = Binning::Simple(50, -1, 1);
  const Binning bins9b = Binning::Simple(300, -0.2, 0.2);

  Spectrum sMuonPCounter("", bins1, loader, varMuonPCounter, kNoSpillCut, cutNuMu);
  Spectrum sRangeP_good("rangep", bins1, loader, varRangeP_good, kNoSpillCut, cutNuMu);
  Spectrum sRangeP_bad("rangep", bins1, loader, varRangeP_bad, kNoSpillCut, cutNuMu && cutBadMuonReco);
  Spectrum sMCSP_good("MCSP", bins1, loader, varMCSP_good, kNoSpillCut, cutNuMu);
  Spectrum sMCSP_bad("MCSP", bins1, loader, varMCSP_bad, kNoSpillCut, cutNuMu && cutBadMuonReco);
  Spectrum sTrueP_good1("true p", bins1, loader, varTrueP_good1, kNoSpillCut, cutNuMu);
  Spectrum sTrueP_bad1("true p", bins1, loader, varTrueP_bad1, kNoSpillCut, cutNuMu && cutBadMuonReco);
  Spectrum sTrueP_good0("true p", bins1, loader, varTrueP_good0, kNoSpillCut, cutNuMu);
  Spectrum sTrueP_bad0("true p", bins1, loader, varTrueP_bad0, kNoSpillCut, cutNuMu && cutBadMuonReco);
  Spectrum sPres_good1("residual p", bins1b, loader, varPres_good1, kNoSpillCut, cutNuMu); 
  Spectrum sPres_bad1("residual p", bins1c, loader, varPres_bad1, kNoSpillCut, cutNuMu && cutBadMuonReco);
  Spectrum sPres_good0("residual p", bins1b, loader, varPres_good0, kNoSpillCut, cutNuMu);
  Spectrum sPres_bad0("residual p", bins1c, loader, varPres_bad0, kNoSpillCut, cutNuMu && cutBadMuonReco);

  Spectrum sEmuon_good("Reconstructed Muon Energy", bins2, loader, varEmuon, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sTrueEmuon_good("True Muon Energy", bins2, loader, varTrueEmuon, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sEresMuon_good("Residual Muon Energy", bins3, loader, varEresMuon, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaCosTheta_good("reco cos(theta) - true cos(theta)", bins9b, loader, varDeltaCosTheta, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaStartx_good("reco start.x - true start.x", bins5b, loader, varDeltaStartx, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaStarty_good("reco start.y - true start.y", bins6b, loader, varDeltaStarty, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaStartz_good("reco start.z - true start.z", bins7b, loader, varDeltaStartz, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaEndx_good("reco end.x - true end.x", bins5b, loader, varDeltaEndx, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaEndy_good("reco end.y - true end.y", bins6b, loader, varDeltaEndy, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaEndz_good("reco end.z - true end.z", bins7b, loader, varDeltaEndz, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaEndStartx_good("reco end.x - true start.x", bins5c, loader, varDeltaEndStartx, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaEndStarty_good("reco end.y - true start.y", bins6c, loader, varDeltaEndStarty, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sDeltaEndStartz_good("reco end.z - true start.z", bins7c, loader, varDeltaEndStartz, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axCosTheta_good("cos(theta)", bins9, varCosTheta);
  const HistAxis axTrueCosTheta_good("cos(theta)", bins9, varTrueCosTheta);
  Spectrum sCosTheta_2D_good(loader, axCosTheta_good, axTrueCosTheta_good, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axStartx_good("start.x", bins5, varStartx);
  const HistAxis axTrueStartx_good("start.x", bins5, varTrueStartx);
  Spectrum sStartx_2D_good(loader, axStartx_good, axTrueStartx_good, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axStarty_good("start.y", bins6, varStarty);
  const HistAxis axTrueStarty_good("start.y", bins6, varTrueStarty);
  Spectrum sStarty_2D_good(loader, axStarty_good, axTrueStarty_good, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axStartz_good("start.z", bins7, varStartz);
  const HistAxis axTrueStartz_good("start.z", bins7, varTrueStartz);
  Spectrum sStartz_2D_good(loader, axStartz_good, axTrueStartz_good, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axEndx_good("end.x", bins5, varEndx);
  const HistAxis axTrueEndx_good("end.x", bins5, varTrueEndx);
  Spectrum sEndx_2D_good(loader, axEndx_good, axTrueEndx_good, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axEndy_good("end.y", bins6, varEndy);
  const HistAxis axTrueEndy_good("end.y", bins6, varTrueEndy);
  Spectrum sEndy_2D_good(loader, axEndy_good, axTrueEndy_good, kNoSpillCut, cutNuMu && cutGoodReco);
  const HistAxis axEndz_good("end.z", bins7, varEndz);
  const HistAxis axTrueEndz_good("end.z", bins7, varTrueEndz);
  Spectrum sEndz_2D_good(loader, axEndz_good, axTrueEndz_good, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sTrueStartRecoEndx_2D_good(loader, axEndx_good, axTrueStartx_good, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sTrueStartRecoEndy_2D_good(loader, axEndy_good, axTrueStarty_good, kNoSpillCut, cutNuMu && cutGoodReco);
  Spectrum sTrueStartRecoEndz_2D_good(loader, axEndz_good, axTrueStartz_good, kNoSpillCut, cutNuMu && cutGoodReco);

  Spectrum sEmuon_bad("Reconstructed Muon Energy", bins2, loader, varEmuon, kNoSpillCut, cutNuMu && cutBadReco);
  Spectrum sTrueEmuon_bad("True Muon Energy", bins2, loader, varTrueEmuon, kNoSpillCut, cutNuMu && cutBadReco);
  Spectrum sEresMuon_bad("Residual Muon Energy", bins3, loader, varEresMuon, kNoSpillCut, cutNuMu && cutBadReco);
  Spectrum sDeltaCosTheta_bad("reco cos(theta) - true cos(theta)", bins9b, loader, varDeltaCosTheta, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaStartx_bad("reco start.x - true start.x", bins5b, loader, varDeltaStartx, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaStarty_bad("reco start.y - true start.y", bins6b, loader, varDeltaStarty, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaStartz_bad("reco start.z - true start.z", bins7b, loader, varDeltaStartz, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaEndx_bad("reco end.x - true end.x", bins5b, loader, varDeltaEndx, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaEndy_bad("reco end.y - true end.y", bins6b, loader, varDeltaEndy, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaEndz_bad("reco end.z - true end.z", bins7b, loader, varDeltaEndz, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaEndStartx_bad("reco end.x - true start.x", bins5c, loader, varDeltaEndStartx, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaEndStarty_bad("reco end.y - true start.y", bins6c, loader, varDeltaEndStarty, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sDeltaEndStartz_bad("reco end.z - true start.z", bins7c, loader, varDeltaEndStartz, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axCosTheta_bad("cos(theta)", bins9, varCosTheta);
  const HistAxis axTrueCosTheta_bad("cos(theta)", bins9, varTrueCosTheta);
  Spectrum sCosTheta_2D_bad(loader, axCosTheta_bad, axTrueCosTheta_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axStartx_bad("start.x", bins5, varStartx);
  const HistAxis axTrueStartx_bad("start.x", bins5, varTrueStartx);
  Spectrum sStartx_2D_bad(loader, axStartx_bad, axTrueStartx_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axStarty_bad("start.y", bins6, varStarty);
  const HistAxis axTrueStarty_bad("start.y", bins6, varTrueStarty);
  Spectrum sStarty_2D_bad(loader, axStarty_bad, axTrueStarty_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axStartz_bad("start.z", bins7, varStartz);
  const HistAxis axTrueStartz_bad("start.z", bins7, varTrueStartz);
  Spectrum sStartz_2D_bad(loader, axStartz_bad, axTrueStartz_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axEndx_bad("end.x", bins5, varEndx);
  const HistAxis axTrueEndx_bad("end.x", bins5, varTrueEndx);
  Spectrum sEndx_2D_bad(loader, axEndx_bad, axTrueEndx_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axEndy_bad("end.y", bins6, varEndy);
  const HistAxis axTrueEndy_bad("end.y", bins6, varTrueEndy);
  Spectrum sEndy_2D_bad(loader, axEndy_bad, axTrueEndy_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  const HistAxis axEndz_bad("end.z", bins7, varEndz);
  const HistAxis axTrueEndz_bad("end.z", bins7, varTrueEndz);
  Spectrum sEndz_2D_bad(loader, axEndz_bad, axTrueEndz_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sTrueStartRecoEndx_2D_bad(loader, axEndx_bad, axTrueStartx_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sTrueStartRecoEndy_2D_bad(loader, axEndy_bad, axTrueStarty_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);
  Spectrum sTrueStartRecoEndz_2D_bad(loader, axEndz_bad, axTrueStartz_bad, kNoSpillCut, cutNuMu && cutBadReco && cutBadMuonReco);

  loader.Go();

  std::cout<<"count_good="<<count_good<<"\n";
  std::cout<<"countMCSP_good="<<countMCSP_good<<"\n";
  std::cout<<"countRangep_good="<<countRangep_good<<"\n";
  std::cout<<"count_bad="<<count_bad<<"\n";
  std::cout<<"countMCSP_bad="<<countMCSP_bad<<"\n";
  std::cout<<"countRangep_bad="<<countRangep_bad<<"\n";

// HISTOGRAMS
  // plot events with residual Enu < -0.5 GeV and residual Emu < -0.25 GeV 
  // plot reconstructed and true muon energy
  TCanvas* c1 = new TCanvas();
  TH1* hEmuon_bad = sEmuon_bad.ToTH1(6.6e20);
  hEmuon_bad->SetTitle("residual Enu < -0.5 GeV");
  hEmuon_bad->GetXaxis()->SetTitle("Energy (GeV)");
  hEmuon_bad->SetLineColor(kTeal);
  hEmuon_bad->Draw("hist");
  TH1* hTrueEmuon_bad = sTrueEmuon_bad.ToTH1(6.6e20);
  hTrueEmuon_bad->SetLineColor(kAzure);
  hTrueEmuon_bad->Draw("same hist");

  TLegend* leg1 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg1->AddEntry(hEmuon_bad, "Reconstructed Energy", "l");
  leg1->AddEntry(hTrueEmuon_bad, "True Energy", "l");
  leg1->Draw();
  c1->Print("GoodBad/Emuon_bad.png");

  // plot residual muon energy
  TCanvas* c2 = new TCanvas();
  TH1* hEresMuon_bad = sEresMuon_bad.ToTH1(6.6e20);
  hEresMuon_bad->SetTitle("residual Enu < -0.5 GeV");
  hEresMuon_bad->GetXaxis()->SetTitle("Residual Energy (GeV)");
  hEresMuon_bad->SetLineColor(kViolet);
  hEresMuon_bad->Draw("hist");
  c2->Print("GoodBad/EresMuon_bad.png");

  // 2D plot of cos(theta) vs true cos(theta)
  TCanvas* c4 = new TCanvas();
  c4->SetLogz();
  TH2* hCosTheta_2D_bad = sCosTheta_2D_bad.ToTH2(6.6e20);
  hCosTheta_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hCosTheta_2D_bad->GetXaxis()->SetTitle("reco cos(theta)");
  hCosTheta_2D_bad->GetYaxis()->SetTitle("true cos(theta)");
  hCosTheta_2D_bad->Draw("colz");
  c4->Print("GoodBad/CosTheta_2D_bad.png");

  TCanvas* c4b = new TCanvas();
  TH1* hDeltaCosTheta_bad = sDeltaCosTheta_bad.ToTH1(6.6e20);
  hDeltaCosTheta_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaCosTheta_bad->GetXaxis()->SetTitle("reco cos(theta) - true cos(theta)");
  hDeltaCosTheta_bad->SetLineColor(kOrange);
  hDeltaCosTheta_bad->Draw("hist");
  c4b->Print("GoodBad/DeltaCosTheta_bad.png");

  // 2D plot of start.x vs true start.x
  TCanvas* c5 = new TCanvas();
  c5->SetLogz();
  TH2* hStartx_2D_bad = sStartx_2D_bad.ToTH2(6.6e20);
  hStartx_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hStartx_2D_bad->GetXaxis()->SetTitle("reco start.x");
  hStartx_2D_bad->GetYaxis()->SetTitle("true start.x");
  hStartx_2D_bad->Draw("colz");
  c5->Print("GoodBad/Startx_2D_bad.png");

  TCanvas* c5b = new TCanvas();
  TH1* hDeltaStartx_bad = sDeltaStartx_bad.ToTH1(6.6e20);
  hDeltaStartx_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaStartx_bad->GetXaxis()->SetTitle("reco start.x - true start.x");
  hDeltaStartx_bad->SetLineColor(kSpring);
  hDeltaStartx_bad->Draw("hist");
  c5b->Print("GoodBad/DeltaStartx_bad.png");

  // 2D plot of start.y vs true start.y
  TCanvas* c6 = new TCanvas();
  c6->SetLogz();
  TH2* hStarty_2D_bad = sStarty_2D_bad.ToTH2(6.6e20);
  hStarty_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hStarty_2D_bad->GetXaxis()->SetTitle("reco start.y");
  hStarty_2D_bad->GetYaxis()->SetTitle("true start.y");
  hStarty_2D_bad->Draw("colz");
  c6->Print("GoodBad/Starty_2D_bad.png");

  TCanvas* c6b = new TCanvas();
  TH1* hDeltaStarty_bad = sDeltaStarty_bad.ToTH1(6.6e20);
  hDeltaStarty_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaStarty_bad->GetXaxis()->SetTitle("reco start.y - true start.y");
  hDeltaStarty_bad->SetLineColor(kTeal);
  hDeltaStarty_bad->Draw("hist");
  c6b->Print("GoodBad/DeltaStarty_bad.png");

  // 2D plot of start.z vs true start.z
  TCanvas* c7 = new TCanvas();
  c7->SetLogz();
  TH2* hStartz_2D_bad = sStartz_2D_bad.ToTH2(6.6e20);
  hStartz_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hStartz_2D_bad->GetXaxis()->SetTitle("reco start.z");
  hStartz_2D_bad->GetYaxis()->SetTitle("true start.z");
  hStartz_2D_bad->Draw("colz");
  c7->Print("GoodBad/Startz_2D_bad.png");

  TCanvas* c7b = new TCanvas();
  TH1* hDeltaStartz_bad = sDeltaStartz_bad.ToTH1(6.6e20);
  hDeltaStartz_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaStartz_bad->GetXaxis()->SetTitle("reco start.z - true start.z");
  hDeltaStartz_bad->SetLineColor(kTeal);
  hDeltaStartz_bad->Draw("hist");
  c7b->Print("GoodBad/DeltaStartz_bad.png");

  // 2D plot of end.x vs true end.x
  TCanvas* c11 = new TCanvas();
  c11->SetLogz();
  TH2* hEndx_2D_bad = sEndx_2D_bad.ToTH2(6.6e20);
  hEndx_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hEndx_2D_bad->GetXaxis()->SetTitle("reco end.x");
  hEndx_2D_bad->GetYaxis()->SetTitle("true end.x");
  hEndx_2D_bad->Draw("colz");
  c11->Print("GoodBad/Endx_2D_bad.png");

  TCanvas* c11b = new TCanvas();
  TH1* hDeltaEndx_bad = sDeltaEndx_bad.ToTH1(6.6e20);
  hDeltaEndx_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaEndx_bad->GetXaxis()->SetTitle("reco end.x - true end.x");
  hDeltaEndx_bad->SetLineColor(kAzure);
  hDeltaEndx_bad->Draw("hist");
  c11b->Print("GoodBad/DeltaEndx_bad.png");

  // 2D plot of end.y vs true end.y
  TCanvas* c12 = new TCanvas();
  c12->SetLogz();
  TH2* hEndy_2D_bad = sEndy_2D_bad.ToTH2(6.6e20);
  hEndy_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hEndy_2D_bad->GetXaxis()->SetTitle("reco end.y");
  hEndy_2D_bad->GetYaxis()->SetTitle("true end.y");
  hEndy_2D_bad->Draw("colz");
  c12->Print("GoodBad/Endy_2D_bad.png");

  TCanvas* c12b = new TCanvas();
  TH1* hDeltaEndy_bad = sDeltaEndy_bad.ToTH1(6.6e20);
  hDeltaEndy_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaEndy_bad->GetXaxis()->SetTitle("reco end.y - true end.y");
  hDeltaEndy_bad->SetLineColor(kViolet);
  hDeltaEndy_bad->Draw("hist");
  c12b->Print("GoodBad/DeltaEndy_bad.png");

  // 2D plot of end.z vs true end.z
  TCanvas* c13 = new TCanvas();
  c13->SetLogz();
  TH2* hEndz_2D_bad = sEndz_2D_bad.ToTH2(6.6e20);
  hEndz_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hEndz_2D_bad->GetXaxis()->SetTitle("reco end.z");
  hEndz_2D_bad->GetYaxis()->SetTitle("true end.z");
  hEndz_2D_bad->Draw("colz");
  c13->Print("GoodBad/Endz_2D_bad.png");

  TCanvas* c13b = new TCanvas();
  TH1* hDeltaEndz_bad = sDeltaEndz_bad.ToTH1(6.6e20);
  hDeltaEndz_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaEndz_bad->GetXaxis()->SetTitle("reco end.z - true end.z");
  hDeltaEndz_bad->SetLineColor(kPink);
  hDeltaEndz_bad->Draw("hist");
  c13b->Print("GoodBad/DeltaEndz_bad.png");

  // 2D plot of end.x vs true start.x
  TCanvas* c17 = new TCanvas();
  c17->SetLogz();
  TH2* hTrueStartRecoEndx_2D_bad = sTrueStartRecoEndx_2D_bad.ToTH2(6.6e20);
  hTrueStartRecoEndx_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hTrueStartRecoEndx_2D_bad->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2D_bad->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2D_bad->Draw("colz");
  c17->Print("GoodBad/TrueStartRecoEndx_2D_bad.png");

  TCanvas* c17b = new TCanvas();
  TH1* hDeltaEndStartx_bad = sDeltaEndStartx_bad.ToTH1(6.6e20);
  hDeltaEndStartx_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaEndStartx_bad->GetXaxis()->SetTitle("reco end.x - true start.x");
  hDeltaEndStartx_bad->SetLineColor(kOrange);
  hDeltaEndStartx_bad->Draw("hist");
  c17b->Print("GoodBad/DeltaEndStartx_bad.png");

  // 2D plot of end.y vs true start.y
  TCanvas* c18 = new TCanvas();
  c18->SetLogz();
  TH2* hTrueStartRecoEndy_2D_bad = sTrueStartRecoEndy_2D_bad.ToTH2(6.6e20);
  hTrueStartRecoEndy_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hTrueStartRecoEndy_2D_bad->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2D_bad->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2D_bad->Draw("colz");
  c18->Print("GoodBad/TrueStartRecoEndy_2D_bad.png");

  TCanvas* c18b = new TCanvas();
  TH1* hDeltaEndStarty_bad = sDeltaEndStarty_bad.ToTH1(6.6e20);
  hDeltaEndStarty_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaEndStarty_bad->GetXaxis()->SetTitle("reco end.y - true start.y");
  hDeltaEndStarty_bad->SetLineColor(kSpring);
  hDeltaEndStarty_bad->Draw("hist");
  c18b->Print("GoodBad/DeltaEndStarty_bad.png");

  // 2D plot of end.z vs true start.z
  TCanvas* c19 = new TCanvas();
  c19->SetLogz();
  TH2* hTrueStartRecoEndz_2D_bad = sTrueStartRecoEndz_2D_bad.ToTH2(6.6e20);
  hTrueStartRecoEndz_2D_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hTrueStartRecoEndz_2D_bad->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2D_bad->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2D_bad->Draw("colz");
  c19->Print("GoodBad/TrueStartRecoEndz_2D_bad.png");

  TCanvas* c19b = new TCanvas();
  TH1* hDeltaEndStartz_bad = sDeltaEndStartz_bad.ToTH1(6.6e20);
  hDeltaEndStartz_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hDeltaEndStartz_bad->GetXaxis()->SetTitle("reco end.z - true start.z");
  hDeltaEndStartz_bad->SetLineColor(kTeal);
  hDeltaEndStartz_bad->Draw("hist");
  c19b->Print("GoodBad/DeltaEndStartz_bad.png");

  // plot events with Eres > -0.5 GeV
  TCanvas* c23 = new TCanvas();
  TH1* hEmuon_good = sEmuon_good.ToTH1(6.6e20);
  hEmuon_good->SetTitle("Eres > -0.5 GeV");
  hEmuon_good->GetXaxis()->SetTitle("Energy (GeV)");
  hEmuon_good->SetLineColor(kTeal);
  hEmuon_good->Draw("hist");
  TH1* hTrueEmuon_good = sTrueEmuon_good.ToTH1(6.6e20);
  hTrueEmuon_good->SetLineColor(kAzure);
  hTrueEmuon_good->Draw("same hist");

  TLegend* leg23 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg23->AddEntry(hEmuon_good, "Reconstructed Energy", "l");
  leg23->AddEntry(hTrueEmuon_good, "True Energy", "l");
  leg23->Draw();
  c23->Print("GoodBad/Emuon_good.png");

  TCanvas* c24 = new TCanvas();
  TH1* hEresMuon_good = sEresMuon_good.ToTH1(6.6e20);
  hEresMuon_good->SetTitle("Eres > -0.5 GeV");
  hEresMuon_good->GetXaxis()->SetTitle("Residual Energy (GeV)");
  hEresMuon_good->SetLineColor(kViolet);
  hEresMuon_good->Draw("hist");
  c24->Print("GoodBad/EresMuon_good.png");

  TCanvas* c26 = new TCanvas();
  c26->SetLogz();
  TH2* hCosTheta_2D_good = sCosTheta_2D_good.ToTH2(6.6e20);
  hCosTheta_2D_good->GetXaxis()->SetTitle("reco cos(theta)");
  hCosTheta_2D_good->GetYaxis()->SetTitle("true cos(theta)");
  hCosTheta_2D_good->SetTitle("Eres > -0.5 GeV");
  hCosTheta_2D_good->Draw("colz");
  c26->Print("GoodBad/CosTheta_2D_good.png");

  TCanvas* c26b = new TCanvas();
  TH1* hDeltaCosTheta_good = sDeltaCosTheta_good.ToTH1(6.6e20);
  hDeltaCosTheta_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaCosTheta_good->GetXaxis()->SetTitle("reco cos(theta) - true cos(theta)");
  hDeltaCosTheta_good->SetLineColor(kOrange);
  hDeltaCosTheta_good->Draw("hist");
  c26b->Print("GoodBad/DeltaCosTheta_good.png");

  TCanvas* c27 = new TCanvas();
  c27->SetLogz();
  TH2* hStartx_2D_good = sStartx_2D_good.ToTH2(6.6e20);
  hStartx_2D_good->SetTitle("Eres > -0.5 GeV");
  hStartx_2D_good->GetXaxis()->SetTitle("reco start.x");
  hStartx_2D_good->GetYaxis()->SetTitle("true start.x");
  hStartx_2D_good->Draw("colz");
  c27->Print("GoodBad/Startx_2D_good.png");

  TCanvas* c27b = new TCanvas();
  TH1* hDeltaStartx_good = sDeltaStartx_good.ToTH1(6.6e20);
  hDeltaStartx_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaStartx_good->GetXaxis()->SetTitle("reco start.x - true start.x");
  hDeltaStartx_good->SetLineColor(kSpring);
  hDeltaStartx_good->Draw("hist");
  c27b->Print("GoodBad/DeltaStartx_good.png");

  TCanvas* c28 = new TCanvas();
  c28->SetLogz();
  TH2* hStarty_2D_good = sStarty_2D_good.ToTH2(6.6e20);
  hStarty_2D_good->SetTitle("Eres > -0.5 GeV");
  hStarty_2D_good->GetXaxis()->SetTitle("reco start.y");
  hStarty_2D_good->GetYaxis()->SetTitle("true start.y");
  hStarty_2D_good->Draw("colz");
  c28->Print("GoodBad/Starty_2D_good.png");

  TCanvas* c28b = new TCanvas();
  TH1* hDeltaStarty_good = sDeltaStarty_good.ToTH1(6.6e20);
  hDeltaStarty_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaStarty_good->GetXaxis()->SetTitle("reco start.y - true start.y");
  hDeltaStarty_good->SetLineColor(kTeal);
  hDeltaStarty_good->Draw("hist");
  c28b->Print("GoodBad/DeltaStarty_good.png");

  TCanvas* c29 = new TCanvas();
  c29->SetLogz();
  TH2* hStartz_2D_good = sStartz_2D_good.ToTH2(6.6e20);
  hStartz_2D_good->SetTitle("Eres > -0.5 GeV");
  hStartz_2D_good->GetXaxis()->SetTitle("reco start.z");
  hStartz_2D_good->GetYaxis()->SetTitle("true start.z");
  hStartz_2D_good->Draw("colz");
  c29->Print("GoodBad/Startz_2D_good.png");

  TCanvas* c29b = new TCanvas();
  TH1* hDeltaStartz_good = sDeltaStartz_good.ToTH1(6.6e20);
  hDeltaStartz_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaStartz_good->GetXaxis()->SetTitle("reco start.z - true start.z");
  hDeltaStartz_good->SetLineColor(kTeal);
  hDeltaStartz_good->Draw("hist");
  c29b->Print("GoodBad/DeltaStartz_good.png");

  TCanvas* c33 = new TCanvas();
  c33->SetLogz();
  TH2* hEndx_2D_good = sEndx_2D_good.ToTH2(6.6e20);
  hEndx_2D_good->SetTitle("Eres > -0.5 GeV");
  hEndx_2D_good->GetXaxis()->SetTitle("reco end.x");
  hEndx_2D_good->GetYaxis()->SetTitle("true end.x");
  hEndx_2D_good->Draw("colz");
  c33->Print("GoodBad/Endx_2D_good.png");

  TCanvas* c33b = new TCanvas();
  TH1* hDeltaEndx_good = sDeltaEndx_good.ToTH1(6.6e20);
  hDeltaEndx_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaEndx_good->GetXaxis()->SetTitle("reco end.x - true end.x");
  hDeltaEndx_good->SetLineColor(kAzure);
  hDeltaEndx_good->Draw("hist");
  c33b->Print("GoodBad/DeltaEndx_good.png");

  TCanvas* c34 = new TCanvas();
  TH2* hEndy_2D_good = sEndy_2D_good.ToTH2(6.6e20);
  c34->SetLogz();
  hEndy_2D_good->SetTitle("Eres > -0.5 GeV");
  hEndy_2D_good->GetXaxis()->SetTitle("reco end.y");
  hEndy_2D_good->GetYaxis()->SetTitle("true end.y");
  hEndy_2D_good->Draw("colz");
  c34->Print("GoodBad/Endy_2D_good.png");

  TCanvas* c34b = new TCanvas();
  TH1* hDeltaEndy_good = sDeltaEndy_good.ToTH1(6.6e20);
  hDeltaEndy_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaEndy_good->GetXaxis()->SetTitle("reco end.y - true end.y");
  hDeltaEndy_good->SetLineColor(kViolet);
  hDeltaEndy_good->Draw("hist");
  c34b->Print("GoodBad/DeltaEndy_good.png");

  TCanvas* c35 = new TCanvas();
  c35->SetLogz();
  TH2* hEndz_2D_good = sEndz_2D_good.ToTH2(6.6e20);
  hEndz_2D_good->SetTitle("Eres > -0.5 GeV");
  hEndz_2D_good->GetXaxis()->SetTitle("reco end.z");
  hEndz_2D_good->GetYaxis()->SetTitle("true end.z");
  hEndz_2D_good->Draw("colz");
  c35->Print("GoodBad/Endz_2D_good.png");

  TCanvas* c35b = new TCanvas();
  TH1* hDeltaEndz_good = sDeltaEndz_good.ToTH1(6.6e20);
  hDeltaEndz_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaEndz_good->GetXaxis()->SetTitle("reco end.z - true end.z");
  hDeltaEndz_good->SetLineColor(kPink);
  hDeltaEndz_good->Draw("hist");
  c35b->Print("GoodBad/DeltaEndz_good.png");

  TCanvas* c39 = new TCanvas();
  c39->SetLogz();
  TH2* hTrueStartRecoEndx_2D_good = sTrueStartRecoEndx_2D_good.ToTH2(6.6e20);
  hTrueStartRecoEndx_2D_good->SetTitle("Eres > -0.5 GeV");
  hTrueStartRecoEndx_2D_good->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2D_good->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2D_good->Draw("colz");
  c39->Print("GoodBad/TrueStartRecoEndx_2D_good.png");

  TCanvas* c39b = new TCanvas();
  TH1* hDeltaEndStartx_good = sDeltaEndStartx_good.ToTH1(6.6e20);
  hDeltaEndStartx_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaEndStartx_good->GetXaxis()->SetTitle("reco end.x - true start.x");
  hDeltaEndStartx_good->SetLineColor(kOrange);
  hDeltaEndStartx_good->Draw("hist");
  c39b->Print("GoodBad/DeltaEndStartx_good.png");

  TCanvas* c40 = new TCanvas();
  TH2* hTrueStartRecoEndy_2D_good = sTrueStartRecoEndy_2D_good.ToTH2(6.6e20);
  c40->SetLogz();
  hTrueStartRecoEndy_2D_good->SetTitle("Eres > -0.5 GeV");
  hTrueStartRecoEndy_2D_good->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2D_good->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2D_good->Draw("colz");
  c40->Print("GoodBad/TrueStartRecoEndy_2D_good.png");

  TCanvas* c40b = new TCanvas();
  TH1* hDeltaEndStarty_good = sDeltaEndStarty_good.ToTH1(6.6e20);
  hDeltaEndStarty_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaEndStarty_good->GetXaxis()->SetTitle("reco end.y - true start.y");
  hDeltaEndStarty_good->SetLineColor(kSpring);
  hDeltaEndStarty_good->Draw("hist");
  c40b->Print("GoodBad/DeltaEndStarty_good.png");

  TCanvas* c41 = new TCanvas();
  c41->SetLogz();
  TH2* hTrueStartRecoEndz_2D_good = sTrueStartRecoEndz_2D_good.ToTH2(6.6e20);
  hTrueStartRecoEndz_2D_good->SetTitle("Eres > -0.5 GeV");
  hTrueStartRecoEndz_2D_good->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2D_good->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2D_good->Draw("colz");
  c41->Print("GoodBad/TrueStartRecoEndz_2D_good.png");

  TCanvas* c41b = new TCanvas();
  TH1* hDeltaEndStartz_good = sDeltaEndStartz_good.ToTH1(6.6e20);
  hDeltaEndStartz_good->SetTitle("residual Enu > -0.5 GeV");
  hDeltaEndStartz_good->GetXaxis()->SetTitle("reco end.z - true start.z");
  hDeltaEndStartz_good->SetLineColor(kTeal);
  hDeltaEndStartz_good->Draw("hist");
  c41b->Print("GoodBad/DeltaEndStartz_good.png");

  // plots to compare rangeP to MCSP
  TCanvas* c42 = new TCanvas();
  TH1* hRangeP_good = sRangeP_good.ToTH1(6.6e20);
  hRangeP_good->SetTitle("Eres > -0.5 GeV");
  hRangeP_good->GetXaxis()->SetTitle("Momentum (GeV)");
  hRangeP_good->SetLineColor(kViolet);
  hRangeP_good->Draw("hist");
  TH1* hTrueP_good1 = sTrueP_good1.ToTH1(6.6e20);
  hTrueP_good1->SetLineColor(kPink);
  hTrueP_good1->Draw("same hist");

  TLegend* leg42 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg42->AddEntry(hRangeP_good, "rangeP", "l");
  leg42->AddEntry(hTrueP_good1, "true momentum", "l");
  leg42->Draw();
  c42->Print("GoodBad/RangeP_good.png");

  TCanvas* c43 = new TCanvas();
  TH1* hRangeP_bad = sRangeP_bad.ToTH1(6.6e20);
  hRangeP_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hRangeP_bad->GetXaxis()->SetTitle("Momentum (GeV)");
  hRangeP_bad->SetLineColor(kOrange);
  hRangeP_bad->Draw("hist");
  TH1* hTrueP_bad1 = sTrueP_bad1.ToTH1(6.6e20);
  hTrueP_bad1->SetLineColor(kPink);
  hTrueP_bad1->Draw("same hist");

  TLegend* leg43 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg43->AddEntry(hRangeP_bad, "rangeP", "l");
  leg43->AddEntry(hTrueP_bad1, "true momentum", "l");
  leg43->Draw();
  c43->Print("GoodBad/RangeP_bad.png");

  TCanvas* c44 = new TCanvas();
  c44->SetTitle("Eres > -0.5 GeV");
  TH1* hPres_good1 = sPres_good1.ToTH1(6.6e20);
  hPres_good1->SetTitle("Eres > -0.5 GeV");
  hPres_good1->GetXaxis()->SetTitle("rangeP Residual (GeV)");
  hPres_good1->SetLineColor(kViolet);
  hPres_good1->Draw("hist");
  c44->Print("GoodBad/Pres1_good.png");

  TCanvas* c44b = new TCanvas();
  c44b->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  TH1* hPres_bad1 = sPres_bad1.ToTH1(6.6e20);
  hPres_bad1->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hPres_bad1->GetXaxis()->SetTitle("rangeP Residual (GeV)");
  hPres_bad1->SetLineColor(kOrange);
  hPres_bad1->Draw("hist");
  c44b->Print("GoodBad/Pres1_bad.png");

  TCanvas* c45 = new TCanvas();
  TH1* hMCSP_good = sMCSP_good.ToTH1(6.6e20);
  hMCSP_good->SetTitle("Eres > -0.5 GeV");
  hMCSP_good->GetXaxis()->SetTitle("Momentum (GeV)");
  hMCSP_good->SetLineColor(kSpring);
  hMCSP_good->Draw("hist");
  TH1* hTrueP_good0 = sTrueP_good0.ToTH1(6.6e20);
  hTrueP_good0->SetLineColor(kAzure);
  hTrueP_good0->Draw("same hist");

  TLegend* leg45 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg45->AddEntry(hMCSP_good, "MCSP", "l");
  leg45->AddEntry(hTrueP_good0, "true momentum", "l");
  leg45->Draw();
  c45->Print("GoodBad/MCSP_good.png");

  TCanvas* c46 = new TCanvas();
  TH1* hMCSP_bad = sMCSP_bad.ToTH1(6.6e20);
  hMCSP_bad->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hMCSP_bad->GetXaxis()->SetTitle("Momentum (GeV)");
  hMCSP_bad->SetLineColor(kTeal);
  hMCSP_bad->Draw("hist");
  TH1* hTrueP_bad0 = sTrueP_bad0.ToTH1(6.6e20);
  hTrueP_bad0->SetLineColor(kAzure);
  hTrueP_bad0->Draw("same hist");

  TLegend* leg46 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg46->AddEntry(hMCSP_bad, "MCSP", "l");
  leg46->AddEntry(hTrueP_bad0, "true momentum", "l");
  leg46->Draw();
  c46->Print("GoodBad/MCSP_bad.png");

  TCanvas* c47 = new TCanvas();
  c47->SetTitle("Eres > -0.5 GeV");
  TH1* hPres_good0 = sPres_good0.ToTH1(6.6e20);
  hPres_good0->SetTitle("Eres > -0.5 GeV");
  hPres_good0->GetXaxis()->SetTitle("MCSP Residual (GeV)");
  hPres_good0->SetLineColor(kSpring);
  hPres_good0->Draw("hist");
  c47->Print("GoodBad/Pres0_good.png");

  TCanvas* c47b = new TCanvas();
  c47b->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  TH1* hPres_bad0 = sPres_bad0.ToTH1(6.6e20);
  hPres_bad0->SetTitle("residual Enu < -0.5 GeV and residual Emu < -0.25 GeV");
  hPres_bad0->GetXaxis()->SetTitle("MCSP Residual (GeV)");
  hPres_bad0->SetLineColor(kTeal);
  hPres_bad0->Draw("hist");
  c47b->Print("GoodBad/Pres0_bad.png");
}
