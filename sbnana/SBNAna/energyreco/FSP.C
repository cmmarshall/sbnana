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

// count0: number of events with Eres >= -0.5
int count0=0;
// count1: number of events with Eres<-0.5 and 1 trk, 0 shws
int count1=0;
// count2: number of events with Eres<-0.5 and 2 trks, 0 shws
int count2=0;
// count3: number of other events with Eres<-0.5
int count3=0;
// count4: number of events with Eres<-0.5
int count4=0;
// count5: number of events contained in reco, but not contained in truth, with Eres<-0.5 and 1 trk, 0 shws
int count5=0;
// count6: number of events contained in truth, but not contained in reco, with Eres<-0.5 and 1 trk, 0 shws
int count6=0;
// count7: used for debugging
int count7=0;
// count8: number of events contained in reco, but not contained in truth, with Eres<-0.5 and 2 trks, 0 shws
int count8=0;
// count9: number of events contained in truth, but not contained in reco, with Eres<-0.5 and 2 trks, 0 shws
int count9=0;
// count10: used for debugging
int count10=0;
// count11: number of events contained in reco, but not contained in truth, with Eres<-0.5 and 3+ trks and/or 1+ shws
int count11=0;
// count12: number of events contained in truth, but not contained in reco, with Eres<-0.5 and 3+ trks and/or 1+ shws
int count12=0;
// count13: used for debugging
int count13=0;
// count14: total number of events
int count14=0;
// count15: number of events in which RFiducial is true and TFiducial is false
int count15=0;
// count16: number of events in which RFiducial is false and TFiducial is true
int count16=0;
// count17: used for debugging
int count17=0;

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
double EresMax = -0.5;

// used for saving final state particle codes
ofstream outputfile1;
ofstream outputfile2;
ofstream outputfile3;

const MultiVar varFSP1([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> FSP;
  double Eres = varEnergyRes(slc);
  if (Eres < EresMax) {
    int ntrk = slc->reco.ntrk;
    int nshw = slc->reco.nshw;
    if (ntrk==1 and nshw==0) {
      int ndaughters = slc->truth.nprim;
      for (int i=0; i<ndaughters; i++) {
        int daughter = slc->truth.prim[i].pdg; 
        double energy = slc->truth.prim[i].genE;
        FSP.push_back(daughter);
        outputfile1.open("FSP1.txt", ios::app);
        outputfile1<<daughter<<" (E="<<energy<<"), ";
        outputfile1.close();
      }
      outputfile1.open("FSP1.txt", ios::app);
      outputfile1<<"\n";
      outputfile1.close();
    }}
  return FSP;
});

const MultiVar varFSP2([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> FSP;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    double Eres = varEnergyRes(slc);
    if (Eres < EresMax) {
      int ndaughters = slc->truth.nprim;
      for (int i=0; i<ndaughters; i++) {
        int daughter = slc->truth.prim[i].pdg; 
        double energy = slc->truth.prim[i].genE;
        FSP.push_back(daughter);
        outputfile2.open("FSP2.txt", ios::app);
        outputfile2<<daughter<<" (E="<<energy<<"), ";
        outputfile2.close();
      }
      outputfile2.open("FSP2.txt", ios::app);
      outputfile2<<"\n";
      outputfile2.close();
    }}
  return FSP;
});

const MultiVar varFSP3([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> FSP;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if ( (ntrk==1 and nshw>0) or (ntrk==2 and nshw>0) or ntrk>=3 ) {
    double Eres = varEnergyRes(slc);
    if (Eres < EresMax) {
      int ndaughters = slc->truth.nprim;
      for (int i=0; i<ndaughters; i++) {
        int daughter = slc->truth.prim[i].pdg; 
        double energy = slc->truth.prim[i].genE;
        FSP.push_back(daughter);
        outputfile3.open("FSP3.txt", ios::app);
        outputfile3<<daughter<<" (E="<<energy<<"), ";
        outputfile3.close();
      }
      outputfile3.open("FSP3.txt", ios::app);
      outputfile3<<"\n";
      outputfile3.close();
    }}
  return FSP;
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

const Var varStarty([](const caf::SRSliceProxy* slc) -> double {
  double start_y=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_y = slc->reco.trk[muon_index].start.y;
  }
  return start_y;
});

const Var varStartz([](const caf::SRSliceProxy* slc) -> double {
  double start_z=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    start_z = slc->reco.trk[muon_index].start.z;
  }
  return start_z;
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

const Var varEndy([](const caf::SRSliceProxy* slc) -> double {
  double end_y=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_y = slc->reco.trk[muon_index].end.y;
  }
  return end_y;
});

const Var varEndz([](const caf::SRSliceProxy* slc) -> double {
  double end_z=-999;
  int muon_index = varMuonTrackIndex(slc);
  if (muon_index >= 0) {
    end_z = slc->reco.trk[muon_index].end.z;
  }
  return end_z;
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

const Var varCounter([](const caf::SRSliceProxy* slc) -> int {
  double Eres = varEnergyRes(slc);
  if (Eres >= EresMax) {
    count0+=1;
  }
  else {
    count4+=1;
    int ntrk = slc->reco.ntrk;
    int nshw = slc->reco.nshw;
    int RecoContained = varMuonTrackContainedness(slc);
    int TrueContained = varTrueMuonTrackContainedness(slc);
    if (ntrk==1 and nshw==0) {
      count1+=1;
      if (RecoContained==1 and TrueContained==0) {
        count5+=1;
      }
      else if (RecoContained==0 and TrueContained==1) {
        count6+=1;
      }
      else {
        count7+=1;
      }}
    else if (ntrk==2 and nshw==0) {
      count2+=1;
      if (RecoContained==1 and TrueContained==0) {
        count8+=1;
      }
      else if (RecoContained==0 and TrueContained==1) {
        count9+=1;
      }
      else {
        count10+=1;
      }}
    else {
      count3+=1;
      if (RecoContained==1 and TrueContained==0) {
        count11+=1;
      }
      else if (RecoContained==0 and TrueContained==1) {
        count12+=1;
      }
      else {
        count13+=1;
      }}}
  bool RFiducial = kRFiducial(slc);
  bool TFiducial = kTFiducial(slc);
  count14+=1;
  if (RFiducial and !TFiducial) {
    count15+=1;
  }
  else if (!RFiducial and TFiducial) {
    count16+=1;
  }
  else {
    count17+=1;
  }
  return 6; // my favorite number!
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

const Cut cut1trk([](const caf::SRSliceProxy* slc) {
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  return ( ntrk==1 && nshw==0 ); 
});

const Cut cut2trk([](const caf::SRSliceProxy* slc) {
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  return ( ntrk==2 && nshw==0 ); 
});

const Cut cut3trk([](const caf::SRSliceProxy* slc) {
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  return ( (ntrk==1 && nshw>0) or (ntrk==2 && nshw>0) or ntrk>=3 ); 
});

const Cut cutGoodReco([](const caf::SRSliceProxy* slc) {
  double Eres = varEnergyRes(slc);
  return (Eres>=EresMax);
});

const Cut cutBadReco([](const caf::SRSliceProxy* slc) {
  double Eres = varEnergyRes(slc);
  return (Eres<EresMax);
});

const Cut cutTrueContainedNotRecoContained([](const caf::SRSliceProxy* slc) {
  int RecoContained = varMuonTrackContainedness(slc); 
  int TrueContained = varTrueMuonTrackContainedness(slc);
  return (RecoContained==0 && TrueContained==1);
});

const Cut cutNotTrueContainedRecoContained([](const caf::SRSliceProxy* slc) {
  int RecoContained = varMuonTrackContainedness(slc); 
  int TrueContained = varTrueMuonTrackContainedness(slc);
  return (RecoContained==1 && TrueContained==0);
});

// DATA FILES
void FSP(const std::string inputName = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf")

{
  SpectrumLoader loader(inputName);

// SPECTRA
  const Binning bins1 = Binning::Simple(3000, -3000, 3000);
  const Binning bins2 = Binning::Simple(300, 0, 2);
  const Binning bins3 = Binning::Simple(300, -2, 1);
  const Binning bins4 = Binning::Simple(100, -1, 1);
  const Binning bins9 = Binning::Simple(50, -1, 1);
  const Binning bins5 = Binning::Simple(50, -400, 400);
  const Binning bins6 = Binning::Simple(50, -200, 200);
  const Binning bins7 = Binning::Simple(50, -1000, 1000);
  const Binning bins8 = Binning::Simple(200, 0, 1000);

  Spectrum sFSP1("Final State Particles", bins1, loader, varFSP1, kNoSpillCut, cutNuMu);
  Spectrum sFSP2("Final State Particles", bins1, loader, varFSP2, kNoSpillCut, cutNuMu);
  Spectrum sFSP3("Final State Particles", bins1, loader, varFSP3, kNoSpillCut, cutNuMu);
  Spectrum sCounter("", bins1, loader, varCounter, kNoSpillCut, cutNuMu);

  Spectrum sEmuon_1trk("Reconstructed Muon Energy", bins2, loader, varEmuon, kNoSpillCut, cutNuMu && cut1trk && cutBadReco);
  Spectrum sTrueEmuon_1trk("True Muon Energy", bins2, loader, varTrueEmuon, kNoSpillCut, cutNuMu && cut1trk && cutBadReco);
  Spectrum sEresMuon_1trk("Residual Muon Energy", bins3, loader, varEresMuon, kNoSpillCut, cutNuMu && cut1trk && cutBadReco);
  Spectrum sTrueCosTheta_Bad_1trk("cos(theta)", bins4, loader, varTrueCosTheta, kNoSpillCut, cutNuMu && cut1trk && cutBadReco);
  Spectrum sTrueCosTheta_Good_1trk("cos(theta)", bins4, loader, varTrueCosTheta, kNoSpillCut, cutNuMu && cut1trk && cutGoodReco);
  const HistAxis axCosTheta_1trk("cos(theta)", bins9, varCosTheta);
  const HistAxis axTrueCosTheta_1trk("cos(theta)", bins9, varTrueCosTheta);
  Spectrum sCosTheta_2D_1trk(loader, axCosTheta_1trk, axTrueCosTheta_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco);
  const HistAxis axStartx_1trk("start.x", bins5, varStartx);
  const HistAxis axTrueStartx_1trk("start.x", bins5, varTrueStartx);
  Spectrum sStartx_2D_1trk(loader, axStartx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStartx_2Db_1trk(loader, axStartx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axStarty_1trk("start.y", bins6, varStarty);
  const HistAxis axTrueStarty_1trk("start.y", bins6, varTrueStarty);
  Spectrum sStarty_2D_1trk(loader, axStarty_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStarty_2Db_1trk(loader, axStarty_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axStartz_1trk("start.z", bins7, varStartz);
  const HistAxis axTrueStartz_1trk("start.z", bins7, varTrueStartz);
  Spectrum sStartz_2D_1trk(loader, axStartz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStartz_2Db_1trk(loader, axStartz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndx_1trk("end.x", bins5, varEndx);
  const HistAxis axTrueEndx_1trk("end.x", bins5, varTrueEndx);
  Spectrum sEndx_2D_1trk(loader, axEndx_1trk, axTrueEndx_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndx_2Db_1trk(loader, axEndx_1trk, axTrueEndx_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndy_1trk("end.y", bins6, varEndy);
  const HistAxis axTrueEndy_1trk("end.y", bins6, varTrueEndy);
  Spectrum sEndy_2D_1trk(loader, axEndy_1trk, axTrueEndy_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndy_2Db_1trk(loader, axEndy_1trk, axTrueEndy_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndz_1trk("end.z", bins7, varEndz);
  const HistAxis axTrueEndz_1trk("end.z", bins7, varTrueEndz);
  Spectrum sEndz_2D_1trk(loader, axEndz_1trk, axTrueEndz_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndz_2Db_1trk(loader, axEndz_1trk, axTrueEndz_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndx_2D_1trk(loader, axEndx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndx_2Db_1trk(loader, axEndx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndy_2D_1trk(loader, axEndy_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndy_2Db_1trk(loader, axEndy_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndz_2D_1trk(loader, axEndz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndz_2Db_1trk(loader, axEndz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut1trk && cutBadReco && cutNotTrueContainedRecoContained);

  Spectrum sEmuon_2trk("Reconstructed Muon Energy", bins2, loader, varEmuon, kNoSpillCut, cutNuMu && cut2trk && cutBadReco);
  Spectrum sTrueEmuon_2trk("True Muon Energy", bins2, loader, varTrueEmuon, kNoSpillCut, cutNuMu && cut2trk && cutBadReco);
  Spectrum sEresMuon_2trk("Residual Muon Energy", bins3, loader, varEresMuon, kNoSpillCut, cutNuMu && cut2trk && cutBadReco);
  Spectrum sTrueCosTheta_Bad_2trk("cos(theta)", bins4, loader, varTrueCosTheta, kNoSpillCut, cutNuMu && cut2trk && cutBadReco);
  Spectrum sTrueCosTheta_Good_2trk("cos(theta)", bins4, loader, varTrueCosTheta, kNoSpillCut, cutNuMu && cut2trk && cutGoodReco);
  const HistAxis axCosTheta_2trk("cos(theta)", bins9, varCosTheta);
  const HistAxis axTrueCosTheta_2trk("cos(theta)", bins9, varTrueCosTheta);
  Spectrum sCosTheta_2D_2trk(loader, axCosTheta_2trk, axTrueCosTheta_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco);
  const HistAxis axStartx_2trk("start.x", bins5, varStartx);
  const HistAxis axTrueStartx_2trk("start.x", bins5, varTrueStartx);
  Spectrum sStartx_2D_2trk(loader, axStartx_2trk, axTrueStartx_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStartx_2Db_2trk(loader, axStartx_2trk, axTrueStartx_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axStarty_2trk("start.y", bins6, varStarty);
  const HistAxis axTrueStarty_2trk("start.y", bins6, varTrueStarty);
  Spectrum sStarty_2D_2trk(loader, axStarty_2trk, axTrueStarty_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStarty_2Db_2trk(loader, axStarty_2trk, axTrueStarty_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axStartz_2trk("start.z", bins7, varStartz);
  const HistAxis axTrueStartz_2trk("start.z", bins7, varTrueStartz);
  Spectrum sStartz_2D_2trk(loader, axStartz_2trk, axTrueStartz_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStartz_2Db_2trk(loader, axStartz_2trk, axTrueStartz_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndx_2trk("end.x", bins5, varEndx);
  const HistAxis axTrueEndx_2trk("end.x", bins5, varTrueEndx);
  Spectrum sEndx_2D_2trk(loader, axEndx_2trk, axTrueEndx_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndx_2Db_2trk(loader, axEndx_2trk, axTrueEndx_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndy_2trk("end.y", bins6, varEndy);
  const HistAxis axTrueEndy_2trk("end.y", bins6, varTrueEndy);
  Spectrum sEndy_2D_2trk(loader, axEndy_2trk, axTrueEndy_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndy_2Db_2trk(loader, axEndy_2trk, axTrueEndy_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndz_2trk("end.z", bins7, varEndz);
  const HistAxis axTrueEndz_2trk("end.z", bins7, varTrueEndz);
  Spectrum sEndz_2D_2trk(loader, axEndz_2trk, axTrueEndz_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndz_2Db_2trk(loader, axEndz_2trk, axTrueEndz_2trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndx_2D_2trk(loader, axEndx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndx_2Db_2trk(loader, axEndx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndy_2D_2trk(loader, axEndy_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndy_2Db_2trk(loader, axEndy_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndz_2D_2trk(loader, axEndz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndz_2Db_2trk(loader, axEndz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut2trk && cutBadReco && cutNotTrueContainedRecoContained);

  Spectrum sEmuon_3trk("Reconstructed Muon Energy", bins2, loader, varEmuon, kNoSpillCut, cutNuMu && cut3trk && cutBadReco);
  Spectrum sTrueEmuon_3trk("True Muon Energy", bins2, loader, varTrueEmuon, kNoSpillCut, cutNuMu && cut3trk && cutBadReco);
  Spectrum sEresMuon_3trk("Residual Muon Energy", bins3, loader, varEresMuon, kNoSpillCut, cutNuMu && cut3trk && cutBadReco);
  Spectrum sTrueCosTheta_Bad_3trk("cos(theta)", bins4, loader, varTrueCosTheta, kNoSpillCut, cutNuMu && cut3trk && cutBadReco);
  Spectrum sTrueCosTheta_Good_3trk("cos(theta)", bins4, loader, varTrueCosTheta, kNoSpillCut, cutNuMu && cut3trk && cutGoodReco);
  const HistAxis axCosTheta_3trk("cos(theta)", bins9, varCosTheta);
  const HistAxis axTrueCosTheta_3trk("cos(theta)", bins9, varTrueCosTheta);
  Spectrum sCosTheta_2D_3trk(loader, axCosTheta_3trk, axTrueCosTheta_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco);
  const HistAxis axStartx_3trk("start.x", bins5, varStartx);
  const HistAxis axTrueStartx_3trk("start.x", bins5, varTrueStartx);
  Spectrum sStartx_2D_3trk(loader, axStartx_3trk, axTrueStartx_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStartx_2Db_3trk(loader, axStartx_3trk, axTrueStartx_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axStarty_3trk("start.y", bins6, varStarty);
  const HistAxis axTrueStarty_3trk("start.y", bins6, varTrueStarty);
  Spectrum sStarty_2D_3trk(loader, axStarty_3trk, axTrueStarty_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStarty_2Db_3trk(loader, axStarty_3trk, axTrueStarty_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axStartz_3trk("start.z", bins7, varStartz);
  const HistAxis axTrueStartz_3trk("start.z", bins7, varTrueStartz);
  Spectrum sStartz_2D_3trk(loader, axStartz_3trk, axTrueStartz_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sStartz_2Db_3trk(loader, axStartz_3trk, axTrueStartz_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndx_3trk("end.x", bins5, varEndx);
  const HistAxis axTrueEndx_3trk("end.x", bins5, varTrueEndx);
  Spectrum sEndx_2D_3trk(loader, axEndx_3trk, axTrueEndx_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndx_2Db_3trk(loader, axEndx_3trk, axTrueEndx_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndy_3trk("end.y", bins6, varEndy);
  const HistAxis axTrueEndy_3trk("end.y", bins6, varTrueEndy);
  Spectrum sEndy_2D_3trk(loader, axEndy_3trk, axTrueEndy_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndy_2Db_3trk(loader, axEndy_3trk, axTrueEndy_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  const HistAxis axEndz_3trk("end.z", bins7, varEndz);
  const HistAxis axTrueEndz_3trk("end.z", bins7, varTrueEndz);
  Spectrum sEndz_2D_3trk(loader, axEndz_3trk, axTrueEndz_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sEndz_2Db_3trk(loader, axEndz_3trk, axTrueEndz_3trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndx_2D_3trk(loader, axEndx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndx_2Db_3trk(loader, axEndx_1trk, axTrueStartx_1trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndy_2D_3trk(loader, axEndy_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndy_2Db_3trk(loader, axEndy_1trk, axTrueStarty_1trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);
  Spectrum sTrueStartRecoEndz_2D_3trk(loader, axEndz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutTrueContainedNotRecoContained);
  Spectrum sTrueStartRecoEndz_2Db_3trk(loader, axEndz_1trk, axTrueStartz_1trk, kNoSpillCut, cutNuMu && cut3trk && cutBadReco && cutNotTrueContainedRecoContained);

  // clear files of data written during the last run of this code
  outputfile1.open("FSP1.txt", ios::trunc);
  outputfile1.close();
  outputfile2.open("FSP2.txt", ios::trunc);
  outputfile2.close();
  outputfile3.open("FSP3.txt", ios::trunc);
  outputfile3.close();

  loader.Go();

  // count0: number of events with Eres >= -0.5
  std::cout<<"count0 = "<<count0<<"\n";
  // count1: number of events with Eres<-0.5 and 1 trk, 0 shws
  std::cout<<"count1 = "<<count1<<"\n";
  // count2: number of events with Eres<-0.5 and 2 trks, 0 shws
  std::cout<<"count2 = "<<count2<<"\n";
  // count3: number of other events with Eres<-0.5
  std::cout<<"count3 = "<<count3<<"\n";
  // count4: number of events with Eres<-0.5
  std::cout<<"count4 = "<<count4<<"\n";
  // count5: number of events contained in reco, but not contained in truth, with Eres<-0.5 and 1 trk, 0 shws
  std::cout<<"count5 = "<<count5<<"\n";
  // count6: number of events contained in truth, but not contained in reco, with Eres<-0.5 and 1 trk, 0 shws
  std::cout<<"count6 = "<<count6<<"\n";
  // count7: used for debugging
  std::cout<<"count7 = "<<count7<<"\n";
  // count8: number of events contained in reco, but not contained in truth, with Eres<-0.5 and 2 trks, 0 shws
  std::cout<<"count8 = "<<count8<<"\n";
  // count9: number of events contained in truth, but not contained in reco, with Eres<-0.5 and 2 trks, 0 shws
  std::cout<<"count9 = "<<count9<<"\n";
  // count10: used for debugging
  std::cout<<"count10 = "<<count10<<"\n";
  // count11: number of events contained in reco, but not contained in truth, with Eres<-0.5 and 3+ trks and/or 1+ shws
  std::cout<<"count11 = "<<count11<<"\n";
  // count12: number of events contained in truth, but not contained in reco, with Eres<-0.5 and 3+ trks and/or 1+ shws
  std::cout<<"count12 = "<<count12<<"\n";
  // count13: used for debugging
  std::cout<<"count13 = "<<count13<<"\n";
  // count14: total number of events
  std::cout<<"count14 = "<<count14<<"\n";
  // count15: number of events in which RFiducial is true and TFiducial is false
  std::cout<<"count15 = "<<count15<<"\n";
  // count16: number of events in which RFiducial is false and TFiducial is true
  std::cout<<"count16 = "<<count16<<"\n";
  // count17: used for debugging 
  std::cout<<"count17 = "<<count17<<"\n";

// HISTOGRAMS
  // plot events with 1 trk and 0 shws
  // plot reconstructed and true muon energy
  TCanvas* c1 = new TCanvas();
  TH1* hEmuon_1trk = sEmuon_1trk.ToTH1(6.6e20);
  hEmuon_1trk->GetXaxis()->SetTitle("Energy (GeV)");
  hEmuon_1trk->SetLineColor(kTeal);
  hEmuon_1trk->Draw("hist");
  TH1* hTrueEmuon_1trk = sTrueEmuon_1trk.ToTH1(6.6e20);
  hTrueEmuon_1trk->SetLineColor(kAzure);
  hTrueEmuon_1trk->Draw("same hist");

  TLegend* leg1 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg1->AddEntry(hEmuon_1trk, "Reconstructed Energy", "l");
  leg1->AddEntry(hTrueEmuon_1trk, "True Energy", "l");
  leg1->Draw();
  c1->Print("FSP/Emuon_1trk.png");

  // plot residual muon energy
  TCanvas* c2 = new TCanvas();
  TH1* hEresMuon_1trk = sEresMuon_1trk.ToTH1(6.6e20);
  hEresMuon_1trk->GetXaxis()->SetTitle("Residual Energy (GeV)");
  hEresMuon_1trk->SetLineColor(kViolet);
  hEresMuon_1trk->Draw("hist");
  c2->Print("FSP/EresMuon_1trk.png");

  // plot cos(theta) for muons
  TCanvas* c3 = new TCanvas();
  TH1* hTrueCosTheta_Good_1trk = sTrueCosTheta_Good_1trk.ToTH1(6.6e20);
  hTrueCosTheta_Good_1trk->GetXaxis()->SetTitle("cos(theta)");
  hTrueCosTheta_Good_1trk->SetLineColor(kPink);
  hTrueCosTheta_Good_1trk->Draw("hist");
  TH1* hTrueCosTheta_Bad_1trk = sTrueCosTheta_Bad_1trk.ToTH1(6.6e20);
  hTrueCosTheta_Bad_1trk->SetLineColor(kOrange);
  hTrueCosTheta_Bad_1trk->Draw("same hist");

  TLegend* leg3 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg3->AddEntry(hTrueCosTheta_Good_1trk, "true cos(theta) for muons w/ Eres>-0.5", "l");
  leg3->AddEntry(hTrueCosTheta_Bad_1trk, "true cos(theta) for muons w/Eres<-0.5", "l");
  leg3->Draw();
  c3->Print("FSP/CosTheta_1trk.png");

  // 2D plot of cos(theta) vs true cos(theta)
  TCanvas* c4 = new TCanvas();
  TH2* hCosTheta_2D_1trk = sCosTheta_2D_1trk.ToTH2(6.6e20);
  hCosTheta_2D_1trk->GetXaxis()->SetTitle("reco cos(theta)");
  hCosTheta_2D_1trk->GetYaxis()->SetTitle("true cos(theta)");
  hCosTheta_2D_1trk->Draw("colz");
  c4->Print("FSP/CosTheta_2D_1trk.png");

  // 2D plot of start.x vs true start.x
  // contained in truth, not in reco
  TCanvas* c5 = new TCanvas();
  TH2* hStartx_2D_1trk = sStartx_2D_1trk.ToTH2(6.6e20);
  hStartx_2D_1trk->SetTitle("contained in truth, not in reco");
  hStartx_2D_1trk->GetXaxis()->SetTitle("reco start.x");
  hStartx_2D_1trk->GetYaxis()->SetTitle("true start.x");
  hStartx_2D_1trk->Draw("colz");
  c5->Print("FSP/Startx_2D_1trk.png");

  // 2D plot of start.y vs true start.y
  // contained in truth, not in reco
  TCanvas* c6 = new TCanvas();
  TH2* hStarty_2D_1trk = sStarty_2D_1trk.ToTH2(6.6e20);
  hStarty_2D_1trk->SetTitle("contained in truth, not in reco");
  hStarty_2D_1trk->GetXaxis()->SetTitle("reco start.y");
  hStarty_2D_1trk->GetYaxis()->SetTitle("true start.y");
  hStarty_2D_1trk->Draw("colz");
  c6->Print("FSP/Starty_2D_1trk.png");

  // 2D plot of start.z vs true start.z
  // contained in truth, not in reco
  TCanvas* c7 = new TCanvas();
  TH2* hStartz_2D_1trk = sStartz_2D_1trk.ToTH2(6.6e20);
  hStartz_2D_1trk->SetTitle("contained in truth, not in reco");
  hStartz_2D_1trk->GetXaxis()->SetTitle("reco start.z");
  hStartz_2D_1trk->GetYaxis()->SetTitle("true start.z");
  hStartz_2D_1trk->Draw("colz");
  c7->Print("FSP/Startz_2D_1trk.png");

  // 2D plot of start.x vs true start.x
  // contained in reco, not in truth
  TCanvas* c8 = new TCanvas();
  TH2* hStartx_2Db_1trk = sStartx_2Db_1trk.ToTH2(6.6e20);
  hStartx_2Db_1trk->SetTitle("contained in reco, not in truth");
  hStartx_2Db_1trk->GetXaxis()->SetTitle("reco start.x");
  hStartx_2Db_1trk->GetYaxis()->SetTitle("true start.x");
  hStartx_2Db_1trk->Draw("colz");
  c8->Print("FSP/Startx_2Db_1trk.png");

  // 2D plot of start.y vs true start.y
  // contained in reco, not in truth
  TCanvas* c9 = new TCanvas();
  TH2* hStarty_2Db_1trk = sStarty_2Db_1trk.ToTH2(6.6e20);
  hStarty_2Db_1trk->SetTitle("contained in reco, not in truth");
  hStarty_2Db_1trk->GetXaxis()->SetTitle("reco start.y");
  hStarty_2Db_1trk->GetYaxis()->SetTitle("true start.y");
  hStarty_2Db_1trk->Draw("colz");
  c9->Print("FSP/Starty_2Db_1trk.png");

  // 2D plot of start.z vs true start.z
  // contained in reco, not in truth
  TCanvas* c10 = new TCanvas();
  TH2* hStartz_2Db_1trk = sStartz_2Db_1trk.ToTH2(6.6e20);
  hStartz_2Db_1trk->SetTitle("contained in reco, not in truth");
  hStartz_2Db_1trk->GetXaxis()->SetTitle("reco start.z");
  hStartz_2Db_1trk->GetYaxis()->SetTitle("true start.z");
  hStartz_2Db_1trk->Draw("colz");
  c10->Print("FSP/Startz_2Db_1trk.png");

  // 2D plot of end.x vs true end.x
  // contained in truth, not in reco
  TCanvas* c11 = new TCanvas();
  TH2* hEndx_2D_1trk = sEndx_2D_1trk.ToTH2(6.6e20);
  hEndx_2D_1trk->SetTitle("contained in truth, not in reco");
  hEndx_2D_1trk->GetXaxis()->SetTitle("reco end.x");
  hEndx_2D_1trk->GetYaxis()->SetTitle("true end.x");
  hEndx_2D_1trk->Draw("colz");
  c11->Print("FSP/Endx_2D_1trk.png");

  // 2D plot of end.y vs true end.y
  // contained in truth, not in reco
  TCanvas* c12 = new TCanvas();
  TH2* hEndy_2D_1trk = sEndy_2D_1trk.ToTH2(6.6e20);
  hEndy_2D_1trk->SetTitle("contained in truth, not in reco");
  hEndy_2D_1trk->GetXaxis()->SetTitle("reco end.y");
  hEndy_2D_1trk->GetYaxis()->SetTitle("true end.y");
  hEndy_2D_1trk->Draw("colz");
  c12->Print("FSP/Endy_2D_1trk.png");

  // 2D plot of end.z vs true end.z
  // contained in truth, not in reco
  TCanvas* c13 = new TCanvas();
  TH2* hEndz_2D_1trk = sEndz_2D_1trk.ToTH2(6.6e20);
  hEndz_2D_1trk->SetTitle("contained in truth, not in reco");
  hEndz_2D_1trk->GetXaxis()->SetTitle("reco end.z");
  hEndz_2D_1trk->GetYaxis()->SetTitle("true end.z");
  hEndz_2D_1trk->Draw("colz");
  c13->Print("FSP/Endz_2D_1trk.png");

  // 2D plot of end.x vs true end.x
  // contained in reco, not in truth
  TCanvas* c14 = new TCanvas();
  TH2* hEndx_2Db_1trk = sEndx_2Db_1trk.ToTH2(6.6e20);
  hEndx_2Db_1trk->SetTitle("contained in reco, not in truth");
  hEndx_2Db_1trk->GetXaxis()->SetTitle("reco end.x");
  hEndx_2Db_1trk->GetYaxis()->SetTitle("true end.x");
  hEndx_2Db_1trk->Draw("colz");
  c14->Print("FSP/Endx_2Db_1trk.png");

  // 2D plot of end.y vs true end.y
  // contained in reco, not in truth
  TCanvas* c15 = new TCanvas();
  TH2* hEndy_2Db_1trk = sEndy_2Db_1trk.ToTH2(6.6e20);
  hEndy_2Db_1trk->SetTitle("contained in reco, not in truth");
  hEndy_2Db_1trk->GetXaxis()->SetTitle("reco end.y");
  hEndy_2Db_1trk->GetYaxis()->SetTitle("true end.y");
  hEndy_2Db_1trk->Draw("colz");
  c15->Print("FSP/Endy_2Db_1trk.png");

  // 2D plot of end.z vs true end.z
  // contained in reco, not in truth
  TCanvas* c16 = new TCanvas();
  TH2* hEndz_2Db_1trk = sEndz_2Db_1trk.ToTH2(6.6e20);
  hEndz_2Db_1trk->SetTitle("contained in reco, not in truth");
  hEndz_2Db_1trk->GetXaxis()->SetTitle("reco end.z");
  hEndz_2Db_1trk->GetYaxis()->SetTitle("true end.z");
  hEndz_2Db_1trk->Draw("colz");
  c16->Print("FSP/Endz_2Db_1trk.png");

  // 2D plot of end.x vs true start.x
  // contained in truth, not in reco
  TCanvas* c17 = new TCanvas();
  TH2* hTrueStartRecoEndx_2D_1trk = sTrueStartRecoEndx_2D_1trk.ToTH2(6.6e20);
  hTrueStartRecoEndx_2D_1trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndx_2D_1trk->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2D_1trk->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2D_1trk->Draw("colz");
  c17->Print("FSP/TrueStartRecoEndx_2D_1trk.png");

  // 2D plot of end.y vs true start.y
  // contained in truth, not in reco
  TCanvas* c18 = new TCanvas();
  TH2* hTrueStartRecoEndy_2D_1trk = sTrueStartRecoEndy_2D_1trk.ToTH2(6.6e20);
  hTrueStartRecoEndy_2D_1trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndy_2D_1trk->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2D_1trk->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2D_1trk->Draw("colz");
  c18->Print("FSP/TrueStartRecoEndy_2D_1trk.png");

  // 2D plot of end.z vs true start.z
  // contained in truth, not in reco
  TCanvas* c19 = new TCanvas();
  TH2* hTrueStartRecoEndz_2D_1trk = sTrueStartRecoEndz_2D_1trk.ToTH2(6.6e20);
  hTrueStartRecoEndz_2D_1trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndz_2D_1trk->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2D_1trk->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2D_1trk->Draw("colz");
  c19->Print("FSP/TrueStartRecoEndz_2D_1trk.png");

  // 2D plot of end.x vs true start.x
  // contained in reco, not in truth
  TCanvas* c20 = new TCanvas();
  TH2* hTrueStartRecoEndx_2Db_1trk = sTrueStartRecoEndx_2Db_1trk.ToTH2(6.6e20);
  hTrueStartRecoEndx_2Db_1trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndx_2Db_1trk->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2Db_1trk->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2Db_1trk->Draw("colz");
  c20->Print("FSP/TrueStartRecoEndx_2Db_1trk.png");

  // 2D plot of end.y vs true start.y
  // contained in reco, not in truth
  TCanvas* c21 = new TCanvas();
  TH2* hTrueStartRecoEndy_2Db_1trk = sTrueStartRecoEndy_2Db_1trk.ToTH2(6.6e20);
  hTrueStartRecoEndy_2Db_1trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndy_2Db_1trk->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2Db_1trk->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2Db_1trk->Draw("colz");
  c21->Print("FSP/TrueStartRecoEndy_2Db_1trk.png");

  // 2D plot of end.z vs true start.z
  // contained in reco, not in truth
  TCanvas* c22 = new TCanvas();
  TH2* hTrueStartRecoEndz_2Db_1trk = sTrueStartRecoEndz_2Db_1trk.ToTH2(6.6e20);
  hTrueStartRecoEndz_2Db_1trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndz_2Db_1trk->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2Db_1trk->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2Db_1trk->Draw("colz");
  c22->Print("FSP/TrueStartRecoEndz_2Db_1trk.png");

// plot events with 2 trks and 0 shws
  TCanvas* c23 = new TCanvas();
  TH1* hEmuon_2trk = sEmuon_2trk.ToTH1(6.6e20);
  hEmuon_2trk->GetXaxis()->SetTitle("Energy (GeV)");
  hEmuon_2trk->SetLineColor(kTeal);
  hEmuon_2trk->Draw("hist");
  TH1* hTrueEmuon_2trk = sTrueEmuon_2trk.ToTH1(6.6e20);
  hTrueEmuon_2trk->SetLineColor(kAzure);
  hTrueEmuon_2trk->Draw("same hist");

  TLegend* leg23 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg23->AddEntry(hEmuon_2trk, "Reconstructed Energy", "l");
  leg23->AddEntry(hTrueEmuon_2trk, "True Energy", "l");
  leg23->Draw();
  c23->Print("FSP/Emuon_2trk.png");

  TCanvas* c24 = new TCanvas();
  TH1* hEresMuon_2trk = sEresMuon_2trk.ToTH1(6.6e20);
  hEresMuon_2trk->GetXaxis()->SetTitle("Residual Energy (GeV)");
  hEresMuon_2trk->SetLineColor(kViolet);
  hEresMuon_2trk->Draw("hist");
  c24->Print("FSP/EresMuon_2trk.png");

  TCanvas* c25 = new TCanvas();
  TH1* hTrueCosTheta_Good_2trk = sTrueCosTheta_Good_2trk.ToTH1(6.6e20);
  hTrueCosTheta_Good_2trk->GetXaxis()->SetTitle("cos(theta)");
  hTrueCosTheta_Good_2trk->SetLineColor(kPink);
  hTrueCosTheta_Good_2trk->Draw("hist");
  TH1* hTrueCosTheta_Bad_2trk = sTrueCosTheta_Bad_2trk.ToTH1(6.6e20);
  hTrueCosTheta_Bad_2trk->SetLineColor(kOrange);
  hTrueCosTheta_Bad_2trk->Draw("same hist");

  TLegend* leg25 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg25->AddEntry(hTrueCosTheta_Good_2trk, "true cos(theta) for muons w/ Eres>-0.5", "l");
  leg25->AddEntry(hTrueCosTheta_Bad_2trk, "true cos(theta) for muons w/Eres<-0.5", "l");
  leg25->Draw();
  c25->Print("FSP/CosTheta_2trk.png");

  TCanvas* c26 = new TCanvas();
  TH2* hCosTheta_2D_2trk = sCosTheta_2D_2trk.ToTH2(6.6e20);
  hCosTheta_2D_2trk->GetXaxis()->SetTitle("reco cos(theta)");
  hCosTheta_2D_2trk->GetYaxis()->SetTitle("true cos(theta)");
  hCosTheta_2D_2trk->Draw("colz");
  c26->Print("FSP/CosTheta_2D_2trk.png");

  TCanvas* c27 = new TCanvas();
  TH2* hStartx_2D_2trk = sStartx_2D_2trk.ToTH2(6.6e20);
  hStartx_2D_2trk->SetTitle("contained in truth, not in reco");
  hStartx_2D_2trk->GetXaxis()->SetTitle("reco start.x");
  hStartx_2D_2trk->GetYaxis()->SetTitle("true start.x");
  hStartx_2D_2trk->Draw("colz");
  c27->Print("FSP/Startx_2D_2trk.png");

  TCanvas* c28 = new TCanvas();
  TH2* hStarty_2D_2trk = sStarty_2D_2trk.ToTH2(6.6e20);
  hStarty_2D_2trk->SetTitle("contained in truth, not in reco");
  hStarty_2D_2trk->GetXaxis()->SetTitle("reco start.y");
  hStarty_2D_2trk->GetYaxis()->SetTitle("true start.y");
  hStarty_2D_2trk->Draw("colz");
  c28->Print("FSP/Starty_2D_2trk.png");

  TCanvas* c29 = new TCanvas();
  TH2* hStartz_2D_2trk = sStartz_2D_2trk.ToTH2(6.6e20);
  hStartz_2D_2trk->SetTitle("contained in truth, not in reco");
  hStartz_2D_2trk->GetXaxis()->SetTitle("reco start.z");
  hStartz_2D_2trk->GetYaxis()->SetTitle("true start.z");
  hStartz_2D_2trk->Draw("colz");
  c29->Print("FSP/Startz_2D_2trk.png");

  TCanvas* c30 = new TCanvas();
  TH2* hStartx_2Db_2trk = sStartx_2Db_2trk.ToTH2(6.6e20);
  hStartx_2Db_2trk->SetTitle("contained in reco, not in truth");
  hStartx_2Db_2trk->GetXaxis()->SetTitle("reco start.x");
  hStartx_2Db_2trk->GetYaxis()->SetTitle("true start.x");
  hStartx_2Db_2trk->Draw("colz");
  c30->Print("FSP/Startx_2Db_2trk.png");

  TCanvas* c31 = new TCanvas();
  TH2* hStarty_2Db_2trk = sStarty_2Db_2trk.ToTH2(6.6e20);
  hStarty_2Db_2trk->SetTitle("contained in reco, not in truth");
  hStarty_2Db_2trk->GetXaxis()->SetTitle("reco start.y");
  hStarty_2Db_2trk->GetYaxis()->SetTitle("true start.y");
  hStarty_2Db_2trk->Draw("colz");
  c31->Print("FSP/Starty_2Db_2trk.png");

  TCanvas* c32 = new TCanvas();
  TH2* hStartz_2Db_2trk = sStartz_2Db_2trk.ToTH2(6.6e20);
  hStartz_2Db_2trk->SetTitle("contained in reco, not in truth");
  hStartz_2Db_2trk->GetXaxis()->SetTitle("reco start.z");
  hStartz_2Db_2trk->GetYaxis()->SetTitle("true start.z");
  hStartz_2Db_2trk->Draw("colz");
  c32->Print("FSP/Startz_2Db_2trk.png");

  TCanvas* c33 = new TCanvas();
  TH2* hEndx_2D_2trk = sEndx_2D_2trk.ToTH2(6.6e20);
  hEndx_2D_2trk->SetTitle("contained in truth, not in reco");
  hEndx_2D_2trk->GetXaxis()->SetTitle("reco end.x");
  hEndx_2D_2trk->GetYaxis()->SetTitle("true end.x");
  hEndx_2D_2trk->Draw("colz");
  c33->Print("FSP/Endx_2D_2trk.png");

  TCanvas* c34 = new TCanvas();
  TH2* hEndy_2D_2trk = sEndy_2D_2trk.ToTH2(6.6e20);
  hEndy_2D_2trk->SetTitle("contained in truth, not in reco");
  hEndy_2D_2trk->GetXaxis()->SetTitle("reco end.y");
  hEndy_2D_2trk->GetYaxis()->SetTitle("true end.y");
  hEndy_2D_2trk->Draw("colz");
  c34->Print("FSP/Endy_2D_2trk.png");

  TCanvas* c35 = new TCanvas();
  TH2* hEndz_2D_2trk = sEndz_2D_2trk.ToTH2(6.6e20);
  hEndz_2D_2trk->SetTitle("contained in truth, not in reco");
  hEndz_2D_2trk->GetXaxis()->SetTitle("reco end.z");
  hEndz_2D_2trk->GetYaxis()->SetTitle("true end.z");
  hEndz_2D_2trk->Draw("colz");
  c35->Print("FSP/Endz_2D_2trk.png");

  TCanvas* c36 = new TCanvas();
  TH2* hEndx_2Db_2trk = sEndx_2Db_2trk.ToTH2(6.6e20);
  hEndx_2Db_2trk->SetTitle("contained in reco, not in truth");
  hEndx_2Db_2trk->GetXaxis()->SetTitle("reco end.x");
  hEndx_2Db_2trk->GetYaxis()->SetTitle("true end.x");
  hEndx_2Db_2trk->Draw("colz");
  c36->Print("FSP/Endx_2Db_2trk.png");

  TCanvas* c37 = new TCanvas();
  TH2* hEndy_2Db_2trk = sEndy_2Db_2trk.ToTH2(6.6e20);
  hEndy_2Db_2trk->SetTitle("contained in reco, not in truth");
  hEndy_2Db_2trk->GetXaxis()->SetTitle("reco end.y");
  hEndy_2Db_2trk->GetYaxis()->SetTitle("true end.y");
  hEndy_2Db_2trk->Draw("colz");
  c37->Print("FSP/Endy_2Db_2trk.png");

  TCanvas* c38 = new TCanvas();
  TH2* hEndz_2Db_2trk = sEndz_2Db_2trk.ToTH2(6.6e20);
  hEndz_2Db_2trk->SetTitle("contained in reco, not in truth");
  hEndz_2Db_2trk->GetXaxis()->SetTitle("reco end.z");
  hEndz_2Db_2trk->GetYaxis()->SetTitle("true end.z");
  hEndz_2Db_2trk->Draw("colz");
  c38->Print("FSP/Endz_2Db_2trk.png");

  TCanvas* c39 = new TCanvas();
  TH2* hTrueStartRecoEndx_2D_2trk = sTrueStartRecoEndx_2D_2trk.ToTH2(6.6e20);
  hTrueStartRecoEndx_2D_2trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndx_2D_2trk->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2D_2trk->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2D_2trk->Draw("colz");
  c39->Print("FSP/TrueStartRecoEndx_2D_2trk.png");

  TCanvas* c40 = new TCanvas();
  TH2* hTrueStartRecoEndy_2D_2trk = sTrueStartRecoEndy_2D_2trk.ToTH2(6.6e20);
  hTrueStartRecoEndy_2D_2trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndy_2D_2trk->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2D_2trk->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2D_2trk->Draw("colz");
  c40->Print("FSP/TrueStartRecoEndy_2D_2trk.png");

  TCanvas* c41 = new TCanvas();
  TH2* hTrueStartRecoEndz_2D_2trk = sTrueStartRecoEndz_2D_2trk.ToTH2(6.6e20);
  hTrueStartRecoEndz_2D_2trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndz_2D_2trk->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2D_2trk->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2D_2trk->Draw("colz");
  c41->Print("FSP/TrueStartRecoEndz_2D_2trk.png");

  TCanvas* c42 = new TCanvas();
  TH2* hTrueStartRecoEndx_2Db_2trk = sTrueStartRecoEndx_2Db_2trk.ToTH2(6.6e20);
  hTrueStartRecoEndx_2Db_2trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndx_2Db_2trk->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2Db_2trk->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2Db_2trk->Draw("colz");
  c42->Print("FSP/TrueStartRecoEndx_2Db_2trk.png");

  TCanvas* c43 = new TCanvas();
  TH2* hTrueStartRecoEndy_2Db_2trk = sTrueStartRecoEndy_2Db_2trk.ToTH2(6.6e20);
  hTrueStartRecoEndy_2Db_2trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndy_2Db_2trk->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2Db_2trk->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2Db_2trk->Draw("colz");
  c43->Print("FSP/TrueStartRecoEndy_2Db_2trk.png");

  TCanvas* c44 = new TCanvas();
  TH2* hTrueStartRecoEndz_2Db_2trk = sTrueStartRecoEndz_2Db_2trk.ToTH2(6.6e20);
  hTrueStartRecoEndz_2Db_2trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndz_2Db_2trk->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2Db_2trk->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2Db_2trk->Draw("colz");
  c44->Print("FSP/TrueStartRecoEndz_2Db_2trk.png");

// plot events with 3+ trks and/or 1+ shws
  TCanvas* c45 = new TCanvas();
  TH1* hEmuon_3trk = sEmuon_3trk.ToTH1(6.6e20);
  hEmuon_3trk->GetXaxis()->SetTitle("Energy (GeV)");
  hEmuon_3trk->SetLineColor(kTeal);
  hEmuon_3trk->Draw("hist");
  TH1* hTrueEmuon_3trk = sTrueEmuon_3trk.ToTH1(6.6e20);
  hTrueEmuon_3trk->SetLineColor(kAzure);
  hTrueEmuon_3trk->Draw("same hist");

  TLegend* leg45 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg45->AddEntry(hEmuon_3trk, "Reconstructed Energy", "l");
  leg45->AddEntry(hTrueEmuon_3trk, "True Energy", "l");
  leg45->Draw();
  c45->Print("FSP/Emuon_3trk.png");

  TCanvas* c46 = new TCanvas();
  TH1* hEresMuon_3trk = sEresMuon_3trk.ToTH1(6.6e20);
  hEresMuon_3trk->GetXaxis()->SetTitle("Residual Energy (GeV)");
  hEresMuon_3trk->SetLineColor(kViolet);
  hEresMuon_3trk->Draw("hist");
  c46->Print("FSP/EresMuon_3trk.png");

  TCanvas* c47 = new TCanvas();
  TH1* hTrueCosTheta_Good_3trk = sTrueCosTheta_Good_3trk.ToTH1(6.6e20);
  hTrueCosTheta_Good_3trk->GetXaxis()->SetTitle("cos(theta)");
  hTrueCosTheta_Good_3trk->SetLineColor(kPink);
  hTrueCosTheta_Good_3trk->Draw("hist");
  TH1* hTrueCosTheta_Bad_3trk = sTrueCosTheta_Bad_3trk.ToTH1(6.6e20);
  hTrueCosTheta_Bad_3trk->SetLineColor(kOrange);
  hTrueCosTheta_Bad_3trk->Draw("same hist");

  TLegend* leg47 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg47->AddEntry(hTrueCosTheta_Good_3trk, "true cos(theta) for muons w/ Eres>-0.5", "l");
  leg47->AddEntry(hTrueCosTheta_Bad_3trk, "true cos(theta) for muons w/Eres<-0.5", "l");
  leg47->Draw();
  c47->Print("FSP/CosTheta_3trk.png");

  TCanvas* c48 = new TCanvas();
  TH2* hCosTheta_2D_3trk = sCosTheta_2D_3trk.ToTH2(6.6e20);
  hCosTheta_2D_3trk->GetXaxis()->SetTitle("reco cos(theta)");
  hCosTheta_2D_3trk->GetYaxis()->SetTitle("true cos(theta)");
  hCosTheta_2D_3trk->Draw("colz");
  c48->Print("FSP/CosTheta_2D_3trk.png");

  TCanvas* c49 = new TCanvas();
  TH2* hStartx_2D_3trk = sStartx_2D_3trk.ToTH2(6.6e20);
  hStartx_2D_3trk->SetTitle("contained in truth, not in reco");
  hStartx_2D_3trk->GetXaxis()->SetTitle("reco start.x");
  hStartx_2D_3trk->GetYaxis()->SetTitle("true start.x");
  hStartx_2D_3trk->Draw("colz");
  c49->Print("FSP/Startx_2D_3trk.png");

  TCanvas* c50 = new TCanvas();
  TH2* hStarty_2D_3trk = sStarty_2D_3trk.ToTH2(6.6e20);
  hStarty_2D_3trk->SetTitle("contained in truth, not in reco");
  hStarty_2D_3trk->GetXaxis()->SetTitle("reco start.y");
  hStarty_2D_3trk->GetYaxis()->SetTitle("true start.y");
  hStarty_2D_3trk->Draw("colz");
  c50->Print("FSP/Starty_2D_3trk.png");

  TCanvas* c51 = new TCanvas();
  TH2* hStartz_2D_3trk = sStartz_2D_3trk.ToTH2(6.6e20);
  hStartz_2D_3trk->SetTitle("contained in truth, not in reco");
  hStartz_2D_3trk->GetXaxis()->SetTitle("reco start.z");
  hStartz_2D_3trk->GetYaxis()->SetTitle("true start.z");
  hStartz_2D_3trk->Draw("colz");
  c51->Print("FSP/Startz_2D_3trk.png");

  TCanvas* c52 = new TCanvas();
  TH2* hStartx_2Db_3trk = sStartx_2Db_3trk.ToTH2(6.6e20);
  hStartx_2Db_3trk->SetTitle("contained in reco, not in truth");
  hStartx_2Db_3trk->GetXaxis()->SetTitle("reco start.x");
  hStartx_2Db_3trk->GetYaxis()->SetTitle("true start.x");
  hStartx_2Db_3trk->Draw("colz");
  c52->Print("FSP/Startx_2Db_3trk.png");

  TCanvas* c53 = new TCanvas();
  TH2* hStarty_2Db_3trk = sStarty_2Db_3trk.ToTH2(6.6e20);
  hStarty_2Db_3trk->SetTitle("contained in reco, not in truth");
  hStarty_2Db_3trk->GetXaxis()->SetTitle("reco start.y");
  hStarty_2Db_3trk->GetYaxis()->SetTitle("true start.y");
  hStarty_2Db_3trk->Draw("colz");
  c53->Print("FSP/Starty_2Db_3trk.png");

  TCanvas* c54 = new TCanvas();
  TH2* hStartz_2Db_3trk = sStartz_2Db_3trk.ToTH2(6.6e20);
  hStartz_2Db_3trk->SetTitle("contained in reco, not in truth");
  hStartz_2Db_3trk->GetXaxis()->SetTitle("reco start.z");
  hStartz_2Db_3trk->GetYaxis()->SetTitle("true start.z");
  hStartz_2Db_3trk->Draw("colz");
  c54->Print("FSP/Startz_2Db_3trk.png");

  TCanvas* c55 = new TCanvas();
  TH2* hEndx_2D_3trk = sEndx_2D_3trk.ToTH2(6.6e20);
  hEndx_2D_3trk->SetTitle("contained in truth, not in reco");
  hEndx_2D_3trk->GetXaxis()->SetTitle("reco end.x");
  hEndx_2D_3trk->GetYaxis()->SetTitle("true end.x");
  hEndx_2D_3trk->Draw("colz");
  c55->Print("FSP/Endx_2D_3trk.png");

  TCanvas* c56 = new TCanvas();
  TH2* hEndy_2D_3trk = sEndy_2D_3trk.ToTH2(6.6e20);
  hEndy_2D_3trk->SetTitle("contained in truth, not in reco");
  hEndy_2D_3trk->GetXaxis()->SetTitle("reco end.y");
  hEndy_2D_3trk->GetYaxis()->SetTitle("true end.y");
  hEndy_2D_3trk->Draw("colz");
  c56->Print("FSP/Endy_2D_3trk.png");

  TCanvas* c57 = new TCanvas();
  TH2* hEndz_2D_3trk = sEndz_2D_3trk.ToTH2(6.6e20);
  hEndz_2D_3trk->SetTitle("contained in truth, not in reco");
  hEndz_2D_3trk->GetXaxis()->SetTitle("reco end.z");
  hEndz_2D_3trk->GetYaxis()->SetTitle("true end.z");
  hEndz_2D_3trk->Draw("colz");
  c57->Print("FSP/Endz_2D_3trk.png");

  TCanvas* c58 = new TCanvas();
  TH2* hEndx_2Db_3trk = sEndx_2Db_3trk.ToTH2(6.6e20);
  hEndx_2Db_3trk->SetTitle("contained in reco, not in truth");
  hEndx_2Db_3trk->GetXaxis()->SetTitle("reco end.x");
  hEndx_2Db_3trk->GetYaxis()->SetTitle("true end.x");
  hEndx_2Db_3trk->Draw("colz");
  c58->Print("FSP/Endx_2Db_3trk.png");

  TCanvas* c59 = new TCanvas();
  TH2* hEndy_2Db_3trk = sEndy_2Db_3trk.ToTH2(6.6e20);
  hEndy_2Db_3trk->SetTitle("contained in reco, not in truth");
  hEndy_2Db_3trk->GetXaxis()->SetTitle("reco end.y");
  hEndy_2Db_3trk->GetYaxis()->SetTitle("true end.y");
  hEndy_2Db_3trk->Draw("colz");
  c59->Print("FSP/Endy_2Db_3trk.png");

  TCanvas* c60 = new TCanvas();
  TH2* hEndz_2Db_3trk = sEndz_2Db_3trk.ToTH2(6.6e20);
  hEndz_2Db_3trk->SetTitle("contained in reco, not in truth");
  hEndz_2Db_3trk->GetXaxis()->SetTitle("reco end.z");
  hEndz_2Db_3trk->GetYaxis()->SetTitle("true end.z");
  hEndz_2Db_3trk->Draw("colz");
  c60->Print("FSP/Endz_2Db_3trk.png");

  TCanvas* c61 = new TCanvas();
  TH2* hTrueStartRecoEndx_2D_3trk = sTrueStartRecoEndx_2D_3trk.ToTH2(6.6e20);
  hTrueStartRecoEndx_2D_3trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndx_2D_3trk->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2D_3trk->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2D_3trk->Draw("colz");
  c61->Print("FSP/TrueStartRecoEndx_2D_3trk.png");

  TCanvas* c62 = new TCanvas();
  TH2* hTrueStartRecoEndy_2D_3trk = sTrueStartRecoEndy_2D_3trk.ToTH2(6.6e20);
  hTrueStartRecoEndy_2D_3trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndy_2D_3trk->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2D_3trk->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2D_3trk->Draw("colz");
  c62->Print("FSP/TrueStartRecoEndy_2D_3trk.png");

  TCanvas* c63 = new TCanvas();
  TH2* hTrueStartRecoEndz_2D_3trk = sTrueStartRecoEndz_2D_3trk.ToTH2(6.6e20);
  hTrueStartRecoEndz_2D_3trk->SetTitle("contained in truth, not in reco");
  hTrueStartRecoEndz_2D_3trk->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2D_3trk->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2D_3trk->Draw("colz");
  c63->Print("FSP/TrueStartRecoEndz_2D_3trk.png");

  TCanvas* c64 = new TCanvas();
  TH2* hTrueStartRecoEndx_2Db_3trk = sTrueStartRecoEndx_2Db_3trk.ToTH2(6.6e20);
  hTrueStartRecoEndx_2Db_3trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndx_2Db_3trk->GetXaxis()->SetTitle("reco end.x");
  hTrueStartRecoEndx_2Db_3trk->GetYaxis()->SetTitle("true start.x");
  hTrueStartRecoEndx_2Db_3trk->Draw("colz");
  c64->Print("FSP/TrueStartRecoEndx_2Db_3trk.png");

  TCanvas* c65 = new TCanvas();
  TH2* hTrueStartRecoEndy_2Db_3trk = sTrueStartRecoEndy_2Db_3trk.ToTH2(6.6e20);
  hTrueStartRecoEndy_2Db_3trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndy_2Db_3trk->GetXaxis()->SetTitle("reco end.y");
  hTrueStartRecoEndy_2Db_3trk->GetYaxis()->SetTitle("true start.y");
  hTrueStartRecoEndy_2Db_3trk->Draw("colz");
  c65->Print("FSP/TrueStartRecoEndy_2Db_3trk.png");

  TCanvas* c66 = new TCanvas();
  TH2* hTrueStartRecoEndz_2Db_3trk = sTrueStartRecoEndz_2Db_3trk.ToTH2(6.6e20);
  hTrueStartRecoEndz_2Db_3trk->SetTitle("contained in reco, not in truth");
  hTrueStartRecoEndz_2Db_3trk->GetXaxis()->SetTitle("reco end.z");
  hTrueStartRecoEndz_2Db_3trk->GetYaxis()->SetTitle("true start.z");
  hTrueStartRecoEndz_2Db_3trk->Draw("colz");
  c66->Print("FSP/TrueStartRecoEndz_2Db_3trk.png");
}
