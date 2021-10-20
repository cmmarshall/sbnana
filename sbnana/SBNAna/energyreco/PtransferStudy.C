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

// VARS
const Var varEreco([](const caf::SRSliceProxy* slc) -> double {
  // Muon energy
  double P_Muon = varMuonTrackCombinedP(slc);
  double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );
  // Track kinetic energy
  std::vector<double> KE_trks = varTrkKE(slc);
  double Total_KE_trks = std::accumulate(KE_trks.begin(), KE_trks.end(), 0.0);
  // Shower energy
  double E_shws = varShwE(slc);
  // Neutrino energy
  double E_nu = E_Muon + Total_KE_trks + E_shws;
  return E_nu;
});

const Var varEmuon([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon = varMuonTrackCombinedP(slc);
  double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );
  return E_Muon;
});

ofstream outputfile;
int counter=0;

const Var varEtransfer([](const caf::SRSliceProxy* slc) -> double {
  double E_nu = varEreco(slc);
  double E_mu = varEmuon(slc);
  double delta_E = E_nu - E_mu;
  if (delta_E==0) {
    counter+=1;
    int ndaughters = slc->truth.nprim;
    for (int i=0; i<ndaughters; i++) {
      int daughter = slc->truth.prim[i].pdg;
      double energy = slc->truth.prim[i].genE;
      outputfile.open("FSP_NoEtransfer.txt", ios::app);
      outputfile<<daughter<<" (E="<<energy<<"), ";
      outputfile.close();
    }
    outputfile.open("FSP_NoEtransfer.txt", ios::app);
    outputfile<<"\n";
    outputfile.close();
  }
  return delta_E;
});

ofstream outputfile2;
int counter2=0;

const Var varTrueEtransfer([](const caf::SRSliceProxy* slc) -> double {
  double E_nu = varNeutrinoTruthE(slc); 
  double E_mu = varMuonTrackMatchedTruthE(slc); 
  double delta_E = E_nu - E_mu;
  if (delta_E==0) {
    counter2+=1;
    int ndaughters = slc->truth.nprim;
    for (int i=0; i<ndaughters; i++) {
      int daughter = slc->truth.prim[i].pdg;
      double energy = slc->truth.prim[i].genE;
      outputfile2.open("FSP_NoTrueEtransfer.txt", ios::app);
      outputfile2<<daughter<<" (E="<<energy<<"), ";
      outputfile2.close();
    }
    outputfile2.open("FSP_NoTrueEtransfer.txt", ios::app);
    outputfile2<<"\n";
    outputfile2.close();
  }
  return delta_E;
});

const Var varPtransfer([](const caf::SRSliceProxy* slc) -> double {
  std::vector<double> P_Muon = varMuonP(slc);
  double px_mu = P_Muon[0];
  double py_mu = P_Muon[1];
  double pz_mu = P_Muon[2];
  double px_nu = 0;
  double py_nu = 0;
  double pz_nu = varEreco(slc); 
  double delta_px = px_mu - px_nu;
  double delta_py = py_mu - py_nu;
  double delta_pz = pz_mu - pz_nu;
  double delta_p = sqrt( delta_px*delta_px + delta_py*delta_py + delta_pz*delta_pz );
  return delta_p;
});

const Var varTruePtransfer([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  double px_mu = slc->reco.trk[muonTrackIndex].truth.p.genp.x;
  double py_mu = slc->reco.trk[muonTrackIndex].truth.p.genp.y;
  double pz_mu = slc->reco.trk[muonTrackIndex].truth.p.genp.z;
  double px_nu = slc->truth.momentum.x;
  double py_nu = slc->truth.momentum.y;
  double pz_nu = slc->truth.momentum.z; 
  double delta_px = px_mu - px_nu;
  double delta_py = py_mu - py_nu;
  double delta_pz = pz_mu - pz_nu;
  double delta_p = sqrt( delta_px*delta_px + delta_py*delta_py + delta_pz*delta_pz );
  return delta_p;
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

const Cut cutMultipleTrks([](const caf::SRSliceProxy* slc) {
  return varEtransfer(slc)>0;
});

const Cut cutOneTrk([](const caf::SRSliceProxy* slc) {
  return varEtransfer(slc)==0;
});

const Cut cutTrueMultipleTrks([](const caf::SRSliceProxy* slc) {
  return varTrueEtransfer(slc)>0;
});

const Cut cutTrueOneTrk([](const caf::SRSliceProxy* slc) {
  return varTrueEtransfer(slc)==0;
});

// DATA FILES
void PtransferStudy(const std::string inputName = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf")

{
  SpectrumLoader loader(inputName);

// Binning
  const Binning bins1 = Binning::Simple(100,0,1);
  const Binning bins2 = Binning::Simple(100,0,1);

  const HistAxis axEtransfer("Transfer Energy", bins1, varEtransfer);
  const HistAxis axPtransfer("Transfer Momentum", bins2, varPtransfer);
  Spectrum sPEtransfer(loader, axPtransfer, axEtransfer, kNoSpillCut, cutNuMu);
  Spectrum sPtransfer("Transfer Momentum", bins2, loader, varPtransfer, kNoSpillCut, cutNuMu && cutOneTrk);
  Spectrum sPEtransfer2(loader, axPtransfer, axEtransfer, kNoSpillCut, cutNuMu && cutMultipleTrks);

  const HistAxis axTrueEtransfer("Transfer Energy", bins1, varTrueEtransfer);
  const HistAxis axTruePtransfer("Transfer Momentum", bins2, varTruePtransfer);
  Spectrum sTruePEtransfer(loader, axTruePtransfer, axTrueEtransfer, kNoSpillCut, cutNuMu);
  Spectrum sTruePtransfer("Transfer Momentum", bins2, loader, varTruePtransfer, kNoSpillCut, cutNuMu && cutTrueOneTrk);
  Spectrum sTruePEtransfer2(loader, axTruePtransfer, axTrueEtransfer, kNoSpillCut, cutNuMu && cutTrueMultipleTrks);

  outputfile.open("FSP_NoEtransfer.txt", ios::trunc);
  outputfile.close();
  outputfile2.open("FSP_NoTrueEtransfer.txt", ios::trunc);
  outputfile2.close();

  loader.Go();

  outputfile.close();
  std::cout<<"counter = "<<counter<<"\n";
  outputfile2.close();
  std::cout<<"counter2 = "<<counter2<<"\n";

  // 2D histogram of Etransfer vs Ptransfer
  TCanvas* c1 = new TCanvas();
  c1->SetLogz();
  TH2* hPEtransfer = sPEtransfer.ToTH2(6.6e20);
  hPEtransfer->GetXaxis()->SetTitle("Momentum Transfer (GeV)");
  hPEtransfer->GetYaxis()->SetTitle("Energy Transfer (GeV)");
  hPEtransfer->Draw("colz");
  c1->Print("PtransferStudy/PEtransfer.png");

  // 1D histogram of Ptransfer in events with Etransfer=0 (1 trk, 0 shws events)
  TCanvas* c2 = new TCanvas();
  TH1* hPtransfer = sPtransfer.ToTH1(6.6e20);
  hPtransfer->GetXaxis()->SetTitle("Momentum Transfer (GeV)");
  hPtransfer->SetLineColor(kCyan);
  hPtransfer->Draw("hist");
  c2->Print("PtransferStudy/Ptransfer.png");

  // 2D histogram of Etransfer vs Ptransfer in events with Etransfer>0 (events with multiple trks and shws)
  TCanvas* c3 = new TCanvas();
  c3->SetLogz();
  TH2* hPEtransfer2 = sPEtransfer2.ToTH2(6.6e20);
  hPEtransfer2->GetXaxis()->SetTitle("Momentum Transfer (GeV)");
  hPEtransfer2->GetYaxis()->SetTitle("Energy Transfer (GeV)");
  hPEtransfer2->Draw("colz");
  c3->Print("PtransferStudy/PEtransfer2.png");

  // 2D hisogram of true Etransfer vs true Ptransfer
  TCanvas* c4 = new TCanvas();
  c4->SetLogz();
  TH2* hTruePEtransfer = sTruePEtransfer.ToTH2(6.6e20);
  hTruePEtransfer->GetXaxis()->SetTitle("Momentum Transfer (GeV)");
  hTruePEtransfer->GetYaxis()->SetTitle("Energy Transfer (GeV)");
  hTruePEtransfer->Draw("colz");
  c4->Print("PtransferStudy/TruePEtransfer.png");

  // 1D histogram of true Ptransfer in events with true Etransfer=0 (1 trk, 0 shws events)
  TCanvas* c5 = new TCanvas();
  TH1* hTruePtransfer = sTruePtransfer.ToTH1(6.6e20);
  hTruePtransfer->GetXaxis()->SetTitle("Momentum Transfer (GeV)");
  hTruePtransfer->SetLineColor(kCyan);
  hTruePtransfer->Draw("hist");
  c5->Print("PtransferStudy/TruePtransfer.png");

  // 2D histogram of true Etransfer vs true Ptransfer in events with true Etransfer>0 (events with multiple trks and shws)
  TCanvas* c6 = new TCanvas();
  c6->SetLogz();
  TH2* hTruePEtransfer2 = sTruePEtransfer2.ToTH2(6.6e20);
  hTruePEtransfer2->GetXaxis()->SetTitle("Momentum Transfer (GeV)");
  hTruePEtransfer2->GetYaxis()->SetTitle("Energy Transfer (GeV)");
  hTruePEtransfer2->Draw("colz");
  c6->Print("PtransferStudy/TruePEtransfer2.png");
}
