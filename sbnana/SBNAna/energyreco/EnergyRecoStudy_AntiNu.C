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
  double E_AntiNu = E_Muon + Total_KE_Proton + Total_E_Pion + E_shws;
  return E_AntiNu;
});

const Var varEnergyRecoB([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = varEnergyReco(slc) + varErecNeutron(slc);
  return E_AntiNu;
});

const Var varEnergyReco1([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw==0) {
    E_AntiNu = varEnergyReco(slc);
  }
  return E_AntiNu;
});

const Var varEnergyReco1b([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw==0) {
    E_AntiNu = varEnergyRecoB(slc);
  }
  return E_AntiNu;
});

const Var varEnergyReco2([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    E_AntiNu = varEnergyReco(slc);
  }
  return E_AntiNu;
});

const Var varEnergyReco2b([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    E_AntiNu = varEnergyRecoB(slc);
  }
  return E_AntiNu;
});

const Var varEnergyReco3([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw>0) {
    E_AntiNu = varEnergyReco(slc);
  }
  else if (ntrk==2 and nshw>0) {
    E_AntiNu = varEnergyReco(slc);
  }
  else if (ntrk>=3) {
    E_AntiNu = varEnergyReco(slc);
  }
  return E_AntiNu;
});

const Var varEnergyReco3b([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw>0) {
    E_AntiNu = varEnergyRecoB(slc);
  }
  else if (ntrk==2 and nshw>0) {
    E_AntiNu = varEnergyRecoB(slc);
  }
  else if (ntrk>=3) {
    E_AntiNu = varEnergyRecoB(slc);
  }
  return E_AntiNu;
});

const Var varNeutrinoTruthE1([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw==0) {
    E_AntiNu = varNeutrinoTruthE(slc);
  }
  return E_AntiNu;
});

const Var varNeutrinoTruthE2([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    E_AntiNu = varNeutrinoTruthE(slc);
  }
  return E_AntiNu;
});

const Var varNeutrinoTruthE3([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw>0) {
    E_AntiNu = varNeutrinoTruthE(slc);
  }
  else if (ntrk==2 and nshw>0) {
    E_AntiNu = varNeutrinoTruthE(slc);
  }
  else if (ntrk>=3) {
    E_AntiNu = varNeutrinoTruthE(slc);
  }
  return E_AntiNu;
});

const Var varEnergyRes([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu_Reco = varEnergyReco(slc);
  double E_AntiNu_True = varNeutrinoTruthE(slc);
  double E_AntiNu_Res = E_AntiNu_Reco - E_AntiNu_True;
  return E_AntiNu_Res;
});

const Var varEnergyResB([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu_Reco = varEnergyRecoB(slc);
  double E_AntiNu_True = varNeutrinoTruthE(slc);
  double E_AntiNu_Res = E_AntiNu_Reco - E_AntiNu_True;
  return E_AntiNu_Res;
});

const Var varEnergyRes1([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw==0) {
    E_AntiNu = varEnergyRes(slc);
  }
  return E_AntiNu;
});

const Var varEnergyRes1b([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw==0) {
    E_AntiNu = varEnergyResB(slc);
  }
  return E_AntiNu;
});

const Var varEnergyRes2([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    E_AntiNu = varEnergyRes(slc);
  }
  return E_AntiNu;
});

const Var varEnergyRes2b([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    E_AntiNu = varEnergyResB(slc);
  }
  return E_AntiNu;
});

const Var varEnergyRes3([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw>0) {
    E_AntiNu = varEnergyRes(slc);
  }
  else if (ntrk==2 and nshw>0) {
    E_AntiNu = varEnergyRes(slc);
  }
  else if (ntrk>=3) {
    E_AntiNu = varEnergyRes(slc);
  }
  return E_AntiNu;
});

const Var varEnergyRes3b([](const caf::SRSliceProxy* slc) -> double {
  double E_AntiNu = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw>0) {
    E_AntiNu = varEnergyResB(slc);
  }
  else if (ntrk==2 and nshw>0) {
    E_AntiNu = varEnergyResB(slc);
  }
  else if (ntrk>=3) {
    E_AntiNu = varEnergyResB(slc);
  }
  return E_AntiNu;
});

const Var varPtrans1([](const caf::SRSliceProxy* slc) -> double {
  double Ptrans = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw==0) {
    Ptrans = varPtrans(slc);
  }
  return Ptrans;
});

const Var varPtrans2([](const caf::SRSliceProxy* slc) -> double {
  double Ptrans = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==2 and nshw==0) {
    Ptrans = varPtrans(slc);
  }
  return Ptrans;
});

const Var varPtrans3([](const caf::SRSliceProxy* slc) -> double {
  double Ptrans = -999;
  int ntrk = slc->reco.ntrk;
  int nshw = slc->reco.nshw;
  if (ntrk==1 and nshw>0) {
    Ptrans = varPtrans(slc);
  }
  else if (ntrk==2 and nshw>0) {
    Ptrans = varPtrans(slc);
  }
  else if (ntrk>=3) {
    Ptrans = varPtrans(slc);
  }
  return Ptrans;
});

// CUTS
const Cut mycutIsNuMuCC([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.iscc && slc->truth.pdg == 14 );
});

const Cut mycutIsAntiNuMuCC([](const caf::SRSliceProxy* slc) {
  return ( kIsNuSlice(slc) && slc->truth.iscc && slc->truth.pdg == -14 );
});

Cut cutNuMu = mycutIsNuMuCC && kRFiducial && cutHasMuonTrack && cutIsMuonTrackLong && cutExitingHadrons;
Cut cutAntiNuMu = mycutIsAntiNuMuCC && kRFiducial  && cutHasMuonTrack && cutIsMuonTrackLong && cutExitingHadrons;

// DATA FILES
// newest file
void EnergyRecoStudy_AntiNu(const std::string inputName = "IcarusProd2021B_BNB_Nu_Cosmics_v09_28_01_01_01_caf")

// rangeP-fixed CAFS
//void EnergyRecoStudy_AntiNu(const std::string inputName = "/pnfs/icarus/scratch/users/jskim/mc/SBNworkshopApril2021__corsika_numu_BNB/CAF/sbncode__v09_27_00_02_FixRangeProton/48388331_*/gen*.root")

// one file
//void EnergyRecoStudy_AntiNu(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/caf-01608068-6435-4a36-93b5-29ead574d963.root")

// twenty-one files
//void EnergyRecoStudy_AntiNu(const std::string inputName = "/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/*.root")

{
  SpectrumLoader loader(inputName);

// Binning
  const Binning bins2 = Binning::Simple(500,0,5);
  int numbins = 50;
  double Emax = 2;
  const Binning bins6 = Binning::Simple(200,0,Emax);
  const Binning bins4 = Binning::Simple(numbins,0,Emax);
  const Binning bins3 = Binning::Simple(500,-3,2);
  const Binning bins7 = Binning::Simple(200,-1,0.5);
  const Binning bins1 = Binning::Simple(numbins,0,0.75);

  Spectrum sEnergyReco("Neutrino Energy (GeV)", bins2, loader, varEnergyReco, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyReco1("Neutrino Energy (GeV)", bins2, loader, varEnergyReco1, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyReco2("Neutrino Energy (GeV)", bins2, loader, varEnergyReco2, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyReco3("Neutrino Energy (GeV)", bins2, loader, varEnergyReco3, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRecoB("Neutrino Energy (GeV)", bins2, loader, varEnergyRecoB, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyReco1b("Neutrino Energy (GeV)", bins2, loader, varEnergyReco1b, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyReco2b("Neutrino Energy (GeV)", bins2, loader, varEnergyReco2b, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyReco3b("Neutrino Energy (GeV)", bins2, loader, varEnergyReco3b, kNoSpillCut, cutAntiNuMu);

  Spectrum sEnergyTrue("Etrue, all events", bins2, loader, varNeutrinoTruthE, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyTrue1("Etrue, 1 trk and 0 shws", bins2, loader, varNeutrinoTruthE1, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyTrue2("Etrue, 2 trks and 0 shws", bins2, loader, varNeutrinoTruthE2, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyTrue3("Etrue, all other events", bins2, loader, varNeutrinoTruthE3, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEreco("Reconstructed Energy", bins6, varEnergyReco);
  const HistAxis axErecoB("Reconstructed Energy", bins6, varEnergyRecoB);
  const HistAxis axEtrue("True Energy (GeV)", bins4, varNeutrinoTruthE);
  Spectrum sRecoTrue2D(loader, axEtrue, axEreco, kNoSpillCut, cutAntiNuMu);
  Spectrum sRecoTrue2Db(loader, axEtrue, axErecoB, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEreco1("Reconstructed Energy", bins6, varEnergyReco1);
  const HistAxis axEreco1b("Reconstructed Energy", bins6, varEnergyReco1b);
  const HistAxis axEtrue1("True Energy (GeV)", bins4, varNeutrinoTruthE1);
  Spectrum sRecoTrue2D_1(loader, axEtrue1, axEreco1, kNoSpillCut, cutAntiNuMu);
  Spectrum sRecoTrue2D_1b(loader, axEtrue1, axEreco1b, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEreco2("Reconstructed Energy", bins6, varEnergyReco2);
  const HistAxis axEreco2b("Reconstructed Energy", bins6, varEnergyReco2b);
  const HistAxis axEtrue2("True Energy (GeV)", bins4, varNeutrinoTruthE2);
  Spectrum sRecoTrue2D_2(loader, axEtrue2, axEreco2, kNoSpillCut, cutAntiNuMu);
  Spectrum sRecoTrue2D_2b(loader, axEtrue2, axEreco2b, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEreco3("Reconstructed Energy", bins6, varEnergyReco3);
  const HistAxis axEreco3b("Reconstructed Energy", bins6, varEnergyReco3b);
  const HistAxis axEtrue3("True Energy (GeV)", bins4, varNeutrinoTruthE3);
  Spectrum sRecoTrue2D_3(loader, axEtrue3, axEreco3, kNoSpillCut, cutAntiNuMu);
  Spectrum sRecoTrue2D_3b(loader, axEtrue3, axEreco3b, kNoSpillCut, cutAntiNuMu);

  Spectrum sEnergyRes("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRes1("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes1, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRes2("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes2, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRes3("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes3, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyResB("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyResB, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRes1b("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes1b, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRes2b("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes2b, kNoSpillCut, cutAntiNuMu);
  Spectrum sEnergyRes3b("Residual Neutrino Energy (GeV)", bins3, loader, varEnergyRes3b, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEres("Residual Energy", bins7, varEnergyRes);
  const HistAxis axEresB("Residual Energy", bins7, varEnergyResB);
  Spectrum sResTrue2D(loader, axEtrue, axEres, kNoSpillCut, cutAntiNuMu);
  Spectrum sResTrue2Db(loader, axEtrue, axEresB, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEres1("Residual Energy", bins7, varEnergyRes1);
  const HistAxis axEres1b("Residual Energy", bins7, varEnergyRes1b);
  Spectrum sResTrue2D_1(loader, axEtrue1, axEres1, kNoSpillCut, cutAntiNuMu);
  Spectrum sResTrue2D_1b(loader, axEtrue1, axEres1b, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEres2("Residual Energy", bins7, varEnergyRes2);
  const HistAxis axEres2b("Residual Energy", bins7, varEnergyRes2b);
  Spectrum sResTrue2D_2(loader, axEtrue2, axEres2, kNoSpillCut, cutAntiNuMu);
  Spectrum sResTrue2D_2b(loader, axEtrue2, axEres2b, kNoSpillCut, cutAntiNuMu);

  const HistAxis axEres3("Residual Energy", bins7, varEnergyRes3);
  const HistAxis axEres3b("Residual Energy", bins7, varEnergyRes3b);
  Spectrum sResTrue2D_3(loader, axEtrue3, axEres3, kNoSpillCut, cutAntiNuMu);
  Spectrum sResTrue2D_3b(loader, axEtrue3, axEres3b, kNoSpillCut, cutAntiNuMu);

  const HistAxis axPtrans("Transverse Momentum", bins1, varPtrans);
  const HistAxis axPtrans1("Transverse Momentum", bins1, varPtrans1);
  const HistAxis axPtrans2("Transverse Momentum", bins1, varPtrans2);
  const HistAxis axPtrans3("Transverse Momentum", bins1, varPtrans3);
  Spectrum sPtransEres2D(loader, axPtrans, axEres, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2Db(loader, axPtrans, axEresB, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2D_1(loader, axPtrans1, axEres1, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2D_1b(loader, axPtrans1, axEres1b, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2D_2(loader, axPtrans2, axEres2, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2D_2b(loader, axPtrans2, axEres2b, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2D_3(loader, axPtrans3, axEres3, kNoSpillCut, cutAntiNuMu);
  Spectrum sPtransEres2D_3b(loader, axPtrans3, axEres3b, kNoSpillCut, cutAntiNuMu);

  loader.Go();

  // all events
  // 1D histogram of Erec 
  TCanvas* c1 = new TCanvas("c1", "c1");
  TH1* hEnergyReco = sEnergyReco.ToTH1(6.6e20);
  hEnergyReco->SetLineColor(kViolet);
  TH1* hEnergyRecoB = sEnergyRecoB.ToTH1(6.6e20);
  hEnergyRecoB->SetLineColor(kAzure);
  TH1* hEnergyTrue = sEnergyTrue.ToTH1(6.6e20);
  hEnergyTrue->SetLineColor(kTeal);
  hEnergyReco->Draw("hist");
  hEnergyRecoB->Draw("same hist");
  hEnergyTrue->Draw("same hist");

  TLegend* leg1 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg1->AddEntry(hEnergyReco, "Ereco, all events", "l");
  leg1->AddEntry(hEnergyRecoB, "Ereco + neutron adjustment", "l");
  leg1->AddEntry(hEnergyTrue, "Etrue", "l");
  leg1->Draw();
  c1->Print("EnergyRecoStudy_AntiNu/Ereco.png");

  // 2D histogram of Etrue vs Erec
  TCanvas* c2 = new TCanvas("c2", "c2");
  TH2* hRecoTrue2D = sRecoTrue2D.ToTH2(6.6e20);
  hRecoTrue2D->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D->GetYaxis()->SetTitle("Reconstructed Energy (GeV)");
  hRecoTrue2D->Draw("colz");
  c2->Print("EnergyRecoStudy_AntiNu/RecoTrue2D.png");

  TCanvas* c14 = new TCanvas("c14", "c14");
  TH2* hRecoTrue2Db = sRecoTrue2Db.ToTH2(6.6e20);
  hRecoTrue2Db->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2Db->GetYaxis()->SetTitle("Adjusted Reconstructed Energy (GeV)");
  hRecoTrue2Db->Draw("colz");
  c14->Print("EnergyRecoStudy_AntiNu/RecoTrue2Db.png");

  // 1D histogram of Eres 
  TCanvas* c3 = new TCanvas("c3", "c3");
  TH1* hEnergyRes = sEnergyRes.ToTH1(6.6e20);
  hEnergyRes->SetLineColor(kViolet);
  hEnergyRes->Draw("hist");
  TH1* hEnergyResB = sEnergyResB.ToTH1(6.6e20);
  hEnergyResB->SetLineColor(kAzure);
  hEnergyResB->Draw("same hist");

  TLegend* leg3 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg3->AddEntry(hEnergyRes, "Eres, all events", "l");
  leg3->AddEntry(hEnergyResB, "Eres + neutron adjustment", "l");
  leg3->Draw();
  c3->Print("EnergyRecoStudy_AntiNu/Eres.png");

  // 2D histogram of Etrue vs Eres
  TCanvas* c13 = new TCanvas("c13", "c13");
  TH2* hResTrue2D = sResTrue2D.ToTH2(6.6e20);
  hResTrue2D->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hResTrue2D->Draw("colz");
  c13->Print("EnergyRecoStudy_AntiNu/ResTrue2D.png");

  TCanvas* c15 = new TCanvas("c15", "c15");
  TH2* hResTrue2Db = sResTrue2Db.ToTH2(6.6e20);
  hResTrue2Db->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2Db->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hResTrue2Db->Draw("colz");
  c15->Print("EnergyRecoStudy_AntiNu/ResTrue2Db.png");

  // scatterplot of Etrue vs mean Eres
  TCanvas* c33 = new TCanvas("c33", "c33");
  TGraphErrors* gr1 = new TGraphErrors();
  TH2* hResTrue2D_copy = sResTrue2D.ToTH2(6.6e20);
  double x, y, ey;
  int index=0;
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_copy->GetMean(2);
    ey = hResTrue2D_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr1->SetPoint(index, x, y);
      gr1->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr1->SetLineColor(kViolet);
  gr1->SetLineWidth(1);
  gr1->SetMarkerColor(kViolet);
  gr1->SetMarkerSize(1);
  gr1->SetMarkerStyle(20);
  gr1->GetXaxis()->SetTitle("True Energy (GeV)");
  gr1->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  gr1->Draw("AP");
  c33->Print("EnergyRecoStudy_AntiNu/AvgEres.png");

  TCanvas* c34 = new TCanvas("c34", "c34");
  TGraphErrors* gr2 = new TGraphErrors();
  TH2* hResTrue2Db_copy = sResTrue2Db.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2Db_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2Db_copy->GetMean(2);
    ey = hResTrue2Db_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr2->SetPoint(index, x, y);
      gr2->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr2->SetLineColor(kAzure);
  gr2->SetLineWidth(1);
  gr2->SetMarkerColor(kAzure);
  gr2->SetMarkerSize(1);
  gr2->SetMarkerStyle(20);
  gr2->GetXaxis()->SetTitle("True Energy (GeV)");
  gr2->GetYaxis()->SetTitle("Avg Adjusted Energy Residual (GeV)");
  gr2->Draw("AP");
  c34->Print("EnergyRecoStudy_AntiNu/AvgEresb.png");

  TCanvas* c41 = new TCanvas("c41", "c41");
  TMultiGraph* mg1 = new TMultiGraph();
  mg1->Add(gr1);
  mg1->Add(gr2);
  mg1->GetXaxis()->SetTitle("True Energy (GeV)");
  mg1->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  mg1->Draw("AP");
  TLegend* leg41 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg41->AddEntry(gr1, "No Adjustment", "l");
  leg41->AddEntry(gr2, "Transverse Momentum Adjustment", "l");
  leg41->Draw();
  c41->Print("EnergyRecoStudy_AntiNu/Scatter1.png");

  // 2D histogram of Ptrans vs Eres
  TCanvas* c25 = new TCanvas("c25", "c25");
  TH2* hPtransEres2D = sPtransEres2D.ToTH2(6.6e20);
  hPtransEres2D->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hPtransEres2D->Draw("colz");
  c25->Print("EnergyRecoStudy_AntiNu/PtransEres2D.png");

  TCanvas* c26 = new TCanvas("c26", "c26");
  TH2* hPtransEres2Db = sPtransEres2Db.ToTH2(6.6e20);
  hPtransEres2Db->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2Db->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hPtransEres2Db->Draw("colz");
  c26->Print("EnergyRecoStudy_AntiNu/PtransEres2Db.png");

  // 1 trk and 0 shws
  TCanvas* c4 = new TCanvas("c4", "c4");
  TH1* hEnergyReco1 = sEnergyReco1.ToTH1(6.6e20);
  hEnergyReco1->SetLineColor(kViolet);
  TH1* hEnergyReco1b = sEnergyReco1b.ToTH1(6.6e20);
  hEnergyReco1b->SetLineColor(kAzure);
  TH1* hEnergyTrue1 = sEnergyTrue1.ToTH1(6.6e20);
  hEnergyTrue1->SetLineColor(kTeal);
  hEnergyReco1->Draw("hist");
  hEnergyReco1b->Draw("same hist");
  hEnergyTrue1->Draw("same hist");

  TLegend* leg4 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg4->AddEntry(hEnergyReco1, "Ereco, 1 trk and 0 shws", "l");
  leg4->AddEntry(hEnergyReco1b, "Ereco + neutron adjustment", "l");
  leg4->AddEntry(hEnergyTrue1, "Etrue", "l");
  leg4->Draw();
  c4->Print("EnergyRecoStudy_AntiNu/Ereco1.png");

  TCanvas* c5 = new TCanvas("c5", "c5");
  TH2* hRecoTrue2D_1 = sRecoTrue2D_1.ToTH2(6.6e20);
  hRecoTrue2D_1->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D_1->GetYaxis()->SetTitle("Reconstructed Energy (GeV)");
  hRecoTrue2D_1->Draw("colz");
  c5->Print("EnergyRecoStudy_AntiNu/RecoTrue2D_1.png");

  TCanvas* c17 = new TCanvas("c17", "c17");
  TH2* hRecoTrue2D_1b = sRecoTrue2D_1b.ToTH2(6.6e20);
  hRecoTrue2D_1b->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D_1b->GetYaxis()->SetTitle("Adjusted Reconstructed Energy (GeV)");
  hRecoTrue2D_1b->Draw("colz");
  c17->Print("EnergyRecoStudy_AntiNu/RecoTrue2D_1b.png");

  TCanvas* c6 = new TCanvas("c6", "c6");
  TH1* hEnergyRes1 = sEnergyRes1.ToTH1(6.6e20);
  hEnergyRes1->SetLineColor(kViolet);
  hEnergyRes1->Draw("hist");
  TH1* hEnergyRes1b = sEnergyRes1b.ToTH1(6.6e20);
  hEnergyRes1b->SetLineColor(kAzure);
  hEnergyRes1b->Draw("same hist");

  TLegend* leg6 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg6->AddEntry(hEnergyRes1, "Eres, 1 trk and 0 shws", "l");
  leg6->AddEntry(hEnergyRes1b, "Eres + neutron adjustment", "l");
  leg6->Draw();
  c6->Print("EnergyRecoStudy_AntiNu/Eres1.png");

  TCanvas* c16 = new TCanvas("c16", "c16");
  TH2* hResTrue2D_1 = sResTrue2D_1.ToTH2(6.6e20);
  hResTrue2D_1->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D_1->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hResTrue2D_1->Draw("colz");
  c16->Print("EnergyRecoStudy_AntiNu/ResTrue2D_1.png");

  TCanvas* c18 = new TCanvas("c18", "c18");
  TH2* hResTrue2D_1b = sResTrue2D_1b.ToTH2(6.6e20);
  hResTrue2D_1b->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D_1b->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hResTrue2D_1b->Draw("colz");
  c18->Print("EnergyRecoStudy_AntiNu/ResTrue2D_1b.png");

  TCanvas* c35 = new TCanvas("c35", "c35");
  TGraphErrors* gr3 = new TGraphErrors();
  TH2* hResTrue2D_1_copy = sResTrue2D_1.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_1_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_1_copy->GetMean(2);
    ey = hResTrue2D_1_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr3->SetPoint(index, x, y);
      gr3->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr3->SetLineColor(kViolet);
  gr3->SetLineWidth(1);
  gr3->SetMarkerColor(kViolet);
  gr3->SetMarkerSize(1);
  gr3->SetMarkerStyle(20);
  gr3->GetXaxis()->SetTitle("True Energy (GeV)");
  gr3->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  gr3->Draw("AP");
  c35->Print("EnergyRecoStudy_AntiNu/AvgEres1.png");

  TCanvas* c36 = new TCanvas("c36", "c36");
  TGraphErrors* gr4 = new TGraphErrors();
  TH2* hResTrue2D_1b_copy = sResTrue2D_1b.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_1b_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_1b_copy->GetMean(2);
    ey = hResTrue2D_1b_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr4->SetPoint(index, x, y);
      gr4->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr4->SetLineColor(kAzure);
  gr4->SetLineWidth(1);
  gr4->SetMarkerColor(kAzure);
  gr4->SetMarkerSize(1);
  gr4->SetMarkerStyle(20);
  gr4->GetXaxis()->SetTitle("True Energy (GeV)");
  gr4->GetYaxis()->SetTitle("Avg Adjusted Energy Residual (GeV)");
  gr4->Draw("AP");
  c36->Print("EnergyRecoStudy_AntiNu/AvgEres1b.png");

  TCanvas* c42 = new TCanvas("c42", "c42");
  TMultiGraph* mg2 = new TMultiGraph();
  mg2->Add(gr3);
  mg2->Add(gr4);
  mg2->GetXaxis()->SetTitle("True Energy (GeV)");
  mg2->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  mg2->Draw("AP");
  TLegend* leg42 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg42->AddEntry(gr3, "No Adjustment", "l");
  leg42->AddEntry(gr4, "Transverse Momentum Adjustment", "l");
  leg42->Draw();
  c42->Print("EnergyRecoStudy_AntiNu/Scatter2.png");

  TCanvas* c27 = new TCanvas("c27", "c27");
  TH2* hPtransEres2D_1 = sPtransEres2D_1.ToTH2(6.6e20);
  hPtransEres2D_1->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D_1->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hPtransEres2D_1->Draw("colz");
  c27->Print("EnergyRecoStudy_AntiNu/PtransEres2D_1.png");

  TCanvas* c28 = new TCanvas("c28", "c28");
  TH2* hPtransEres2D_1b = sPtransEres2D_1b.ToTH2(6.6e20);
  hPtransEres2D_1b->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D_1b->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hPtransEres2D_1b->Draw("colz");
  c28->Print("EnergyRecoStudy_AntiNu/PtransEres2D_1b.png");

  // 2 trks and 0 shws
  TCanvas* c7 = new TCanvas("c7", "c7");
  TH1* hEnergyReco2 = sEnergyReco2.ToTH1(6.6e20);
  hEnergyReco2->SetLineColor(kViolet);
  TH1* hEnergyReco2b = sEnergyReco2b.ToTH1(6.6e20);
  hEnergyReco2b->SetLineColor(kAzure);
  TH1* hEnergyTrue2 = sEnergyTrue2.ToTH1(6.6e20);
  hEnergyTrue2->SetLineColor(kTeal);
  hEnergyReco2->Draw("hist");
  hEnergyReco2b->Draw("same hist");
  hEnergyTrue2->Draw("same hist");

  TLegend* leg7 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg7->AddEntry(hEnergyReco2, "Ereco, 2 trks and 0 shws", "l");
  leg7->AddEntry(hEnergyReco2b, "Ereco + neutron adjustment", "l");
  leg7->AddEntry(hEnergyTrue2, "Etrue", "l");
  leg7->Draw();
  c7->Print("EnergyRecoStudy_AntiNu/Ereco2.png");

  TCanvas* c8 = new TCanvas("c8", "c8");
  TH2* hRecoTrue2D_2 = sRecoTrue2D_2.ToTH2(6.6e20);
  hRecoTrue2D_2->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D_2->GetYaxis()->SetTitle("Reconstructed Energy (GeV)");
  hRecoTrue2D_2->Draw("colz");
  c8->Print("EnergyRecoStudy_AntiNu/RecoTrue2D_2.png");

  TCanvas* c20 = new TCanvas("c20", "c20");
  TH2* hRecoTrue2D_2b = sRecoTrue2D_2b.ToTH2(6.6e20);
  hRecoTrue2D_2b->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D_2b->GetYaxis()->SetTitle("Adjusted Reconstructed Energy (GeV)");
  hRecoTrue2D_2b->Draw("colz");
  c20->Print("EnergyRecoStudy_AntiNu/RecoTrue2D_2b.png");

  TCanvas* c9 = new TCanvas("c9", "c9");
  TH1* hEnergyRes2 = sEnergyRes2.ToTH1(6.6e20);
  hEnergyRes2->SetLineColor(kViolet);
  hEnergyRes2->Draw("hist");
  TH1* hEnergyRes2b = sEnergyRes2b.ToTH1(6.6e20);
  hEnergyRes2b->SetLineColor(kAzure);
  hEnergyRes2b->Draw("same hist");

  TLegend* leg9 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg9->AddEntry(hEnergyRes2, "Eres, 2 trks and 0 shws", "l");
  leg9->AddEntry(hEnergyRes2b, "Eres + neutron adjustment", "l");
  leg9->Draw();
  c9->Print("EnergyRecoStudy_AntiNu/Eres2.png");

  TCanvas* c19 = new TCanvas("c19", "c19");
  TH2* hResTrue2D_2 = sResTrue2D_2.ToTH2(6.6e20);
  hResTrue2D_2->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D_2->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hResTrue2D_2->Draw("colz");
  c19->Print("EnergyRecoStudy_AntiNu/ResTrue2D_2.png");

  TCanvas* c21 = new TCanvas("c21", "c21");
  TH2* hResTrue2D_2b = sResTrue2D_2b.ToTH2(6.6e20);
  hResTrue2D_2b->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D_2b->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hResTrue2D_2b->Draw("colz");
  c21->Print("EnergyRecoStudy_AntiNu/ResTrue2D_2b.png");

  TCanvas* c37 = new TCanvas("c37", "c37");
  TGraphErrors* gr5 = new TGraphErrors();
  TH2* hResTrue2D_2_copy = sResTrue2D_2.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_2_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_2_copy->GetMean(2);
    ey = hResTrue2D_2_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr5->SetPoint(index, x, y);
      gr5->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr5->SetLineColor(kViolet);
  gr5->SetLineWidth(1);
  gr5->SetMarkerColor(kViolet);
  gr5->SetMarkerSize(1);
  gr5->SetMarkerStyle(20);
  gr5->GetXaxis()->SetTitle("True Energy (GeV)");
  gr5->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  gr5->Draw("AP");
  c37->Print("EnergyRecoStudy_AntiNu/AvgEres2.png");

  TCanvas* c38 = new TCanvas("c38", "c38");
  TGraphErrors* gr6 = new TGraphErrors();
  TH2* hResTrue2D_2b_copy = sResTrue2D_2b.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_2b_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_2b_copy->GetMean(2);
    ey = hResTrue2D_2b_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr6->SetPoint(index, x, y);
      gr6->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr6->SetLineColor(kAzure);
  gr6->SetLineWidth(1);
  gr6->SetMarkerColor(kAzure);
  gr6->SetMarkerSize(1);
  gr6->SetMarkerStyle(20);
  gr6->GetXaxis()->SetTitle("True Energy (GeV)");
  gr6->GetYaxis()->SetTitle("Avg Adjusted Energy Residual (GeV)");
  gr6->Draw("AP");
  c38->Print("EnergyRecoStudy_AntiNu/AvgEres2b.png");

  TCanvas* c43 = new TCanvas("c43", "c43");
  TMultiGraph* mg3 = new TMultiGraph();
  mg3->Add(gr5);
  mg3->Add(gr6);
  mg3->GetXaxis()->SetTitle("True Energy (GeV)");
  mg3->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  mg3->Draw("AP");
  TLegend* leg43 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg43->AddEntry(gr5, "No Adjustment", "l");
  leg43->AddEntry(gr6, "Transverse Momentum Adjustment", "l");
  leg43->Draw();
  c43->Print("EnergyRecoStudy_AntiNu/Scatter3.png");

  TCanvas* c29 = new TCanvas("c29", "c29");
  TH2* hPtransEres2D_2 = sPtransEres2D_2.ToTH2(6.6e20);
  hPtransEres2D_2->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D_2->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hPtransEres2D_2->Draw("colz");
  c29->Print("EnergyRecoStudy_AntiNu/PtransEres2D_2.png");

  TCanvas* c30 = new TCanvas("c30", "c30");
  TH2* hPtransEres2D_2b = sPtransEres2D_2b.ToTH2(6.6e20);
  hPtransEres2D_2b->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D_2b->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hPtransEres2D_2b->Draw("colz");
  c30->Print("EnergyRecoStudy_AntiNu/PtransEres2D_2b.png");

  // all other events
  TCanvas* c10 = new TCanvas("c10", "c10");
  TH1* hEnergyReco3 = sEnergyReco3.ToTH1(6.6e20);
  hEnergyReco3->SetLineColor(kViolet);
  TH1* hEnergyReco3b = sEnergyReco3b.ToTH1(6.6e20);
  hEnergyReco3b->SetLineColor(kAzure);
  TH1* hEnergyTrue3 = sEnergyTrue3.ToTH1(6.6e20);
  hEnergyTrue3->SetLineColor(kTeal);
  hEnergyReco3->Draw("hist");
  hEnergyReco3b->Draw("same hist");
  hEnergyTrue3->Draw("same hist");

  TLegend* leg10 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg10->AddEntry(hEnergyReco3, "Ereco, all other events", "l");
  leg10->AddEntry(hEnergyReco3b, "Ereco + neutron adjustment", "l");
  leg10->AddEntry(hEnergyTrue3, "Etrue", "l");
  leg10->Draw();
  c10->Print("EnergyRecoStudy_AntiNu/Ereco3.png");

  TCanvas* c11 = new TCanvas("c11", "c11");
  TH2* hRecoTrue2D_3 = sRecoTrue2D_3.ToTH2(6.6e20);
  hRecoTrue2D_3->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D_3->GetYaxis()->SetTitle("Reconstructed Energy (GeV)");
  hRecoTrue2D_3->Draw("colz");
  c11->Print("EnergyRecoStudy_AntiNu/RecoTrue2D_3.png");

  TCanvas* c23 = new TCanvas("c23", "c23");
  TH2* hRecoTrue2D_3b = sRecoTrue2D_3b.ToTH2(6.6e20);
  hRecoTrue2D_3b->GetXaxis()->SetTitle("True Energy (GeV)");
  hRecoTrue2D_3b->GetYaxis()->SetTitle("Adjusted Reconstructed Energy (GeV)");
  hRecoTrue2D_3b->Draw("colz");
  c23->Print("EnergyRecoStudy_AntiNu/RecoTrue2D_3b.png");

  TCanvas* c12 = new TCanvas("c12", "c12");
  TH1* hEnergyRes3 = sEnergyRes3.ToTH1(6.6e20);
  hEnergyRes3->SetLineColor(kViolet);
  hEnergyRes3->Draw("hist");
  TH1* hEnergyRes3b = sEnergyRes3b.ToTH1(6.6e20);
  hEnergyRes3b->SetLineColor(kAzure);
  hEnergyRes3b->Draw("same hist");

  TLegend* leg12 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg12->AddEntry(hEnergyRes3, "Eres, all other events", "l");
  leg12->AddEntry(hEnergyRes3b, "Eres + neutron adjustment", "l");
  leg12->Draw();
  c12->Print("EnergyRecoStudy_AntiNu/Eres3.png");

  TCanvas* c22 = new TCanvas("c22", "c22");
  TH2* hResTrue2D_3 = sResTrue2D_3.ToTH2(6.6e20);
  hResTrue2D_3->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D_3->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hResTrue2D_3->Draw("colz");
  c22->Print("EnergyRecoStudy_AntiNu/ResTrue2D_3.png");

  TCanvas* c24 = new TCanvas("c24", "c24");
  TH2* hResTrue2D_3b = sResTrue2D_3b.ToTH2(6.6e20);
  hResTrue2D_3b->GetXaxis()->SetTitle("True Energy (GeV)");
  hResTrue2D_3b->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hResTrue2D_3b->Draw("colz");
  c24->Print("EnergyRecoStudy_AntiNu/ResTrue2D_3b.png");

  TCanvas* c39 = new TCanvas("c39", "c39");
  TGraphErrors* gr7 = new TGraphErrors();
  TH2* hResTrue2D_3_copy = sResTrue2D_3.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_3_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_3_copy->GetMean(2);
    ey = hResTrue2D_3_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr7->SetPoint(index, x, y);
      gr7->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr7->SetLineColor(kViolet);
  gr7->SetLineWidth(1);
  gr7->SetMarkerColor(kViolet);
  gr7->SetMarkerSize(1);
  gr7->SetMarkerStyle(20);
  gr7->GetXaxis()->SetTitle("True Energy (GeV)");
  gr7->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  gr7->Draw("AP");
  c39->Print("EnergyRecoStudy_AntiNu/AvgEres3.png");

  TCanvas* c40 = new TCanvas("c40", "c40");
  TGraphErrors* gr8 = new TGraphErrors();
  TH2* hResTrue2D_3b_copy = sResTrue2D_3b.ToTH2(6.6e20);
  for (int i=0; i<numbins; i++) {
    int binmin = i;
    int binmax = i+1;
    x = binmin*(Emax/numbins);
    hResTrue2D_3b_copy->GetXaxis()->SetRange(binmin,binmax);
    y = hResTrue2D_3b_copy->GetMean(2);
    ey = hResTrue2D_3b_copy->GetMean(12);
    if (!(y==0 and ey==0)) {
      gr8->SetPoint(index, x, y);
      gr8->SetPointError(index, 0, ey);
      index+=1;
    }}
  index=0;
  gr8->SetLineColor(kAzure);
  gr8->SetLineWidth(1);
  gr8->SetMarkerColor(kAzure);
  gr8->SetMarkerSize(1);
  gr8->SetMarkerStyle(20);
  gr8->GetXaxis()->SetTitle("True Energy (GeV)");
  gr8->GetYaxis()->SetTitle("Avg Adjusted Energy Residual (GeV)");
  gr8->Draw("AP");
  c40->Print("EnergyRecoStudy_AntiNu/AvgEres3b.png");

  TCanvas* c44 = new TCanvas("c44", "c44");
  TMultiGraph* mg4 = new TMultiGraph();
  mg4->Add(gr7);
  mg4->Add(gr8);
  mg4->GetXaxis()->SetTitle("True Energy (GeV)");
  mg4->GetYaxis()->SetTitle("Average Energy Residual (GeV)");
  mg4->Draw("AP");
  TLegend* leg44 = new TLegend(0.6, 0.7, 0.85, 0.8);
  leg44->AddEntry(gr7, "No Adjustment", "l");
  leg44->AddEntry(gr8, "Transverse Momentum Adjustment", "l");
  leg44->Draw();
  c44->Print("EnergyRecoStudy_AntiNu/Scatter4.png");

  TCanvas* c31 = new TCanvas("c31", "c31");
  TH2* hPtransEres2D_3 = sPtransEres2D_3.ToTH2(6.6e20);
  hPtransEres2D_3->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D_3->GetYaxis()->SetTitle("Residual Energy (GeV)");
  hPtransEres2D_3->Draw("colz");
  c31->Print("EnergyRecoStudy_AntiNu/PtransEres2D_3.png");

  TCanvas* c32 = new TCanvas("c32", "c32");
  TH2* hPtransEres2D_3b = sPtransEres2D_3b.ToTH2(6.6e20);
  hPtransEres2D_3b->GetXaxis()->SetTitle("Transverse Momentum (GeV)");
  hPtransEres2D_3b->GetYaxis()->SetTitle("Adjusted Residual Energy (GeV)");
  hPtransEres2D_3b->Draw("colz");
  c32->Print("EnergyRecoStudy_AntiNu/PtransEres2D_3b.png");
}
