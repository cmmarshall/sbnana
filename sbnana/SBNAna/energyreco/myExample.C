#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"
#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"
#include "sbnana/CAFAna/Core/EnsembleRatio.h"
#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

using namespace ana;

#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include "canvas_margin.h"

#include "myMuonSelection.h"
#include "myProtonSelection.h"
#include "myEstimator.h"
#include "myTruth.h"
#include "myEventSelection.h"

double TargetPOT = 6.6e20;
TString str_TargetPOT = "6.6e20 POT";

void myExample(){


  //==== my personal ROOT plotting style. Check canvas_margin.h
  setTDRStyle();

  SpectrumLoader loader("/pnfs/icarus/persistent/users/icaruspro/SBNworkshopApril2021/CAF/corsika_nue_BNB/0/0/caf-01608068-6435-4a36-93b5-29ead574d963.root");

  //==== Baseline event selection defined in myEventSelection.h
  Cut cutMyCut = kIsNuMuCC && cutHasMuonTrack && cutIsMuonTrackLong && cutHasProtonTrack;

  //==== Set output directory
  TString outputDir = "./output/myAnalyzer/BindingE40MeV/";
  gSystem->mkdir(outputDir, kTRUE);

  //==== output ROOT file to save some histograms or canvases
  TString outputName = "output.root";
  TFile *outputfile = new TFile(outputDir+outputName,"RECREATE");

  //==== Binnings
  double xEnergyMin = 0.; // GeV
  double xEnergyMax = 10.;
  double dxEnergy = 0.1;
  const Binning binsEnergy = Binning::Simple( int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax );

  //==== HistAxis to be converted to Spectrum
  //==== Muon
  //====   Reco
  const HistAxis axMuonTrackRangeP("MuonTrackRangeP", binsEnergy, varMuonTrackRangeP);
  const HistAxis axMuonTrackMCSP("MuonTrackMCSP", binsEnergy, varMuonTrackMCSP);
  const HistAxis axMuonTrackCombinedP("MuonTrackCombinedP", binsEnergy, varMuonTrackCombinedP);
  const HistAxis axMuonTrackCaloP("MuonTrackCaloP", binsEnergy, varMuonTrackCaloP);
  //====   Truth
  const HistAxis axMuonTrackMatchedTruthP("MuonTrackMatchedTruthP", binsEnergy, varMuonTrackMatchedTruthP);
  
  //==== Spectrum
  //==== Muon
  //====   Reco
  Spectrum sMuonTrackRangeP(loader, axMuonTrackRangeP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackMCSP(loader, axMuonTrackMCSP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCombinedP(loader, axMuonTrackCombinedP, kNoSpillCut, cutMyCut);
  Spectrum sMuonTrackCaloP(loader, axMuonTrackCaloP, kNoSpillCut, cutMyCut);
  //====   Truth
  Spectrum sMuonTrackMatchedTruthP(loader, axMuonTrackMatchedTruthP, kNoSpillCut, cutMyCut);

  //==== This is the call that actually fills in the spectrum
  loader.Go();

  //==== Drawing

  //====   For auto-y
  double y_max = -999.;
  double y_min = 999999999;

  //====   POT on the top left corner

  TLatex latex_POT;
  latex_POT.SetNDC();
  latex_POT.SetTextSize(0.035);

  //====   Muon momentum (P)

  y_max = -999.;
  TCanvas *c_MuonP = new TCanvas("c_MuonP", "", 800, 800);
  canvas_margin(c_MuonP);
  c_MuonP->cd();

  TH1D *hist_dummy_MuonP = new TH1D("hist_dummy_MuonP", "", int( (xEnergyMax-xEnergyMin)/dxEnergy ), xEnergyMin, xEnergyMax);
  hist_axis(hist_dummy_MuonP);
  hist_dummy_MuonP->GetXaxis()->SetTitle("Muon momentum (GeV)");
  hist_dummy_MuonP->GetXaxis()->SetRangeUser(0.,5.0);
  hist_dummy_MuonP->Draw("hist");

  TH1 *hMuonTrackRangeP = sMuonTrackRangeP.ToTH1(TargetPOT);
  hMuonTrackRangeP->SetName("hMuonTrackRangeP");
  hMuonTrackRangeP->SetLineColor(kGreen);
  hMuonTrackRangeP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackRangeP->GetMaximum() );

  TH1 *hMuonTrackMCSP = sMuonTrackMCSP.ToTH1(TargetPOT);
  hMuonTrackMCSP->SetName("hMuonTrackMCSP");
  hMuonTrackMCSP->SetLineColor(kBlue);
  hMuonTrackMCSP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackMCSP->GetMaximum() );

  TH1 *hMuonTrackCombinedP = sMuonTrackCombinedP.ToTH1(TargetPOT);
  hMuonTrackCombinedP->SetName("hMuonTrackCombinedP");
  hMuonTrackCombinedP->SetLineColor(kBlack);
  hMuonTrackCombinedP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCombinedP->GetMaximum() );

  TH1 *hMuonTrackCaloP = sMuonTrackCaloP.ToTH1(TargetPOT);
  hMuonTrackCaloP->SetName("hMuonTrackCaloP");
  hMuonTrackCaloP->SetLineColor(kViolet);
  hMuonTrackCaloP->SetLineWidth(2);
  y_max = max( y_max, hMuonTrackCaloP->GetMaximum() );

  TH1 *hMuonTrackMatchedTruthP = sMuonTrackMatchedTruthP.ToTH1(TargetPOT);
  hMuonTrackMatchedTruthP->SetName("hMuonTrackMatchedTruthP");
  hMuonTrackMatchedTruthP->SetLineColor(kRed);
  hMuonTrackMatchedTruthP->SetLineWidth(2);
  hMuonTrackMatchedTruthP->SetLineStyle(3);
  y_max = max( y_max, hMuonTrackMatchedTruthP->GetMaximum() );

  hMuonTrackRangeP->Draw("histsame");
  hMuonTrackMCSP->Draw("histsame");
  hMuonTrackCombinedP->Draw("histsame");
  hMuonTrackCaloP->Draw("histsame");
  hMuonTrackMatchedTruthP->Draw("histsame");

  TLegend *lg_MuonP = new TLegend(0.35, 0.70, 0.90, 0.90);
  lg_MuonP->SetBorderSize(0);
  lg_MuonP->SetFillStyle(0);
  lg_MuonP->AddEntry(hMuonTrackRangeP, "Range", "l");
  lg_MuonP->AddEntry(hMuonTrackMCSP, "MCS", "l");
  lg_MuonP->AddEntry(hMuonTrackCombinedP, "Range (Contained) + MCS (Exiting)", "l");
  lg_MuonP->AddEntry(hMuonTrackCaloP, "Bestplane calo", "l");
  lg_MuonP->AddEntry(hMuonTrackMatchedTruthP, "Truth", "l");

  lg_MuonP->Draw();

  latex_POT.DrawLatex(0.20, 0.96, str_TargetPOT);

  hist_dummy_MuonP->GetYaxis()->SetRangeUser(0., y_max*1.1);
  c_MuonP->SaveAs(outputDir+"MuonP.pdf");

  outputfile->cd();
  c_MuonP->Write();

  c_MuonP->Close();

}

