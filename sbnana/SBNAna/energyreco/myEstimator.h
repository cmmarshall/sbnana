#ifndef myEstimator_h
#define myEstimator_h

#include "myMuonSelection.h"
#include "myProtonSelection.h"

#define bindingE 0.040

//==== Muon estimator

//====   Reco

const Var varMuonTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    return trk_Muon.rangeP.p_muon;
  }
  else{
    return -999;
  }
});

const Var varMuonTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    return trk_Muon.mcsP.fwdP_muon;
  }   
  else{
    return -999;
  }   
}); 

const Var varMuonTrackCombinedP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    int mtc = varMuonTrackContainedness(slc); // 0 : exiting, 1 : contained
    if(mtc==0) return varMuonTrackMCSP(slc);
    else if(mtc==1) return varMuonTrackRangeP(slc);
    else{
      cout << "[varMuonTrackCombinedP] wtf?" << endl;
      exit(EXIT_FAILURE);
      return varMuonTrackMCSP(slc); //
    }
  }   
  else{
    return -999;
  }   
}); 

const Var varMuonTrackCaloP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double best_calo = 0.;
    if(trk_Muon.bestplane == 0){
      best_calo = trk_Muon.calo0.ke;
    }
    else if(trk_Muon.bestplane == 1){
      best_calo = trk_Muon.calo1.ke;
    }
    else{
      best_calo = trk_Muon.calo2.ke;
    }

    double KE_Muon = best_calo/1000.;
    double E_Muon = KE_Muon+M_MUON;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;

  }
  else{
    return -999;
  }
});

const Var varMuonTrackCaloKE([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double best_calo = 0.;
    if(trk_Muon.bestplane == 0){
      best_calo = trk_Muon.calo0.ke;
    }
    else if(trk_Muon.bestplane == 1){
      best_calo = trk_Muon.calo1.ke;
    }
    else{
      best_calo = trk_Muon.calo2.ke;
    }

    double KE_Muon = best_calo/1000.;
    return KE_Muon;

  }
  else{
    return -999;
  }
});

//====   Truth
//====     track.p has the information of the best match SRTrueParticle
//====     https://github.com/SBNSoftware/sbnanaobj/blob/733b185944bf4565ddb451c9aee3ade8715d4114/sbnanaobj/StandardRecord/SRTrackTruth.h#L33

const Var varMuonTrackMatchedTruthE([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    double E_Muon = trk_Muon.truth.p.genE;
    return E_Muon;
  }
  else{
    return -999.;
  }
});

const Var varMuonTrackMatchedTruthP([](const caf::SRSliceProxy* slc) -> double {
  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){ 
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);
    double E_Muon = trk_Muon.truth.p.genE;
    double P_Muon = sqrt( E_Muon*E_Muon - M_MUON*M_MUON );
    return P_Muon;
  }
  else{
    return -999.;
  }
});

//====   Residual study

const Var varMuonTrackRangePResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackRangeP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackMCSPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackMCSP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackCombinedPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCombinedP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackCaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff;
});

const Var varMuonTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackRangeP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackMCSPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackMCSP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackCombinedPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCombinedP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});

const Var varMuonTrackCaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Muon_Reco = varMuonTrackCaloP(slc);
  double P_Muon_Truth = varMuonTrackMatchedTruthP(slc);
  double P_Muon_Diff = P_Muon_Reco-P_Muon_Truth;
  return P_Muon_Diff/P_Muon_Truth;
});


//==== Proton

//====   Reco

const Var varProtonTrackRangeP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    return trk_Proton.rangeP.p_proton;
  }
  else{
    return -999;
  }
});

const Var varProtonTrackMCSP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    return trk_Proton.mcsP.fwdP_proton;
  }
  else{
    return -999;
  }
});

const Var varProtonTrackCombinedP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    int mtc = varProtonTrackContainedness(slc); // 0 : exiting, 1 : contained
    if(mtc==0) return varProtonTrackMCSP(slc);
    else if(mtc==1) return varProtonTrackRangeP(slc);
    else{
      cout << "[varProtonTrackCombinedP] wtf?" << endl;
      exit(EXIT_FAILURE);
      return varProtonTrackMCSP(slc); //
    }
  }
  else{
    return -999;
  }
});

const Var varProtonTrackCaloP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = 0.;
    if(trk_Proton.bestplane == 0){
      best_calo = trk_Proton.calo0.ke;
    }
    else if(trk_Proton.bestplane == 1){
      best_calo = trk_Proton.calo1.ke;
    }
    else{
      best_calo = trk_Proton.calo2.ke;
    }

    double KE_Proton = best_calo/1000.;
    double E_Proton = KE_Proton+M_PROTON;
    double P_Proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return P_Proton;

  }
  else{
    return -999;
  }
});

const Var varProtonTrackCaloKE([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double best_calo = 0.;
    if(trk_Proton.bestplane == 0){
      best_calo = trk_Proton.calo0.ke;
    }
    else if(trk_Proton.bestplane == 1){
      best_calo = trk_Proton.calo1.ke;
    }
    else{
      best_calo = trk_Proton.calo2.ke;
    }

    double KE_Proton = best_calo/1000.;
    return KE_Proton;

  }
  else{
    return -999;
  }
});

//====   Truth
//====     track.p has the information of the best match SRTrueParticle
//====     https://github.com/SBNSoftware/sbnanaobj/blob/733b185944bf4565ddb451c9aee3ade8715d4114/sbnanaobj/StandardRecord/SRTrackTruth.h#L33

const Var varProtonTrackMatchedTruthE([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    double E_Proton = trk_Proton.truth.p.genE;
    return E_Proton;
  }
  else{
    return -999.;
  }
});

const Var varProtonTrackMatchedTruthP([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varProtonTrackIndex(slc);
  if(protonTrackIndex>=0){
    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);
    double E_Proton = trk_Proton.truth.p.genE;
    double P_Proton = sqrt( E_Proton*E_Proton - M_PROTON*M_PROTON );
    return P_Proton;
  }
  else{
    return -999.;
  }
});

//====   Residual study

const Var varProtonTrackRangePResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackRangeP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff;
});

const Var varProtonTrackMCSPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackMCSP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff;
});

const Var varProtonTrackCombinedPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackCombinedP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff;
});

const Var varProtonTrackCaloPResidual([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackCaloP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff;
});

const Var varProtonTrackRangePResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackRangeP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff/P_Proton_Truth;
});

const Var varProtonTrackMCSPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackMCSP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff/P_Proton_Truth;
});

const Var varProtonTrackCombinedPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackCombinedP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff/P_Proton_Truth;
});

const Var varProtonTrackCaloPResidualFraction([](const caf::SRSliceProxy* slc) -> double {
  double P_Proton_Reco = varProtonTrackCaloP(slc);
  double P_Proton_Truth = varProtonTrackMatchedTruthP(slc);
  double P_Proton_Diff = P_Proton_Reco-P_Proton_Truth;
  return P_Proton_Diff/P_Proton_Truth;
});


//==== Neutrino

//====   Reco

const Var varNeutrinoCaloEnergy([](const caf::SRSliceProxy* slc) -> double {

  double KE_Muon = varMuonTrackCaloKE(slc); // calo-based
  //==== calo-based
  double E_Muon = KE_Muon+M_MUON;

  double KE_Proton = varProtonTrackCaloKE(slc);

  double Eb = bindingE; // TODO binding energy

  double E_Nu = E_Muon + KE_Proton + Eb;

  return E_Nu;

});

const Var varNeutrinoCombinedEnergy([](const caf::SRSliceProxy* slc) -> double {

  double P_Muon = varMuonTrackCombinedP(slc); // comdinbed
  //==== combined
  double E_Muon = sqrt( P_Muon*P_Muon + M_MUON*M_MUON );

  double KE_Proton = varProtonTrackCaloKE(slc);

  double Eb = bindingE; // TODO binding energy will be added later

  double E_Nu = E_Muon + KE_Proton + Eb;

  return E_Nu;

});

//==== https://s3.cern.ch/inspire-prod-files-9/93642a13c46438d97680971700e2013c
const Var varNeutrinoQE([](const caf::SRSliceProxy* slc) -> double {

  int muonTrackIndex = varMuonTrackIndex(slc);
  if(muonTrackIndex>=0){
    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double Eb = bindingE; // TODO binding energy will be added later
    double mu_costh = trk_Muon.costh;
    double mu_p = varMuonTrackCombinedP(slc);
    double mu_E = sqrt(mu_p*mu_p+M_MUON*M_MUON);

    double EQE_num = M_NEUTRON*M_NEUTRON - (M_PROTON-Eb)*(M_PROTON-Eb) - M_MUON*M_MUON + 2.*(M_PROTON-Eb)*mu_E;
    double EQE_den = 2.*(M_PROTON-Eb-mu_E+mu_p*mu_costh);

    return EQE_num/EQE_den;

  }
  else{
    return -999.;
  }

});

//====   Fake-Reco, using truth

const Var varNeutrinoFakeRecoEnergy([](const caf::SRSliceProxy* slc) -> double {

  double E_Muon = varMuonTrackMatchedTruthE(slc);

  double E_Proton = varProtonTrackMatchedTruthE(slc);
  double KE_Proton = E_Proton - M_PROTON;

  double Eb = bindingE; // TODO binding energy will be added later

  double E_Nu = E_Muon + KE_Proton + Eb;

  return E_Nu;

});

//====   Truth

const Var varNeutrinoTruthE = SIMPLEVAR(truth.E);

//====   Residual

const Var varNeutrinoCaloEnergyResidual([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoCaloEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff;

});

const Var varNeutrinoCombinedEnergyResidual([](const caf::SRSliceProxy* slc) -> double {
  
  double E_Neutrino_Reco = varNeutrinoCombinedEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff;

});

const Var varNeutrinoQEResidual([](const caf::SRSliceProxy* slc) -> double {
  
  double E_Neutrino_Reco = varNeutrinoQE(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff;

});

const Var varNeutrinoCaloEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoCaloEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff/E_Neutrino_Truth;

});

const Var varNeutrinoCombinedEnergyResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoCombinedEnergy(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff/E_Neutrino_Truth;

});

const Var varNeutrinoQEResidualFraction([](const caf::SRSliceProxy* slc) -> double {

  double E_Neutrino_Reco = varNeutrinoQE(slc);
  double E_Neutrino_Truth = varNeutrinoTruthE(slc);
  double E_Neutrino_Diff = E_Neutrino_Reco-E_Neutrino_Truth;
  return E_Neutrino_Diff/E_Neutrino_Truth;

});

const Var varMuonTruthP([](const caf::SRSliceProxy* slc) -> double {

  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 13 ){
      double this_E = slc->truth.prim.at(i).genE;
      if(this_E>max_E){
        max_E = this_E;
      }
    }
  }
  double max_P = sqrt(max_E*max_E-M_MUON*M_MUON);
  return max_P;

});
const Var varProtonTruthP([](const caf::SRSliceProxy* slc) -> double {

  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
      double this_E = slc->truth.prim.at(i).genE;
      if(this_E>max_E){
        max_E = this_E;
      }
    }
  }
  double max_P = sqrt(max_E*max_E-M_PROTON*M_PROTON);
  return max_P;

});


//==== TEST

const Var varTruthMuonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

  //cout << "[varGetMatchedTrackForProton] called" << endl;

  int truth_index(-1);
  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 13 ){
      double this_E = slc->truth.prim.at(i).genE;
      if(this_E>max_E){
        max_E = this_E;
        truth_index = i;
      }
    }
  }

  if(truth_index<0) return -999.;

  //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
  double truth_E = slc->truth.prim.at(truth_index).genE;

  int PTrackInd(-1);
  for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
    auto const& trk = slc->reco.trk.at(i);
    double this_matched_E = trk.truth.p.genE;
    if( abs(trk.truth.p.pdg)==13 && fabs(truth_E-this_matched_E)/truth_E < 0.01 ){
      PTrackInd = i;
    }
  }
  if(PTrackInd<0) return -999.;
  //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
  //cout << "[varGetMatchedTrackForProton] Found track's matched genE = " << slc->reco.trk.at(PTrackInd).truth.p.genE << endl;

  return PTrackInd;

});

const Var varTruthMuonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){

    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Proton, Chi2Muon;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    return Chi2Proton;

  }
  else{
    return 999999;
  }
});

const Var varTruthMuonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthMuonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){
    
    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Proton, Chi2Muon;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    return Chi2Muon;

  }
  else{
    return 999999;
  }

});

const Var varTruthProtonMatchedTrackIndex([](const caf::SRSliceProxy* slc) -> double {

  //cout << "[varGetMatchedTrackForProton] called" << endl;

  int truth_index(-1);
  double max_E(-999);
  for(std::size_t i(0); i < slc->truth.prim.size(); ++i){
    if( abs(slc->truth.prim.at(i).pdg)== 2212 ){
      double this_E = slc->truth.prim.at(i).genE;
      if(this_E>max_E){
        max_E = this_E;
        truth_index = i;
      }
    }
  }

  if(truth_index<0) return -999.;

  //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
  double truth_E = slc->truth.prim.at(truth_index).genE;

  int PTrackInd(-1);
  for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
    auto const& trk = slc->reco.trk.at(i);
    double this_matched_E = trk.truth.p.genE;
    if( abs(trk.truth.p.pdg)==2212 && fabs(truth_E-this_matched_E)/truth_E < 0.01 ){
      PTrackInd = i;
    }
  }
  if(PTrackInd<0) return -999.;
  //cout << "[varGetMatchedTrackForProton] Truth proton genE = " << slc->truth.prim.at(truth_index).genE << endl;
  //cout << "[varGetMatchedTrackForProton] Found track's matched genE = " << slc->reco.trk.at(PTrackInd).truth.p.genE << endl;

  return PTrackInd;

});

const Var varTruthProtonMatchedTrackChi2Proton([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){

    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Proton, Chi2Muon;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    return Chi2Proton;

  }
  else{
    return 999999;
  }
});

const Var varTruthProtonMatchedTrackChi2Muon([](const caf::SRSliceProxy* slc) -> double {
  int protonTrackIndex = varTruthProtonMatchedTrackIndex(slc);
  if(protonTrackIndex>=0){
    
    auto const& trk = slc->reco.trk.at(protonTrackIndex);
    float Chi2Proton, Chi2Muon;
    if(trk.bestplane == 0){
      Chi2Proton = trk.chi2pid0.chi2_proton;
      Chi2Muon = trk.chi2pid0.chi2_muon;
    }
    else if(trk.bestplane == 1){
      Chi2Proton = trk.chi2pid1.chi2_proton;
      Chi2Muon = trk.chi2pid1.chi2_muon;
    }
    else{
      Chi2Proton = trk.chi2pid2.chi2_proton;
      Chi2Muon = trk.chi2pid2.chi2_muon;
    }

    return Chi2Muon;

  }
  else{
    return 999999;
  }

});

#endif
