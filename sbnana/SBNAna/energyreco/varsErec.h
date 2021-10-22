#ifndef varsErec_h
#define varsErec_h

const MultiVar varProtonRangeKE([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> KE;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  // if the event involved a muon
  if (muonTrackIndex >= 0) {
    int i=0;
    // loop over the tracks
    for (auto const& trk : slc->reco.trk) {
      if (i != muonTrackIndex) {
        double Chi2Proton, Chi2Pion;
        if (trk.bestplane == 0) {
          Chi2Proton = trk.chi2pid0.chi2_proton;
          Chi2Pion = trk.chi2pid0.chi2_pion;
        }
        else if (trk.bestplane == 1) {
          Chi2Proton = trk.chi2pid1.chi2_proton;
          Chi2Pion = trk.chi2pid1.chi2_pion;
        }
        else {
          Chi2Proton = trk.chi2pid2.chi2_proton;
          Chi2Pion = trk.chi2pid2.chi2_pion;
        }
        // if the track was left by a proton
        if (Chi2Proton < Chi2Pion or (Chi2Proton-Chi2Pion) < 75) {
          double p=trk.rangeP.p_proton;
          double E = sqrt( p*p + M_PROTON*M_PROTON );
          KE.push_back( E - M_PROTON );
        }}
      i += 1;
    }}
  return KE;
});

const MultiVar varTrkKE([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> KE;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  // if the event involved a muon
  if (muonTrackIndex >= 0) {
    int i=0;
    // loop over the tracks
    for (auto const& trk : slc->reco.trk) {
      if (i != muonTrackIndex) {
        double p=trk.rangeP.p_proton;
        double E = sqrt( p*p + M_PROTON*M_PROTON );
        KE.push_back( E - M_PROTON );
      }
      i += 1;
    }}
  return KE;
});

const MultiVar varPionRangeE([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> E;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  // if the event involved a muon
  if (muonTrackIndex >= 0) {
    int i=0;
    // loop over the tracks
    for (auto const& trk : slc->reco.trk) {
      if (i != muonTrackIndex) {
        double Chi2Proton, Chi2Pion;
        if (trk.bestplane == 0) {
          Chi2Proton = trk.chi2pid0.chi2_proton;
          Chi2Pion = trk.chi2pid0.chi2_pion;
        }
        else if (trk.bestplane == 1) {
          Chi2Proton = trk.chi2pid1.chi2_proton;
          Chi2Pion = trk.chi2pid1.chi2_pion;
        }
        else {
          Chi2Proton = trk.chi2pid2.chi2_proton;
          Chi2Pion = trk.chi2pid2.chi2_pion;
        }
        // if the track was left by a proton
        if (Chi2Proton > Chi2Pion and (Chi2Proton-Chi2Pion) > 75) {
          double p=trk.rangeP.p_pion;
          E.push_back( sqrt(p*p + M_PION*M_PION) );
        }}
      i += 1;
    }}
  return E;
});

const MultiVar varProtonCombinedKE([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> KE;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  // if the event involved a muon
  if (muonTrackIndex >= 0) {
    int i=0;
    // loop over the tracks
    for (auto const& trk : slc->reco.trk) {
      if (i != muonTrackIndex) {
        double Chi2Proton, Chi2Pion;
        if (trk.bestplane == 0) {
          Chi2Proton = trk.chi2pid0.chi2_proton;
          Chi2Pion = trk.chi2pid0.chi2_pion;
        }
        else if (trk.bestplane == 1) {
          Chi2Proton = trk.chi2pid1.chi2_proton;
          Chi2Pion = trk.chi2pid1.chi2_pion;
        }
        else {
          Chi2Proton = trk.chi2pid2.chi2_proton;
          Chi2Pion = trk.chi2pid2.chi2_pion;
        }
        // if the track was left by a proton
        if (Chi2Proton < Chi2Pion or (Chi2Proton-Chi2Pion) < 75) {
          // check if the particle escaped from the detector
          double XMargin = 25.;
          double YMargin = 25.;
          double ZMarginUp = 30.;
          double ZMarginDown = 50.;
          bool isContained;
          // cryo0
          if ( trk.end.x < 0 ) {
            isContained = ( !isnan(trk.end.x) &&
                          ( trk.end.x < -71.1 - XMargin && trk.end.x > -369.33 + XMargin ) &&
                            !isnan(trk.end.y) &&
                          ( trk.end.y > -181.7 + YMargin && trk.end.y < 134.8 - YMargin ) &&
                            !isnan(trk.end.z) &&
                          ( trk.end.z > -895.95 + ZMarginUp && trk.end.z < 895.95 - ZMarginDown ) );
          }
          // cryo1
          else {
            isContained = ( !isnan(trk.end.x) &&
                          ( trk.end.x > 71.1 + XMargin && trk.end.x < 369.33 - XMargin ) &&
                            !isnan(trk.end.y) &&
                          ( trk.end.y > -181.7 + YMargin && trk.end.y < 134.8 - YMargin ) &&
                            !isnan(trk.end.z) &&
                          ( trk.end.z > -895.95 + ZMarginUp && trk.end.z < 895.95 - ZMarginDown ) );
          }
          double p=0.;
          if (isContained) {
            p = trk.rangeP.p_proton;
          }
          else {
            p = trk.mcsP.fwdP_proton;
          }
          double E = sqrt( p*p + M_PROTON*M_PROTON );
          KE.push_back( E - M_PROTON );
        }}
      i += 1;
    }}
  return KE;
});

const MultiVar varPionCombinedE([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> E;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  // if the event involved a muon
  if (muonTrackIndex >= 0) {
    int i=0;
    // loop over the tracks
    for (auto const& trk : slc->reco.trk) {
      if (i != muonTrackIndex) {
        double Chi2Proton, Chi2Pion;
        if (trk.bestplane == 0) {
          Chi2Proton = trk.chi2pid0.chi2_proton;
          Chi2Pion = trk.chi2pid0.chi2_pion;
        }
        else if (trk.bestplane == 1) {
          Chi2Proton = trk.chi2pid1.chi2_proton;
          Chi2Pion = trk.chi2pid1.chi2_pion;
        }
        else {
          Chi2Proton = trk.chi2pid2.chi2_proton;
          Chi2Pion = trk.chi2pid2.chi2_pion;
        }
        // if the track was left by a proton
        if (Chi2Proton > Chi2Pion and (Chi2Proton-Chi2Pion) > 75) {
          // check if the particle escaped from the detector
          double XMargin = 25.;
          double YMargin = 25.;
          double ZMarginUp = 30.;
          double ZMarginDown = 50.;
          bool isContained;
          // cryo0
          if ( trk.end.x < 0 ) {
            isContained = ( !isnan(trk.end.x) &&
                          ( trk.end.x < -71.1 - XMargin && trk.end.x > -369.33 + XMargin ) &&
                            !isnan(trk.end.y) &&
                          ( trk.end.y > -181.7 + YMargin && trk.end.y < 134.8 - YMargin ) &&
                            !isnan(trk.end.z) &&
                          ( trk.end.z > -895.95 + ZMarginUp && trk.end.z < 895.95 - ZMarginDown ) );
          }
          // cryo1
          else {
            isContained = ( !isnan(trk.end.x) &&
                          ( trk.end.x > 71.1 + XMargin && trk.end.x < 369.33 - XMargin ) &&
                            !isnan(trk.end.y) &&
                          ( trk.end.y > -181.7 + YMargin && trk.end.y < 134.8 - YMargin ) &&
                            !isnan(trk.end.z) &&
                          ( trk.end.z > -895.95 + ZMarginUp && trk.end.z < 895.95 - ZMarginDown ) );
          }
          double p=0.;
          if (isContained) {
            p = trk.rangeP.p_pion;
          }
          else {
            p = trk.mcsP.fwdP_pion;
          }
          E.push_back( sqrt(p*p + M_PION*M_PION) );
        }}
      i += 1;
    }}
  return E;
});

const Var varShwE([](const caf::SRSliceProxy* slc) -> double {
  double total = 0.;
  for (auto const& shw : slc->reco.shw) {
    total += (shw.bestplane_energy/1000.);
  }
  return total;
});

const Cut cutExitingHadrons([](const caf::SRSliceProxy* slc) {
  int count=0;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  if (muonTrackIndex >= 0) {
    int i=0;
    for (auto const& trk : slc->reco.trk) {
      if (i != muonTrackIndex) {
        double XMargin = 25.;
        double YMargin = 25.;
        double ZMarginUp = 30.;
        double ZMarginDown = 50.;
        bool isContained;
        if ( trk.end.x < 0 ) {
          isContained = ( !isnan(trk.end.x) &&
                        ( trk.end.x < -71.1 - XMargin && trk.end.x > -369.33 + XMargin ) &&
                          !isnan(trk.end.y) &&
                        ( trk.end.y > -181.7 + YMargin && trk.end.y < 134.8 - YMargin ) &&
                          !isnan(trk.end.z) &&
                        ( trk.end.z > -895.95 + ZMarginUp && trk.end.z < 895.95 - ZMarginDown ) );
        }
        else {
          isContained = ( !isnan(trk.end.x) &&
                        ( trk.end.x > 71.1 + XMargin && trk.end.x < 369.33 - XMargin ) &&
                          !isnan(trk.end.y) &&
                        ( trk.end.y > -181.7 + YMargin && trk.end.y < 134.8 - YMargin ) &&
                          !isnan(trk.end.z) &&
                        ( trk.end.z > -895.95 + ZMarginUp && trk.end.z < 895.95 - ZMarginDown ) );
        }
        if (!isContained) {
          count += 1;
        }}
      i+=1;}}
  return count==0;
});

#endif
