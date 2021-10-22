#ifndef varsTransverseP_h
#define varsTransverseP_h

const MultiVar varMuonP([](const caf::SRSliceProxy* slc) -> std::vector<double> {
  std::vector<double> momentum;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  if (muonTrackIndex >= 0) {
    int i=0;
    for (auto const& trk : slc->reco.trk) {
      if (i == muonTrackIndex) {
        double p = varMuonTrackCombinedP(slc);
        double theta = acos(trk.costh);
        double phi = trk.phi;
        double px = p*sin(theta)*cos(phi);
        double py = p*sin(theta)*sin(phi);
        double pz = p*cos(theta);
        momentum.push_back(px);
        momentum.push_back(py);
        momentum.push_back(pz);
      }
      i+=1;
    }}
  return momentum;
});     

const Var varPin([](const caf::SRSliceProxy* slc) -> double {
  double Pin=0.;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  if (muonTrackIndex >= 0) {
    int i=0;
    std::vector<double> p_muon = varMuonP(slc);
    double px_muon = p_muon[0];
    double py_muon = p_muon[1];
    double pxy_muon = sqrt( px_muon*px_muon + py_muon*py_muon ); 
    Pin += pxy_muon;
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
        double p=0;
        if (Chi2Proton < Chi2Pion or (Chi2Proton-Chi2Pion) < 75) {
          if (isContained) {
            p = trk.rangeP.p_proton;
          }
          else {
            p = trk.mcsP.fwdP_proton;
          }}
        if (Chi2Proton > Chi2Pion and (Chi2Proton-Chi2Pion) > 75) {
          if (isContained) {
            p = trk.rangeP.p_pion;
          }
          else {
            p = trk.mcsP.fwdP_pion;
          }}
        double theta = acos(trk.costh);
        double phi = trk.phi;
        double px = p*sin(theta)*cos(phi);
        double py = p*sin(theta)*sin(phi);
        Pin += px*(px_muon/pxy_muon)+py*(py_muon/pxy_muon);
      } 
    i+=1;
  }}
  return Pin;
});
        
const Var varPout([](const caf::SRSliceProxy* slc) -> double {
  double Pout=0.;
  int  muonTrackIndex = varMuonTrackIndex(slc);
  if (muonTrackIndex >= 0) {
    int i=0;
    std::vector<double> p_muon = varMuonP(slc);
    double px_muon = p_muon[0];
    double py_muon = p_muon[1];
    double pxy_muon = sqrt( px_muon*px_muon + py_muon*py_muon ); 
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
        double p=0;
        if (Chi2Proton < Chi2Pion or (Chi2Proton-Chi2Pion) < 75) {
          if (isContained) {
            p = trk.rangeP.p_proton;
          }
          else {
            p = trk.mcsP.fwdP_proton;
          }}
        if (Chi2Proton > Chi2Pion and (Chi2Proton-Chi2Pion) > 75) {
          if (isContained) {
            p = trk.rangeP.p_pion;
          }
          else {
            p = trk.mcsP.fwdP_pion;
          }}
        double theta = acos(trk.costh);
        double phi = trk.phi;
        double px = p*sin(theta)*cos(phi);
        double py = p*sin(theta)*sin(phi);
        Pout += px*(py_muon/pxy_muon)-py*(px_muon/pxy_muon);
      } 
    i+=1;
  }}
  return Pout;
});

const Var varPtrans([](const caf::SRSliceProxy* slc) -> double {
  double Pin = varPin(slc);
  double Pout = varPout(slc);
  double Ptrans = sqrt( Pin*Pin + Pout*Pout );
  return Ptrans;
});

double Mneutron = 0.939565413; //GeV
double Ebinding = 0.030; //GeV

const Var varErecNeutron([](const caf::SRSliceProxy* slc) -> double {
  double Ptrans = varPtrans(slc);
  double Eneutron = sqrt( Ptrans*Ptrans + Mneutron*Mneutron );
  double Tneutron = Eneutron - Mneutron;
  double ErecNeutron = Tneutron + Ebinding;
  return ErecNeutron;
});

#endif
