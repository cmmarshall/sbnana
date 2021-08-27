#ifndef myProtonSelection_h
#define myProtonSelection_h

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "myConstants.h"

//==== Longest of the (proton-chi2) tracks
//==== - loop over tracks, and check its chi2
//==== - if multiple tracks are found, pick the longest

const Var varProtonTrackIndex([](const caf::SRSliceProxy* slc) -> int {
  int PTrackInd(-1);
  double max_TrackLength = -999;
  for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
    auto const& trk = slc->reco.trk.at(i);

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

    //if( Chi2Muon > 60 && Chi2Proton < 30 ){
    if( Chi2Proton<60 ){

			if(trk.len > max_TrackLength){
				max_TrackLength = trk.len;
				PTrackInd = i;
			}

    }

  }
  return PTrackInd;
});

const Var varProtonTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
  int protonTrackIndex = varProtonTrackIndex(slc);
  //==== -1 : No track found
  //==== 0 : Exiting
  //==== 1 : Contained
  if(protonTrackIndex>=0){

    auto const& trk_Proton = slc->reco.trk.at(protonTrackIndex);

    double XMargin = 25.;
    double YMargin = 25.;
    double ZMarginUp = 30.;
    double ZMarginDown = 50.;

    bool isContained;
    //==== Cryo0
    if( trk_Proton.end.x < 0 ){
      isContained = ( !isnan(trk_Proton.end.x) &&
                    ( trk_Proton.end.x < -71.1 - XMargin && trk_Proton.end.x > -369.33 + XMargin ) &&
                    !isnan(trk_Proton.end.y) &&
                    ( trk_Proton.end.y > -181.7 + YMargin && trk_Proton.end.y < 134.8 - YMargin ) &&
                    !isnan(trk_Proton.end.z) &&
                    ( trk_Proton.end.z > -895.95 + ZMarginUp && trk_Proton.end.z < 895.95 - ZMarginDown ) );
    }
    //==== Cryo1
    else{
      isContained = ( !isnan(trk_Proton.end.x) &&
                    ( trk_Proton.end.x > 71.1 + XMargin && trk_Proton.end.x < 369.33 - XMargin ) &&
                    !isnan(trk_Proton.end.y) &&
                    ( trk_Proton.end.y > -181.7 + YMargin && trk_Proton.end.y < 134.8 - YMargin ) &&
                    !isnan(trk_Proton.end.z) &&
                    ( trk_Proton.end.z > -895.95 + ZMarginUp && trk_Proton.end.z < 895.95 - ZMarginDown ) );
    }
    if(isContained) return 1;
    else return 0;
  }
  else{
    return -1;
  }

});


#endif
