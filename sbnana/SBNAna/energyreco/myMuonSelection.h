#ifndef myMuonSelection_h
#define myMuonSelection_h

#include "sbnana/SBNAna/Cuts/NumuCutsIcarus202106.h"
#include "sbnana/SBNAna/Cuts/Cuts.h"
#include "myConstants.h"

//==== Longest picked, and then check chi2
//==== If the longest track does not satisfy chi2, then no track is assigned for this slice

const Var varLongestTrackIndex([](const caf::SRSliceProxy* slc) -> int {
  int PTrackInd(-1);
  double max_TrackLength = -999;
  for(std::size_t i(0); i < slc->reco.trk.size(); ++i){
    auto const& trk = slc->reco.trk.at(i);
    if(trk.len > max_TrackLength){
      max_TrackLength = trk.len;
      PTrackInd = i;
    }
  }
  return PTrackInd;
});

//==== Longest of the (muon-chi2) tracks
//==== - loop over tracks, and check its chi2
//==== - if multiple tracks are found, pick the longest

const Var varMuonTrackIndex([](const caf::SRSliceProxy* slc) -> int {
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

    if( Chi2Proton > 60 && Chi2Muon < 30 ){

			if(trk.len > max_TrackLength){
				max_TrackLength = trk.len;
				PTrackInd = i;
			}

    }

  }
  return PTrackInd;
});

const Var varisMuonLongest([](const caf::SRSliceProxy* slc) -> int {

  int longestTrackIndex = varLongestTrackIndex(slc);
  int muonTrackIndex = varMuonTrackIndex(slc);


  //==== if muonTrackIndex<0, longestTrackIndex is also negative (i.e., not found)
  if(muonTrackIndex<0){
    return -2;
  }
  //==== though muonTrackIndex is found, longestTrackIndex can be negative.
  //==== e.g., the longest track failed the chi2, but the trailing track(s) passes
  else if(muonTrackIndex>0 && longestTrackIndex<0){
    return -1;
  }
  //==== now, both are found, but safe-gaurd
  else if(muonTrackIndex>=0 && longestTrackIndex>=0){
    return 1;
  }
  else{
    //cout << "[varisMuonLongest] longestTrackIndex = " << longestTrackIndex << ", muonTrackIndex = " << muonTrackIndex << endl;
    return 0;
  }

  //==== if muonTrackIndex found (i.e., >0), longestTrackIndex is always found


});

const Var varMuonTrackContainedness([](const caf::SRSliceProxy* slc) -> int {
  int muonTrackIndex = varMuonTrackIndex(slc);
  //==== -1 : No track found
  //==== 0 : Exiting
  //==== 1 : Contained
  if(muonTrackIndex>=0){

    auto const& trk_Muon = slc->reco.trk.at(muonTrackIndex);

    double XMargin = 25.;
    double YMargin = 25.;
    double ZMarginUp = 30.;
    double ZMarginDown = 50.;

    bool isContained;
    //==== Cryo0
    if( trk_Muon.end.x < 0 ){
      isContained = ( !isnan(trk_Muon.end.x) &&
                    ( trk_Muon.end.x < -71.1 - XMargin && trk_Muon.end.x > -369.33 + XMargin ) &&
                    !isnan(trk_Muon.end.y) &&
                    ( trk_Muon.end.y > -181.7 + YMargin && trk_Muon.end.y < 134.8 - YMargin ) &&
                    !isnan(trk_Muon.end.z) &&
                    ( trk_Muon.end.z > -895.95 + ZMarginUp && trk_Muon.end.z < 895.95 - ZMarginDown ) );
    }
    //==== Cryo1
    else{
      isContained = ( !isnan(trk_Muon.end.x) &&
                    ( trk_Muon.end.x > 71.1 + XMargin && trk_Muon.end.x < 369.33 - XMargin ) &&
                    !isnan(trk_Muon.end.y) &&
                    ( trk_Muon.end.y > -181.7 + YMargin && trk_Muon.end.y < 134.8 - YMargin ) &&
                    !isnan(trk_Muon.end.z) &&
                    ( trk_Muon.end.z > -895.95 + ZMarginUp && trk_Muon.end.z < 895.95 - ZMarginDown ) );
    }
    if(isContained) return 1;
    else return 0;
  }
  else{
    return -1;
  }

});

#endif
