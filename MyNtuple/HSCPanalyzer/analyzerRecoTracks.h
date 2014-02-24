#ifndef ANALYZERRECOTRACKS_H
#define ANALYZERRECOTRACKS_H
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "histograms.h"
#include <iostream>
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------


namespace reco
{

  std::vector<Track_s> trackCollection;

  void getCandidateTrackCollection(){

    trackCollection.clear();
    int isolation = 0;
    int all=0;


    for(int i=0; i<nTrack; i++){

      if(Track[i].pt<50.)                                            continue;
      //if(std::abs(Track[i].eta)>2.1)                                 continue;
      //if(std::abs(Track[i].eta)>1.42 && std::abs(Track[i].eta<1.65)) continue;
      //if(std::abs(Track[i].eta)>0.15 && std::abs(Track[i].eta<0.35)) continue;
      //if(std::abs(Track[i].eta)>1.55 && std::abs(Track[i].eta<1.85)) continue;
      //if(Track[i].d0>0.2)                                            continue;
      //if(Track[i].dz>5)                                              continue;
      //if(Track[i].numberOfValidHits<5)                               continue;
      //if(Track[i].numberOfLostHits>0)                                continue;
      //if(Track[i].trackerExpectedHitsInner_numberOfLostHits>0)       continue;
      //if(Track[i].trackerExpectedHitsOuter_numberOfLostHits<3)       continue;


      // Isolation
      double sumPt = 0;
      for(int j=0; j<nTrack; j++){

	if(i==j) continue;
	double dPhi = std::abs(Track[j].phi - Track[i].phi);
	double dEta = std::abs(Track[j].eta - Track[i].eta);
	double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2));

	if(dR<0.3)  sumPt += Track[j].pt;
      }
      
      if(sumPt/Track[i].pt>0.05)                                     continue;
      
      trackCollection.push_back(Track[i]);
    }
  }

}
#endif

