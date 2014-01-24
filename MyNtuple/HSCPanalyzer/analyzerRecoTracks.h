#ifndef ANALYZERRECOTRACKS_H
#define ANALYZERRECOTRACKS_H
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include <iostream>
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------


namespace reco
{

  struct GenParticle_s chipGenParticle;
  struct GenParticle_s chimGenParticle;

  struct Track_s chipSimTrack;
  struct Track_s chimSimTrack;


  int nChi0 = 0;

  /*****************************************
 Find chargino in GenParticle collection 
  *****************************************/
  void findChipmInGenParticleCollection(){

    bool zeroChip = true;
    bool zeroChim = true;

    for(int i=0; i<nGenParticle; i++){

      if(abs(GenParticle[i].pdgId)==1000024){

	if(GenParticle[i].pdgId>0 && zeroChip){

	  chipGenParticle = GenParticle[i];
	  zeroChip = false;
	}
	else if(GenParticle[i].pdgId<0 && zeroChim){

	  chimGenParticle = GenParticle[i];
	  zeroChim = false;
	}	      
      }
      if(!zeroChip && !zeroChim) break;
    }

    if(zeroChip || zeroChim) cout<<"To few charginos in GenParticle collection!"<<endl;
  }

}
#endif

