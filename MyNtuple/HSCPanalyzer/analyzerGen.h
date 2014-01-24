#ifndef ANALYZERGEN_H
#define ANALYZERGEN_H
//-----------------------------------------------------------------------------
#include "analyzer.h"
#include <iostream>
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------

struct GenParticle_s chipGenParticle;
struct GenParticle_s chimGenParticle;

struct SimTrack_s chipSimTrack;
struct SimTrack_s chimSimTrack;

struct SimTrack_s chi0SimTrack[10];
struct SimTrack_s piSimTrack[1000];

struct SimVertex_s chipOriginSimVertex;
struct SimVertex_s chimOriginSimVertex;

struct SimVertex_s chipDecaySimVertex;
struct SimVertex_s chimDecaySimVertex;

bool decayedChip = false;
bool decayedChim = false;

int nChi0 = 0;
int nPi   = 0;

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

/*****************************************
 Find chargino in SimTrack collection 
*****************************************/
void findChipmInSimTrackCollection(){

  bool zeroChip = true;
  bool zeroChim = true;

  int nChipm = 0;

  for(int i=0; i<nSimTrack; i++){

    if(abs(SimTrack[i].type)==1000024){
	      
      if(SimTrack[i].type>0 && zeroChip){
	chipSimTrack = SimTrack[i];
	zeroChip = false;
	nChipm += 1;
      }
      else if(SimTrack[i].type<0 && zeroChim){
	chimSimTrack = SimTrack[i];
	zeroChim = false;
	nChipm += 1;
      }	      
    }
    if(!zeroChip && !zeroChim) break;
  }
  if(zeroChip || zeroChim) cout<<"Not all chargino tracks were found!"<<endl;
}


/*****************************************
 Find neutralino in SimTrack collection 
*****************************************/
void findChi0InSimTrackCollection(){

  nChi0 = 0;
  for(int i=0; i<nSimTrack; i++){
    if(abs(SimTrack[i].type)==1000022){
      chi0SimTrack[nChi0] = SimTrack[i];
      nChi0 += 1;
    }
  }
}

/*****************************************
 Find neutralino in SimTrack collection 
*****************************************/
void findPiInSimTrackCollection(){

  nPi = 0;
  for(int i=0; i<nSimTrack; i++){
    if(abs(SimTrack[i].type)==211){
      piSimTrack[nChi0] = SimTrack[i];
      nPi += 1;
    }
  }
}


/****************************************************
 Find chargino origin vertex in SimVertex collection 
****************************************************/
void findChipmOriginInSimVertexCollection(){

  bool zeroChip = true;
  bool zeroChim = true;

  for(int i=0; i<nSimVertex; i++){

    if(SimVertex[i].vertexId==chipSimTrack.vertIndex){

      chipOriginSimVertex = SimVertex[i];
      zeroChip=false;
    }
    if(SimVertex[i].vertexId==chimSimTrack.vertIndex){

      chimOriginSimVertex = SimVertex[i];
      zeroChim=false;
    }
    if(!zeroChip && !zeroChim) break;
  }

  if(zeroChip || zeroChim) cout<<"Not all chargino origin vertices were found!"<<endl;
}


/****************************************************
 Find chargino decay vertex in SimVertex collection 
****************************************************/
void findChipmDecayInSimVertexCollection(){

  bool zeroChip = true;
  bool zeroChim = true;

  bool decayed   = false;
  bool decayedPi = false;

  decayedChip = false;
  decayedChim = false;

  for(int i=0; i<nSimVertex; i++){

    decayed = false;

    if(SimVertex[i].parentIndex==chipSimTrack.trackId || SimVertex[i].parentIndex==chimSimTrack.trackId){

      for(int j=0; j<nChi0; j++){
	if(SimVertex[i].vertexId == chi0SimTrack[j].vertIndex) decayed = true;
      }
      for(int j=0; j<nPi; j++){
	if(SimVertex[i].vertexId == piSimTrack[j].vertIndex) decayedPi = true;
      }
      if(!decayed) continue;

      if(!decayedPi){
	//continue;	
      }
 
      if(SimVertex[i].parentIndex==chipSimTrack.trackId){
	
	chipDecaySimVertex = SimVertex[i];
	zeroChip    = false;
	decayedChip = true;
      }
      else if(SimVertex[i].parentIndex==chimSimTrack.trackId){
	
	chimDecaySimVertex = SimVertex[i];
	zeroChim    = false;
	decayedChim = true;
      }
    }
    if(!zeroChip && !zeroChim) break;
  }
}

/****************************************************
 Find decay particles of chargino decay
****************************************************/
void findDecayParticles(){


  for(int i=0; i<nSimTrack; i++){

    if(decayedChip){
      if(SimTrack[i].vertIndex == chipDecaySimVertex.vertexId){
	//cout<<"Chip decay particle: "<<SimTrack[i].type<<endl;
      }
    }

    if(decayedChim){
      
      if(SimTrack[i].vertIndex == chimDecaySimVertex.vertexId){
	//cout<<"Chim decay particle: "<<SimTrack[i].type<<endl;
      }
    }
  }
}



#endif

