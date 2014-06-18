#ifndef ANALYZERFUNCTIONSGENSIM_H
#define ANALYZERFUNCTIONSGENSIM_H
//-----------------------------------------------------------------------------
#include "analyzerGenSim.h"
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

std::vector<SimTrack_s> AllPions;
std::vector<SimTrack_s> PionsFromDecay;
std::vector<SimTrack_s> Chi0FromDecay;

bool decayedChip = false;
bool decayedChim = false;

int nChi0 = 0;
int nPi   = 0;

int decayed = 0;
int decayVertexFound=0;
int decayVertexFound2=0;
int decayVertexFound3=0;
int decayVertexFound4=0;

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

  if(zeroChip || zeroChim){
    cout<<"To few charginos in GenParticle collection!"<<endl;
    cout<<"zeroChip = "<<zeroChip<<endl;
    cout<<"zeroChim = "<<zeroChim<<endl;
  }
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
  if(zeroChip || zeroChim){
    cout<<"Not all chargino tracks were found!"<<endl;
    cout<<"zeroChip = "<<zeroChip<<endl;
    cout<<"zeroChim = "<<zeroChim<<endl;
  }
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
      Chi0FromDecay.push_back(SimTrack[i]);
    }
  }
}

/*****************************************
 Find all pions in SimTrack collection 
*****************************************/
void findPiInSimTrackCollection(){

  nPi = 0;
  for(int i=0; i<nSimTrack; i++){
    if(abs(SimTrack[i].type)==211){
      piSimTrack[nPi] = SimTrack[i];
      AllPions.push_back(SimTrack[i]);
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

  //decayedChip = false;
  //decayedChim = false;

  for(int i=0; i<nSimVertex; i++){

    decayed = false;

    if(SimVertex[i].parentIndex==chipSimTrack.trackId || SimVertex[i].parentIndex==chimSimTrack.trackId){

      decayVertexFound += 1;

      for(unsigned int j=0; j<nChi0; j++){
	if(SimVertex[i].vertexId == chi0SimTrack[j].vertIndex) decayed = true;
      }
      for(unsigned int j=0; j<nPi; j++){
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
	decayVertexFound2 += 1;
      }
      if(SimVertex[i].parentIndex==chimSimTrack.trackId){
	
	chimDecaySimVertex = SimVertex[i];
	zeroChim    = false;
	decayedChim = true;
	decayVertexFound2 += 1;
      }
    }
    if(!zeroChip && !zeroChim) break;
  }
}

/****************************************************
 Find decay particles of chargino decay
****************************************************/
void findDecayParticles(){


  for(unsigned int i=0; i<nSimTrack; i++){

    if(decayedChip){
      if(SimTrack[i].vertIndex == chipDecaySimVertex.vertexId){

	cout<<"chipSimTrack.momentum_pt = "<<chipSimTrack.momentum_pt<<endl;
	cout<<"SimTrack[i].type = "<<SimTrack[i].type<<endl;
	cout<<"SimTrack[i].pt = "<<SimTrack[i].momentum_pt<<endl;

	if(abs(SimTrack[i].type)==211) PionsFromDecay.push_back(SimTrack[i]);
	else decayVertexFound4 += 1;
	if(abs(SimTrack[i].type) !=211 && abs(SimTrack[i].type) !=1000022)
	  {
	    cout<<"Chip decay particle: "<<SimTrack[i].type<<endl;
	    cout<<"Chip decay particle: "<<SimTrack[i].momentum_pt<<endl;
	  }
      }
    }

    if(decayedChim){

     
      
      if(SimTrack[i].vertIndex == chimDecaySimVertex.vertexId){

	cout<<"chimSimTrack.momentum_pt = "<<chimSimTrack.momentum_pt<<endl;
	cout<<"SimTrack[i].type = "<<SimTrack[i].type<<endl;
	cout<<"SimTrack[i].pt = "<<SimTrack[i].momentum_pt<<endl;
	
	if(abs(SimTrack[i].type)==211) PionsFromDecay.push_back(SimTrack[i]);
	else decayVertexFound4 += 1;
	if(abs(SimTrack[i].type) !=211 && abs(SimTrack[i].type) !=1000022)
	  {
	    cout<<"Chim decay particle: "<<SimTrack[i].type<<endl;
	    cout<<"Chim decay particle: "<<SimTrack[i].momentum_pt<<endl;
	  }
      }
    }
  }
}



#endif

