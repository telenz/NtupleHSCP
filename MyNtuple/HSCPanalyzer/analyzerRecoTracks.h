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

  struct PFJet_s leadingJet;
  std::vector<PFJet_s> jetCollection;
  
  std::vector<Track_s> trackCollection;

  bool zeroChip = true;
  bool zeroChim = true;

  int nChi0 = 0;

  /*****************************************
 Find chargino in GenParticle collection 
  *****************************************/
  void findChipmInGenParticleCollection(){

    zeroChip = true;
    zeroChim = true;

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

    //if(zeroChip || zeroChim) cout<<"To few charginos in GenParticle collection!"<<endl;
  }


  void getSelectedJetCollection(){

    jetCollection.clear();

    for(int i=0; i<nPFJet; i++){

      bool isCharginoCandidate = false;

      if(PFJet[i].pt<30.)                          continue;
      if(std::abs(PFJet[i].eta)>4.5)               continue;
      if(PFJet[i].neutralHadronEnergyFraction>0.7) continue;
      if(PFJet[i].chargedEmEnergyFraction>0.5)     continue;

      // Clean collection from canditate tracks (=charginos)
      for(int j=0; j<reco::trackCollection.size(); j++){

	double dPhi = std::abs(reco::trackCollection[j].phi - PFJet[i].phi);
	double dEta = std::abs(reco::trackCollection[j].eta - PFJet[i].eta);
	double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2));

	if(dR<0.5){
	  //cout<<"reco::trackCollection[j].pt = "<<reco::trackCollection[j].pt<<endl;
	  //cout<<"PFJet[i].pt = "<<PFJet[i].pt<<endl;
	  isCharginoCandidate = true;
	  break;
	}
      }

      if(isCharginoCandidate)                     continue;
      jetCollection.push_back(PFJet[i]);
    }
  }

  bool leadingJetRequirementsFullfilled(){

    if(jetCollection.size()==0)                          return false;
    if(jetCollection[0].pt<110.)                         return false;
    if(std::abs(jetCollection[0].eta)>2.4)               return false;
    if(jetCollection[0].chargedHadronEnergyFraction<0.2) return false;
    if(jetCollection[0].neutralEmEnergyFraction>0.7)     return false;
      
    leadingJet=jetCollection[0];
    return true;
  }

  
  bool isGoodVertex(){
    
    if(Vertex[0].z > 24.)                                                return false;
    if(std::sqrt(std::pow(Vertex[0].x,2) + std::pow(Vertex[0].y,2)) > 2) return false;
    if(Vertex[0].ndof < 4)                                               return false;
    return true;
    
  }


  void getCandidateTrackCollection(){

    trackCollection.clear();
    int isolation = 0;
    int all=0;

    for(int i=0; i<nTrack; i++){

      if(Track[i].pt<50.)                                            continue;
      if(std::abs(Track[i].eta)>2.1)                                 continue;
      if(std::abs(Track[i].eta)>1.42 && std::abs(Track[i].eta<1.65)) continue;
      if(std::abs(Track[i].eta)>0.15 && std::abs(Track[i].eta<0.35)) continue;
      if(std::abs(Track[i].eta)>1.55 && std::abs(Track[i].eta<1.85)) continue;
      if(Track[i].d0>0.2)                                            continue;
      if(Track[i].dz>5)                                              continue;
      if(Track[i].numberOfValidHits<5)                               continue;
      if(Track[i].numberOfLostHits>0)                                continue;
      if(Track[i].trackerExpectedHitsInner_numberOfLostHits>0)       continue;
      if(Track[i].trackerExpectedHitsOuter_numberOfLostHits<3)       continue;


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

