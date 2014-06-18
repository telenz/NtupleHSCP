#ifndef ANALYZERFUNCTIONS_H
#define ANALYZERFUNCTIONS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------



struct GenParticle_s chipGenParticle;
struct GenParticle_s chimGenParticle;

std::vector<GenParticle_s> chipmGenParticle;
  

bool zeroChip = true;
bool zeroChim = true;
bool zeroJet  = true;

int nChi0 = 0;
int nJet=0;

/*****************************************
  Find chargino in GenParticle collection 
*****************************************/
void findChipmInGenParticleCollection(){

  zeroChip = true;
  zeroChim = true;

  chipmGenParticle.clear();

  for(int i=0; i<nGenParticle; i++){

    if(abs(GenParticle[i].pdgId)==1000024){

      if(GenParticle[i].pdgId>0 && zeroChip){

	chipGenParticle = GenParticle[i];
	chipmGenParticle.push_back(GenParticle[i]);
	zeroChip = false;
      }
      else if(GenParticle[i].pdgId<0 && zeroChim){

	chimGenParticle = GenParticle[i];
	chipmGenParticle.push_back(GenParticle[i]);
	zeroChim = false;
      }	      
    }
    if(!zeroChip && !zeroChim) break;
  }

  //if(zeroChip || zeroChim) cout<<"To few charginos in GenParticle collection!"<<endl;
}

//--------------------------------------------------------------------------------------------------
   
struct GenParticle_s  findLeadingJetInGenParticleCollectionWithPtGt30(){

  struct GenParticle_s leadingJetGenParticle;
  zeroJet = true;
  for(int i=0; i<nGenParticle; i++){

    if(abs(GenParticle[i].pdgId)<=6){
      if(GenParticle[i].pt<30) continue;
      leadingJetGenParticle = GenParticle[i];
      nJet+=1;
      zeroJet = false;
      break;
    }
	       
  }
  return leadingJetGenParticle;
}


//--------------------------------------------------------------------------------------------------

std::vector<PFJet_s>  getSelectedJetCollection(){

  std::vector<PFJet_s> jetCollection;
  jetCollection.clear();
  for(int i=0; i<nPFJet; i++){

    //bool isCharginoCandidate = false;

    if(PFJet[i].pt<30.)                          continue;
    if(std::abs(PFJet[i].eta)>4.5)               continue;
    if(PFJet[i].neutralHadronEnergyFraction>0.7) continue;
    if(PFJet[i].chargedEmEnergyFraction>0.5)     continue;

    jetCollection.push_back(PFJet[i]);
  }

  return jetCollection;
}

//--------------------------------------------------------------------------------------------------

bool areTwoJetsBackToBack(std::vector<PFJet_s> jetColl){


  for(unsigned int i=0; i<jetColl.size(); i++){

    for(unsigned int j=i+1; j<jetColl.size(); j++){

      double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-jetColl[j].phi));
      if(dPhi>2.5) return true;     
 
    }
    
  }

  return false;


}


//--------------------------------------------------------------------------------------------------

bool leadingJetRequirementsFullfilled(struct PFJet_s* leadingJet, TH1D* countsEventCuts){

  if(leadingJet==0)                               return false;
  countsEventCuts->Fill("nonEmptyJetColl", 1);
  if(leadingJet->pt<110.)                         return false;
  countsEventCuts->Fill("leadingJetPtGt110GeV", 1);
  if(std::abs(leadingJet->eta)>2.4)               return false;
  countsEventCuts->Fill("absLeadJetEtaLt2p4", 1);
  if(leadingJet->chargedHadronEnergyFraction<0.2) return false;
  countsEventCuts->Fill("CHEFgt0p2", 1);
  if(leadingJet->neutralEmEnergyFraction>0.7)     return false;
  countsEventCuts->Fill("NEmEFle0p7", 1);

  return true;
}

//--------------------------------------------------------------------------------------------------
  
bool isGoodVertex(){
    
  if(Vertex[0].z > 24.)                                                return false;
  if(std::sqrt(std::pow(Vertex[0].x,2) + std::pow(Vertex[0].y,2)) > 2) return false;
  if(Vertex[0].ndof < 4)                                               return false;
  return true;
}

//--------------------------------------------------------------------------------------------------

bool triggerFired(struct PFJet_s* leadingJet){

  if(leadingJet==0)               return false;
  if(leadingJet->pt<80.)           return false;
  if(PFMET[0].et<100.)            return false;
  return true;
}

//--------------------------------------------------------------------------------------------------

bool isTrackReconstructedJet(struct Track_s track){

  for(unsigned int i=0; i<PFJet.size(); i++){
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(track.phi-PFJet[i].phi));
    double dEta = std::abs(track.eta - PFJet[i].eta);
    double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2)); 

    if(dR<0.5) return true;
  }

  return false;
}


//--------------------------------------------------------------------------------------------------

double trackEnergyIn0p3Cone(struct Track_s *track, std::vector<Track_s>& trkColl,unsigned int position){

  double sumPt = 0;
  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;
  for(unsigned int j=0; j<trkColl.size(); j++){
      
    if(j==position) continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(trkColl[j].phi - track->phi));
    dEta = std::abs(trkColl[j].eta - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.3)  sumPt += trkColl[j].pt;
    }

  return sumPt;
}

//--------------------------------------------------------------------------------------------------

#endif

