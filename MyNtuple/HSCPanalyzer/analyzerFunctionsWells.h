#ifndef ANALYZERFUNCTIONSWELLS_H
#define ANALYZERFUNCTIONSWELLS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
using namespace std;
//-----------------------------------------------------------------------------



struct evt::GenParticle_s chipGenParticle;
struct evt::GenParticle_s chimGenParticle;

std::vector<evt::GenParticle_s> chipmGenParticle;


std::vector<double> etaCSC, phiCSC, etaEcal, phiEcal;  

bool zeroChip = true;
bool zeroChim = true;
bool zeroJet  = true;

int nChi0 = 0;

/*****************************************
  Find chargino in GenParticle collection 
*****************************************/
void findChipmInGenParticleCollection(){

  zeroChip = true;
  zeroChim = true;

  chipmGenParticle.clear();

  for(int i=0; i<evt::nGenParticle; i++){

    if(abs(evt::GenParticle[i].pdgId)==1000024){

      if(evt::GenParticle[i].pdgId>0 && zeroChip){

	chipGenParticle = evt::GenParticle[i];
	chipmGenParticle.push_back(evt::GenParticle[i]);
	zeroChip = false;
      }
      else if(evt::GenParticle[i].pdgId<0 && zeroChim){

	chimGenParticle = evt::GenParticle[i];
	chipmGenParticle.push_back(evt::GenParticle[i]);
	zeroChim = false;
      }	      
    }
    if(!zeroChip && !zeroChim) break;
  }

  //if(zeroChip || zeroChim) cout<<"To few charginos in GenParticle collection!"<<endl;
}

//--------------------------------------------------------------------------------------------------
   
struct evt::GenParticle_s  findLeadingJetInGenParticleCollectionWithPtGt30(){

  struct evt::GenParticle_s leadingJetGenParticle;
  zeroJet = true;
  for(int i=0; i<evt::nGenParticle; i++){

    if(abs(evt::GenParticle[i].pdgId)<=6){
      if(evt::GenParticle[i].pt<30) continue;
      leadingJetGenParticle = evt::GenParticle[i];
      zeroJet = false;
      break;
    }
	       
  }
  return leadingJetGenParticle;
}


//--------------------------------------------------------------------------------------------------

std::vector<evt::Jet_s>  getSelectedJetCollection(){

  std::vector<evt::Jet_s> jetCollection;
  jetCollection.clear();
  for(unsigned int i=0; i<evt::Jet.size(); i++){

    //bool isCharginoCandidate = false;

    if(evt::Jet[i].pt<30.)                          continue;
    if(std::abs(evt::Jet[i].eta)>4.5)               continue;
    //if(Jet[i].neutralHadronEnergyFraction>0.7) continue;
    //if(Jet[i].chargedEmEnergyFraction>0.5)     continue;

    jetCollection.push_back(evt::Jet[i]);
  }

  return jetCollection;
}

//--------------------------------------------------------------------------------------------------

bool areTwoJetsBackToBack(std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    for(unsigned int j=i+1; j<jetColl.size(); j++){

      double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-jetColl[j].phi));
      if(dPhi>=2.5) return true;     
      
    }
  }

  return false;

}

//--------------------------------------------------------------------------------------------------

bool isMetInJetDirection(std::vector<evt::Jet_s> jetColl, struct evt::MET_s *met){

  for(unsigned int i=0; i<jetColl.size(); i++){

    if(i>1) break;
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-met->phi));
    if(dPhi<0.5) return true;     
    
  }

  return false;

}


//--------------------------------------------------------------------------------------------------

bool leadingJetRequirementsFullfilled(struct evt::Jet_s* leadingJet, TH1D* countsEventCuts){

  if(leadingJet==0)                               return false;
  countsEventCuts->Fill("nonEmptyJetColl", 1);
  if(leadingJet->pt<110.)                         return false;
  countsEventCuts->Fill("leadingJetPtGt110GeV", 1);
  if(std::abs(leadingJet->eta)>2.4)               return false;
  countsEventCuts->Fill("absLeadJetEtaLt2p4", 1);
  if(leadingJet->chargedHadronEnergyFraction<0.2) return false;
  countsEventCuts->Fill("CHEFgt0p2", 1);
  if(leadingJet->chargedEmEnergyFraction>0.5)     return false;
  countsEventCuts->Fill("CHEmEFle0p5", 1);
  if(leadingJet->neutralHadronEnergyFraction>0.7) return false;
  countsEventCuts->Fill("NHEFle0p7", 1);
  if(leadingJet->neutralEmEnergyFraction>0.7)     return false;
  countsEventCuts->Fill("NEmEFle0p7", 1);

  return true;
}

//--------------------------------------------------------------------------------------------------
  
bool isGoodVertex(){
    
  if(evt::Vertex[0].z > 24.)                                                     return false;
  if(std::sqrt(std::pow(evt::Vertex[0].x,2) + std::pow(evt::Vertex[0].y,2)) > 2) return false;
  if(evt::Vertex[0].ndof < 4)                                                    return false;
  return true;
}

//--------------------------------------------------------------------------------------------------

bool triggerFired(struct evt::Jet_s* leadingJet){

  if(leadingJet==0)               return false;
  if(leadingJet->pt<80.)          return false;
  if(evt::MET[0].et<100.)         return false;
  return true;
}

//--------------------------------------------------------------------------------------------------

bool isTrackReconstructedJet(struct evt::Track_s track, std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(track.phi-jetColl[i].phi));
    double dEta = std::abs(track.eta - jetColl[i].eta);
    double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2)); 

    if(dR<0.5) return true;
  }

  return false;
}


//--------------------------------------------------------------------------------------------------

double trackEnergyIn0p3Cone(struct evt::Track_s *track, std::vector<evt::Track_s>& trkColl){

  double sumPt = 0;
  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;
  int matchedTracks = 0;
  
  for(unsigned int j=0; j<trkColl.size(); j++){
        
    dPhi = std::abs(TVector2::Phi_mpi_pi(trkColl[j].phi - track->phi));
    dEta = std::abs(trkColl[j].eta - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.0000001){
      matchedTracks += 1;
      continue;

    }

    if(dR<0.3)  sumPt += trkColl[j].pt;
  }

  if(matchedTracks>1) cout<<"More than one match!!!!!!!!!!!!!!"<<endl;
  else if(matchedTracks==0) cout<<"No matched track!!!!!!!!!!!!!"<<endl;

  return sumPt;
}

//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
bool getTrkIsMatchedDeadEcal(struct evt::Track_s *track){

  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;  

  for(unsigned int i=0; i<etaEcal.size(); i++){
        
    dPhi = std::abs(TVector2::Phi_mpi_pi(phiEcal[i] - track->phi));
    
    dEta = std::abs(etaEcal[i] - track->eta);
    
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.05) return true;
  }  
  return false;
}
//--------------------------------------------------------------------------------------------------
bool getTrkIsMatchedBadCSC(struct evt::Track_s *track){

  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;  

  for(unsigned int i=0; i<etaCSC.size(); i++){
        
    dPhi = std::abs(TVector2::Phi_mpi_pi(phiCSC[i] - track->phi));
    
    dEta = std::abs(etaCSC[i] - track->eta);
    
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.25) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------

#endif

