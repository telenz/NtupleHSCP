#ifndef ANALYZERFUNCTIONSWELLS_H
#define ANALYZERFUNCTIONSWELLS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
#include "analyzerWells.h"
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

    //if(evt::Jet[i].pt<30.)                          continue;
    //if(std::abs(evt::Jet[i].eta)>4.5)               continue;
    //if(Jet[i].neutralHadronEnergyFraction>0.7) continue;
    //if(Jet[i].chargedEmEnergyFraction>0.5)     continue;

    jetCollection.push_back(evt::Jet[i]);
  }

  return jetCollection;
}

//--------------------------------------------------------------------------------------------------

std::vector<evt::Jet_s>  getSubleadingJetCollection(){

  std::vector<evt::Jet_s> jetCollection;
  jetCollection.clear();
  for(unsigned int i=0; i<evt::Jet.size(); i++){

    //bool isCharginoCandidate = false;

    if(evt::Jet[i].pt<=30.)                          continue;
    if(std::abs(evt::Jet[i].eta)>=4.5)               continue;
    if(evt::Jet[i].neutralHadronEnergyFraction>=0.7) continue;
    if(evt::Jet[i].chargedEmEnergyFraction>=0.5)     continue;

    jetCollection.push_back(evt::Jet[i]);
  }

  return jetCollection;
}

//--------------------------------------------------------------------------------------------------

bool areTwoJetsBackToBack(std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    for(unsigned int j=i+1; j<jetColl.size(); j++){

      double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-jetColl[j].phi));
 
      if(abs(dPhi)>=2.5) return true;     
      
    }
  }

  return false;

}

//--------------------------------------------------------------------------------------------------

bool isMetInJetDirection(std::vector<evt::Jet_s> jetColl, struct evt::MET_s *met){

  for(unsigned int i=0; i<jetColl.size(); i++){

    if(i>1) break;
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-met->phi));
    if(dPhi<=0.5) return true;     
    
  }

  return false;

}


//--------------------------------------------------------------------------------------------------

bool leadingJetRequirementsFullfilled(struct evt::Jet_s* leadingJet, TH1D* countsEventCuts){

  if(leadingJet==0)                               return false;
  //countsEventCuts->Fill("nonEmptyJetColl", evt::weight);
  if(leadingJet->pt<=110.)                         return false;
  countsEventCuts->Fill("leadingJetPtGt110GeV", evt::weight);
  if(std::abs(leadingJet->eta)>=2.4)               return false;
  countsEventCuts->Fill("absLeadJetEtaLt2p4", evt::weight);
  if(leadingJet->chargedHadronEnergyFraction<=0.2) return false;
  countsEventCuts->Fill("CHEFgt0p2", evt::weight);
  if(leadingJet->chargedEmEnergyFraction>=0.5)     return false;
  countsEventCuts->Fill("CHEmEFle0p5", evt::weight);
  if(leadingJet->neutralHadronEnergyFraction>=0.7) return false;
  countsEventCuts->Fill("NHEFle0p7", evt::weight);
  if(leadingJet->neutralEmEnergyFraction>=0.7)     return false;
  countsEventCuts->Fill("NEmEFle0p7", evt::weight);

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
 
  for(unsigned int j=0; j<trkColl.size(); j++){
        
    dPhi = std::abs(TVector2::Phi_mpi_pi(trkColl[j].phi - track->phi));
    dEta = std::abs(trkColl[j].eta - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);
    
    if(dR<0.3)  sumPt += trkColl[j].pt;
  }

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
bool isWithinIntermoduleGapsOfECAL(struct evt::Track_s *track){

  if(track->eta<-1.14018   && track->eta>-1.1439)     return true;
  if(track->eta<-0.791884  && track->eta>-0.796051)   return true;
  if(track->eta<-0.44356   && track->eta>-0.447911)   return true;
  if(track->eta<0.00238527 && track->eta>-0.00330793) return true;
  if(track->eta<0.446183   && track->eta>0.441949)    return true;
  if(track->eta<0.793955   && track->eta>0.789963)    return true;
  if(track->eta<1.14164    && track->eta>1.13812)     return true;
  
  return false;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCuts(std::vector<evt::Track_s> inputColl, bool keepCriteria)
{

  std::vector<evt::Track_s> outputColl;

  for(unsigned int i=0; i<inputColl.size(); i++){
    if(!keepCriteria) continue;
    outputColl.push_back(inputColl[i]);
  }

  return outputColl;

}
//--------------------------------------------------------------------------------------------------



#endif

