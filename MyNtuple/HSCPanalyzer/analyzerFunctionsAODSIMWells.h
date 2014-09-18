#ifndef ANALYZERFUNCTIONSAODSIMWELLS_H
#define ANALYZERFUNCTIONSAODSIMWELLS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------

bool isTrackCaloIsolated(struct Track_s *track){


  double caloTowerIso05RhoCorr = 
    std::max(0.,track->caloHadDeltaRp5 + track->caloEMDeltaRp5 - sdouble_value*TMath::Pi()*0.5*0.5); 
  //std::max(0.,track->caloHadDeltaRp5 + track->caloEMDeltaRp5 + sdouble_value*TMath::Pi()*0.5*0.5); 
 
  if(caloTowerIso05RhoCorr<10) return true;
  else{

    //cout<<"track->caloHadDeltaRp5 = "<<track->caloHadDeltaRp5<<endl;
    //cout<<"track->caloEMDeltaRp5 = "<<track->caloEMDeltaRp5<<endl;
    //cout<<"sdouble_value*TMath::Pi()*0.5*0.5 = "<<sdouble_value*TMath::Pi()*0.5*0.5<<endl;
  }
  return false;

}

//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedTau(struct Track_s* track){

  double dPhi = 0;
  double dEta = 0;
  double dR   = 0; 


  for(unsigned int i=0; i<Tau.size(); i++){

    if(Tau[i].pt<30)                                    continue;
    if(abs(Tau[i].eta)>2.3)                             continue;
    if(Tau[i].byLooseCombinedIsolationDeltaBetaCorr<=0) continue;
    if(Tau[i].decayModeFinding<=0)                      continue;
    if(Tau[i].againstElectronLoose<=0)                  continue;
    if(Tau[i].againstMuonTight<=0)                      continue;

    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Tau[i].phi));
    dEta = std::abs(track->eta - Tau[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedElectron(struct Track_s* track){

  double dPhi = 0;
  double dEta = 0;
  double dR   = 0; 

  for(unsigned int i=0; i<Electron.size(); i++){

    if(Electron[i].mvaNonTrigV0<0) continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Electron[i].phi));
    dEta = std::abs(track->eta - Electron[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedMuon(struct Track_s* track){

  double dPhi = 0;
  double dEta = 0;
  double dR   = 0; 

  for(unsigned int i=0; i<Muon.size(); i++){
    
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Muon[i].phi));
    dEta = std::abs(track->eta - Muon[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedLepton(struct Track_s* track){

  double dPhi = 0;
  double dEta = 0;
  double dR   = 0; 

  for(unsigned int i=0; i<Electron.size(); i++){

    if(Electron[i].mvaNonTrigV0<0) continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Electron[i].phi));
    dEta = std::abs(track->eta - Electron[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  for(unsigned int i=0; i<Muon.size(); i++){
    
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Muon[i].phi));
    dEta = std::abs(track->eta - Muon[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------

#endif

