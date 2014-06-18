#ifndef ANALYZERFUNCTIONSAODSIM_H
#define ANALYZERFUNCTIONSAODSIM_H
//-----------------------------------------------------------------------------
//#include "analyzer_AODSIM.h"
#include <iostream>
#include "TVector2.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------


//--------------------------------------------------------------------------------------------------

bool isTrackCaloIsolated(struct Track_s *track){

  double energySum = 0;
  double dPhi = 0;
  double dEta = 0;
  double dR   = 0;
  
  for(unsigned int i=0; i<CaloCluster1.size(); i++){

    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-CaloCluster1[i].phi));
    dEta = std::abs(track->eta - CaloCluster1[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if( dR<0.5 ) energySum += CaloCluster[i].energy;   
  }

  if(energySum<10) return true;

  return false;

}

//--------------------------------------------------------------------------------------------------

bool isTrackReconstructedLepton(struct Track_s* track){

  double dPhi = 0;
  double dEta = 0;
  double dR   = 0; 

  for(unsigned int i=0; i<GsfElectron.size(); i++){

    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-GsfElectron[i].phi));
    dEta = std::abs(track->eta - GsfElectron[i].eta);
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

