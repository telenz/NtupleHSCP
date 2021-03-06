#ifndef ANALYZERSELECTIONAODSIM_H
#define ANALYZERSELECTIONAODSIM_H
//-----------------------------------------------------------------------------
#include "analyzerAODSIM.h"
#include "analyzerHistograms.h"
#include "analyzerFunctions.h"
#include <iostream>
#include <vector>
using namespace std;
//-----------------------------------------------------------------------------

std::vector<Track_s> getChipmInTrackCollection(std::vector<Track_s>& inputCollection);
std::vector<Track_s> getFakeTracksInTrackCollection(const std::vector<Track_s>& inputCollection);
std::vector<Track_s> getCandidateTrackCollection(std::vector<Track_s>& Track, TH1D* countsTrackCriteria);
std::vector<Track_s> getCandidateTrackCollection_SoftCuts(std::vector<Track_s>& Track, TH1D* countsTrackCriteria);
void matchTrackToGenParticle(std::vector<Track_s>& inputCollection);

int mismatchedGenChipmToTrack = 0;
int matchedGenChipmToTrack    = 0;

class Event
{

 public:
  std::vector<Track_s> Track;
  std::vector<PFJet_s> PFJet;
  std::vector<PFMET_s> PFMET;
  std::vector<Vertex_s> Vertex;
  struct GenParticle_s leadJetGenParticle;
  TH1D *countsTrackCriteria;
  TH1D *countsEventCuts;
  bool triggerCut;
  bool preselectionJetCut;
  bool metCut;
  bool leadingJetCut;
  bool deltaPhiCut;
  bool vertexCut;
  bool trackCandidateCut;
  bool trackCandidateSoftCut;
  bool onlyChipm;
  bool noChi;
  Hist hist;
  double mass;
 
  
 public:Event(TString histName, outputFile ofile_):
  hist(histName, ofile_)
    { 

      countsTrackCriteria = new TH1D("countsTrackCriteria","countsTrackCriteria",1,0,1);
      countsTrackCriteria->SetBit(TH1::kCanRebin);
      countsTrackCriteria->SetStats(0);
      countsTrackCriteria->Fill("beforeTrackCriteria", 0);
      countsTrackCriteria->Fill("PtGreater50GeV", 0);
      countsTrackCriteria->Fill("EtaLess2p1", 0);
      countsTrackCriteria->Fill("EtaLess1p42Gt1p65", 0);
      countsTrackCriteria->Fill("EtaLess0p15Gt0p35", 0);
      countsTrackCriteria->Fill("EtaLess1p55Gt1p85", 0);
      countsTrackCriteria->Fill("d0Less0p2mm", 0);
      countsTrackCriteria->Fill("dZLess5mm", 0);
      countsTrackCriteria->Fill("NOfValidHitsGe7", 0);
      countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", 0);
      countsTrackCriteria->Fill("NOfLostHitsInnerEq0", 0);
      countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", 0);
      countsTrackCriteria->Fill("InJetCollectionR0p5", 0);
      countsTrackCriteria->Fill("InLeptonCollectionR0p15", 0);
      countsTrackCriteria->Fill("deltaPtless0p25", 0);
      countsTrackCriteria->Fill("CaloIsolation0p5", 0);
      countsTrackCriteria->Fill("NOfLostHitsOuterGe3", 0);
    
      
      countsEventCuts = new TH1D("countsEventCuts","countsEventCuts",1,0,1);
      countsEventCuts->SetBit(TH1::kCanRebin);
      countsEventCuts->SetStats(0);
      
      triggerCut            = true;
      preselectionJetCut    = true;
      metCut                = true;
      leadingJetCut         = true;
      deltaPhiCut           = true;
      vertexCut             = false;
      trackCandidateCut     = true;
      trackCandidateSoftCut = false;
      onlyChipm             = false;
      noChi                 = false;
    }

  int Selection()
  {
    
    countsEventCuts->Fill("noCuts", 1);
    Track.clear();
    PFJet.clear();
    PFMET.clear();
    Vertex.clear();
    Track=evt::Track;
    PFJet=evt::PFJet;
    PFMET=evt::PFMET;
    Vertex=evt::Vertex;
    
    // Trigger Cut
    if(triggerCut)
      {
	if(edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5 == 0) return 0;
      }
    countsEventCuts->Fill("triggerCut", 1);

    // Vertex cut
    hist.FillVertexHistograms(Vertex);
    if(vertexCut)
      {
	if(!isGoodVertex()) return 0;
      }
    countsEventCuts->Fill("vertexCuts", 1);

    // MET cut
    if(metCut)
      {
    	if(PFMET.size()==0) return 0;
	bool oneMetGt100 = false;
	for(int i=0; i<PFMET.size(); i++){
	  if(PFMET[i].pt>100){
	    oneMetGt100=true;
	  }
	}
	hist.hMet->Fill(PFMET[0].pt);
	if(!oneMetGt100) return 0;
      }
    countsEventCuts->Fill("metCut", 1);

    // Fill new jet collection
    if(preselectionJetCut)
      {
	PFJet = getSelectedJetCollection();
      }

    // Leading Jet Cut
    hist.hnPFJet->Fill(PFJet.size());
    if(PFJet.size()!=0) hist.h1stjetpt->Fill(PFJet[0].pt);
    if(leadingJetCut)
      {
	if(PFJet.size()==0) return 0;
	if(!leadingJetRequirementsFullfilled(&PFJet[0], countsEventCuts)) return 0;
      }
    countsEventCuts->Fill("1stJetCuts", 1);

    // DeltaPhi Cut
    if(deltaPhiCut)
      {
    	for(unsigned int i=0; i<PFJet.size(); i++){
    
	  for(unsigned int j=i+1; j<PFJet.size(); j++){
    
	    double dPhi = std::abs(TVector2::Phi_mpi_pi(PFJet[i].phi-PFJet[j].phi));
	    hist.hDeltaPhi->Fill(dPhi);     
	  }
	  
	}
	if(areTwoJetsBackToBack(PFJet)) return 0;
      }
    countsEventCuts->Fill("DeltaPhiCut", 1);

    if(onlyChipm)
      {
    	Track.clear();
    	Track = getChipmInTrackCollection(evt::Track);
      }
    else if(noChi)
      {
    	Track.clear();
    	Track = getFakeTracksInTrackCollection(evt::Track);
      }
        
    // Track Candidate Cut
    if(trackCandidateCut)
      {
	Track = getCandidateTrackCollection(Track,countsTrackCriteria);
	if(Track.size()==0) return 0;
      }
    
    if(trackCandidateSoftCut)
      {
    	Track = getCandidateTrackCollection_SoftCuts(Track,countsTrackCriteria);
	if(Track.size()==0) return 0;
      }
       
    matchTrackToGenParticle(Track);
    hist.FillTrackVariables_AODSIM(Track);
    
   
    //-----------------------------------------------

    if(chipmGenParticle.size()>0) hist.FillGenParticleHistograms(chipmGenParticle);
    
    return 0;

    }
    
    //-----------------------------------------------

    
    //-----------------------------------------------
};


//--------------------------------------------------------------------------------------------------

std::vector<Track_s> getCandidateTrackCollection(std::vector<Track_s>& Track, TH1D* countsTrackCriteria){

  std::vector<Track_s> trackCollection;
  trackCollection.clear();

  for(unsigned int i=0; i<Track.size(); i++){

    countsTrackCriteria->Fill("beforeTrackCriteria", 1);

    if(Track[i].pt<50.)                                                      continue;
    countsTrackCriteria->Fill("PtGreater50GeV", 1);
    //.................................................................................//
    if(std::abs(Track[i].eta)>2.1)                                           continue;
    countsTrackCriteria->Fill("EtaLess2p1", 1);
    if(std::abs(Track[i].eta)>1.42 && std::abs(Track[i].eta<1.65))           continue;
    countsTrackCriteria->Fill("EtaLess1p42Gt1p65", 1);
    if(std::abs(Track[i].eta)>0.15 && std::abs(Track[i].eta<0.35))           continue;
    countsTrackCriteria->Fill("EtaLess0p15Gt0p35", 1);
    if(std::abs(Track[i].eta)>1.55 && std::abs(Track[i].eta<1.85))           continue;
    countsTrackCriteria->Fill("EtaLess1p55Gt1p85", 1);
    //.................................................................................//
    if(Track[i].d0>0.2)                                                      continue;
    countsTrackCriteria->Fill("d0Less0p2mm", 1);
    //.................................................................................//
    if(Track[i].dz>5)                                                        continue;
    countsTrackCriteria->Fill("dZLess5mm", 1);
    //.................................................................................//
    if(Track[i].numberOfValidHits<7)                                         continue;
    countsTrackCriteria->Fill("NOfValidHitsGe7", 1);
    //.................................................................................//
    if(Track[i].numberOfLostHits>0)                                          continue;
    countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", 1);
    //.................................................................................//
    if(Track[i].trackerExpectedHitsInner_numberOfLostHits>0)                 continue;
    countsTrackCriteria->Fill("NOfLostHitsInnerEq0", 1);
    //.................................................................................//
    double sumPt = trackEnergyIn0p3Cone(&Track[i],Track, i);
    if(sumPt/Track[i].pt>0.05)                                     continue;
    countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", 1);
    //.................................................................................//
    if(isTrackReconstructedJet(Track[i])) continue;
    countsTrackCriteria->Fill("InJetCollectionR0p5", 1);
    //.................................................................................//
    if(isTrackReconstructedLepton(&Track[i])) continue;
    countsTrackCriteria->Fill("InLeptonCollectionR0p15", 1);
    //.................................................................................//
    if(Track[i].ptError/Track[i].pt<0.25)                          continue;
    countsTrackCriteria->Fill("deltaPtless0p25", 1);
    //.................................................................................//
    if(!isTrackCaloIsolated(&Track[i])) continue;
    countsTrackCriteria->Fill("CaloIsolation0p5", 1);
    //.................................................................................//
    if(Track[i].trackerExpectedHitsOuter_numberOfLostHits<3)       continue;
    countsTrackCriteria->Fill("NOfLostHitsOuterGe3", 1);
    //.................................................................................//
    
    trackCollection.push_back(Track[i]);
  }

  return trackCollection;
}

//--------------------------------------------------------------------------------------------------
std::vector<Track_s> getCandidateTrackCollection_SoftCuts(std::vector<Track_s>& Track, TH1D* countsTrackCriteria){

  std::vector<Track_s> trackCollection;
  trackCollection.clear();

  for(unsigned int i=0; i<Track.size(); i++){

    countsTrackCriteria->Fill("beforeTrackCriteria", 1);
    if(Track[i].pt<10.)                                                      continue;
    countsTrackCriteria->Fill("PtGreater10GeV", 1);
    countsTrackCriteria->Fill("EtaLess2p1", 1);
    //if(Track[i].d0>0.2)                                                      continue;
    countsTrackCriteria->Fill("d0Less0p2mm", 1);
    //if(Track[i].dz>5)                                                        continue;
    countsTrackCriteria->Fill("dZLess5mm", 1);
    //if(Track[i].numberOfValidHits<5)                                         continue;
    countsTrackCriteria->Fill("NOfValidHitsGe5", 1);
    if(Track[i].numberOfLostHits>0)                                          continue;
    countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", 1);
    if(Track[i].trackerExpectedHitsInner_numberOfLostHits>0)                 continue;
    countsTrackCriteria->Fill("NOfLostHitsInnerEq0", 1);
    //if(sumPt/Track[i].pt>0.05)                                     continue;
    countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", 1);
    //if(Track[i].dEdxHitsNPHarm2_1000<3)                            continue;
    //countsTrackCriteria->Fill("dEdxHarm2Less3", 1);
    //if(Track[i].ptError/Track[i].pt<0.25)                            continue;
    countsTrackCriteria->Fill("deltaPtless0p25", 1);
    //if(Track[i].trackerExpectedHitsOuter_numberOfLostHits<3)                 continue;
    countsTrackCriteria->Fill("NOfLostHitsOuterGe3", 1);

    trackCollection.push_back(Track[i]);
  }

  return trackCollection;
}


//---------------------------------------------------------------------------------------------------------
// Get only the Chipm from the full Track Collection

std::vector<Track_s> getChipmInTrackCollection(std::vector<Track_s>& inputCollection){

  std::vector<Track_s> chipmTrackCollection;
  chipmTrackCollection.clear();


  double dPhichip = 0;
  double dEtachip = 0;
  double dRchip   = 0;
  double dPhichim = 0;
  double dEtachim = 0;
  double dRchim   = 0;


  for(unsigned int i=0; i<inputCollection.size(); i++){

    if(!zeroChip){
      dPhichip = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chipGenParticle.phi));
      dEtachip = std::abs(inputCollection[i].eta - chipGenParticle.eta);
      dRchip   = std::sqrt( dPhichip*dPhichip + dEtachip*dEtachip );
      if(dRchip<0.001){
	chipmTrackCollection.push_back(inputCollection[i]);
	matchedGenChipmToTrack += 1;
      }
    }
    if(!zeroChim){
      dPhichim = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chimGenParticle.phi));
      dEtachim = std::abs(inputCollection[i].eta - chimGenParticle.eta);
      dRchim   = std::sqrt(dPhichim*dPhichim + dEtachim*dEtachim );
      if(dRchim<0.001){
	chipmTrackCollection.push_back(inputCollection[i]);
	matchedGenChipmToTrack += 1;
      }
    }
  }
  return chipmTrackCollection;

}


//---------------------------------------------------------------------------------------------------------
// Match Tracks to a generator particle in the GenCollection

void matchTrackToGenParticle(std::vector<Track_s>& inputCollection){

  double dPhi = 0;
  double dEta = 0;
  double dR   = 0;

  for(unsigned int i=0; i<inputCollection.size(); i++){

    inputCollection[i].pdgId=0;
    inputCollection[i].beta=10;
    double dRsaved = 0.5;

    for(unsigned int j=0; j<GenParticle.size(); j++){

      dEta = std::abs(inputCollection[i].eta - GenParticle[j].eta);
      if(dEta>dRsaved) continue;
      dPhi = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - GenParticle[j].phi));
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );
      if(dR<dRsaved){
	  inputCollection[i].pdgId=GenParticle[j].pdgId;
	  inputCollection[i].beta=GenParticle[j].p/GenParticle[j].energy;
	  dRsaved = dR;
      }
    }
  }
  
}
//---------------------------------------------------------------------------------------------------------
std::vector<Track_s> getFakeTracksInTrackCollection(const std::vector<Track_s>& inputCollection){

  std::vector<Track_s> outputCollection;
  outputCollection.clear();
  double dRchip = 0;
  double dRchim = 0;

  for(unsigned int i=0; i<inputCollection.size(); i++){

    dRchip=10000;
    dRchim=10000;

    if(!zeroChip){
      double dPhichip = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chipGenParticle.phi));
      double dEtachip = std::abs(inputCollection[i].eta - chipGenParticle.eta);
      dRchip   = std::sqrt(pow(dPhichip,2) + pow(dEtachip,2));
    }
    if(!zeroChim){
      double dPhichim = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chimGenParticle.phi));
      double dEtachim = std::abs(inputCollection[i].eta - chimGenParticle.eta);
      dRchim   = std::sqrt(pow(dPhichim,2) + pow(dEtachim,2));
    }
    if(dRchim>0.01 && dRchip>0.01) outputCollection.push_back(inputCollection[i]);
  }
  return outputCollection;

}



#endif

