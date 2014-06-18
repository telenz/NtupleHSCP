#ifndef ANALYZERRECOTRACKS_H
#define ANALYZERRECOTRACKS_H
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "histograms.h"
#include "analyzerRECOClasses.h"
#include <iostream>
#include <vector>
using namespace std;
//-----------------------------------------------------------------------------

std::vector<Track_s> getChipmInTrackCollection(std::vector<Track_s> inputCollection);
std::vector<Track_s> matchTrackToGenParticle(std::vector<Track_s> inputCollection);
std::vector<Track_s> getFakeTracksInTrackCollection(std::vector<Track_s> inputCollection);
std::vector<Track_s> getCandidateTrackCollection(std::vector<Track_s> Track, TH1D* countsTrackCriteria);
std::vector<Track_s> getCandidateTrackCollection_SoftCuts(std::vector<Track_s> Track, TH1D* countsTrackCriteria);

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
  bool preselectionJetCut;
  bool metCut;
  bool leadingJetCut;
  bool vertexCut;
  bool trackCandidateCut;
  bool trackCandidateSoftCut;
  bool onlyChipm;
  bool noChi;
  Hist hist;
  SetOfHistogramsDeDx DeDxHist;
  double mass;
 
  
 public:Event(TString histName, outputFile ofile_):
  hist(histName, ofile_),
    DeDxHist(histName)
    { 
      countsTrackCriteria = new TH1D("countsTrackCriteria","countsTrackCriteria",1,0,1);
      countsTrackCriteria->SetBit(TH1::kCanRebin);
      countsTrackCriteria->SetStats(0);
      countsEventCuts = new TH1D("countsEventCuts","countsEventCuts",1,0,1);
      countsEventCuts->SetBit(TH1::kCanRebin);
      countsEventCuts->SetStats(0);
      
      preselectionJetCut    = true;
      metCut                = false;
      leadingJetCut         = false;
      vertexCut             = true;
      trackCandidateCut     = false;
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
    

    if(onlyChipm)
      {
	Track.clear();
	Track = getChipmInTrackCollection(evt::Track);
      }
    if(noChi)
      {
	Track.clear();
	Track = getFakeTracksInTrackCollection(evt::Track);
      }
        
    if(preselectionJetCut)
      {
	PFJet = getSelectedJetCollection();
      }

    // Vertex cut
    hist.FillVertexHistograms(Vertex);
    if(vertexCut)
      {
	if(!isGoodVertex()) return 0;
      }
    countsEventCuts->Fill("vertexCuts", 1);

    // Leading Jet Cut
    hist.hnPFJet->Fill(PFJet.size());
    if(PFJet.size()==0) hist.h1stjetpt->Fill(PFJet[0].pt);
    if(leadingJetCut)
      {
	if(PFJet.size()==0) return 0;
	if(!leadingJetRequirementsFullfilled(&PFJet[0], countsEventCuts)) return 0;
      }
    countsEventCuts->Fill("1stJetCut", 1);

    // MET cut
    if(metCut)
      {
	if(PFMET.size()==0) return 0;
	if(PFMET[0].et<100) return 0;
      }
    hist.hMet->Fill(PFMET[0].et);
    countsEventCuts->Fill("metCut", 1);
 
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
    Track = matchTrackToGenParticle(Track);
    hist.FillTrackVariables(Track);
    DeDxHist.FillDeDxHistograms(Track);
    
   
    //-----------------------------------------------

    if(chipmGenParticle.size()>0) hist.FillGenParticleHistograms(chipmGenParticle);
    return 0;
    }
    
    //-----------------------------------------------

    
    //-----------------------------------------------
};




std::vector<Track_s> getCandidateTrackCollection(std::vector<Track_s> Track, TH1D* countsTrackCriteria){

  std::vector<Track_s> trackCollection;
  trackCollection.clear();

  for(unsigned int i=0; i<Track.size(); i++){

    countsTrackCriteria->Fill("beforeTrackCriteria", 1);
    if(Track[i].pt<50.)                                                      continue;
    countsTrackCriteria->Fill("PtGreater50GeV", 1);
    if(std::abs(Track[i].eta)>2.1)                                           continue;
    //if(std::abs(Track[i].eta)>1.42 && std::abs(Track[i].eta<1.65)) continue;
    //if(std::abs(Track[i].eta)>0.15 && std::abs(Track[i].eta<0.35)) continue;
    //if(std::abs(Track[i].eta)>1.55 && std::abs(Track[i].eta<1.85)) continue;
    countsTrackCriteria->Fill("EtaLess2p1", 1);
    //if(Track[i].d0>0.2)                                                      continue;
    countsTrackCriteria->Fill("d0Less0p2mm", 1);
    //if(Track[i].dz>5)                                                        continue;
    countsTrackCriteria->Fill("dZLess5mm", 1);
    if(Track[i].numberOfValidHits<5)                                         continue;
    countsTrackCriteria->Fill("NOfValidHitsGe5", 1);
    if(Track[i].numberOfLostHits>0)                                          continue;
    countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", 1);
    if(Track[i].trackerExpectedHitsInner_numberOfLostHits>0)                 continue;
    countsTrackCriteria->Fill("NOfLostHitsInnerEq0", 1);
    
    // Isolation
    double sumPt = 0;
    for(unsigned int j=0; j<Track.size(); j++){
      
      if(i==j) continue;
      double dPhi = std::abs(Track[j].phi - Track[i].phi);
      double dEta = std::abs(Track[j].eta - Track[i].eta);
      double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2));

      if(dR<0.3)  sumPt += Track[j].pt;
    }
      
    if(sumPt/Track[i].pt>0.05)                                     continue;
    countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", 1);
    if(Track[i].dEdxHitsNPHarm2_1000<3)                            continue;
    countsTrackCriteria->Fill("dEdxHarm2Less3", 1);
    if(Track[i].ptError/Track[i].pt<0.25)                          continue;
    countsTrackCriteria->Fill("deltaPtless0p25", 1);
    cout<<"Track[i].trackerExpectedHitsOuter_numberOfLostHits = "<<Track[i].trackerExpectedHitsOuter_numberOfLostHits<<endl;
    if(Track[i].trackerExpectedHitsOuter_numberOfLostHits<3)       continue;
    countsTrackCriteria->Fill("NOfLostHitsOuterGe3", 1);

    trackCollection.push_back(Track[i]);
  }

  return trackCollection;
}


std::vector<Track_s> getCandidateTrackCollection_SoftCuts(std::vector<Track_s> Track, TH1D* countsTrackCriteria){

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
    
    // Isolation
    double sumPt = 0;
    for(unsigned int j=0; j<Track.size(); j++){
      
      if(i==j) continue;
      double dPhi = std::abs(Track[j].phi - Track[i].phi);
      double dEta = std::abs(Track[j].eta - Track[i].eta);
      double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2));

      if(dR<0.3)  sumPt += Track[j].pt;
    }
      
    //if(sumPt/Track[i].pt>0.05)                                     continue;
    countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", 1);
    if(Track[i].dEdxHitsNPHarm2_1000<3)                            continue;
    countsTrackCriteria->Fill("dEdxHarm2Less3", 1);
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

std::vector<Track_s> getChipmInTrackCollection(std::vector<Track_s> inputCollection){

  std::vector<Track_s> chipmTrackCollection;
  chipmTrackCollection.clear();

  for(unsigned int i=0; i<inputCollection.size(); i++){

    if(!zeroChip){
      double dPhichip = std::abs(inputCollection[i].phi - chipGenParticle.phi);
      double dEtachip = std::abs(inputCollection[i].eta - chipGenParticle.eta);
      double dRchip   = std::sqrt(pow(dPhichip,2) + pow(dEtachip,2));
      if(dRchip<0.001){
	chipmTrackCollection.push_back(inputCollection[i]);
	matchedGenChipmToTrack += 1;
      }
    }
    if(!zeroChim){
      double dPhichim = std::abs(inputCollection[i].phi - chimGenParticle.phi);
      double dEtachim = std::abs(inputCollection[i].eta - chimGenParticle.eta);
      double dRchim   = std::sqrt(pow(dPhichim,2) + pow(dEtachim,2));
      if(dRchim<0.001){
	chipmTrackCollection.push_back(inputCollection[i]);
	matchedGenChipmToTrack += 1;
      }
    }
  }
  return chipmTrackCollection;

}


//---------------------------------------------------------------------------------------------------------
// Get only the Chipm from the full Track Collection

std::vector<Track_s> matchTrackToGenParticle(std::vector<Track_s> inputCollection){

  for(unsigned int i=0; i<inputCollection.size(); i++){

    inputCollection[i].pdgId=0;
    inputCollection[i].beta=10;
    double dRsaved = 100.0;

    for(unsigned int j=0; j<GenParticle.size(); j++){
      
      double dPhi = std::abs(inputCollection[i].phi - GenParticle[j].phi);
      double dEta = std::abs(inputCollection[i].eta - GenParticle[j].eta);
      double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2));
      if(dR<0.5){
	if(dRsaved>dR){
	  inputCollection[i].pdgId=GenParticle[j].pdgId;
	  inputCollection[i].beta=GenParticle[j].p/GenParticle[j].energy;
	  dRsaved = dR;
	}
      }
    }
    //cout<<"dRsaved = "<<dRsaved<<endl;
  }

  return inputCollection;
}
//---------------------------------------------------------------------------------------------------------
std::vector<Track_s> getFakeTracksInTrackCollection(std::vector<Track_s> inputCollection){

  std::vector<Track_s> outputCollection;
  outputCollection.clear();
  double dRchip = 0;
  double dRchim = 0;

  for(unsigned int i=0; i<inputCollection.size(); i++){

    dRchip=10000;
    dRchim=10000;

    if(!zeroChip){
      double dPhichip = std::abs(inputCollection[i].phi - chipGenParticle.phi);
      double dEtachip = std::abs(inputCollection[i].eta - chipGenParticle.eta);
      dRchip   = std::sqrt(pow(dPhichip,2) + pow(dEtachip,2));
    }
    if(!zeroChim){
      double dPhichim = std::abs(inputCollection[i].phi - chimGenParticle.phi);
      double dEtachim = std::abs(inputCollection[i].eta - chimGenParticle.eta);
      dRchim   = std::sqrt(pow(dPhichim,2) + pow(dEtachim,2));
    }
    if(dRchim>0.01 && dRchip>0.01) outputCollection.push_back(inputCollection[i]);
  }
  return outputCollection;

}



#endif

