#ifndef ANALYZERRECOCLASSES_H
#define ANALYZERRECOCLASSES_H
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "TVector3.h"
#include <iostream>
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------

class SetOfDeDxHistograms
{

 public:
  TH1D* hDeDx; 
  TH1D* hDeDx7Hits;  
  TH1D* hDeDx3Hits;  

  TH2D* hPDeDx;
  TH2D* hPDeDx7Hits;
  TH2D* hPDeDx3Hits;

  TH2D* hNValidDeDx;
  TH2D* hNValidDeDx7Hits;
  TH2D* hNValidDeDx3Hits;
  
  
  public:SetOfDeDxHistograms(TString histName)
  {
    hDeDx      = new TH1D("hDeDx" + histName , "hDeDx2" + histName ,400,0,40);  
    hDeDx7Hits = new TH1D("hDeDx7Hits" + histName , "hDeDx7Hits" + histName ,400,0,40);  
    hDeDx3Hits = new TH1D("hDeDx3Hits" + histName , "hDeDx3Hits" + histName ,400,0,40);  

    hPDeDx       = new TH2D("hPDeDx" + histName , "hnPDeDx" + histName ,150,0,1500,400,0,40);
    hPDeDx7Hits  = new TH2D("hPDeDx7Hits" + histName , "hPDeDx7Hits" + histName ,150,0,1500,400,0,40);
    hPDeDx3Hits  = new TH2D("hPDeDx3Hits" + histName , "hPDeDx3Hits" + histName ,150,0,1500,400,0,40);
    
    hNValidDeDx       = new TH2D("hNValidDeDx" + histName , "hnNValidDeDx" + histName ,50,0,50,400,0,40);
    hNValidDeDx7Hits  = new TH2D("hNValidDeDx7Hits" + histName , "hNValidDeDx7Hits" + histName ,50,0,50,400,0,40);
    hNValidDeDx3Hits  = new TH2D("hNValidDeDx3Hits" + histName , "hNValidDeDx3Hits" + histName ,50,0,50,400,0,40);
   
    }
};

class SetOfHistogramsDeDx
{

 public:
  class SetOfDeDxHistograms harm2;
  class SetOfDeDxHistograms trun40;
  class SetOfDeDxHistograms median;

 public:SetOfHistogramsDeDx(TString histName): 
    harm2("Harm2"),
    trun40("Trun40"),
    median("Median")
      {};

  void FillDeDxHistograms(std::vector<Track_s> trackCollection)
  {

    for(unsigned int i=0; i<trackCollection.size(); i++){
      TVector3 trackVec;
      trackVec.SetPtEtaPhi(trackCollection[i].pt,trackCollection[i].eta,trackCollection[i].phi);
      double p = std::sqrt(std::pow(trackCollection[i].pt,2) + std::pow(trackVec.Pz(),2));

      harm2.hDeDx        -> Fill(trackCollection[i].dEdxHitsNPHarm2_1000);
      harm2.hDeDx7Hits   -> Fill(trackCollection[i].dEdxHitsNPHarm2_7);
      harm2.hDeDx3Hits   -> Fill(trackCollection[i].dEdxHitsNPHarm2_3);
      trun40.hDeDx       -> Fill(trackCollection[i].dEdxHitsNPTrun40_1000);
      trun40.hDeDx7Hits  -> Fill(trackCollection[i].dEdxHitsNPTrun40_7);
      trun40.hDeDx3Hits  -> Fill(trackCollection[i].dEdxHitsNPTrun40_3);
      median.hDeDx       -> Fill(trackCollection[i].dEdxHitsNPMedian_1000);
      median.hDeDx7Hits  -> Fill(trackCollection[i].dEdxHitsNPMedian_7);
      median.hDeDx3Hits  -> Fill(trackCollection[i].dEdxHitsNPMedian_3);

      harm2.hPDeDx       -> Fill(p, trackCollection[i].dEdxHitsNPHarm2_1000);
      harm2.hPDeDx7Hits  -> Fill(p, trackCollection[i].dEdxHitsNPHarm2_7);
      harm2.hPDeDx3Hits  -> Fill(p, trackCollection[i].dEdxHitsNPHarm2_3);
      trun40.hPDeDx      -> Fill(p, trackCollection[i].dEdxHitsNPTrun40_1000);
      trun40.hPDeDx7Hits -> Fill(p, trackCollection[i].dEdxHitsNPTrun40_7);
      trun40.hPDeDx3Hits -> Fill(p, trackCollection[i].dEdxHitsNPTrun40_3);
      median.hPDeDx      -> Fill(p, trackCollection[i].dEdxHitsNPMedian_1000);
      median.hPDeDx7Hits -> Fill(p, trackCollection[i].dEdxHitsNPMedian_7);
      median.hPDeDx3Hits -> Fill(p, trackCollection[i].dEdxHitsNPMedian_3);

      harm2.hNValidDeDx       -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPHarm2_1000);
      harm2.hNValidDeDx7Hits  -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPHarm2_7);
      harm2.hNValidDeDx3Hits  -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPHarm2_3);
      trun40.hNValidDeDx      -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPTrun40_1000);
      trun40.hNValidDeDx7Hits -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPTrun40_7);
      trun40.hNValidDeDx3Hits -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPTrun40_3);
      median.hNValidDeDx      -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPMedian_1000);
      median.hNValidDeDx7Hits -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPMedian_7);
      median.hNValidDeDx3Hits -> Fill(trackCollection[i].numberOfValidHits, trackCollection[i].dEdxHitsNPMedian_3);
    }
  }

};

#endif

