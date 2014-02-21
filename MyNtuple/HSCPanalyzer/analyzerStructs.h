#ifndef ANALYZERSTRUCTS_H
#define ANALYZERSTRUCTS_H
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "TVector3.h"
#include <iostream>
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------

struct SetOfTrackVariables
{

  SetOfTrackVariables()
  {

    double dedx = 0;
    double dedx7Hits = 0;
    double dedx5Hits = 0;
    double dedx3Hits = 0;

  }
};


class SetOfVertexHistograms
{

 public:
  TH1D* hVertexNdof;
  TH1D* hVertexZ;
  TH1D* hVertexRho;


 public:SetOfVertexHistograms()
  {
    hVertexNdof = new TH1D("hVertexNdof", "hVertexNdof",500,0,500);  
    hVertexZ    = new TH1D("hVertexZ", "hVertexZ",160,-40,40);  
    hVertexRho  = new TH1D("hVertexRho", "hVertexRho",20,0,10);
 
  };
    void FillVertexHistograms(int i)
  {
    hVertexNdof -> Fill(Vertex[i].ndof);
    hVertexZ    -> Fill(Vertex[i].z);
    hVertexRho  -> Fill(std::sqrt(std::pow(Vertex[i].x,2) + std::pow(Vertex[i].y,2)));
  }


};

class SetOfHistograms
{

 public:
  TH1D* hDeDx; 
  TH1D* hDeDx7Hits;  
  TH1D* hDeDx3Hits;  

  TH2D* hPDeDx;
  TH2D* hPDeDx7Hits;
  TH2D* hPDeDx3Hits;
  
  
  public:SetOfHistograms(TString histName)
  {
    hDeDx      = new TH1D("hDeDx" + histName, "hDeDx2" + histName,100,0,100);  
    hDeDx7Hits = new TH1D("hDeDx7Hits" + histName, "hDeDx7Hits" + histName,100,0,100);  
    hDeDx3Hits = new TH1D("hDeDx3Hits" + histName, "hDeDx3Hits" + histName,100,0,100);  

    hPDeDx       = new TH2D("hPDeDx" + histName, "hnPDeDx" + histName,150,0,1500,400,0,40);
    hPDeDx7Hits  = new TH2D("hPDeDx7Hits" + histName, "hPDeDx7Hits" + histName,150,0,1500,400,0,40);
    hPDeDx3Hits  = new TH2D("hPDeDx3Hits" + histName, "hPDeDx3Hits" + histName,150,0,1500,400,0,40);
   
    }
};

class SetOfHistogramsChipm
{

 public:
  TH1D *htrackpt;
  TH1D *htracketa;
  TH1D *htrackd0;
  TH1D *htrackNOfValidHits; 
  TH1D *htrackNOfLostHits; 
  TH1D *htrackdz;
  TH1D *hsumPtDeltaR0p3;
  class SetOfHistograms harm2;
  class SetOfHistograms trun40;
  class SetOfHistograms median;

 public:SetOfHistogramsChipm(TString histName): 
    harm2("Harm2"+histName),
    trun40("Trun40"+histName),
    median("Median"+histName)
      {
	htrackpt   = new TH1D("htrackpt" + histName,"htrackpt" + histName,150,0,1500);
	htracketa  = new TH1D("htracketa" + histName,"htracketa" + histName,100,-5,5);
	htrackd0   = new TH1D("htrackd0" + histName,"htrackd0" + histName,100,-2,2);
	htrackdz   = new TH1D("htrackdz" + histName,"htrackdz" + histName,100,-50,50);
	htrackNOfValidHits = new TH1D("htrackNOfValidHits" + histName,"htrackNOfValidHits" + histName,50,0,50);
	htrackNOfLostHits = new TH1D("htrackNOfLostHits" + histName,"htrackNOfLostHits" + histName,50,0,50);
	hsumPtDeltaR0p3 = new TH1D("hsumPtDeltaR0p3" + histName,"hsumPtDeltaR0p3" + histName,100,0,0.5);
	
      };

  void FillTrackHistograms(int i)
    {
      htrackpt->Fill(reco::trackCollection[i].pt);
      htracketa->Fill(reco::trackCollection[i].eta);
      htrackd0->Fill(reco::trackCollection[i].d0);
      htrackdz->Fill(reco::trackCollection[i].dz);
      htrackNOfValidHits->Fill(reco::trackCollection[i].numberOfValidHits);
      htrackNOfLostHits->Fill(reco::trackCollection[i].numberOfLostHits);
      
    }

  void FillDeDxHistograms(int i)
    {
      TVector3 trackVec;
      trackVec.SetPtEtaPhi(reco::trackCollection[i].pt,reco::trackCollection[i].eta,reco::trackCollection[i].phi);
      double p = std::sqrt(std::pow(reco::trackCollection[i].pt,2) + std::pow(trackVec.Pz(),2));

      harm2.hDeDx        -> Fill(reco::trackCollection[i].dEdxHitsNPHarm2_1000);
      harm2.hDeDx7Hits   -> Fill(reco::trackCollection[i].dEdxHitsNPHarm2_7);
      harm2.hDeDx3Hits   -> Fill(reco::trackCollection[i].dEdxHitsNPHarm2_3);
      trun40.hDeDx       -> Fill(reco::trackCollection[i].dEdxHitsNPTrun40_1000);
      trun40.hDeDx7Hits  -> Fill(reco::trackCollection[i].dEdxHitsNPTrun40_7);
      trun40.hDeDx3Hits  -> Fill(reco::trackCollection[i].dEdxHitsNPTrun40_3);
      median.hDeDx       -> Fill(reco::trackCollection[i].dEdxHitsNPMedian_1000);
      median.hDeDx7Hits  -> Fill(reco::trackCollection[i].dEdxHitsNPMedian_7);
      median.hDeDx3Hits  -> Fill(reco::trackCollection[i].dEdxHitsNPMedian_3);

      harm2.hPDeDx       -> Fill(p, reco::trackCollection[i].dEdxHitsNPHarm2_1000);
      harm2.hPDeDx7Hits  -> Fill(p, reco::trackCollection[i].dEdxHitsNPHarm2_7);
      harm2.hPDeDx3Hits  -> Fill(p, reco::trackCollection[i].dEdxHitsNPHarm2_3);
      trun40.hPDeDx      -> Fill(p, reco::trackCollection[i].dEdxHitsNPTrun40_1000);
      trun40.hPDeDx7Hits -> Fill(p, reco::trackCollection[i].dEdxHitsNPTrun40_7);
      trun40.hPDeDx3Hits -> Fill(p, reco::trackCollection[i].dEdxHitsNPTrun40_3);
      median.hPDeDx      -> Fill(p, reco::trackCollection[i].dEdxHitsNPMedian_1000);
      median.hPDeDx7Hits -> Fill(p, reco::trackCollection[i].dEdxHitsNPMedian_7);
      median.hPDeDx3Hits -> Fill(p, reco::trackCollection[i].dEdxHitsNPMedian_3);
    }

};

#endif

