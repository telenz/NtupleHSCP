#ifndef ANALYZERHISTOGRAMS_H
#define ANALYZERHISTOGRAMS_H

#include "analyzerAODSIM.h"
#include "TH1.h"
#include "TVector3.h"

using namespace std;

const double K = 2.529;
const double C = 2.772;

class Hist
{

 public:
  TH1D *htrackPt;
  TH1D *htrackP;
  TH1D *htrackEta;
  TH1D *htrackd0;
  TH1D *htrackdz;
  TH1D *htrackNValid;
  TH1D *htrackNLostMid;
  TH1D *htrackNLostInner;
  TH1D *htrackNLostOuter;
  TH1D *hRelativePtError;
  TH1D *htrackPdgId;
  TH1D *htrackgenParticle;
  //TH2D *htrackPdgIdMass;
  //TH2D *htrackPMass;
  //TH2D *htrackgenBetaMass;
  //TH2D *htrackgenBetaGammaMass;

  //TH1D *hMass;
  //TH1D *hMass_1Hits;
  //TH1D *hMass_3Hits;
  //TH1D *hMass_7Hits;

  TH1D* h1stjetpt;
  TH1D* hMet;
  TH1D *hnPFMet;   
  TH1D *hnPFJet;
  
  TH1D *hDeltaPhi;

  TH1D* hVertexNdof;
  TH1D* hVertexZ;
  TH1D* hVertexRho;
  TH1D *hnVertex;

  TH1D* hgenPtChi;
  TH1D* hgenEtaChi;
  TH1D* hgenPhiChi;
  TH1D* hgenBetaChi;
  TH1D *hgenPChi;

  


 public:Hist(TString histName, outputFile ofile_)
    {
      ofile_.file_->mkdir(histName);
      ofile_.file_->cd(histName);
      htrackPt = new TH1D("htrackPt","htrackPt",150,0,1500);
      htrackP = new TH1D("htrackP","htrackP",150,0,1500);
      htrackEta = new TH1D("htrackEta","htrackEta",100,-5,5);
      htrackd0 = new TH1D("htrackd0","htrackd0",100,-5,5);
      htrackdz = new TH1D("htrackdz","htrackdz",100,-50,50);
      htrackNValid = new TH1D("htrackNValid","htrackNValid",40,0,40);
      htrackNLostMid = new TH1D("htrackNLostMid","htrackNLostMid",10,0,10);
      htrackNLostInner = new TH1D("htrackNLostInner","htrackNLostInner",10,0,10);
      htrackNLostOuter = new TH1D("htrackNLostOuter","htrackNLostOuter",20,0,20);
      hRelativePtError = new TH1D("hRelativePtError" ,"hRelativePtError",300,0,3);

      htrackPdgId = new TH1D("htrackPdgId","htrackPdgId",500,0,500);
      htrackgenParticle = new TH1D("htrackgenParticle","htrackgenParticle",1,0,1);
      htrackgenParticle->Fill("unmatched", 0);
      htrackgenParticle->Fill("d", 0);
      htrackgenParticle->Fill("u", 0);
      htrackgenParticle->Fill("s", 0);
      htrackgenParticle->Fill("c", 0);
      htrackgenParticle->Fill("b", 0);
      htrackgenParticle->Fill("t", 0);
      htrackgenParticle->Fill("e", 0);
      htrackgenParticle->Fill("mu", 0);
      htrackgenParticle->Fill("#tau", 0);
      htrackgenParticle->Fill("g", 0);
      htrackgenParticle->Fill("#gamma", 0);
      htrackgenParticle->Fill("pi", 0);
      htrackgenParticle->Fill("mesons", 0);
      htrackgenParticle->Fill("baryons", 0);
      htrackgenParticle->Fill("others", 0);
      //htrackPdgIdMass = new TH2D("htrackPdgIdMass","htrackPdgIdMass",500,0,500,75,0,1500);
      //htrackPMass = new TH2D("htrackPMass","htrackPMass",150,0,1500,75,0,1500);
      //htrackgenBetaMass = new TH2D("htrackgenBetaMass","htrackgenBetaMass",100,0,1,75,0,1500);
      //htrackgenBetaGammaMass = new TH2D("htrackgenBetaGammaMass","htrackgenBetaGammaMass",30000,0,3000,75,0,1500);

      //hMass = new TH1D("hMass","hMass",75,0,1500);
      //hMass_7Hits = new TH1D("hMass_7Hits","hMass_7Hits",75,0,1500);
      //hMass_3Hits = new TH1D("hMass_3Hits","hMass_3Hits",75,0,1500);
      //hMass_1Hits = new TH1D("hMass_1Hits","hMass_1Hits",75,0,1500);

      hnPFMet    = new TH1D("hnPFMet","hnPFMet",200,0,200);
      hnPFJet    = new TH1D("hnPFJet","hnPFJet",20,0,20);
      hDeltaPhi  = new TH1D("hDeltaPhi","hDeltaPhi",32,0,3.2);
      h1stjetpt  = new TH1D("h1stjetpt","h1stjetpt",150,0,150);
      hMet       = new TH1D("hMet","hMe",150,0,150);

      hnVertex    = new TH1D("hnVertex","hnVertex",200,0,200);
      hVertexNdof = new TH1D("hVertexNdof", "hVertexNdof",500,0,500);  
      hVertexZ    = new TH1D("hVertexZ", "hVertexZ",160,-40,40);  
      hVertexRho  = new TH1D("hVertexRho", "hVertexRho",20,0,10);

      hgenPtChi   = new TH1D("hgenPtChi","hgenPtChi",150,0,1500);
      hgenPChi    = new TH1D("hgenPChi","hgenPChi",150,0,1500);
      hgenEtaChi  = new TH1D("hgenEtaChi","hgenEtaChi",200,-5,5);
      hgenPhiChi  = new TH1D("hgenPhiChi","hgenPhiChi",100,0,3.142);
      hgenBetaChi = new TH1D("hgenBetaChi","hgenBetaChi",100,0,1);

      
     
    };

  void FillTrackVariables_AODSIM(std::vector<Track_s> inputCollection)
  {  

    for(unsigned int i=0; i<inputCollection.size(); i++){
      TVector3 trackVector;
      trackVector.SetPtEtaPhi(inputCollection[i].pt,inputCollection[i].eta,inputCollection[i].phi);
      double p = std::sqrt(std::pow(inputCollection[i].pt,2) + std::pow(trackVector.Pz(),2));
      //double mass =  sqrt((inputCollection[i].dEdxHitsNPHarm2_1000 - C)*pow(p,2)/K);
      //double mass_7Hits =  sqrt((inputCollection[i].dEdxHitsNPHarm2_7 - C)*pow(p,2)/K);
      //double mass_3Hits =  sqrt((inputCollection[i].dEdxHitsNPHarm2_3 - C)*pow(p,2)/K);
      //double mass_1Hits =  sqrt((inputCollection[i].dEdxHitsNPHarm2_1 - C)*pow(p,2)/K);
      //hMass->Fill(mass);
      //hMass_7Hits->Fill(mass_7Hits);
      //hMass_3Hits->Fill(mass_3Hits);
      //hMass_1Hits->Fill(mass_1Hits);
      htrackP->Fill(p);
      htrackPt->Fill(inputCollection[i].pt);
      htrackEta->Fill(inputCollection[i].eta);
      htrackd0->Fill(inputCollection[i].d0);
      htrackdz->Fill(inputCollection[i].dz);
      htrackNValid->Fill(inputCollection[i].numberOfValidHits);
      htrackNLostMid->Fill(inputCollection[i].numberOfLostHits);
      htrackNLostInner->Fill(inputCollection[i].trackerExpectedHitsInner_numberOfLostHits);
      htrackNLostOuter->Fill(inputCollection[i].trackerExpectedHitsOuter_numberOfLostHits);
      hRelativePtError->Fill(inputCollection[i].ptError/inputCollection[i].pt);
      //htrackPMass-> Fill(p,mass);
      
      //if(inputCollection[i].beta*1/sqrt(1.-pow(inputCollection[i].beta,2))>3000){
     
      if(abs(inputCollection[i].pdgId==0))        htrackgenParticle->Fill("unmatched", 1);
      else if(abs(inputCollection[i].pdgId)==1)   htrackgenParticle->Fill("d", 1);
      else if(abs(inputCollection[i].pdgId)==2)   htrackgenParticle->Fill("u", 1);
      else if(abs(inputCollection[i].pdgId)==3)   htrackgenParticle->Fill("s", 1);
      else if(abs(inputCollection[i].pdgId)==4)   htrackgenParticle->Fill("c", 1);
      else if(abs(inputCollection[i].pdgId)==5)   htrackgenParticle->Fill("b", 1);
      else if(abs(inputCollection[i].pdgId)==6)   htrackgenParticle->Fill("t", 1);
      else if(abs(inputCollection[i].pdgId)==11)  htrackgenParticle->Fill("e", 1);
      else if(abs(inputCollection[i].pdgId)==13)  htrackgenParticle->Fill("mu", 1);
      else if(abs(inputCollection[i].pdgId)==15)  htrackgenParticle->Fill("#tau", 1);
      else if(abs(inputCollection[i].pdgId)==21)  htrackgenParticle->Fill("g", 1);
      else if(abs(inputCollection[i].pdgId)==22)  htrackgenParticle->Fill("#gamma", 1);
      else if(abs(inputCollection[i].pdgId)==211) htrackgenParticle->Fill("pi", 1);
      else if(abs(inputCollection[i].pdgId)<1000 &&abs(inputCollection[i].pdgId)!=211 && abs(inputCollection[i].pdgId)>23)   htrackgenParticle->Fill("mesons", 1);
      else if(abs(inputCollection[i].pdgId)>1000 && abs(inputCollection[i].pdgId)<10000)   htrackgenParticle->Fill("baryons", 1);
      else                                        htrackgenParticle->Fill("others", 1);

      if(inputCollection[i].beta<10){
	htrackPdgId->Fill(abs(inputCollection[i].pdgId));
	//htrackPdgIdMass->Fill(abs(inputCollection[i].pdgId),mass);
	//htrackgenBetaMass-> Fill(inputCollection[i].beta,mass);
	//if(inputCollection[i].beta>0.9999999){
	//  htrackgenBetaGammaMass-> Fill(2500,mass);
	//}
	//else{
	//  htrackgenBetaGammaMass-> Fill(inputCollection[i].beta*1/sqrt(1.-pow(inputCollection[i].beta,2)),mass);
	//	}
      }
    }
  };

  void FillVertexHistograms(std::vector<Vertex_s> inputCollection)
  {
    hnVertex ->Fill(inputCollection.size());
    for(unsigned int i=0; i<inputCollection.size(); i++){
      hVertexNdof -> Fill(inputCollection[i].ndof);
      hVertexZ    -> Fill(inputCollection[i].z);
      hVertexRho  -> Fill(std::sqrt(std::pow(inputCollection[i].x,2) + std::pow(inputCollection[i].y,2)));
    }
  }

  void FillGenParticleHistograms(std::vector<GenParticle_s> inputCollection)
  {
    for(unsigned int i=0; i<inputCollection.size(); i++){

      hgenBetaChi-> Fill(inputCollection[i].p/inputCollection[i].energy);
      hgenPtChi  -> Fill(inputCollection[i].pt);
      hgenPChi  -> Fill(inputCollection[i].p);
      hgenEtaChi -> Fill(inputCollection[i].eta);
      hgenPhiChi -> Fill(inputCollection[i].phi);
    }
  }

};


#endif
