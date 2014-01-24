#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLine.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TObject.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TColorWheel.h"
#include "TLatex.h"
#include "TROOT.h"
#include "plotStyle.h"

TH1D* GetTH1D(TString filename, TString objectName);
TH2D* GetTH2D(TString filename, TString objectName);
int plotHistograms1D(TString histoName,TString fileName, TString saveName, TString titleName, TString xTitle, double yMax, bool log);
int plotHistograms2D(TString histoName,TString fileName, TString saveName, TString titleName, TString xTitle, TString yTitle, double max);


int plotAllHistogram(){


  bool gensim = false;
  bool reco   = true;

  TString title, sample, sourceFile;
  //TString targetName = "plotsHSCPdEdx/";
  TString targetName = "plotsHSCPTest/";
  
 
  //double hbar = 6.582119*pow(10,-25);   // in [GeV*s]
  //double c0   = 299792458;              // in [m/s]


  double ctau[15] = {0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0, 8.0, 10., 15., 20., 30.};
  
  std::vector<int> mass;
  mass.push_back(100);
  //mass.push_back(800);
  std::vector<int> width;
  //width.push_back(0);
  width.push_back(14);


  if((TString) targetName == "plotsHSCPTest/"){

    sourceFile = "analyzerTest_histograms" + sample + ".root";
    plotHistograms2D("hZVsRho",sourceFile,targetName+"hZVsRhoChipm" + sample + ".pdf",title,"z","#rho",-1);



    return 0;
  }

  for(unsigned int m=0; m<mass.size(); m++){
    for(unsigned int i=0; i<width.size(); i++){

      sample = (TString) "_m" + (long) mass[m] + (TString) "_width" + (long) width[i];
      //sample = "";
      sourceFile = "analyzerdEdx_histograms" + sample + ".root";
      cout<<"Histograms taken from following file: "<<sourceFile<<endl;

      if(reco){

	title = Form("c#tau = %2.2f m, mass = %i GeV",ctau[width[i]],mass[m]);
	//2D 800 GeV
	/*
	plotHistograms2D("hPDeDx",sourceFile,targetName+"hPDeDx" + sample + ".pdf",title,"p [GeV]","dE/dx",1400);
	plotHistograms2D("hPDeDxChipm",sourceFile,targetName+"hPDeDxChipm" + sample + ".pdf",title,"p [GeV]","dE/dx",30);
	plotHistograms2D("hPDeDxNoChipm",sourceFile,targetName+"hPDeDxNoChipm" + sample + ".pdf",title,"p [GeV]","dE/dx",1000);
	plotHistograms2D("hPDeDx7Hits",sourceFile,targetName+"hPDeDx7Hits" + sample + ".pdf",title,"p [GeV]","dE/dx",1400);
	plotHistograms2D("hPDeDxChipm7Hits",sourceFile,targetName+"hPDeDxChipm7Hits" + sample + ".pdf",title,"p [GeV]","dE/dx",30);
	plotHistograms2D("hPDeDxNoChipm7Hits",sourceFile,targetName+"hPDeDxNoChipm7Hits" + sample + ".pdf",title,"p [GeV]","dE/dx",1000);
	*/
	//2D
	plotHistograms2D("hPDeDx",sourceFile,targetName+"hPDeDx" + sample + ".pdf",title,"p [GeV]","dE/dx",1400);
	plotHistograms2D("hPDeDxChipm",sourceFile,targetName+"hPDeDxChipm" + sample + ".pdf",title,"p [GeV]","dE/dx",100);
	plotHistograms2D("hPDeDxNoChipm",sourceFile,targetName+"hPDeDxNoChipm" + sample + ".pdf",title,"p [GeV]","dE/dx",1400);
	plotHistograms2D("hPDeDx7Hits",sourceFile,targetName+"hPDeDx7Hits" + sample + ".pdf",title,"p [GeV]","dE/dx",1400);
	plotHistograms2D("hPDeDxChipm7Hits",sourceFile,targetName+"hPDeDxChipm7Hits" + sample + ".pdf",title,"p [GeV]","dE/dx",100);
	plotHistograms2D("hPDeDxNoChipm7Hits",sourceFile,targetName+"hPDeDxNoChipm7Hits" + sample + ".pdf",title,"p [GeV]","dE/dx",1400);
	//1D
	plotHistograms1D("hgenPChipm",sourceFile, targetName+"hgenPChipm" + sample + ".pdf",title,"p^{#Chi^{#pm}} [GeV]",-1,1);
	plotHistograms1D("hDeDx",sourceFile, targetName+"hDeDx" + sample + ".pdf",title,"dE/dx [?]",-1,0);
	plotHistograms1D("hDeDx7Hits",sourceFile, targetName+"hDeDx7Hits" + sample + ".pdf",title,"dE/dx [?]",-1,0);


      }

      /*
      if(gensim){


	if(i==0) title = "#Gamma = 9 * 10^{-14}";
	//if(i==0) title = "#Gamma = 6.0 * 10^{-16}";
	else if(i==1) title = "#Gamma = 5 * 10^{-15}";
	else if(i==2) title = "#Gamma = 1 * 10^{-15}";
	else if(i==3) title = "#Gamma = 5 * 10^{-16}";
	else if(i==4) title = "#Gamma = 7 * 10^{-20}";

	title = Form("c#tau = %2.2f m, mass = %i GeV",ctau[width[i]],mass[m]);
	//plotHistograms1D("trkLengthChargino",Form("analyzer_histograms_Width%i.root",i),Form("plots/ChipmTrkLength_Width%i.pdf",i),title,"s_{Track}");
	//plotHistograms1D("ptSimTrackChipm",Form("analyzer_histograms_Width%i.root",i),Form("plots/ChipmPtSimTrack_Width%i.pdf",i),title,"p_{T}^{#Chi^{#pm}}");
	//plotHistograms2D("ptVsDecayedChipm",Form("analyzer_histograms_Width%i.root",i),Form("plots/ChipmPtVsDecayed_Width%i.pdf",i),title,"p_{T}^{#Chi^{#pm}}","decayed");

	//GenParticles plots
	plotHistograms1D("hgenPtChipm",sourceFile, targetName+"hgenPtChipm_" + sample + ".pdf",title,"p_{T}^{#Chi^{#pm}} [GeV]",-1,1);
	plotHistograms1D("hgenEtaChipm",sourceFile,targetName+"hgenEtaChipm_" + sample + ".pdf",title,"#eta^{#Chi^{#pm}}",-1,0);
	plotHistograms1D("hgenPhiChipm",sourceFile,targetName+"hgenPhiChipm_" + sample + ".pdf",title,"#Phi^{#Chi^{#pm}}",-1,0);
	plotHistograms1D("hZChipm",sourceFile,targetName+"hZChipm_" + sample + ".pdf",title,"z",-1,0);
	plotHistograms1D("htrkLengthDeltaRayChipm",sourceFile, targetName+"htrkLengthDeltaRayChipm_" + sample + ".pdf",title,"s^{Delta ray}_{Track} [cm]",-1,1);
	plotHistograms1D("hDeltaRayPDGIdChipm",sourceFile, targetName+"hDeltaRayPDGIdChipm_" + sample + ".pdf",title,"pdg Id",-1,0);
	plotHistograms1D("hgenEtaNotDecayedChipm",sourceFile, targetName+"hgenEtaNotDecayedChipm_" + sample + ".pdf",title,"#eta^{#Chi^{#pm}}_{not decayed}",-1,0);
	plotHistograms1D("hTimeOfFlightChipm",sourceFile, targetName+"hTimeOfFlightChipm_" + sample + ".pdf",title,"TOF^{#Chi^{#pm}} [s]",-1,1);


	plotHistograms1D("htrkLengthChipm",sourceFile,targetName+"hChipmTrkLength_" + sample + ".pdf",title,"s_{Track} [cm]",900,1);
	plotHistograms1D("htrkLengthChip",sourceFile,targetName+"hChipTrkLength_" + sample + ".pdf",title,"s_{Track}^{#Chi^{+}} [cm]",120,1);
	plotHistograms1D("htrkLengthChim",sourceFile,targetName+"hChimTrkLength_" + sample + ".pdf",title,"s_{Track}^{#Chi^{-}} [cm]",120,1);
	//plotHistograms1D("ptSimTrackChipm",Form("analyzer_histograms.root"),Form("plots/ChipmPtSimTrack.pdf"),title,"p_{T}^{#Chi^{#pm}}");
	//plotHistograms2D("ptVsDecayedChipm",Form("analyzer_histograms.root"),Form("plots/ChipmPtVsDecayed.pdf"),title,"p_{T}^{#Chi^{#pm}}","decayed");
	plotHistograms2D("hRhoVsZChipm",sourceFile,targetName+"hRhoVsZChipm_" + sample + ".pdf",title,"#rho","z");
	plotHistograms2D("hZVsRhoChipm",sourceFile,targetName+"hZVsRhoChipm_" + sample + ".pdf",title,"z","#rho");
      }
      */
    }
  }

  return 0;

}


int plotHistograms1D(TString histoName,TString fileName, TString saveName, TString titleName, TString xTitle, double yMax, bool log){

  TeresaPlottingStyle::init();

  TH1D *histo;
  TCanvas* canvas1;
  //TString fileName;
  //TLatex* info;
      

  canvas1 = new TCanvas("canvas1","canvas",0,0,500,500);
  canvas1 ->cd();
  

  //fileName = "analyzer_histograms_Width0.root";
  histo = GetTH1D(fileName,histoName);
  histo->SetTitle(titleName);
  histo->GetXaxis()->SetTitle(xTitle);

  if(yMax>0) histo->SetMaximum(yMax);
  if(histo->FindLastBinAbove(12)==0) histo->GetXaxis()->SetRangeUser(0,15);

  if(log) canvas1 -> SetLogy();
  else histo->SetMinimum(0);

  histo ->Draw("E");


  canvas1->SaveAs(saveName);

  return 0;
  
}

int plotHistograms2D(TString histoName,TString fileName, TString saveName, TString titleName, TString xTitle, TString yTitle, double max){

  TeresaPlottingStyle::init();

  TH2D *histo;
  TCanvas* canvas1;

  canvas1 = new TCanvas("canvas1","canvas",0,0,500,500);
  canvas1 ->cd();
  //canvas1 -> SetLogy();
  //canvas1 -> SetLogx();
  //canvas1 -> SetLogy();

  histo = GetTH2D(fileName,histoName);
  histo->SetTitle(titleName);
  //histo->GetYaxis()->SetRangeUser(0,20.);
  if(max>0) histo->SetMaximum(max);
  histo ->Draw("COLZ");
  histo->GetXaxis()->SetTitle(xTitle);
  histo->GetYaxis()->SetTitle(yTitle);

  canvas1->SaveAs(saveName);

  return 0;
  
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Some Functions


TH1D* GetTH1D(TString filename, TString objectName){

  TH1D* object = 0;

  TFile* file =  TFile::Open(filename);
  if(file != 0) file -> GetObject(objectName,object);
  object -> SetDirectory(0);
  delete file;

  return object;
}

TH2D* GetTH2D(TString filename, TString objectName){

  TH2D* object = 0;

  TFile* file =  TFile::Open(filename);
  if(file != 0) file -> GetObject(objectName,object);
  object -> SetDirectory(0);
  delete file;

  return object;
}
