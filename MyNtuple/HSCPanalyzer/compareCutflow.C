#include <iostream>
#include <fstream>
#include <iomanip>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

int compareCutflow(TString fileName="/afs/desy.de/user/t/tlenz/HSCPworkdir/analyzer/CMSSW_5_3_2_patch4/src/Ntuples/MyNtuple/HSCPanalyzer/analyzerWells_histograms.root")
{


 TFile *file = TFile::Open(fileName,"READ");
 TH1D* hCutflow; 
 file->GetObject("AllSelectedTracks/countsEventCuts",hCutflow);

 //hCutflow->Draw();

 ofstream cutflow;
 cutflow.open("cutflow.txt");

 cutflow<<"CUTFLOW:"<<endl<<endl;
 cout<<"CUTFLOW:"<<endl<<endl; 

 for(int i=1; i<hCutflow->GetNbinsX(); i++){

   if(strcmp(hCutflow->GetXaxis()->GetBinLabel(i), "")==0) break;

   cutflow<<left<<setw(6)<<i<<left<<":   ";
   cutflow<<left<<setw(35)<<hCutflow->GetXaxis()->GetBinLabel(i);
   cutflow<<std::setprecision(6)<<hCutflow->GetBinContent(i)<<endl;

   cout<<left<<setw(6)<<i<<left<<":   ";
   cout<<left<<setw(35)<<hCutflow->GetXaxis()->GetBinLabel(i);
   cout<<std::setprecision(6)<<hCutflow->GetBinContent(i)<<endl;
 }
 

 cutflow.close();

 return 0;


}
