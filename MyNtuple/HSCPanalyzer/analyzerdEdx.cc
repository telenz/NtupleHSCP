//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "analyzerRecoTracks.h"
#include "analyzerStructs.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <time.h>
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line
  clock_t t;
  t = clock();
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples

  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  //cout << "Number of events: " << nevents << endl;

  // Select variables to be read

  selectVariables(stream);


  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the
  // following line

  // TApplication app("analyzer", &argc, argv);

  /**
     Notes 1
     -------
     1. Use
     ofile = outputFile(cmdline.outputfile, stream)

     to skim events to output file in addition to writing out histograms.

     2. Use
     ofile.addEvent(event-weight)

     to specify that the current event is to be added to the output file.
     If omitted, the event-weight is defaulted to 1.

     3. Use
     ofile.count(cut-name, event-weight)

     to keep track, in the count histogram, of the number of events
     passing a given cut. If omitted, the event-weight is taken to be 1.
     If you want the counts in the count histogram to appear in a given
     order, specify the order, before entering the event loop, as in
     the example below

     ofile.count("NoCuts", 0)
     ofile.count("GoodEvent", 0)
     ofile.count("Vertex", 0)
     ofile.count("MET", 0)

     Notes 2
     -------
     By default all variables are saved. Before the event loop, you can use
  
     select(objectname)
	  
     e.g.,
	
     select("GenParticle")
  
     to declare that you intend to select objects of this type. The
     selection is done using

     select(objectname, index)
	  
     e.g.,
	  
     select("GenParticle", 3),
  
     which is called within the event loop. Call saveSelectedObjects()
     before a call to addEvent if you wish to save the selected objects.
     All other objects are saved by default.
	 
     NB: If you declare your intention to select objects of a given type
     by calling select(objectname), but subsequently fail to select
     them using select(objectname, index) then none will be saved!
  */
  

  outputFile ofile(cmdline.outputfilename);

  //---------------------------------------------------------------------------
  // Declare histograms
  //---------------------------------------------------------------------------


  class SetOfHistogramsChipm NoChipm("NoChipm");
  class SetOfHistogramsChipm Chipm("Chipm");
  class SetOfHistogramsChipm All("All");

  class SetOfVertexHistograms vertices;

  TH1D *hnPFMet    = new TH1D("hnPFMet","hnPFMet",200,0,200);
  TH1D *hnVertex   = new TH1D("hnVertex","hnVertex",200,0,200);
  TH1D *hnPFJet    = new TH1D("hnPFJet","hnPFJet",20,0,20);
  TH1D *hgenPChipm = new TH1D("hgenPChipm","hgenPChipm",150,0,1500);
  TH1D *hMet       = new TH1D("hMet","hMet",150,0,1500);
  TH1D *h1stjetpt  = new TH1D("h1stjetpt","h1stjetpt",150,0,1500);
  TH1D *hDeltaPhi  = new TH1D("hDeltaPhi","hDeltaPhi",32,0,3.2);
  

  //---------------------------------------------------------------------------
  // Declaration of Variables
  //---------------------------------------------------------------------------

  double dPhi    = 0;
  double dEta    = 0;
  double dR      = 0;
  double dRchip  = 0;
  double dRchim  = 0;
  double dPtchip = 0;
  double dPtchim = 0;
  int Nchip = 0;
  int Nchim = 0;
  int chipmFound = 0;
  int n = 10000 ;
  TVector3 chip;
  TVector3 chim;
  
  int count = 0;
  int beforeVertexCut = 0;
  int afterVertexCut  = 0;
  int beforeMETCut    = 0;
  int afterMETCut     = 0;
  int before1stJetCut = 0;
  int after1stJetCut  = 0;
  //---------------------------------------------------------------------------
  bool vertexCut     = true;
  bool metCut        = true;
  bool leadingJetCut = true;
  bool candidateCut  = true;

  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry <nevents; ++entry)
    {
	  
      // Read event into memory
      stream.read(entry);
      // NB: call to clear object selection map (indexmap)
      initialize();
	  
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file analyzer.h to find out what structs
      // are available. Alternatively, you can call individual fill functions.
      fillObjects();


      /******************************************************************************************************************************
       ******************************************************************************************************************************
       ******************************************************************************************************************************
       *****************************************************************************************************************************/

      ofile.count("NoCuts", 1);
      //-------------------------------------------------------------- Cuts ---------------------------------------------------------

      //-------------- Vertex Cut ------------------
      beforeVertexCut +=1;
      if(!reco::isGoodVertex() && vertexCut){
	continue;
      }
      afterVertexCut +=1;
      //-------------- Vertex Cut ------------------
      ofile.count("vertexCuts", 1);



      reco::getCandidateTrackCollection();
      //-------------- Candidate cut ------------------
      if(candidateCut){
	if(reco::trackCollection.size()==0) continue;
      }
      //-----------------------------------------------
      ofile.count("trackCuts", 1);

      hnPFMet->Fill(nPFMET);
      hnVertex->Fill(nVertex);
      for(int i=0; i<nVertex; i++){
	vertices.FillVertexHistograms(i);
      }
	  
      //----------------- MET Cut ------------------
      beforeMETCut += 1;
      hMet->Fill(PFMET[0].et);
      if(metCut){
	if(nPFMET==1){
	  if(PFMET[0].et<220.) continue;
	}
	else{
	  continue;
	}
      }
      afterMETCut += 1;
      //---------------- MET Cut ------------------
      ofile.count("metCut", 1);

      reco::findChipmInGenParticleCollection();
	  
      
      reco::getSelectedJetCollection();
      hnPFJet->Fill(reco::jetCollection.size());
	 
      //----------------- 1st Jet Cuts ------------------
      before1stJetCut += 1;
      if(reco::jetCollection.size()>0)  h1stjetpt->Fill(reco::jetCollection[0].pt);
      if(reco::jetCollection.size()>1){
	hDeltaPhi->Fill(std::abs(TVector2::Phi_mpi_pi(reco::jetCollection[0].phi-reco::jetCollection[1].phi)));
      }
      if(leadingJetCut && !reco::leadingJetRequirementsFullfilled()) continue;
      after1stJetCut += 1;
      ofile.count("1stJetCut", 1);
      //----------------- 1st Jet Cuts ------------------
	  

      // ---------------------
      // -- event selection --
      // ---------------------
	  
      if(!reco::zeroChip){
	chip.SetPtEtaPhi(reco::chipGenParticle.pt,reco::chipGenParticle.eta,reco::chipGenParticle.phi);
	hgenPChipm -> Fill(sqrt(pow(reco::chipGenParticle.pt,2)+pow(chip.Pz(),2)));
      }
      if(!reco::zeroChim){
	chim.SetPtEtaPhi(reco::chimGenParticle.pt,reco::chimGenParticle.eta,reco::chimGenParticle.phi);
	hgenPChipm -> Fill(sqrt(pow(reco::chimGenParticle.pt,2)+pow(chim.Pz(),2)));
      }


      // Loop over track collection
      Nchip=0;
      Nchim=0;
      for(int i=0; i<reco::trackCollection.size(); i++){

	All.FillTrackHistograms(i);
	All.FillDeDxHistograms(i);
	    
	dPhi    = std::abs(reco::chipGenParticle.phi - reco::trackCollection[i].phi);
	dEta    = std::abs(reco::chipGenParticle.eta - reco::trackCollection[i].eta);
	dRchip  = std::sqrt(pow(dPhi,2) + pow(dEta,2));
	dPtchip = std::abs(reco::trackCollection[i].pt - reco::chipGenParticle.pt)/reco::chipGenParticle.pt;
	dPhi    = std::abs(reco::chimGenParticle.phi - reco::trackCollection[i].phi);
	dEta    = std::abs(reco::chimGenParticle.eta - reco::trackCollection[i].eta);
	dRchim  = std::sqrt(pow(dPhi,2) + pow(dEta,2));
	dPtchim = std::abs(reco::trackCollection[i].pt - reco::chimGenParticle.pt)/reco::chimGenParticle.pt;

	    
	if((!reco::zeroChip && dRchip<0.1 && dPtchip<0.3 && Nchip==0) || (!reco::zeroChim && dRchim<0.1 && dPtchim<0.3 && Nchim == 0)){
	  chipmFound += 1;

	  if(!reco::zeroChip && dRchip<0.1 && dPtchip<0.3) Nchip += 1;
	  if(!reco::zeroChim && dRchim<0.1 && dPtchim<0.3) Nchim += 1;

	  Chipm.FillTrackHistograms(i);
	  Chipm.FillDeDxHistograms(i);
	  //Chipm.hsumPtDeltaR0p3->Fill(sumPtDeltaR0p3);
	}
	else{

	  NoChipm.FillTrackHistograms(i);
	  NoChipm.FillDeDxHistograms(i);
	  //NoChipm.hsumPtDeltaR0p3->Fill(sumPtDeltaR0p3);
	      
	}
	    
      }
	  
      // ---------------------
      // -- fill histograms --
      // ---------------------


      //------------------------------------------------------------------------------------------------------------------------------------
      //-------- GenParticle Helper  Collection --------------------------------------------------------------------------------------------
      //------------------------------------------------------------------------------------------------------------------------------------

      count +=1;
    }//end of loop over events

  cout<<count<<" events passed the selection."<<endl;  
  cout<<endl<<"chipmFound = "<<chipmFound<<endl<<endl;
  cout<<"matching efficiency = "<<chipmFound/(2.*n)<<endl<<endl;
  cout<<endl<<"It took "<<(clock()-t)/CLOCKS_PER_SEC<<" sec"<<endl;


  cout<<"Vertex cut efficiency  = "<<1.*afterVertexCut/beforeVertexCut<<endl;
  cout<<"MET cut efficiency     = "<<1.*afterMETCut/beforeMETCut<<endl;
  cout<<"before1stJetCut        = "<<before1stJetCut<<endl;
  cout<<"1st Jet cut efficiency = "<<1.*after1stJetCut/before1stJetCut<<endl;

  stream.close();
  ofile.close();
  
  cout<<endl;
  return 0;
}
