
//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "analyzerFunctionsGenSim.h"
#include "analyzerGenSim.h"
#include "TLorentzVector.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Get file list and histogram filename from command line

  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples

  vector<string> filenames = getFilenames(cmdline.filelist);
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  cout << "Number of events: " << nevents << endl;

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
  //SimTrack Collextion

  TH2D *hZVsRhoChipm  = new TH2D("hZVsRho","hZVsRho",140,0,1400,8000,0,800);
  TH1D *htrackpt      = new TH1D("htrackpt","htrackpt",150,0,1500);
  TH1D *ptOfPions     = new TH1D("ptOfPions","ptOfPions",100,0,1.5);
  TH1D *pOfPions      = new TH1D("pOfPions","pOfPions",100,0,4);
  TH1D *ptOfChi0      = new TH1D("ptOfChi0","ptOfChi0",100,0,1000);
  TH1D *pOfChi0       = new TH1D("pOfChi0","pOfChi0",100,0,1000);
  TH1D *ptOfChipm     = new TH1D("ptOfChipm","ptOfChipm",100,0,1000);
  TH1D *pOfChipm      = new TH1D("pOfChipm","pOfChipm",100,0,1000);
  TH1D *ptOfAllPions  = new TH1D("ptOfAllPions","ptOfAllPions",100,0,1000);
  TH1D *pOfAllPions   = new TH1D("pOfAllPions","pOfAllPions",100,0,1000);
  
  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry <nevents; ++entry)
	{

	  
	  /*
	  for(int i=0; i<nSimTrack; i++){
	    if(abs(SimTrack[i].type)==1000024){   
	      for(int j=0; j<nSimVertex; j++){
		if(SimVertex[j].parentIndex==SimTrack[i].trackId){
		  double Rho = sqrt(pow(SimVertex[j].position_x,2)+pow(SimVertex[j].position_y,2));
		  hZVsRhoChipm -> Fill(abs(SimVertex[j].position_z),Rho);
		}
	      }
	    }
	  }
	  */
	  // Read event into memory
	  stream.read(entry);
	  // NB: call to clear object selection map (indexmap)
	  initialize();
	  
	  // Uncomment the following line if you wish to copy variables into
	  // structs. See the header file analyzer.h to find out what structs
	  // are available. Alternatively, you can call individual fill functions.
	  fillObjects();


	  decayedChip = false;
	  decayedChim = false;
	  AllPions.clear();
	  PionsFromDecay.clear();
	  Chi0FromDecay.clear();

	  for(int i=0; i<SimTrack.size(); i++)
	    {
	    htrackpt->Fill(SimTrack[i].momentum_pt);
	    }

	  //-----------------------------------------------------------------------------------------
	  //-------- Find Charginos and Neutralinos in collection  ----------------------------------
	  //-----------------------------------------------------------------------------------------
	  findChipmInGenParticleCollection();
	  findChipmInSimTrackCollection();
	  findChi0InSimTrackCollection();
	  findChipmDecayInSimVertexCollection();
	     
	  //----------------------------------------------------------------------------------------
	  //--------- Calculate decay position of chargino -----------------------------------------
	  //----------------------------------------------------------------------------------------

	  double Rho;
	  if(decayedChip){
	    Rho = sqrt(pow(chipDecaySimVertex.position_x,2)+pow(chipDecaySimVertex.position_y,2));
	    hZVsRhoChipm -> Fill(abs(chipDecaySimVertex.position_z),Rho);
	    decayVertexFound3+=1;
	  }

	  if(decayedChim){
	    Rho = sqrt(pow(chimDecaySimVertex.position_x,2)+pow(chimDecaySimVertex.position_y,2));
	    hZVsRhoChipm -> Fill(abs(chimDecaySimVertex.position_z),Rho);
	    decayVertexFound3+=1;
	  }

	  //----------------------------------------------------------------------------------------
	  //--------- Find the pions in the SIMTrack Collection ------------------------------------
	  //----------------------------------------------------------------------------------------
	  findPiInSimTrackCollection();
	  findDecayParticles();
	  if(PionsFromDecay.size() > 2) cout<<"PionsFromDecay.size() = "<<PionsFromDecay.size()<<endl;
	  for(unsigned int i=0; i<PionsFromDecay.size(); i++) 
	    {

	      TLorentzVector pi;
	      pi.SetPtEtaPhiE(PionsFromDecay[i].momentum_pt,PionsFromDecay[i].momentum_eta,PionsFromDecay[i].momentum_phi,PionsFromDecay[i].momentum_energy);
	      ptOfPions->Fill(PionsFromDecay[i].momentum_pt);
	      pOfPions->Fill(pi.P());
	    }

	  for(unsigned int i=0; i<Chi0FromDecay.size(); i++) 
	    {

	      TLorentzVector chi0;
	      chi0.SetPtEtaPhiE(Chi0FromDecay[i].momentum_pt,Chi0FromDecay[i].momentum_eta,Chi0FromDecay[i].momentum_phi,Chi0FromDecay[i].momentum_energy);
	      ptOfChi0->Fill(Chi0FromDecay[i].momentum_pt);
	      pOfChi0->Fill(chi0.P());
	    }


	  for(unsigned int i=0; i<AllPions.size(); i++) 
	    {

	      TLorentzVector pi;
	      pi.SetPtEtaPhiE(AllPions[i].momentum_pt,AllPions[i].momentum_eta,AllPions[i].momentum_phi,AllPions[i].momentum_energy);
	      ptOfAllPions->Fill(AllPions[i].momentum_pt);
	      pOfAllPions->Fill(pi.P());
	    }

  	   	    	
	}//end of loop over events
  
  cout<<"decayVertexFound  = "<<decayVertexFound<<endl;
  cout<<"decayVertexFound2 = "<<decayVertexFound2<<endl;
  cout<<"decayVertexFound3 = "<<decayVertexFound3<<endl;
  cout<<"decayVertexFound4 = "<<decayVertexFound4<<endl;

  stream.close();
  ofile.close();
  
  return 0;
}
