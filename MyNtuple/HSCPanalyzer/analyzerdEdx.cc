//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "analyzerRecoTracks.h"
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

  TH2D *hPDeDx, *hPDeDxChipm, *hPDeDxNoChipm;
  TH2D *hPDeDx7Hits, *hPDeDxChipm7Hits, *hPDeDxNoChipm7Hits;
  TH1D *hDeDx, *hDeDx7Hits, *hgenPChipm;

  hgenPChipm         = new TH1D("hgenPChipm","hgenPChipm",150,0,1500);
  hDeDx              = new TH1D("hDeDx","hDeDx",100,0,100);
  hPDeDx             = new TH2D("hPDeDx","hnPDeDx",150,0,1500,400,0,40);
  hPDeDxChipm        = new TH2D("hPDeDxChipm","hPDeDxChipm",150,0,1500,400,0,40);
  hPDeDxNoChipm      = new TH2D("hPDeDxNoChipm","hPDeDxNoChipm",150,0,1500,400,0,40);

  hDeDx7Hits         = new TH1D("hDeDx7Hits","hDeDx7Hits",100,0,100);
  hPDeDx7Hits        = new TH2D("hPDeDx7Hits","hPDeDx7Hits",150,0,1500,400,0,40);
  hPDeDxChipm7Hits   = new TH2D("hPDeDxChipm7Hits","hPDeDxChipm7Hits",150,0,1500,400,0,40);
  hPDeDxNoChipm7Hits = new TH2D("hPDeDxNoChipm7Hits","hPDeDxNoChipm7Hits",150,0,1500,400,0,40);

  //---------------------------------------------------------------------------
  // Declaration of Variables
  //---------------------------------------------------------------------------

  double dPhi    = 0;
  double dEta    = 0;
  double dRchip  = 0;
  double dRchim  = 0;
  int chipmFound = 0;
  int n = 10000 ;
  
  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry <n; ++entry)
	{
	  //cout<<endl<<"entry = "<<entry<<endl;
	  // Read event into memory
	  stream.read(entry);
	  // NB: call to clear object selection map (indexmap)
	  initialize();
	  
	  // Uncomment the following line if you wish to copy variables into
	  // structs. See the header file analyzer.h to find out what structs
	  // are available. Alternatively, you can call individual fill functions.
	  fillObjects();
	  

	  // ---------------------
	  // -- event selection --
	  // ---------------------
	  reco::findChipmInGenParticleCollection();

	  hgenPChipm -> Fill(std::abs(std::pow(reco::chipGenParticle.pt,2)+std::pow(reco::chipGenParticle.pz,2)));
	  hgenPChipm -> Fill(std::abs(std::pow(reco::chimGenParticle.pt,2)+std::pow(reco::chimGenParticle.pz,2)));

	  for(int i=0; i<nTrack; i++){

	    dPhi   = std::abs(reco::chipGenParticle.phi - Track[i].phi);
	    dEta   = std::abs(reco::chipGenParticle.eta - Track[i].eta);
	    dRchip = std::sqrt(pow(dPhi,2) + pow(dEta,2));
	    dPhi   = std::abs(reco::chimGenParticle.phi - Track[i].phi);
	    dEta   = std::abs(reco::chimGenParticle.eta - Track[i].eta);
	    dRchim = std::sqrt(pow(dPhi,2) + pow(dEta,2));

	    hDeDx      -> Fill(Track[i].dEdxNPHarm2);
	    hDeDx7Hits -> Fill(Track[i].dEdxHitsNPHarm2_7);

	    double p = std::sqrt(std::pow(Track[i].pt,2) + std::pow(Track[i].pz,2));
	    
	    if( (dRchip<0.1) | (dRchim <0.1) ){
	      chipmFound += 1;
	      hPDeDxChipm      -> Fill(p, Track[i].dEdxNPHarm2);
	      hPDeDxChipm7Hits -> Fill(p, Track[i].dEdxHitsNPHarm2_7);
	    }
	    else if((dRchip>0.5) && (dRchim>0.5)){
	      hPDeDxNoChipm      -> Fill(p, Track[i].dEdxNPHarm2);
	      hPDeDxNoChipm7Hits -> Fill(p, Track[i].dEdxHitsNPHarm2_7);
	    }
	  }
	  
	  
	  
	  // ---------------------
	  // -- fill histograms --
	  // ---------------------

	  for(int i=0; i<nTrack; i++){
	    double p = std::sqrt(std::pow(Track[i].pt,2) + std::pow(Track[i].pz,2));
	    hPDeDx      -> Fill(p, Track[i].dEdxNPHarm2);
	    hPDeDx7Hits -> Fill(p, Track[i].dEdxHitsNPHarm2_7);
	  }

	  //------------------------------------------------------------------------------------------------------------------------------------
	  //-------- GenParticle Helper  Collection --------------------------------------------------------------------------------------------
	  //------------------------------------------------------------------------------------------------------------------------------------

 	
	}//end of loop over events

  
  cout<<endl<<"chipmFound = "<<chipmFound<<endl<<endl;
  cout<<"matching efficiency = "<<chipmFound/(2.*n)<<endl<<endl;


  stream.close();
  ofile.close();
  
  return 0;
}
