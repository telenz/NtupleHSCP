//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "analyzerGen.h"
#include "analyzer.h"
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
  TString ChipmType[3];  //1. p+m, 1. p, 3. m
  ChipmType[0] = "Chipm";
  ChipmType[1] = "Chip";
  ChipmType[2] = "Chim";
 
  //GenParticle Collection
  TH1D *hgenPtChipm[3], *hgenEtaChipm[3], *hgenPhiChipm[3], *hgenMChipm[3];
  TH1D *hNumQuarks, *hNumLeptons, *hNumBosons;

  //SimTrack Collextion
  TH1D *htrkPtChipm[3], *htrkEtaChipm[3], *htrkPhiChipm[3], *htrkLengthChipm[3], *htrkLengthDeltaRay[3], *hDeltaRayPDGId[3], *htrkChargeChipm[3], *htrkPtOverPNotDecayed[3], *htrkPtOverPDecayed[3], *hTimeOfFlight[3];
  TH1D *hXChipm[3], *hYChipm[3], *hZChipm[3], *hBetaDecayed[3], *htrkEtaNotDecayed[3], *hgenEtaNotDecayed[3];
  TH2D *hRhoVsZChipm[3],*hZVsRhoChipm[3], *hZVsRhoChipmDeltaRay[3];
  TH2D *hgenEtaDecayed[3], *hgenPtDecayed[3], *hgenBetatDecayed[3],  *hgenBetaDecayed[3];


  for(int i=0; i<3; i++){

    //GenParticle Collection
    hgenPtChipm[i] = new TH1D("hgenPt"+ChipmType[i],"hgenPt"+ChipmType[i],150,0,1500);
    hgenEtaChipm[i] = new TH1D("hgenEta"+ChipmType[i],"hgenEta"+ChipmType[i],100,-5,5);
    hgenPhiChipm[i] = new TH1D("hgenPhi"+ChipmType[i],"hgenPhi"+ChipmType[i],50,0,TMath::Pi());
    hgenMChipm[i] = new TH1D("hgenM"+ChipmType[i],"hgenM"+ChipmType[i],1000,0,1000);
    hTimeOfFlight[i] = new TH1D("hTimeOfFlight"+ChipmType[i],"hTimeOfFlight"+ChipmType[i],1000,pow(10,-19),pow(10,-15));

    //SimTrack Collextion
    htrkPtChipm[i] = new TH1D("htrkPt"+ChipmType[i],"htrkPt"+ChipmType[i],150,0,1500);
    htrkEtaChipm[i] = new TH1D("htrkEta"+ChipmType[i],"htrkEta"+ChipmType[i],100,-5,5);
    htrkPhiChipm[i] = new TH1D("htrkPhi"+ChipmType[i],"htrkPhi"+ChipmType[i],50,0,TMath::Pi());
    htrkLengthChipm[i]    = new TH1D("htrkLength"+ChipmType[i],"htrkLength"+ChipmType[i],150,0,1500);
    htrkLengthDeltaRay[i]    = new TH1D("htrkLengthDeltaRay"+ChipmType[i],"htrkLengthDeltaRay"+ChipmType[i],150,0,1500);
    hDeltaRayPDGId[i]    = new TH1D("hDeltaRayPDGId"+ChipmType[i],"hDeltaRayPDGId"+ChipmType[i],1000,-500,500);
    htrkChargeChipm[i]    = new TH1D("htrkCharge"+ChipmType[i],"htrkCharge"+ChipmType[i],3,-1,2);
    hXChipm[i] = new TH1D("hX"+ChipmType[i],"hX"+ChipmType[i],140,0,1400);
    hYChipm[i] = new TH1D("hY"+ChipmType[i],"hY"+ChipmType[i],140,0,1400);
    hZChipm[i] = new TH1D("hZ"+ChipmType[i],"hZ"+ChipmType[i],140,0,1400);
    
    hRhoVsZChipm[i] = new TH2D("hRhoVsZ"+ChipmType[i],"hRhoVsZ"+ChipmType[i],140,0,1400,140,0,1400);
    hZVsRhoChipm[i] = new TH2D("hZVsRho"+ChipmType[i],"hZVsRho"+ChipmType[i],140,0,1400,80,0,800);
    hZVsRhoChipmDeltaRay[i] = new TH2D("hZVsRhoDeltaRay"+ChipmType[i],"hZVsRhoDeltaRay"+ChipmType[i],140,0,1400,80,0,800);
    hgenEtaDecayed[i]  = new TH2D("hgenEtaDecayed"+ChipmType[i],"hgenEtaDecayed"+ChipmType[i],50,-5,5,2,0,2);
    hgenPtDecayed[i]  = new TH2D("hgenPtDecayed"+ChipmType[i],"hgenPtDecayed"+ChipmType[i],150,0,1500,2,0,2);
    hBetaDecayed[i]     = new TH1D("hBetaDecayed"+ChipmType[i],"hBetaDecayed"+ChipmType[i],100,0,1);
    htrkEtaNotDecayed[i]     = new TH1D("htrkEtaNotDecayed"+ChipmType[i],"htrkEtaNotDecayed"+ChipmType[i],100,-5,5);
    htrkPtOverPNotDecayed[i]     = new TH1D("htrkPtOverPNotDecayed"+ChipmType[i],"htrkPtOverPNotDecayed"+ChipmType[i],100,0,1);
    htrkPtOverPDecayed[i]     = new TH1D("htrkPtOverPDecayed"+ChipmType[i],"htrkPtOverPDecayed"+ChipmType[i],100,0,1);
    hgenEtaNotDecayed[i]     = new TH1D("hgenEtaNotDecayed"+ChipmType[i],"hgenEtaNotDecayed"+ChipmType[i],100,-5,5);
    hgenBetatDecayed[i] = new TH2D("hgenBetatDecayed"+ChipmType[i],"hgenBetatDecayed"+ChipmType[i],100,0,1,2,0,2);
    hgenBetaDecayed[i]  = new TH2D("hgenBetaDecayed"+ChipmType[i],"hgenBetaDecayed"+ChipmType[i],100,0,1,2,0,2);
  }

  hNumQuarks = new TH1D("hNumQuarks","hNumQuarks",1000,0,1000);
  hNumLeptons = new TH1D("hNumLeptons","hNumLeptons",10,0,10);
  hNumBosons = new TH1D("hNumBosons","hNumBosons",10,0,10);

  //TH2D* hptVsDecayedChip   = new TH2D("hptVsDecayedChip","hptVsDecayedChip",100,0,1000,2,0,2);
  //TH2D* hptVsDecayedChim   = new TH2D("hptVsDecayedChim","hptVsDecayedChim",100,0,1000,2,0,2);
  //TH1D* hsurvivalProbabilityTrackChipm = new TH1D("hsurvivalProbabilityTrackChipm","hsurvivalProbabilityTrackChipm",100,0,1);

  //---------------------------------------------------------------------------
  // Declaration of Variables
  //---------------------------------------------------------------------------
  TLorentzVector trkChipVector, trkChimVector, genChipVector, genChimVector;

  int nQuarks  = 0;
  int nLeptons = 0;
  int nBosons  = 0;

  //cout<<"maximal Flight length = "<<sqrt(pow(rDetectorMax,2) + pow(lDetector,2))<<" m."<<endl;

  
  //---------------------------------------------------------------------------
  // Loop over events
  //---------------------------------------------------------------------------

  for(int entry=0; entry <nevents; ++entry)
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

	  
	  
	  // ---------------------
	  // -- fill histograms --
	  // ---------------------

	  nQuarks  = 0;
	  nLeptons = 0;
	  nBosons  = 0;

	  //------------------------------------------------------------------------------------------------------------------------------------
	  //-------- GenParticle Helper  Collection --------------------------------------------------------------------------------------------
	  //------------------------------------------------------------------------------------------------------------------------------------

	  findChipmInGenParticleCollection();
	  

	  hgenPtChipm[0]  -> Fill(chipGenParticle.pt);
	  hgenEtaChipm[0] -> Fill(chipGenParticle.eta);
	  hgenPhiChipm[0] -> Fill(chipGenParticle.phi);
	  hgenMChipm[0]   -> Fill(chipGenParticle.mass);
	  hgenPtChipm[0]  -> Fill(chimGenParticle.pt);
	  hgenEtaChipm[0] -> Fill(chimGenParticle.eta);
	  hgenPhiChipm[0] -> Fill(chimGenParticle.phi);
	  hgenMChipm[0]   -> Fill(chimGenParticle.mass);

	  genChipVector.SetPtEtaPhiM(chipGenParticle.pt,chipGenParticle.eta,chipGenParticle.phi,chipGenParticle.mass);		
	  hgenPtChipm[1]  -> Fill(chipGenParticle.pt);
	  hgenEtaChipm[1] -> Fill(chipGenParticle.eta);
	  hgenPhiChipm[1] -> Fill(chipGenParticle.phi);
	  hgenMChipm[1]   -> Fill(chipGenParticle.mass);
	  
	  genChimVector.SetPtEtaPhiM(chimGenParticle.pt,chimGenParticle.eta,chimGenParticle.phi,chimGenParticle.mass);	
	  hgenPtChipm[2]  -> Fill(chimGenParticle.pt);
	  hgenEtaChipm[2] -> Fill(chimGenParticle.eta);
	  hgenPhiChipm[2] -> Fill(chimGenParticle.phi);
	  hgenMChipm[2]   -> Fill(chimGenParticle.mass);
	  
	  //if(abs(GenParticle_pdgId[i])<=6) nQuarks += 1;
	  //else if(abs(GenParticle_pdgId[i])>10 && abs(GenParticle_pdgId[i])<=16) nLeptons += 1;
	  //else if(abs(GenParticle_pdgId[i])>20 && abs(GenParticle_pdgId[i])<=38) nBosons += 1;
	  
	  //hNumQuarks->Fill(nQuarks);
	  //hNumLeptons->Fill(nLeptons);
	  //hNumBosons->Fill(nBosons);
	  
  
	  //--------------------------------------------------------------------------------------------------------------------------
	  //-------- SimTrack Collection ---------------------------------------------------------------------------------------------
	  //--------------------------------------------------------------------------------------------------------------------------
	  findChipmInSimTrackCollection();
	  findChi0InSimTrackCollection();
	  
	  trkChipVector.SetPtEtaPhiE(chipSimTrack.momentum_pt,chipSimTrack.momentum_eta,chipSimTrack.momentum_phi,chipSimTrack.momentum_energy);
	  trkChimVector.SetPtEtaPhiE(chimSimTrack.momentum_pt,chimSimTrack.momentum_eta,chimSimTrack.momentum_phi,chimSimTrack.momentum_energy);
	  htrkPtChipm[0]     -> Fill(chipSimTrack.momentum_pt);
	  htrkEtaChipm[0]    -> Fill(chipSimTrack.momentum_eta);
	  htrkPhiChipm[0]    -> Fill(chipSimTrack.momentum_phi);
	  htrkChargeChipm[0] -> Fill(chipSimTrack.charge);
	  htrkPtChipm[0]     -> Fill(chimSimTrack.momentum_pt);
	  htrkEtaChipm[0]    -> Fill(chimSimTrack.momentum_eta);
	  htrkPhiChipm[0]    -> Fill(chimSimTrack.momentum_phi);
	  htrkChargeChipm[0] -> Fill(chimSimTrack.charge);
	      
	  if(chipSimTrack.noVertex==1) cout<<"no origin vertex found for chip !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	  if(chimSimTrack.noVertex==1) cout<<"no origin vertex found for chim !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

	  htrkPtChipm[1]     -> Fill(chipSimTrack.momentum_pt);
	  htrkEtaChipm[1]    -> Fill(chipSimTrack.momentum_eta);
	  htrkPhiChipm[1]    -> Fill(chipSimTrack.momentum_phi);
	  htrkChargeChipm[1] -> Fill(chipSimTrack.charge);
		
	  htrkPtChipm[2]     -> Fill(chimSimTrack.momentum_pt);
	  htrkEtaChipm[2]    -> Fill(chimSimTrack.momentum_eta);
	  htrkPhiChipm[2]    -> Fill(chimSimTrack.momentum_phi);
	  htrkChargeChipm[2] -> Fill(chimSimTrack.charge);
	  	  
	  //-----------------------------------------------------------------------------------------------------
	  //--------- SimVertex Collection ----------------------------------------------------------------------
	  //-----------------------------------------------------------------------------------------------------

	  findChipmOriginInSimVertexCollection();
	  findChipmDecayInSimVertexCollection();

	  double Rho;
	  double TransverseTrackLength;
	  double TrackLength, TrackLengthChip, TrackLengthChim;  

	  if(decayedChip){

	    Rho                   = sqrt(pow(chipDecaySimVertex.position_x,2)+pow(chipDecaySimVertex.position_y,2));
	    TransverseTrackLength = sqrt(pow(chipDecaySimVertex.position_x-chipOriginSimVertex.position_x,2)+pow(chipDecaySimVertex.position_y-chipOriginSimVertex.position_y,2));
	    TrackLength           = sqrt(pow(TransverseTrackLength,2)+pow(chipDecaySimVertex.position_z-chipOriginSimVertex.position_z,2));

	    htrkLengthChipm[1] -> Fill(TrackLength);
	    hRhoVsZChipm[1]    -> Fill(Rho,abs(chipDecaySimVertex.position_z));
	    hZVsRhoChipm[1]    -> Fill(abs(chipDecaySimVertex.position_z),Rho);
	    hXChipm[1]         -> Fill(abs(chipDecaySimVertex.position_x));
	    hYChipm[1]         -> Fill(abs(chipDecaySimVertex.position_y));
	    hZChipm[1]         -> Fill(abs(chipDecaySimVertex.position_z));
	    htrkLengthChipm[0] -> Fill(TrackLength);
	    hRhoVsZChipm[0]    -> Fill(Rho,abs(chipDecaySimVertex.position_z));
	    hZVsRhoChipm[0]    -> Fill(abs(chipDecaySimVertex.position_z),Rho);
	    hXChipm[0]         -> Fill(abs(chipDecaySimVertex.position_x));
	    hYChipm[0]         -> Fill(abs(chipDecaySimVertex.position_y));
	    hZChipm[0]         -> Fill(abs(chipDecaySimVertex.position_z));

	    TrackLengthChip = TrackLength;
	  }

	  if(decayedChim){
	    Rho                   = sqrt(pow(chimDecaySimVertex.position_x,2)+pow(chimDecaySimVertex.position_y,2));
	    TransverseTrackLength = sqrt(pow(chimDecaySimVertex.position_x-chimOriginSimVertex.position_x,2)+pow(chimDecaySimVertex.position_y-chimOriginSimVertex.position_y,2));
	    TrackLength           = sqrt(pow(TransverseTrackLength,2)+pow(chimDecaySimVertex.position_z-chimOriginSimVertex.position_z,2));

	    htrkLengthChipm[2] -> Fill(TrackLength);
	    hRhoVsZChipm[2]    -> Fill(Rho,abs(chimDecaySimVertex.position_z));
	    hZVsRhoChipm[2]    -> Fill(abs(chimDecaySimVertex.position_z),Rho);
	    hXChipm[2]         -> Fill(abs(chimDecaySimVertex.position_x));
	    hYChipm[2]         -> Fill(abs(chimDecaySimVertex.position_y));
	    hZChipm[2]         -> Fill(abs(chimDecaySimVertex.position_z));
	    htrkLengthChipm[0] -> Fill(TrackLength);
	    hRhoVsZChipm[0]    -> Fill(Rho,abs(chimDecaySimVertex.position_z));
	    hZVsRhoChipm[0]    -> Fill(abs(chimDecaySimVertex.position_z),Rho);
	    hXChipm[0]         -> Fill(abs(chimDecaySimVertex.position_x));
	    hYChipm[0]         -> Fill(abs(chimDecaySimVertex.position_y));
	    hZChipm[0]         -> Fill(abs(chimDecaySimVertex.position_z));

	    TrackLengthChim = TrackLength;
	  }
	  	   	  
	  //-----------------------------------------------------------------------------------------------------
	  //--------- SimTrack Collection ----------------------------------------------------------------------
	  //-----------------------------------------------------------------------------------------------------


	  double flightTime =   trkChipVector.T();
	  flightTime = TrackLength*trkChipVector.E()/trkChipVector.P();
	  flightTime = 0.0;
	  flightTime = TrackLengthChip/100*trkChipVector.E()/(trkChipVector.P()*pow(TMath::C(),2));

	  if(!decayedChip) {
	    htrkEtaNotDecayed[0]     -> Fill(chipSimTrack.momentum_eta);
	    htrkEtaNotDecayed[1]     -> Fill(chipSimTrack.momentum_eta);
	    htrkPtOverPNotDecayed[0] -> Fill(chipSimTrack.momentum_pt/trkChipVector.P());
	    htrkPtOverPNotDecayed[1] -> Fill(chipSimTrack.momentum_pt/trkChipVector.P());
	    hgenEtaNotDecayed[0]     -> Fill(chipGenParticle.eta);
	    hgenEtaNotDecayed[1]     -> Fill(chipGenParticle.eta);
	  }
	  else{
	    htrkPtOverPDecayed[0]    -> Fill(chipSimTrack.momentum_pt/trkChipVector.P());
	    htrkPtOverPDecayed[1]    -> Fill(chipSimTrack.momentum_pt/trkChipVector.P());
	    hTimeOfFlight[0]         -> Fill(flightTime);
	    hTimeOfFlight[1]         -> Fill(flightTime);
	    hBetaDecayed[0]          -> Fill(chipGenParticle.pt/genChipVector.E());
	    hBetaDecayed[1]          -> Fill(chipGenParticle.pt/genChipVector.E());
	  }

	  
	  flightTime = TrackLengthChim/100*trkChimVector.E()/(trkChimVector.P()*pow(TMath::C(),2));
	  if(!decayedChim) {
	    htrkEtaNotDecayed[0]     -> Fill(chimSimTrack.momentum_eta);
	    htrkEtaNotDecayed[2]     -> Fill(chimSimTrack.momentum_eta);
	    htrkPtOverPNotDecayed[0] -> Fill(chimSimTrack.momentum_pt/trkChimVector.P());
	    htrkPtOverPNotDecayed[2] -> Fill(chimSimTrack.momentum_pt/trkChimVector.P());
	    hgenEtaNotDecayed[0]     -> Fill(chimGenParticle.eta);
	    hgenEtaNotDecayed[2]     -> Fill(chimGenParticle.eta);
	  }
	  else{
	    htrkPtOverPDecayed[0]    -> Fill(chimSimTrack.momentum_pt/trkChimVector.P());
	    htrkPtOverPDecayed[2]    -> Fill(chimSimTrack.momentum_pt/trkChimVector.P());
	    hTimeOfFlight[0]         -> Fill(flightTime);
	    hTimeOfFlight[2]         -> Fill(flightTime);
	    hBetaDecayed[0]          -> Fill(chimGenParticle.pt/genChimVector.E());
	    hBetaDecayed[2]          -> Fill(chimGenParticle.pt/genChimVector.E());
	  }



	  for(int i=0; i<evt::nSimVertex; i++){

	    if(evt::SimVertex[i].parentIndex==chipSimTrack.trackId){

	      if(evt::SimVertex[i].vertexId == chi0SimTrack[0].vertIndex || evt::SimVertex[i].vertexId == chi0SimTrack[1].vertIndex) continue;
	   
	      Rho                   = sqrt(pow(SimVertex[i].position_x,2)+pow(SimVertex[i].position_y,2));
	      TransverseTrackLength = sqrt(pow(SimVertex[i].position_x-chipOriginSimVertex.position_x,2)+pow(SimVertex[i].position_y-chipOriginSimVertex.position_y,2));
	      TrackLength           = sqrt(pow(TransverseTrackLength,2)+pow(SimVertex[i].position_z-chipOriginSimVertex.position_z,2));
	      
	      htrkLengthDeltaRay[0] -> Fill(TrackLength);
	      hDeltaRayPDGId[0]     -> Fill(SimTrack[i].type);
	      htrkLengthDeltaRay[1] -> Fill(TrackLength);
	      hDeltaRayPDGId[1]     -> Fill(SimTrack[i].type);
	    }
	    else if(evt::SimVertex[i].parentIndex==chimSimTrack.trackId){

	      if(evt::SimVertex[i].vertexId == chi0SimTrack[0].vertIndex || evt::SimVertex[i].vertexId == chi0SimTrack[1].vertIndex) continue;
	   
	      Rho                   = sqrt(pow(SimVertex[i].position_x,2)+pow(SimVertex[i].position_y,2));
	      TransverseTrackLength = sqrt(pow(SimVertex[i].position_x-chimOriginSimVertex.position_x,2)+pow(SimVertex[i].position_y-chimOriginSimVertex.position_y,2));
	      TrackLength           = sqrt(pow(TransverseTrackLength,2)+pow(SimVertex[i].position_z-chimOriginSimVertex.position_z,2));
	      
	      htrkLengthDeltaRay[0] -> Fill(TrackLength);
	      hDeltaRayPDGId[0]     -> Fill(SimTrack[i].type);
	      htrkLengthDeltaRay[2] -> Fill(TrackLength);
	      hDeltaRayPDGId[2]     -> Fill(SimTrack[i].type);
	    }
	  }	  

	  hgenEtaDecayed[1] -> Fill(chipGenParticle.eta,decayedChip);
	  hgenEtaDecayed[0] -> Fill(chipGenParticle.eta,decayedChip);
	  hgenPtDecayed[1]  -> Fill(chipGenParticle.pt,decayedChip);
	  hgenPtDecayed[0]  -> Fill(chipGenParticle.pt,decayedChip);
	  hgenEtaDecayed[2] -> Fill(chimGenParticle.eta,decayedChim);
	  hgenEtaDecayed[0] -> Fill(chimGenParticle.eta,decayedChim);
	  hgenPtDecayed[2]  -> Fill(chimGenParticle.pt,decayedChim);
	  hgenPtDecayed[0]  -> Fill(chimGenParticle.pt,decayedChim);
  	
	}//end of loop over events

  cout<<endl<<endl<<htrkLengthChipm[1]->Integral()<<" charginos+ decayed."<<endl;
  cout<<htrkLengthChipm[2]->Integral()<<" charginos- decayed."<<endl;
  cout<<hgenMChipm[0]->Integral()<<" charginos existed."<<endl;
  cout<<htrkPtChipm[0]->Integral()<<" chargino tracks were found."<<endl;
  cout<<(htrkLengthChipm[1]->Integral()+htrkLengthChipm[2]->Integral())/hgenMChipm[0]->Integral()*100.<<"% of all charginos decayed."<<endl;
  cout<<(hgenMChipm[0]->Integral()-htrkPtChipm[0]->Integral())/hgenMChipm[0]->Integral()*100.<<"% of generated charginos couldn't be associated with a track."<<endl;
  cout<<"mean decay length Chip = "<<htrkLengthChipm[1]->GetMean()<<" cm."<<endl;
  cout<<"mean decay length Chim = "<<htrkLengthChipm[2]->GetMean()<<" cm."<<endl<<endl<<endl;
  cout<<endl<<endl;
  
  stream.close();
  ofile.close();
  
  return 0;
}
