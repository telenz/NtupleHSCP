//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
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
  TH1D *htrkPtChipm[3], *htrkEtaChipm[3], *htrkPhiChipm[3], *htrkLengthChipm[3], *htrkLengthDeltaRay[3], *hDeltaRayPDGId[3], *htrkChargeChipm[3];
  TH1D *hZChipm[3], *hBetaDecayed[3], *htrkEtaNotDecayed[3], *hgenEtaNotDecayed[3];
  TH2D *hRhoVsZChipm[3],*hZVsRhoChipm[3];
  TH2D *hgenEtaDecayed[3], *hgenBetatDecayed[3],  *hgenBetaDecayed[3], *hSimHitsDecayed[3];


  for(int i=0; i<3; i++){

    //GenParticle Collection
    hgenPtChipm[i] = new TH1D("hgenPt"+ChipmType[i],"hgenPt"+ChipmType[i],100,0,1000);
    hgenEtaChipm[i] = new TH1D("hgenEta"+ChipmType[i],"hgenEta"+ChipmType[i],100,-5,5);
    hgenPhiChipm[i] = new TH1D("hgenPhi"+ChipmType[i],"hgenPhi"+ChipmType[i],50,0,TMath::Pi());
    hgenMChipm[i] = new TH1D("hgenM"+ChipmType[i],"hgenM"+ChipmType[i],1000,0,1000);

    //SimTrack Collextion
    htrkPtChipm[i] = new TH1D("htrkPt"+ChipmType[i],"htrkPt"+ChipmType[i],100,0,1000);
    htrkEtaChipm[i] = new TH1D("htrkEta"+ChipmType[i],"htrkEta"+ChipmType[i],100,-5,5);
    htrkPhiChipm[i] = new TH1D("htrkPhi"+ChipmType[i],"htrkPhi"+ChipmType[i],50,0,TMath::Pi());
    htrkLengthChipm[i]    = new TH1D("htrkLength"+ChipmType[i],"htrkLength"+ChipmType[i],150,0,1500);
    htrkLengthDeltaRay[i]    = new TH1D("htrkLengthDeltaRay"+ChipmType[i],"htrkLengthDeltaRay"+ChipmType[i],150,0,1500);
    hDeltaRayPDGId[i]    = new TH1D("hDeltaRayPDGId"+ChipmType[i],"hDeltaRayPDGId"+ChipmType[i],1000,-500,500);
    htrkChargeChipm[i]    = new TH1D("htrkCharge"+ChipmType[i],"htrkCharge"+ChipmType[i],3,-1,2);
    hZChipm[i] = new TH1D("hZ"+ChipmType[i],"hZ"+ChipmType[i],140,0,1400);
    
    hRhoVsZChipm[i] = new TH2D("hRhoVsZ"+ChipmType[i],"hRhoVsZ"+ChipmType[i],140,0,1400,140,0,1400);
    hZVsRhoChipm[i] = new TH2D("hZVsRho"+ChipmType[i],"hZVsRho"+ChipmType[i],140,0,1400,80,0,800);
    hgenEtaDecayed[i]  = new TH2D("hgenEtaDecayed"+ChipmType[i],"hgenEtaDecayed"+ChipmType[i],100,-5,5,2,0,2);
    hBetaDecayed[i]     = new TH1D("hBetaDecayed"+ChipmType[i],"hBetaDecayed"+ChipmType[i],100,0,1);
    htrkEtaNotDecayed[i]     = new TH1D("htrkEtaNotDecayed"+ChipmType[i],"htrkEtaNotDecayed"+ChipmType[i],100,-5,5);
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
  double vertex;
  double SumFlightLengthChip = 0;
  double SumFlightLengthChim = 0;
 
  bool decayedChim = false;
  bool decayedChip = false;
 

  TLorentzVector trkChipmVector;

  TLorentzVector genChipmVector[3];

  double rDetectorMax = 1.0;
  double rDetectorMin = 0.045;
  double lDetector = 2.8;
  double teta = 0;
  double sDecay = 0;
  double sDecayMin = 0;
  double sDecayMax = 0;
  double beta = 0;
  double gamma = 0;
  double hbar = 6.58212*pow(10,-24);    //[hbar] = GeV*s
  double decayWidth= 6.*pow(10,-15);    //[decayWidth] = GeV
  
  double maxFlightLength = 0;
  double minFlightLength = 0;
  bool breakBool = true;


  int nQuarks = 0;
  int nLeptons = 0;
  int nBosons = 0;

  int countCharg = 0;
  int countNvertex = 0;
  int vertexChip = 0 ;
  int vertexChim = 0;

  int allright = 0;

  cout<<"maximal Flight length = "<<sqrt(pow(rDetectorMax,2) + pow(lDetector,2))<<" m."<<endl;

  
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

	  int aux = 0;
	  nQuarks = 0;
	  nLeptons = 0;
	  nBosons = 0;
	  decayedChip = false;
	  decayedChim = false;

	  //-----------------------------------------------------------------------------------------------------------------------------
	  //-------- GenParticle Helper  Collection ---------------------------------------------------------------------------------------------
	  
	  for(int i=0; i<nGenParticle; i++){

	    if(abs(GenParticle_pdgId[i])==1000024){

	      hgenPtChipm[0]->Fill(GenParticle_pt[i]);
	      hgenEtaChipm[0]->Fill(GenParticle_eta[i]);
	      hgenPhiChipm[0]->Fill(GenParticle_phi[i]);
	      hgenMChipm[0]->Fill(GenParticle_mass[i]);

	      aux += 1;	      
	      if(GenParticle_pdgId[i]>0){
		
		genChipmVector[1].SetPtEtaPhiM(GenParticle_pt[i],GenParticle_eta[i],GenParticle_phi[i],GenParticle_mass[i]);		
		hgenPtChipm[1]->Fill(GenParticle_pt[i]);
		hgenEtaChipm[1]->Fill(GenParticle_eta[i]);
		hgenPhiChipm[1]->Fill(GenParticle_phi[i]);
		hgenMChipm[1]->Fill(GenParticle_mass[i]);
	      }
	      else if(GenParticle_pdgId[i]<0){
		
		//genChipmVector[2].SetPtEtaPhiM(GenParticle_pt[i],GenParticle_eta[i],GenParticle_phi[i],GenParticle_mass[i]);	
		//hgenPtChipm[2]->Fill(GenParticle_pt[i]);
		//hgenEtaChipm[2]->Fill(GenParticle_eta[i]);
		//hgenPhiChipm[2]->Fill(GenParticle_phi[i]);
		//hgenMChipm[2]->Fill(GenParticle_mass[i]);
				
	      }	      
	    }
	    
	    if(abs(GenParticle_pdgId[i])<=6) nQuarks += 1;
	    else if(abs(GenParticle_pdgId[i])>10 && abs(GenParticle_pdgId[i])<=16) nLeptons += 1;
	    else if(abs(GenParticle_pdgId[i])>20 && abs(GenParticle_pdgId[i])<=38) nBosons += 1;
	    if(aux>4) cout<<"There are more than two charginos in the sample!!!!!"<<endl;
	    
	  }
	  
	  //-----------------------------------------------------------------------------------------------------------------------------
	  //-----------------------------------------------------------------------------------------------------------------------------
	  
	  
	  hNumQuarks->Fill(nQuarks);
	  hNumLeptons->Fill(nLeptons);
	  hNumBosons->Fill(nBosons);

	  int trackidxChipm[3]  = {-1,-1,-1};
	  int trackidxChi0[3]   = {-1,-1,-1};
	  int VertexidxChip[3]  = {-1,-1,-1};
	  int VertexidxChim[3]  = {-1,-1,-1};
	  int ChipmOriginVtx[3] = {-1,-1,-1};
	  int Chi0OriginVtx[3]  = {-1,-1,-1};
	  int ChipmDecayVtx[3]  = {-1,-1,-1};
	  int Chi0DecayVtx[3]   = {-1,-1,-1};

	  int ChipmOriginVtx_X[3] = {-1,-1,-1};
	  int ChipmOriginVtx_Y[3] = {-1,-1,-1};
	  int ChipmOriginVtx_Z[3] = {-1,-1,-1};
	  int ChipmOriginVtx_T[3] = {-1,-1,-1};

	  int ChipmDecayVertex_X[3] = {-1,-1,-1};
	  int ChipmDecayVertex_Y[3] = {-1,-1,-1};
	  int ChipmDecayVertex_Z[3] = {-1,-1,-1};
	  int ChipmDecayVertex_T[3] = {-1,-1,-1};


	  bool firstNeutralino = false;
	  bool secondNeutralino = false;

	  double TrackLength = -100000000000000000;
	  double TransverseTrackLength = -100000000000000000;
	  double Rho = -100000000000000000;

	  int originVtxNotSame = false;

	  aux = 0;
	  
	  //-------- SimTrack Collection ---------------------------------------------------------------------------------------------
	  
	  for(int i=0; i<nSimTrack; i++){

	    //cout<<"SimTrack_type[i] = "<<SimTrack_type[i]<<endl;

	    if(abs(SimTrack_type[i])==1000024){
	      //cout<<"chargino found!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	      trkChipmVector.SetPtEtaPhiE(SimTrack_momentum_pt[i],SimTrack_momentum_eta[i],SimTrack_momentum_phi[i],SimTrack_momentum_energy[i]);
	      htrkPtChipm[0]->Fill(SimTrack_momentum_pt[i]);
	      htrkEtaChipm[0]->Fill(SimTrack_momentum_eta[i]);
	      htrkPhiChipm[0]->Fill(SimTrack_momentum_phi[i]);
	      htrkChargeChipm[0]->Fill(SimTrack_charge[i]);
	      
	      if(SimTrack_noVertex[i]==1) cout<<"no origin vertex found!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

	      if(SimTrack_type[i]>0){

		aux += 1;
		trackidxChipm[1]  = SimTrack_trackId[i];
		ChipmOriginVtx[1] = SimTrack_vertIndex[i];

		htrkPtChipm[1]->Fill(SimTrack_momentum_pt[i]);
		htrkEtaChipm[1]->Fill(SimTrack_momentum_eta[i]);
		htrkPhiChipm[1]->Fill(SimTrack_momentum_phi[i]);
		htrkChargeChipm[1]->Fill(SimTrack_charge[i]);
		
	      }
	      else if(SimTrack_type[i]<0){

		aux += 1;
		trackidxChipm[2]  = SimTrack_trackId[i];
		ChipmOriginVtx[2] = SimTrack_vertIndex[i];
		
		
		htrkPtChipm[2]->Fill(SimTrack_momentum_pt[i]);
		htrkEtaChipm[2]->Fill(SimTrack_momentum_eta[i]);
		htrkPhiChipm[2]->Fill(SimTrack_momentum_phi[i]);
		htrkChargeChipm[2]->Fill(SimTrack_charge[i]);
		
	      }

	      if(aux>2) cout<<"Too many Chipm tracks!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	    }
	    if(abs(SimTrack_type[i])==1000022){
	      //cout<<"neutralino found!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	      if(!firstNeutralino){
		trackidxChi0[1] == SimTrack_trackId[i];
		Chi0OriginVtx[1] = SimTrack_vertIndex[i];
		firstNeutralino = true;
	      }
	      else if(firstNeutralino && !secondNeutralino){
		trackidxChi0[2] == SimTrack_trackId[i];
		Chi0OriginVtx[2] = SimTrack_vertIndex[i];
		secondNeutralino = true;
	      }
	      else cout<<"There are more than 2 neutralinos in the sample!!!!!!!!!!!!!!!!!"<<endl;
	    }
	  }
	  
	  aux=0;
	  //-----------------------------------------------------------------------------------------------------
	  
	  vertexChip = 0;
	  vertexChim = 0;
	  
	  //--------- SimVertex Collection ----------------------------------------------------------------------
	  
	  for(int i=0; i<nSimVertex; i++){

	    // Set x,y,z of ORIGIN vertex
	    
	    if(SimVertex_vertexId[i]==ChipmOriginVtx[1]){
	      //cout<<"origin vertex id + = "<<SimVertex_vertexId[i]<<endl;
	      ChipmOriginVtx_X[1] = SimVertex_position_x[i];
	      ChipmOriginVtx_Y[1] = SimVertex_position_y[i];
	      ChipmOriginVtx_Z[1] = SimVertex_position_z[i];
	      ChipmOriginVtx_T[1] = SimVertex_position_t[i];
	    }
	    if(SimVertex_vertexId[i]==ChipmOriginVtx[2]){
	      //cout<<"origin vertex id - = "<<SimVertex_vertexId[i]<<endl;
	      ChipmOriginVtx_X[2] = SimVertex_position_x[i];
	      ChipmOriginVtx_Y[2] = SimVertex_position_y[i];
	      ChipmOriginVtx_Z[2] = SimVertex_position_z[i];
	      ChipmOriginVtx_T[2] = SimVertex_position_t[i];      
	    }
	     
	    
	    int idxChipm = -1000;
	    
	    if(SimVertex_parentIndex[i]==trackidxChipm[1] || SimVertex_parentIndex[i]==trackidxChipm[2]){
	      
	      if(SimVertex_parentIndex[i]==trackidxChipm[1]) idxChipm = 1;
	      else if(SimVertex_parentIndex[i]==trackidxChipm[2]) idxChipm = 2;
	      else{
		cout<<"Something wrong in trackidx of Chipm!!!!!!!!!!!!"<<endl;
		return 1;
	      }

	      if(trackidxChipm[idxChipm]<0) continue;

	      if(idxChipm==1){
		VertexidxChip[vertexChip] = SimVertex_vertexId[i];
		vertexChip += 1;
		if(SimVertex_vertexId[i] == Chi0OriginVtx[1] || SimVertex_vertexId[i] == Chi0OriginVtx[2]) allright +=1;
		else continue;
		decayedChip = true;
		
	      }
	      else if(idxChipm==2){
		VertexidxChim[vertexChim] = SimVertex_vertexId[i];
		vertexChim += 1;
		if(SimVertex_vertexId[i] == Chi0OriginVtx[1] || SimVertex_vertexId[i] == Chi0OriginVtx[2]) allright +=1;
		else continue;
		decayedChim = true;	
	      }	    
	      ChipmDecayVertex_X[idxChipm] = SimVertex_position_x[i];
	      ChipmDecayVertex_Y[idxChipm] = SimVertex_position_y[i];
	      ChipmDecayVertex_Z[idxChipm] = SimVertex_position_z[i];
	      ChipmDecayVertex_T[idxChipm] = SimVertex_position_t[i];

	      Rho = sqrt(pow(ChipmDecayVertex_X[idxChipm],2)+pow(ChipmDecayVertex_Y[idxChipm],2));
	      TransverseTrackLength = sqrt(pow(ChipmDecayVertex_X[idxChipm]-ChipmOriginVtx_X[idxChipm],2)+pow(ChipmDecayVertex_Y[idxChipm]-ChipmOriginVtx_Y[idxChipm],2));
	      TrackLength = sqrt(pow(TransverseTrackLength,2)+pow(ChipmDecayVertex_Z[idxChipm]-ChipmOriginVtx_Z[idxChipm],2));
	     
	      ChipmDecayVertex_X[idxChipm] = SimVertex_position_x[i];
	      ChipmDecayVertex_Y[idxChipm] = SimVertex_position_y[i];
	      ChipmDecayVertex_Z[idxChipm] = SimVertex_position_z[i];
	      ChipmDecayVertex_T[idxChipm] = SimVertex_position_t[i];
          
	      htrkLengthChipm[idxChipm]->Fill(TrackLength);
	      hRhoVsZChipm[idxChipm]->Fill(Rho,ChipmDecayVertex_Z[idxChipm]);
	      hZVsRhoChipm[idxChipm]->Fill(ChipmDecayVertex_Z[idxChipm],Rho);
	      hZChipm[idxChipm]->Fill(ChipmDecayVertex_Z[idxChipm]);

	      SumFlightLengthChim += TransverseTrackLength;
   
	      htrkLengthChipm[0]->Fill(TrackLength);
	      hRhoVsZChipm[0]->Fill(Rho,ChipmDecayVertex_Z[idxChipm]);
	      hZVsRhoChipm[0]->Fill(ChipmDecayVertex_Z[idxChipm],Rho);
	      hZChipm[0]->Fill(ChipmDecayVertex_Z[idxChipm]);

	      aux+=1;
	    } 	    	
	  } // end of vertex collection loop
	  
	  
	  //cout<<endl<<"vertexChim = "<<vertexChim<<endl;	  
	  if(!decayedChim){
	    //cout<<"entry = "<<entry<<endl;
	    //cout<<"not decayed Chim"<<endl;
	  }
	  if(!decayedChip){
	    //cout<<"entry = "<<entry<<endl;
	    //cout<<"not decayed Chip"<<endl;
	  }
	  
	  
	    for(int i=0; i<nSimTrack; i++){

	    if(SimTrack_trackId[i] == trackidxChipm[1]){
	    if(!decayedChip) {
	    htrkEtaNotDecayed[0]->Fill(SimTrack_momentum_eta[i]);
		htrkEtaNotDecayed[1]->Fill(SimTrack_momentum_eta[i]);
	      }

	    }
	    if(SimTrack_trackId[i] == trackidxChipm[2]){
	      if(!decayedChim) {
		htrkEtaNotDecayed[0]->Fill(SimTrack_momentum_eta[i]);
		htrkEtaNotDecayed[2]->Fill(SimTrack_momentum_eta[i]);
	      }
	    }
	    for(int j=0; j<vertexChim; j++){
	      
	      if(SimTrack_vertIndex[i]==VertexidxChim[j]) {
		
		ChipmDecayVertex_X[2] = SimVertex_position_x[VertexidxChim[j]];
		ChipmDecayVertex_Y[2] = SimVertex_position_y[VertexidxChim[j]];
		ChipmDecayVertex_Z[2] = SimVertex_position_z[VertexidxChim[j]];
		ChipmDecayVertex_T[2] = SimVertex_position_t[VertexidxChim[j]];
		
		Rho = sqrt(pow(ChipmDecayVertex_X[2],2)+pow(ChipmDecayVertex_Y[2],2));
		TransverseTrackLength = sqrt(pow(ChipmDecayVertex_X[2]-ChipmOriginVtx_X[2],2)+pow(ChipmDecayVertex_Y[2]-ChipmOriginVtx_Y[2],2));
		TrackLength = sqrt(pow(TransverseTrackLength,2)+pow(ChipmDecayVertex_Z[2]-ChipmOriginVtx_Z[2],2));
		
		htrkLengthDeltaRay[0]->Fill(TrackLength);
		htrkLengthDeltaRay[2]->Fill(TrackLength);

		hDeltaRayPDGId[0]->Fill(SimTrack_type[i]);
		hDeltaRayPDGId[2]->Fill(SimTrack_type[i]);
	      }
	    }
	    
	    for(int j=0; j<vertexChip; j++){
	      
	      if(SimTrack_vertIndex[i]==VertexidxChip[j] ){
		
		ChipmDecayVertex_X[1] = SimVertex_position_x[VertexidxChip[j]];
		ChipmDecayVertex_Y[1] = SimVertex_position_y[VertexidxChip[j]];
		ChipmDecayVertex_Z[1] = SimVertex_position_z[VertexidxChip[j]];
		ChipmDecayVertex_T[1] = SimVertex_position_t[VertexidxChip[j]];
		
		Rho = sqrt(pow(ChipmDecayVertex_X[1],2)+pow(ChipmDecayVertex_Y[1],2));
		TransverseTrackLength = sqrt(pow(ChipmDecayVertex_X[1]-ChipmOriginVtx_X[1],2)+pow(ChipmDecayVertex_Y[1]-ChipmOriginVtx_Y[1],2));
		TrackLength = sqrt(pow(TransverseTrackLength,2)+pow(ChipmDecayVertex_Z[1]-ChipmOriginVtx_Z[1],2));
		
		htrkLengthDeltaRay[0]->Fill(TrackLength);
		htrkLengthDeltaRay[1]->Fill(TrackLength);

		hDeltaRayPDGId[0]->Fill(SimTrack_type[i]);
		hDeltaRayPDGId[1]->Fill(SimTrack_type[i]);
	      }	       
	    }
	  }
	  

	  
	  for(int i=0; i<nGenParticle; i++){
	    if(abs(GenParticle_pdgId[i])==1000024){
	      if(GenParticle_pdgId[i]>0){
		hgenEtaDecayed[1]->Fill(GenParticle_eta[i],decayedChip);
		hgenEtaDecayed[0]->Fill(GenParticle_eta[i],decayedChip);
	      }
	      else{
		hgenEtaDecayed[2]->Fill(GenParticle_eta[i],decayedChim);
		hgenEtaDecayed[0]->Fill(GenParticle_eta[i],decayedChim);
	      }
	    }
	  }
	  
	  
	  for(int i=0; i<nGenParticle; i++){


	    if(abs(GenParticle_pdgId[i])==1000024){
	      
	      if(GenParticle_pdgId[i]>0){

		if(!decayedChip){
		  hgenEtaNotDecayed[0]->Fill(GenParticle_eta[i]);
		  hgenEtaNotDecayed[1]->Fill(GenParticle_eta[i]);
		}
	
		hgenBetaDecayed[1]->Fill(GenParticle_pt[i]/genChipmVector[1].E(),decayedChip);
		hgenBetaDecayed[0]->Fill(GenParticle_pt[i]/genChipmVector[1].E(),decayedChip);
		if(decayedChip){
		  countCharg += 1;
		  hBetaDecayed[0]->Fill(GenParticle_pt[i]/genChipmVector[1].E());
		  hBetaDecayed[1]->Fill(GenParticle_pt[i]/genChipmVector[1].E());
		}
		
	      }
	      else if(GenParticle_pdgId[i]<0){

		if(!decayedChim){
		  hgenEtaNotDecayed[0]->Fill(GenParticle_eta[i]);
		  hgenEtaNotDecayed[2]->Fill(GenParticle_eta[i]);
		}
	
		hgenBetaDecayed[2]->Fill(GenParticle_pt[i]/genChipmVector[2].E(),decayedChim);
		hgenBetaDecayed[0]->Fill(GenParticle_pt[i]/genChipmVector[2].E(),decayedChim);
		if(decayedChim){
		  countCharg += 1;
		  hBetaDecayed[0]->Fill(GenParticle_pt[i]/genChipmVector[2].E());
		  hBetaDecayed[2]->Fill(GenParticle_pt[i]/genChipmVector[2].E());
		}
			
	      }
	      
	    }
	  }
	  
	  //ptVsDecayedChip -> Fill(SimTrack_momentum_pt[trackidxChip],decayedChip);
	  //ptVsDecayedChim -> Fill(SimTrack_momentum_pt[trackidxChim],decayedChim);	  

	}//end of loop over events
  
  cout<<endl<<endl<<htrkLengthChipm[1]->Integral()<<" charginos+ decayed."<<endl;
  cout<<htrkLengthChipm[2]->Integral()<<" charginos- decayed."<<endl;
  cout<<hgenMChipm[0]->Integral()<<" charginos existed."<<endl;
  cout<<htrkPtChipm[0]->Integral()<<" chargino tracks were found."<<endl;
  cout<<(htrkLengthChipm[1]->Integral()+htrkLengthChipm[2]->Integral())/hgenMChipm[0]->Integral()*100.<<"% of all charginos decayed."<<endl;
  cout<<(hgenMChipm[0]->Integral()-htrkPtChipm[0]->Integral())/hgenMChipm[0]->Integral()*100.<<"% of generated charginos couldn't be associated with a track."<<endl;
  //cout<<(1.-(1.*nChipDecay+1.*nChimDecay)/(1.*nChipm))*100.<<"% of all charginos did NOT decay."<<endl;
  cout<<"mean decay length Chip = "<<htrkLengthChipm[1]->GetMean()<<" cm."<<endl;
  cout<<"mean decay length Chim = "<<htrkLengthChipm[2]->GetMean()<<" cm."<<endl<<endl<<endl;
  cout<<"countCharg = "<<countCharg<<endl;
  cout<<" countNvertex = "<< countNvertex<<endl;
  cout<<"alright = "<<allright<<endl;

  cout<<endl<<endl;
  
  stream.close();
  ofile.close();
  
  return 0;
}
