//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "analyzerdEdx.h"
#include "analyzerFunctions.h"
#include "analyzerRecoTracks.h"
#include "analyzerRECOClasses.h"
#include "analyzerClassesAODSIM.h"
#include "config.h"
#include "histograms.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include "TProfile.h"
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

  TH1D *hMet100GeV = new TH1D("hMet100GeV","hMet100GeV",150,0,1500);
  TH1D *h1stjetpt100GeV  = new TH1D("h1stjetpt100GeV","h1stjetpt100GeV",150,50,1500);
  TH1D *h1stjetpt80GeV   = new TH1D("h1stjetpt80GeV","h1stjetpt80GeV",150,50,1500);
  TH1D *h1stjetpt110GeV  = new TH1D("h1stjetpt110GeV","h1stjetpt110GeV",150,50,1500);
  TH1D *h1stjetpt120GeV  = new TH1D("h1stjetpt120GeV","h1stjetpt120GeV",150,50,1500);


  class Event all("AllTracks",ofile);
  all.onlyChipm=false;
  all.trackCandidateCut=false;
  class Event onlyChi("OnlyChiTracks",ofile);
  onlyChi.onlyChipm=true;
  onlyChi.trackCandidateCut=false;
  class Event noChi("NoChiTracks",ofile);
  noChi.noChi=true;
  noChi.trackCandidateCut=false;
  class Event onlyChiSelectedTracks("OnlyChiSelectedTracks",ofile);
  onlyChiSelectedTracks.onlyChipm=true;
  onlyChiSelectedTracks.trackCandidateCut=true;
  class Event allSelectedTracks("AllSelectedTracks",ofile);
  allSelectedTracks.onlyChipm=false;
  allSelectedTracks.trackCandidateCut=true;
  class Event allSelectedTracksMET1stJET("AllSelectedTracksMET1stJET",ofile);
  allSelectedTracksMET1stJET.onlyChipm=false;
  allSelectedTracksMET1stJET.trackCandidateCut=true;
  allSelectedTracksMET1stJET.metCut=true;
  allSelectedTracksMET1stJET.leadingJetCut=true;
  class Event onlyChiSelectedTracksMET1stJET("OnlyChiSelectedTracksMET1stJET",ofile);
  onlyChiSelectedTracksMET1stJET.onlyChipm=true;
  onlyChiSelectedTracksMET1stJET.trackCandidateCut=true;
  onlyChiSelectedTracksMET1stJET.metCut=true;
  onlyChiSelectedTracksMET1stJET.leadingJetCut=true;
  class Event onlyChiSoftSelectedTracks("OnlyChiSoftSelectedTracks",ofile);
  onlyChiSoftSelectedTracks.onlyChipm=true;
  onlyChiSoftSelectedTracks.trackCandidateSoftCut=true;
  class Event allSoftSelectedTracks("AllSoftSelectedTracks",ofile);
  allSoftSelectedTracks.onlyChipm=false;
  allSoftSelectedTracks.trackCandidateSoftCut=true;


  //class SetOfVertexHistograms vertices;

 
  //TH1D *h1stGenJetPtGt30  = new TH1D("h1stGenJetPtGt30","h1stGenJetPtGt30",150,0,1500);
  //TH2D *h1stGenJetPtGt30ChipmPt  = new TH2D("h1stGenJetPtGt30ChipmPt","h1stGenJetPtGt30ChipmPt",150,0,1500,150,0,1500);
  //TH2D *h1stGenJetPtGt30RecoJetPt  = new TH2D("h1stGenJetPtGt30RecoJetPt","h1stGenJetPtGt30RecoJetPt",150,0,1500,150,0,1500);
  
  
  
  //---------------------------------------------------------------------------
  // Declaration of Variables
  //---------------------------------------------------------------------------

  // Mass Reconstruction
  //double K = 2.529;
  //double C = 2.772;
  /*
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
  */

  //---------------------------------------------------------------------------

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
      ofile.count("NoCuts", 0);
      //-------------------------------------------------------------- Look for special collections  --------------------------------
      findChipmInGenParticleCollection();

      if(nPFMET==1){
	if(PFMET[0].et>100.)   hMet100GeV->Fill(PFMET[0].et);
      }
      std::vector<PFJet_s> jetCollection      = getSelectedJetCollection();
      if(jetCollection.size()>0 && jetCollection[0].pt>100) h1stjetpt100GeV->Fill(jetCollection[0].pt);
      if(jetCollection.size()>0 && jetCollection[0].pt>80)  h1stjetpt80GeV->Fill(jetCollection[0].pt);
      if(jetCollection.size()>0 && jetCollection[0].pt>110) h1stjetpt110GeV->Fill(jetCollection[0].pt);
      if(jetCollection.size()>0 && jetCollection[0].pt>120) h1stjetpt120GeV->Fill(jetCollection[0].pt);
      // onlyChi.PFJet   = getSelectedJetCollection();
      // allTracks.PFJet = getSelectedJetCollection();
      
      //onlyChi.leadJetGenParticle               = findLeadingJetInGenParticleCollectionWithPtGt30();
      //allTracks.leadJetGenParticle = findLeadingJetInGenParticleCollectionWithPtGt30();

      //onlyChi.Track   = getChipmInTrackCollection(Track);
      //allTracks.Track = Track;
      
      //onlyChiSelectedTracks = OnlyChi;
      //allSelectedTracks = all;
      //-------------------------------------------------------------- Cuts ---------------------------------------------------------

      //if(!triggerFired(&jetCollection[0]))  continue;
      //onlyChiSelectedTracks.Track = getCandidateTrackCollection(&onlyChi);
      //all.Track                   = getCandidateTrackCollection(&all);
 
      //all.hist.FillTrackVariables(all.Track);
      //onlyChi.hist.FillTrackVariables(onlyChi.Tracks)
     
      all.Selection();
      allSelectedTracks.Selection();
      onlyChi.Selection();
      noChi.Selection();
      onlyChiSelectedTracks.Selection();
      allSelectedTracksMET1stJET.Selection();
      onlyChiSelectedTracksMET1stJET.Selection();
      onlyChiSoftSelectedTracks.Selection();
      allSoftSelectedTracks.Selection();
      

      // a few jet related histograms 
      //if(!zeroJet)                  h1stGenJetPtGt30->Fill(leadJetGenParticle.pt);
      //if(!zeroJet&&!zeroChip) h1stGenJetPtGt30ChipmPt->Fill(leadJetGenParticle.pt,chipGenParticle.pt);
      //if(!zeroJet&&!zeroChim) h1stGenJetPtGt30ChipmPt->Fill(leadJetGenParticle.pt,chimGenParticle.pt);
      
      //if(!zeroJet && jetCollection.size()>0) h1stGenJetPtGt30RecoJetPt->Fill(leadJetGenParticle.pt,jetCollection[0].pt);
      

      // Loop over track collection
      //Nchip=0;
      //Nchim=0;

      //onlyChi.SetOFHistogramsChipmTracks.FillDeDxHistograms(onlyChi.Track);
      //onlyChiSelectedTracks.SetOFHistogramsChipmTracks.FillDeDxHistograms(onlyChiSelectedTracks.Track);
      //all.SetOFHistogramsChipmTracks.FillDeDxHistograms(all.Track);
      
	  
      // ---------------------
      // -- fill histograms --
      // ---------------------


      //------------------------------------------------------------------------------------------------------------------------------------
      //-------- GenParticle Helper  Collection --------------------------------------------------------------------------------------------
      //------------------------------------------------------------------------------------------------------------------------------------

      
    }//end of loop over events


  

  //TH1D* projectionX = new TH1D("projectionX","projectionX",100,0,100);
  /*
  for(int x=1; x<Chipm.harm2.hPDeDx->GetNbinsX(); x++){

    for(int y=1; y<Chipm.harm2.hPDeDx->GetNbinsY(); y++){


      mean += Chipm.harm2.hPDeDx->GetBinContent(x,y);



    }
  }
  */

  /*
  TProfile* chipm =  Chipm.harm2.hPDeDx->ProfileX("_pfx", 1, -1,"o");
  TProfile* all   =  All.harm2.hPDeDx->ProfileX("_pfx", 1, -1,"o");
  chipm->GetYaxis()->SetRangeUser(0,40);
  all->GetYaxis()->SetRangeUser(0,40);
  
  TF1 *f1 = new TF1("f1","2.529*[0]**2/x**2+2.772",700,1000);
  f1->SetParameter(0,800);
  chipm->Fit("f1","QR");
  cout<<endl<<"Parameter = "<<f1->GetParameter(0)<<" GeV"<<endl;
  cout<<"Error = "<<f1->GetParError(0)<<endl<<endl;
  cout<<"Chi2/ndof = "<<f1->GetChisquare()<<"/"<<f1->GetNDF()<<endl<<endl;
  all->Fit("f1","QR");
  cout<<endl<<"Parameter = "<<f1->GetParameter(0)<<" GeV"<<endl;
  cout<<"Error = "<<f1->GetParError(0)<<endl<<endl;
  cout<<"Chi2/ndof = "<<f1->GetChisquare()<<"/"<<f1->GetNDF()<<endl<<endl<<endl;
  */
  
  /*
  cout<<endl<<"chipmFound = "<<chipmFound<<endl<<endl;
  cout<<"matching efficiency = "<<chipmFound/(2.*n)<<endl<<endl;
  cout<<endl<<"It took "<<(clock()-t)/CLOCKS_PER_SEC<<" sec"<<endl;

  cout<<"nJet = "<<nJet<<endl;


  cout<<"Vertex cut efficiency  = "<<1.*afterVertexCut/beforeVertexCut<<endl;
  cout<<"MET cut efficiency     = "<<1.*afterMETCut/beforeMETCut<<endl;
  cout<<"before1stJetCut        = "<<before1stJetCut<<endl;
  cout<<"1st Jet cut efficiency = "<<1.*after1stJetCut/before1stJetCut<<endl;

  cout<<"MisMatched = "<<mismatchedGenChipmToTrack<<endl<<endl;
  */


  //------------------------------------------------------------------
  ofstream text;
  text.open("importantFindings.txt");
  text<<"Matching of Chipm in GenParticles and in track Collection: "<<endl;
  text<<"Mismatched Chipm (when abs(trackPt/genPt-1)>0.8)="<<mismatchedGenChipmToTrack<<endl;
  text<<"Matched Chipm (when deltaR<0.5) = "<<matchedGenChipmToTrack<<endl<<endl;
  //------------------------------------------------------------------

  stream.close();
  ofile.close();
  
  cout<<endl;
  return 0;
}
