#ifndef ANALYZERCLASSESAODSIMWELLS_H
#define ANALYZERCLASSESAODSIMWELLS_H

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
    hVertexNdof -> Fill(evt::Vertex[i].ndof);
    hVertexZ    -> Fill(evt::Vertex[i].z);
    hVertexRho  -> Fill(std::sqrt(std::pow(evt::Vertex[i].x,2) + std::pow(evt::Vertex[i].y,2)));
  }


};



#endif
