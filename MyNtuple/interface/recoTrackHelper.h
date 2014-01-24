#ifndef RECOTRACKHELPER_H
#define RECOTRACKHELPER_H
//-----------------------------------------------------------------------------
// Subsystem:   Ntuples
// Package:     MyNtuple
// Description: TheNtupleMaker helper class for reco::Track
// Created:     Fri Jan 17 10:01:44 2014
// Author:      Teresa Lenz      
//-----------------------------------------------------------------------------
#include <algorithm>
#include <iostream>
#include <vector>
#include <map>
#include "PhysicsTools/TheNtupleMaker/interface/HelperFor.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPCaloInfo.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPDeDxInfo.h" 
//-----------------------------------------------------------------------------
// Definitions:
//   helper:        object of class TrackHelper
//   helped object: object of class reco::Track
//
//
// The following variables are automatically defined and available to
//       all methods:
//
//         blockname          name of config. buffer object (config block) 
//         buffername         name of buffer in config block
//         labelname          name of label in config block (for getByLabel)
//         parameter          parameter (as key, value pairs)
//                            accessed as in the following example:
//
//                            string param = parameter("label");
//
//         0. hltconfig       pointer to HLTConfigProvider
//                            Note: this will be zero if HLTConfigProvider
//                                  has not been properly initialized
//
//         1. config          pointer to global ParameterSet object
//         2. event           pointer to the current event
//         3. object          pointer to the current helped object
//         4. oindex          index of current helped object
//
//         5. index           index of item(s) returned by helper.
//                            Note 1: an item is associated with all
//                                    helper methods (think of it as an
//                                    extension of the helped object)
//                                  
//                            Note 2: index may differ from oindex if,
//                                    for a given helped object, the count
//                                    variable (see below) differs from 1.
//
//         6. count           number of items per helped object (default=1)
//                            Note:
//                                  count = 0  ==> current helped object is
//                                                 to be skipped
//
//                                  count = 1  ==> current helped object is
//                                                 to be kept
//
//                                  count > 1  ==> current helped object is
//                                                 associated with "count"
//                                                 items, where each item
//                                                 is associated with all the
//                                                 helper methods
//
//       variables 0-6 are initialized by TheNtupleMaker.
//       variables 0-5 should not be changed.
//       variable    6 can be changed by the helper to control whether a
//                     helped object should be kept or generates more items
//-----------------------------------------------------------------------------
namespace reco
{
  /// A helper class for reco::Track.
  class TrackHelper : public HelperFor<reco::Track>
  {
  public:
	///
	TrackHelper();

	virtual ~TrackHelper();

	// Uncomment if this class does some event-level analysis
	virtual void analyzeEvent();
	 
	// Uncomment if this class does some object-level analysis
	virtual void analyzeObject();

	// ---------------------------------------------------------
	// -- User access methods go here
	// ---------------------------------------------------------
	
	//For DeDxNPHarm2
	reco::DeDxData dEdxNPHarm2Track;
	edm::ValueMap<reco::DeDxData> dEdxTrackMap;
	//For DeDxHitsNPHarm2
	susybsm::HSCPDeDxInfo dEdxHitsNPHarm2Track;
	edm::ValueMap<susybsm::HSCPDeDxInfo> dEdxHitsTrackMap;

	double dEdxHits(unsigned int nHits){

	  std::vector<double> vect_charge;
	  for(unsigned int i=0; i<dEdxHitsNPHarm2Track.charge.size(); i++){

	    if(dEdxHitsNPHarm2Track.detIds[i]<3) continue;   // skip pixels
	    if(!dEdxHitsNPHarm2Track.shapetest[i]) continue; // shape test ???
	    double Norm =3.61e-06*265;                       // ??
	    Norm *=10.0;                                     //mm --> cm
	    vect_charge.push_back(Norm*dEdxHitsNPHarm2Track.charge[i]/dEdxHitsNPHarm2Track.pathlength[i]);
	    if(vect_charge.size()==nHits) break;
	  }

	  int size = vect_charge.size();

	  //harmonic 2
	  double expo = -2;
	  double aux = 0;
	  for(int i = 0; i< size; i ++){
	    aux+=pow(vect_charge[i],expo);
	  }
	  return (size>0)?pow(aux/size,1./expo):0.;
	};
	
  private:
	// -- User internals
	edm::Handle<std::vector<reco::Track> > trackCollectionHandle;

  public:
    // ---------------------------------------------------------
    // -- Access Methods
    // ---------------------------------------------------------

	double dEdxNPHarm2() {return dEdxNPHarm2Track.dEdx();};
	unsigned int dEdxNPNoM() {return dEdxNPHarm2Track.numberOfSaturatedMeasurements();};
	
	double dEdxHitsNPHarm2(int nHits) {
	  return dEdxHits(nHits);
	};
	


	// WARNING: some methods may fail to compile because of coding
	//          problems in one of the CMSSW base classes. If so,
	//          just comment out the offending method and try again.
  


    // from reco::TrackBase
    reco::TrackBase::TrackAlgorithm algo() const { return object->algo(); }

    // from reco::TrackBase
    reco::TrackBase::TrackAlgorithm algoByName(std::string name) const
    { return object->algoByName(name); }

    // from reco::TrackBase
    std::string algoName() const { return object->algoName(); }

    // from reco::TrackBase
    /**
    std::string algoName(reco::TrackBase::TrackAlgorithm a) const
    { return object->algoName(a); }
    */

    // from reco::TrackBase
    int charge() const { return object->charge(); }

    // from reco::TrackBase
    double chi2() const { return object->chi2(); }

    // from reco::TrackBase
    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> >
    covariance() const
    { return object->covariance(); }

    // from reco::TrackBase
    double covariance(int i, int j) const { return object->covariance(i, j); }

    // from reco::TrackBase
    unsigned int covIndex(unsigned int i, unsigned int j) const
    { return object->covIndex(i, j); }

    // from reco::TrackBase
    double d0() const { return object->d0(); }

    // from reco::TrackBase
    double d0Error() const { return object->d0Error(); }

    // from reco::TrackBase
    double dsz() const { return object->dsz(); }

    // from reco::TrackBase
    /**
    double dsz(math::XYZPoint myBeamSpot) const
    { return object->dsz(myBeamSpot); }
    */

    // from reco::TrackBase
    double dszError() const { return object->dszError(); }

    // from reco::TrackBase
    double dxy() const { return object->dxy(); }

    // from reco::TrackBase
    /**
    double dxy(math::XYZPoint myBeamSpot) const
    { return object->dxy(myBeamSpot); }
    */

    // from reco::TrackBase
    /**
    double dxy(reco::BeamSpot theBeamSpot) const
    { return object->dxy(theBeamSpot); }
    */

    // from reco::TrackBase
    double dxyError() const { return object->dxyError(); }

    // from reco::TrackBase
    double dz() const { return object->dz(); }

    // from reco::TrackBase
    /**
    double dz(math::XYZPoint myBeamSpot) const
    { return object->dz(myBeamSpot); }
    */

    // from reco::TrackBase
    double dzError() const { return object->dzError(); }

    // from reco::TrackBase
    double error(int i) const { return object->error(i); }

    // from reco::TrackBase
    double eta() const { return object->eta(); }

    // from reco::TrackBase
    double etaError() const { return object->etaError(); }

    // from reco::Track
    const reco::TrackExtraRef extra() const { return object->extra(); }

    // from reco::Track
    unsigned short found() const { return object->found(); }

    // from reco::TrackBase
    const reco::HitPattern hitPattern() const { return object->hitPattern(); }

    // from reco::Track
    unsigned int innerDetId() const { return object->innerDetId(); }

    // from reco::Track
    const math::XYZVector innerMomentum() const
    { return object->innerMomentum(); }

    // from reco::Track
    bool innerOk() const { return object->innerOk(); }

    // from reco::Track
    const math::XYZPoint innerPosition() const
    { return object->innerPosition(); }

    // from reco::Track
    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> >
    innerStateCovariance() const
    { return object->innerStateCovariance(); }

    // from reco::TrackBase
    bool isLooper() const { return object->isLooper(); }

    // from reco::TrackBase
    double lambda() const { return object->lambda(); }

    // from reco::TrackBase
    double lambdaError() const { return object->lambdaError(); }

    // from reco::Track
    unsigned short lost() const { return object->lost(); }

    // from reco::TrackBase
    const math::XYZVector momentum() const { return object->momentum(); }

    // from reco::TrackBase
    double ndof() const { return object->ndof(); }

    // from reco::TrackBase
    signed char nLoops() const { return object->nLoops(); }

    // from reco::TrackBase
    double normalizedChi2() const { return object->normalizedChi2(); }

    // from reco::TrackBase
    unsigned short numberOfLostHits() const
    { return object->numberOfLostHits(); }

    // from reco::TrackBase
    unsigned short numberOfValidHits() const
    { return object->numberOfValidHits(); }

    // from reco::Track
    unsigned int outerDetId() const { return object->outerDetId(); }

    // from reco::Track
    double outerEta() const { return object->outerEta(); }

    // from reco::Track
    const math::XYZVector outerMomentum() const
    { return object->outerMomentum(); }

    // from reco::Track
    bool outerOk() const { return object->outerOk(); }

    // from reco::Track
    double outerP() const { return object->outerP(); }

    // from reco::Track
    double outerPhi() const { return object->outerPhi(); }

    // from reco::Track
    const math::XYZPoint outerPosition() const
    { return object->outerPosition(); }

    // from reco::Track
    double outerPt() const { return object->outerPt(); }

    // from reco::Track
    double outerPx() const { return object->outerPx(); }

    // from reco::Track
    double outerPy() const { return object->outerPy(); }

    // from reco::Track
    double outerPz() const { return object->outerPz(); }

    // from reco::Track
    double outerRadius() const { return object->outerRadius(); }

    // from reco::Track
    ROOT::Math::SMatrix<double,5,5,ROOT::Math::MatRepSym<double,5> >
    outerStateCovariance() const
    { return object->outerStateCovariance(); }

    // from reco::Track
    double outerTheta() const { return object->outerTheta(); }

    // from reco::Track
    double outerX() const { return object->outerX(); }

    // from reco::Track
    double outerY() const { return object->outerY(); }

    // from reco::Track
    double outerZ() const { return object->outerZ(); }

    // from reco::TrackBase
    double p() const { return object->p(); }

    // from reco::TrackBase
    double parameter(int i) const { return object->parameter(i); }

    // from reco::TrackBase
    ROOT::Math::SVector<double,5> parameters() const
    { return object->parameters(); }

    // from reco::TrackBase
    double phi() const { return object->phi(); }

    // from reco::TrackBase
    double phiError() const { return object->phiError(); }

    // from reco::TrackBase
    double pt() const { return object->pt(); }

    // from reco::TrackBase
    double ptError() const { return object->ptError(); }

    // from reco::TrackBase
    double px() const { return object->px(); }

    // from reco::TrackBase
    double py() const { return object->py(); }

    // from reco::TrackBase
    double pz() const { return object->pz(); }

    // from reco::TrackBase
    double qoverp() const { return object->qoverp(); }

    // from reco::TrackBase
    double qoverpError() const { return object->qoverpError(); }

    // from reco::TrackBase
    /**
    bool quality(reco::TrackBase::TrackQuality q) const
    { return object->quality(q); }
    */

    // from reco::TrackBase
    reco::TrackBase::TrackQuality qualityByName(std::string name) const
    { return object->qualityByName(name); }

    // from reco::TrackBase
    int qualityMask() const { return object->qualityMask(); }

    // from reco::TrackBase
    /**
    std::string qualityName(reco::TrackBase::TrackQuality q) const
    { return object->qualityName(q); }
    */

    // from reco::Track
    TrackingRecHitRef recHit(size_t i) const { return object->recHit(i); }

    // from reco::Track
    size_t recHitsSize() const { return object->recHitsSize(); }

    // from reco::TrackBase
    const math::XYZPoint referencePoint() const
    { return object->referencePoint(); }

    // from reco::Track
    const reco::TrackResiduals residuals() const
    { return object->residuals(); }

    // from reco::Track
    double residualX(int position) const
    { return object->residualX(position); }

    // from reco::Track
    double residualY(int position) const
    { return object->residualY(position); }

    // from reco::Track
    PropagationDirection seedDirection() const
    { return object->seedDirection(); }

    // from reco::Track
    edm::RefToBase<TrajectorySeed> seedRef() const
    { return object->seedRef(); }

    // from reco::TrackBase
    double theta() const { return object->theta(); }

    // from reco::TrackBase
    double thetaError() const { return object->thetaError(); }

    // from reco::TrackBase
    const reco::HitPattern trackerExpectedHitsInner() const
    { return object->trackerExpectedHitsInner(); }

    // from reco::TrackBase
    const reco::HitPattern trackerExpectedHitsOuter() const
    { return object->trackerExpectedHitsOuter(); }

    // from reco::TrackBase
    double validFraction() const { return object->validFraction(); }

    // from reco::TrackBase
    const math::XYZPoint vertex() const { return object->vertex(); }

    // from reco::TrackBase
    double vx() const { return object->vx(); }

    // from reco::TrackBase
    double vy() const { return object->vy(); }

    // from reco::TrackBase
    double vz() const { return object->vz(); }
  };
}
#endif
