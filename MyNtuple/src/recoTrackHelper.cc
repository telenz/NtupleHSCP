//-----------------------------------------------------------------------------
// Subsystem:   Ntuples
// Package:     MyNtuple
// Description: TheNtupleMaker helper class for reco::Track
// Created:     Fri Jan 17 10:01:44 2014
// Author:      Teresa Lenz      
//-----------------------------------------------------------------------------
#include "Ntuples/MyNtuple/interface/recoTrackHelper.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include <iostream>

//-----------------------------------------------------------------------------
using namespace std;
using namespace reco;
//-----------------------------------------------------------------------------
// This constructor is called once per job
TrackHelper::TrackHelper()
  : HelperFor<reco::Track>() {}
    
TrackHelper::~TrackHelper() {}

// -- Called once per event
void TrackHelper::analyzeEvent()
{
// For DeDxNPHarm2
edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxNPHarm2TrackHandle;
event->getByLabel("dedxNPHarm2", dEdxNPHarm2TrackHandle);
dEdxTrackMap = *dEdxNPHarm2TrackHandle.product();

event->getByLabel(labelname,trackCollectionHandle);

// For DeDxHitsNPHarm2
edm::Handle<edm::ValueMap<susybsm::HSCPDeDxInfo> > dEdxHitsNPHarm2TrackHandle;
event->getByLabel("dedxHitInfo", dEdxHitsNPHarm2TrackHandle);
dEdxHitsTrackMap = *dEdxHitsNPHarm2TrackHandle.product();

event->getByLabel(labelname,trackCollectionHandle);

}

// -- Called once per object
void TrackHelper::analyzeObject()
{
// For DeDxNPHarm2 
reco::TrackRef track  = reco::TrackRef( trackCollectionHandle, oindex);
dEdxNPHarm2Track = dEdxTrackMap[track];

// For DeDxHitsNPHarm2 
dEdxHitsNPHarm2Track = dEdxHitsTrackMap[track];

}

// -- User access methods
//double TrackHelper::someMethod()  const
//{
//  return  //-- some-value --
//}
