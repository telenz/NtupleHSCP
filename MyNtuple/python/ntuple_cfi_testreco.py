#-------------------------------------------------------------------------
# Created: Mon Feb 10 14:07:10 2014 by mkntuplecfi.py
#-------------------------------------------------------------------------
import FWCore.ParameterSet.Config as cms
demo =\
cms.EDAnalyzer("TheNtupleMaker",
               ntupleName = cms.untracked.string("ntuple.root"),
               analyzerName = cms.untracked.string("analyzer.cc"),


# NOTE: the names listed below will be the prefixes for
#       the associated C++ variables created by mkanalyzer.py
#       and the asscociated C++ structs.

               buffers =
               cms.untracked.
               vstring(
    'TriggerResults',
    'GenParticle',
    'PFJet',
    'PFMET',
    'Track',
    'Vertex'
    ),

               TriggerResults =
               cms.untracked.
               vstring(
    'edmTriggerResults               TriggerResults                    1',
    #---------------------------------------------------------------------
    'bool  accept()',
    'bool  wasrun(unsigned int i)'
    ),
               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    200',
    #---------------------------------------------------------------------
    'double  p()',
    'double  et()',
    'double  px()',
    'double  py()',
    'double  pz()',
    'double  pt()',
    'double  phi()',
    'double  eta()',
    'int  pdgId()'
    ),
               PFJet =
               cms.untracked.
               vstring(
    'recoPFJet                       ak5PFJetsPt15                   200',
    #---------------------------------------------------------------------
    'double  p()',
    'double  energy()',
    'double  et()',
    'double  pt()',
    'double  phi()',
    'double  eta()',
    'float  chargedHadronEnergyFraction()',
    'float  neutralHadronEnergyFraction()',
    'float  neutralEmEnergyFraction()'
    ),
               PFMET =
               cms.untracked.
               vstring(
    'recoPFMET                       pfMet                           200',
    #---------------------------------------------------------------------
    'double  energy()',
    'double  et()',
    'double  pt()',
    'double  phi()',
    'double  eta()'
    ),
               Track =
               cms.untracked.
               vstring(
    'recoTrack                       TrackRefitter                   200',
    #---------------------------------------------------------------------
    'double  p()',
    'double  pt()',
    'double  px()',
    'double  py()',
    'double  pz()',
    'double  phi()',
    'double  eta()',
    'double  d0()',
    'double  dz()',
    'double  ptError()',
    'double  etaError()',
    'double  phiError()',
    'double  d0Error()',
    'double  dzError()',
    'unsigned short  numberOfValidHits()',
    'unsigned short  numberOfLostHits()'
    ),
               Vertex =
               cms.untracked.
               vstring(
    'recoVertex                      offlinePrimaryVertices          200',
    #---------------------------------------------------------------------
    'bool  isValid()',
    'bool  isFake()',
    'size_t  tracksSize()',
    'double  chi2()',
    'double  ndof()',
    'double  normalizedChi2()',
    'double  x()',
    'double  y()',
    'double  z()',
    'double  xError()',
    'double  yError()',
    'double  zError()',
    'double  covariance(int i, int j)',
    'bool  hasRefittedTracks()',
    'unsigned int  nTracks(float minWeight=5.0e-1)'
    )
               )
