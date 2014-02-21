#-------------------------------------------------------------------------
# Created: Fri Jan 17 11:59:05 2014 by mkntuplecfi.py
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
    'edmTriggerResultsHelper',
    'Track',
    'GenParticle',
    'PFJet',
    'PFMET',
    'Vertex'
    ),

               edmTriggerResultsHelper =
               cms.untracked.
               vstring(
    'edmTriggerResultsHelper TriggerResults::HLT 1',
    #---------------------------------------------------------------------
    ' int value("HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5")',
    ' int prescale("HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5")',
    #' int value("HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5")',
    #' int prescale("HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5")',
    ),

               Track =
               cms.untracked.
               vstring(
    'recoTrackHelper                       TrackRefitter                   2000000',
    #---------------------------------------------------------------------
    'double pt()',
    'double eta()',
    'double phi()',
    'unsigned short  numberOfValidHits()',
    'unsigned short  numberOfLostHits()',
    'double  d0()',
    'double  dz()',
    'double  d0Error()',
    'double  dzError()',
    'double  ptError()',
    'double  dEdxNPHarm2()',
    'double  dEdxNPTru40()',
    'unsigned int     dEdxNPNoM()',
    'double  dEdxHitsNPHarm2(1000)',
    'double  dEdxHitsNPHarm2(7)',
    'double  dEdxHitsNPHarm2(5)',
    'double  dEdxHitsNPHarm2(3)',
    'double  dEdxHitsNPHarm2(2)',
    'double  dEdxHitsNPHarm2(1)',
    'double  dEdxHitsNPTrun40(1000)',
    'double  dEdxHitsNPTrun40(7)',
    'double  dEdxHitsNPTrun40(5)',
    'double  dEdxHitsNPTrun40(3)',
    'double  dEdxHitsNPTrun40(2)',
    'double  dEdxHitsNPTrun40(1)',
    'double  dEdxHitsNPMedian(1000)',
    'double  dEdxHitsNPMedian(7)',
    'double  dEdxHitsNPMedian(5)',
    'double  dEdxHitsNPMedian(3)',
    'double  dEdxHitsNPMedian(2)',
    'double  dEdxHitsNPMedian(1)',
    'unsigned short trackerExpectedHitsOuter().numberOfLostHits()',
    'unsigned short trackerExpectedHitsInner().numberOfLostHits()',
    ),
               
               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    1000',
    #---------------------------------------------------------------------
    'int  charge()',
    'double  energy()',
    'double  et()',
    'double  p()',
    'double  pt()',
    'double  phi()',
    'double  eta()',
    'size_t  numberOfDaughters()',
    'size_t  numberOfMothers()',
    'double  mass()',
    'double  vx()',
    'double  vy()',
    'double  vz()',
    'int  pdgId()',
    'int  status()'
    ),
                              
               PFJet =
               cms.untracked.
               vstring(
    'recoPFJet                       ak5PFJetsPt15                   200',
    #---------------------------------------------------------------------
    'double  energy()',
    'double  et()',
    'double  p()',
    'double  pt()',
    'double  phi()',
    'double  eta()',
    'float  chargedHadronEnergyFraction()',
    'float  neutralHadronEnergyFraction()',
    'float  neutralEmEnergyFraction()',
    'float  chargedEmEnergyFraction()'
    ),
               
               PFMET =
               cms.untracked.
               vstring(
    'recoPFMET                       pfMet                           200',
    #---------------------------------------------------------------------
    'double  energy()',
    'double  et()',
    'double  p()',
    'double  pt()',
    'double  phi()',
    'double  eta()'
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
    )
               
               )
