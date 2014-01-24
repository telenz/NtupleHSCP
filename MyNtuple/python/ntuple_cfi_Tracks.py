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
    'Track',
    'GenParticle'
    ),

               Track =
               cms.untracked.
               vstring(
    'recoTrackHelper                       TrackRefitter                   20000000',
    #---------------------------------------------------------------------
    'double  pt()',
    'double  px()',
    'double  py()',
    'double  pz()',
    'double eta()',
    'double phi()',
    'double  dEdxNPHarm2()',
    'unsigned int     dEdxNPNoM()',
    'double  dEdxHitsNPHarm2(1000)',
    'double  dEdxHitsNPHarm2(7)',
    ),
               
               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    1000',
    #---------------------------------------------------------------------
    'int  charge()',
    'double  p()',
    'double  energy()',
    'double  et()',
    'double  px()',
    'double  py()',
    'double  pz()',
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
    )
               )
