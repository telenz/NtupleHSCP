#-------------------------------------------------------------------------
# Created: Tue Dec 10 15:30:30 2013 by mkntuplecfi.py
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
    'GenRunInfoProduct',
    'SimTrack',
    'SimVertex',
    'GenParticle'
    ),

               GenRunInfoProduct =
               cms.untracked.
               vstring(
    'GenRunInfoProduct               generator                         1',
    #---------------------------------------------------------------------
    'double  crossSection()'
    ),
               SimTrack =
               cms.untracked.
               vstring(
    'SimTrack                        g4SimHits                       200',
    #---------------------------------------------------------------------
    'float  charge()',
    'int  vertIndex()',
    'bool  noVertex()',
    'int  genpartIndex()',
    'bool  noGenpart()',
    'int  type()',
    'unsigned int  trackId()'
    ),
               SimVertex =
               cms.untracked.
               vstring(
    'SimVertex                       g4SimHits                       200',
    #---------------------------------------------------------------------
    'int  parentIndex()',
    'bool  noParent()',
    'unsigned int  vertexId()'
    ),
               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    200',
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
