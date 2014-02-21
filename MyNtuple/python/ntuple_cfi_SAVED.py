#-------------------------------------------------------------------------
# Created: Fri Oct 11 13:12:51 2013 by mkntuplecfi.py
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
    'SimTrack',
    'SimVertex',
    'GenParticle',
    'recoGenParticleHelper'
    ),

               SimTrack =
               cms.untracked.
               vstring(
    'SimTrack                        g4SimHits                       10000',
    #---------------------------------------------------------------------
    'float  charge()',
    'int  vertIndex()',
    'bool  noVertex()',
    'int  genpartIndex()',
    'bool  noGenpart()',
    'int  type()',
    'unsigned int  trackId()',
    'double momentum().pt()',
    'double momentum().phi()',
    'double momentum().eta()',
    'double momentum().energy()',
    ),
               SimVertex =
               cms.untracked.
               vstring(
    'SimVertex                       g4SimHits                       10000',
    #---------------------------------------------------------------------
    'int  parentIndex()',
    'bool  noParent()',
    'unsigned int  vertexId()',
    'double position().x()',
    'double position().y()',
    'double position().z()',
    'double position().t()',
    ),
               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    10000',
    #---------------------------------------------------------------------
    'double  vx()',
    'double  vy()',
    'double  vz()',
    'double  pt()',
    'double  et()',
    'double  pt()',
    'double  eta()',
    'double  phi()',
    'int  pdgId()',
    'int  status()',
    ),
               recoGenParticleHelper =
               cms.untracked.
               vstring(
    "recoGenParticleHelper           genParticles                     10000",
    #---------------------------------------------------------------------
#    'double  mass()',
#    'double  vx()',
#    'double  vy()',
#    'double  vz()',
#    'int  pdgId()',
#    'int  motherRef(size_t i=0)->charge()',
#    'double  motherRef(size_t i=0)->mass()',
#    'int  motherRef(size_t i=0)->pdgId()'
    "    int   firstMother()",
    "    int   lastMother()",
    "    int   firstDaughter()",
    "    int   lastDaughter()",
    "    int   charge()",
    "    int   pdgId()",
    "    int   status()",
    " double   pt()",
    " double   eta()",
    " double   phi()",
    " double   mass()"
    )
               )
