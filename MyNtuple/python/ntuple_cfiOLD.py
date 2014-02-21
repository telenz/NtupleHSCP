#-------------------------------------------------------------------------
# Created: Mon Sep 30 14:18:17 2013 by mkntuplecfi.py
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
    'GenParticle',
    'recoGenParticleHelper'
    ),

               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    200',
    #---------------------------------------------------------------------
    'double  mass()',
    'double  vx()',
    'double  vy()',
    'double  vz()',
    #'int  pdgId()',
    #'double  motherRef(size_t i=0)->mass()',
    #'int  motherRef(size_t i=0)->pdgId()'
    ),
               
               recoGenParticleHelper =
               cms.untracked.
               vstring(
    "recoGenParticleHelper           genParticles                     500",
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
    " double   mass()",
    "    int   status()"
    )
               )
