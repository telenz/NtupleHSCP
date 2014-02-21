#-------------------------------------------------------------------------
# Created: Fri Jan 17 14:53:46 2014 by mkntuplecfi.py
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
    'GenParticle'
    ),

               GenParticle =
               cms.untracked.
               vstring(
    'recoGenParticle                 genParticles                    200',
    #---------------------------------------------------------------------
    'double  px()',
    'double  py()',
    'double  pz()',
    'double  pt()'
    )
               )
