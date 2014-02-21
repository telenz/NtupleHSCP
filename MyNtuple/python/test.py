#-------------------------------------------------------------------------
# Created: Fri Feb 21 10:48:21 2014 by mkntuplecfi.py
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
    'GenEventInfoProduct',
    'GenRunInfoProduct'
    ),

               GenEventInfoProduct =
               cms.untracked.
               vstring(
    'GenEventInfoProduct             generator                         1',
    #---------------------------------------------------------------------
    'double  weight()'
    ),
               GenRunInfoProduct =
               cms.untracked.
               vstring(
    'GenRunInfoProduct               generator                         1',
    #---------------------------------------------------------------------
    'double  filterEfficiency()',
    'double  crossSection()'
    )
               )
