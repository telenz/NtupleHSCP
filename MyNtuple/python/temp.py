#-------------------------------------------------------------------------
# Created: Mon Feb 10 18:02:35 2014 by mkntuplecfi.py
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
    'PFJet'
    ),

               PFJet =
               cms.untracked.
               vstring(
    'recoPFJet                       ak5PFJetsPt15                   200',
    #---------------------------------------------------------------------
    'float  chargedEmEnergyFraction()',
    'float  neutralEmEnergyFraction()'
    )
               )
