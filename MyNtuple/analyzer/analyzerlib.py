# -----------------------------------------------------------------------------
#  File:        analyzerlib.py
#  Description: Analyzer for ntuples created by TheNtupleMaker
#  Created:     Fri Jan 24 05:53:13 2014 by mkanalyzer.py
#  Author:      Teresa Lenz
# -----------------------------------------------------------------------------
import os, sys, re
from ROOT import *
from string import atof, atoi, split, strip, replace
from array import array
from math import *
from time import sleep
from sys import exit
# -----------------------------------------------------------------------------
gSystem.AddIncludePath("-Iinclude")
gROOT.ProcessLine(".L src/treestream.cc+")
gROOT.ProcessLine(".L src/pdg.cc+")
# -----------------------------------------------------------------------------
# -- Classes, procedures and functions
# -----------------------------------------------------------------------------
class outputFile:

	def __init__(self, filename, stream=None, savecount=50000):
		if stream != None:
			print "events will be skimmed to file", filename			
			self.tree = stream.tree().CloneTree(0)
			self.weight = Double()
			self.b_weight = self.tree.Branch("eventWeight", self.weight,
											 "eventWeight/D")
			self.SAVECOUNT = savecount
		else:
			self.tree = None
			self.b_weight = None

		self.entry = 0

		self.filename = filename
		self.file = TFile(filename, "recreate")

		self.hist = TH1F("counts", "", 1, 0, 1)
		self.hist.SetBit(TH1.kCanRebin)
		self.hist.SetStats(0)

		self.b_weight = 0

	def addEvent(self, weight=1.0):
		if self.tree == None: return

		self.file = self.tree.GetCurrentFile()
		self.file.cd()
		self.tree.Fill()

		self.entry += 1		
		if self.entry % self.SAVECOUNT == 0:
			self.tree.AutoSave("SaveSelf")

	def count(self, cond, w=1.0):
		self.hist.Fill(cond, w)

	def close(self):
		print "==> histograms saved to file", self.filename
		if self.tree != None:			
			print "==> events skimmed to file", self.filename
			self.file = self.tree.GetCurrentFile()

		self.file.cd()
		#self.file.Write("", TObject.kOverwrite)
		self.file.Write()
		self.file.ls()
		self.file.Close()
# -----------------------------------------------------------------------------
class commandLine:
	def __init__(self):
		pass

def decodeCommandLine():
	argv = sys.argv
	argc = len(argv)

	cl = commandLine()
	cl.progname = split(os.path.basename(argv[0]),'.')[0]

	if argc > 1:
		cl.filelist = argv[1]
	else:
		cl.filelist = "filelist.txt"

	if argc > 2: 
		cl.outputfilename = argv[2] # 2nd (optional) command line argument
	else:
		cl.outputfilename = cl.progname + "_histograms"

	# Make sure extension is ".root"
	if cl.outputfilename[:-5] != ".root": cl.outputfilename += ".root"
	print "==> output to:", cl.outputfilename

	return cl
# -----------------------------------------------------------------------------
def error(message):
	print "** error ** " + message
	sys.exit(0)
# -----------------------------------------------------------------------------
#  Read ntuple filenames from file list

def getFilenames(filelist):
	if not os.path.exists(filelist):
		error("unable to open file: " + filelist)

	# Get ntuple file names
	filenames = filter(lambda x: x != "",
					   map(strip, open(filelist).readlines()))
	v = vector("string")()
	for filename in filenames: v.push_back(filename)
	return v
# -----------------------------------------------------------------------------
TEXTFONT = 42
TEXTSIZE = 0.031
#------------------------------------------------------------------------------
def setStyle():
	gROOT.SetStyle("Pub")
	style = gROOT.GetStyle("Pub")
	style.SetPalette(1)

	# For the canvas
	style.SetCanvasBorderMode(0)
	style.SetCanvasColor(kWhite)
	style.SetCanvasDefH(500)
	style.SetCanvasDefW(500)
	style.SetCanvasDefX(0)
	style.SetCanvasDefY(0)

	# For the pad
	style.SetPadBorderMode(0)
	style.SetPadColor(kWhite)
	style.SetPadGridX(kFALSE)
	style.SetPadGridY(kTRUE)
	style.SetGridColor(kGreen)
	style.SetGridStyle(3)
	style.SetGridWidth(1)

	# For the frame
	style.SetFrameBorderMode(0)
	style.SetFrameBorderSize(1)
	style.SetFrameFillColor(0)
	style.SetFrameFillStyle(0)
	style.SetFrameLineColor(1)
	style.SetFrameLineStyle(1)
	style.SetFrameLineWidth(1)

	# For the histogram
	style.SetHistLineColor(1)
	style.SetHistLineStyle(0)
	style.SetHistLineWidth(1)
	style.SetEndErrorSize(2)
	#style.SetErrorX(0.)
	style.SetMarkerSize(0.1)
	style.SetMarkerStyle(20)

	# For the fit/function
	style.SetOptFit(1)
	style.SetFitFormat("5.4g")
	style.SetFuncColor(2)
	style.SetFuncStyle(1)
	style.SetFuncWidth(1)

	#For the date
	style.SetOptDate(0)

	# For the statistics box
	style.SetOptFile(0)
	style.SetOptStat("")
	# To display the mean and RMS
	#style.SetOptStat("mr") 
	style.SetStatColor(kWhite)
	style.SetStatFont(TEXTFONT)
	style.SetStatFontSize(TEXTSIZE)
	style.SetStatTextColor(1)
	style.SetStatFormat("6.4g")
	style.SetStatBorderSize(1)
	style.SetStatH(0.2)
	style.SetStatW(0.3)

	# Margins
	style.SetPadTopMargin(0.05)
	style.SetPadBottomMargin(0.16)
	style.SetPadLeftMargin(0.16)
	style.SetPadRightMargin(0.16)

	# For the global title
	style.SetOptTitle(0)
	style.SetTitleFont(TEXTFONT)
	style.SetTitleColor(1)
	style.SetTitleTextColor(1)
	style.SetTitleFillColor(10)
	style.SetTitleFontSize(TEXTSIZE*1.1)

	# For the axis titles
	style.SetTitleColor(1, "XYZ")
	style.SetTitleFont(TEXTFONT, "XYZ")
	style.SetTitleSize(TEXTSIZE*1.2, "XYZ") # 0,05
	style.SetTitleXOffset(1.25) # 0.9
	style.SetTitleYOffset(1.25) # 1.25

	# For the axis labels
	style.SetLabelColor(1, "XYZ")
	style.SetLabelFont(TEXTFONT, "XYZ")
	style.SetLabelOffset(0.006, "XYZ")
	style.SetLabelSize(TEXTSIZE*1.2, "XYZ")

	# For the axis
	style.SetAxisColor(1, "XYZ")
	style.SetStripDecimals(kTRUE)
	style.SetTickLength(0.03, "XYZ")
	style.SetNdivisions(505, "XYZ")

	# To get tick marks on the opposite side of the frame
	style.SetPadTickX(1)  
	style.SetPadTickY(1)

	# Change for log plots
	style.SetOptLogx(0)
	style.SetOptLogy(0)
	style.SetOptLogz(0)

	# Postscript options
	style.SetPaperSize(20.,20.)

	style.cd()
# -----------------------------------------------------------------------------
#  Define variables to be read
# -----------------------------------------------------------------------------
GenParticle_charge	= vector("int")(100000)
GenParticle_energy	= vector("double")(100000)
GenParticle_et	= vector("double")(100000)
GenParticle_eta	= vector("double")(100000)
GenParticle_mass	= vector("double")(100000)
GenParticle_numberOfDaughters	= vector("size_t")(100000)
GenParticle_numberOfMothers	= vector("size_t")(100000)
GenParticle_p	= vector("double")(100000)
GenParticle_pdgId	= vector("int")(100000)
GenParticle_phi	= vector("double")(100000)
GenParticle_pt	= vector("double")(100000)
GenParticle_px	= vector("double")(100000)
GenParticle_py	= vector("double")(100000)
GenParticle_pz	= vector("double")(100000)
GenParticle_status	= vector("int")(100000)
GenParticle_vx	= vector("double")(100000)
GenParticle_vy	= vector("double")(100000)
GenParticle_vz	= vector("double")(100000)
GenRunInfoProduct_crossSection	= Double()
SimTrack_charge	= vector("float")(1000000)
SimTrack_genpartIndex	= vector("int")(1000000)
SimTrack_momentum_energy	= vector("double")(1000000)
SimTrack_momentum_eta	= vector("double")(1000000)
SimTrack_momentum_phi	= vector("double")(1000000)
SimTrack_momentum_pt	= vector("double")(1000000)
SimTrack_noGenpart	= vector("int")(1000000)
SimTrack_noVertex	= vector("int")(1000000)
SimTrack_trackId	= vector("unsigned int")(1000000)
SimTrack_type	= vector("int")(1000000)
SimTrack_vertIndex	= vector("int")(1000000)
SimVertex_noParent	= vector("int")(1000000)
SimVertex_parentIndex	= vector("int")(1000000)
SimVertex_position_t	= vector("double")(1000000)
SimVertex_position_x	= vector("double")(1000000)
SimVertex_position_y	= vector("double")(1000000)
SimVertex_position_z	= vector("double")(1000000)
SimVertex_vertexId	= vector("unsigned int")(1000000)
nGenParticle	= Long()
nSimTrack	= Long()
nSimVertex	= Long()

stream.select("recoGenParticle_genParticles.charge", GenParticle_charge)
stream.select("recoGenParticle_genParticles.energy", GenParticle_energy)
stream.select("recoGenParticle_genParticles.et", GenParticle_et)
stream.select("recoGenParticle_genParticles.eta", GenParticle_eta)
stream.select("recoGenParticle_genParticles.mass", GenParticle_mass)
stream.select("recoGenParticle_genParticles.numberOfDaughters", GenParticle_numberOfDaughters)
stream.select("recoGenParticle_genParticles.numberOfMothers", GenParticle_numberOfMothers)
stream.select("recoGenParticle_genParticles.p", GenParticle_p)
stream.select("recoGenParticle_genParticles.pdgId", GenParticle_pdgId)
stream.select("recoGenParticle_genParticles.phi", GenParticle_phi)
stream.select("recoGenParticle_genParticles.pt", GenParticle_pt)
stream.select("recoGenParticle_genParticles.px", GenParticle_px)
stream.select("recoGenParticle_genParticles.py", GenParticle_py)
stream.select("recoGenParticle_genParticles.pz", GenParticle_pz)
stream.select("recoGenParticle_genParticles.status", GenParticle_status)
stream.select("recoGenParticle_genParticles.vx", GenParticle_vx)
stream.select("recoGenParticle_genParticles.vy", GenParticle_vy)
stream.select("recoGenParticle_genParticles.vz", GenParticle_vz)
stream.select("GenRunInfoProduct_generator.crossSection", GenRunInfoProduct_crossSection)
stream.select("SimTrack_g4SimHits.charge", SimTrack_charge)
stream.select("SimTrack_g4SimHits.genpartIndex", SimTrack_genpartIndex)
stream.select("SimTrack_g4SimHits.momentum_energy", SimTrack_momentum_energy)
stream.select("SimTrack_g4SimHits.momentum_eta", SimTrack_momentum_eta)
stream.select("SimTrack_g4SimHits.momentum_phi", SimTrack_momentum_phi)
stream.select("SimTrack_g4SimHits.momentum_pt", SimTrack_momentum_pt)
stream.select("SimTrack_g4SimHits.noGenpart", SimTrack_noGenpart)
stream.select("SimTrack_g4SimHits.noVertex", SimTrack_noVertex)
stream.select("SimTrack_g4SimHits.trackId", SimTrack_trackId)
stream.select("SimTrack_g4SimHits.type", SimTrack_type)
stream.select("SimTrack_g4SimHits.vertIndex", SimTrack_vertIndex)
stream.select("SimVertex_g4SimHits.noParent", SimVertex_noParent)
stream.select("SimVertex_g4SimHits.parentIndex", SimVertex_parentIndex)
stream.select("SimVertex_g4SimHits.position_t", SimVertex_position_t)
stream.select("SimVertex_g4SimHits.position_x", SimVertex_position_x)
stream.select("SimVertex_g4SimHits.position_y", SimVertex_position_y)
stream.select("SimVertex_g4SimHits.position_z", SimVertex_position_z)
stream.select("SimVertex_g4SimHits.vertexId", SimVertex_vertexId)
stream.select("nrecoGenParticle_genParticles", nGenParticle)
stream.select("nSimTrack_g4SimHits", nSimTrack)
stream.select("nSimVertex_g4SimHits", nSimVertex)

