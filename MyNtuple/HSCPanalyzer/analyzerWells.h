#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Thu Jul 17 10:30:40 2014 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
// -- System

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>

#include "analyzerutil.h"
#include "treestream.h"
#include "pdg.h"

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"

namespace evt {
//-----------------------------------------------------------------------------
// --- Declare variables
//-----------------------------------------------------------------------------
std::vector<double>	Electron_energy(200,0);
std::vector<double>	Electron_et(200,0);
std::vector<double>	Electron_eta(200,0);
std::vector<float>	Electron_mvaNonTrigV0(200,0);
std::vector<double>	Electron_phi(200,0);
std::vector<double>	Electron_pt(200,0);
std::vector<double>	Electron_pz(200,0);
std::vector<int>	GenParticle_charge(2000,0);
std::vector<double>	GenParticle_energy(2000,0);
std::vector<double>	GenParticle_et(2000,0);
std::vector<double>	GenParticle_eta(2000,0);
std::vector<double>	GenParticle_mass(2000,0);
std::vector<double>	GenParticle_p(2000,0);
std::vector<int>	GenParticle_pdgId(2000,0);
std::vector<double>	GenParticle_phi(2000,0);
std::vector<double>	GenParticle_pt(2000,0);
std::vector<double>	GenParticle_pz(2000,0);
std::vector<float>	Jet_chargedEmEnergyFraction(200,0);
std::vector<float>	Jet_chargedHadronEnergyFraction(200,0);
std::vector<double>	Jet_energy(200,0);
std::vector<double>	Jet_et(200,0);
std::vector<double>	Jet_eta(200,0);
std::vector<float>	Jet_neutralEmEnergyFraction(200,0);
std::vector<float>	Jet_neutralHadronEnergyFraction(200,0);
std::vector<double>	Jet_phi(200,0);
std::vector<double>	Jet_pt(200,0);
std::vector<double>	Jet_pz(200,0);
std::vector<double>	MET_energy(200,0);
std::vector<double>	MET_et(200,0);
std::vector<double>	MET_eta(200,0);
std::vector<double>	MET_phi(200,0);
std::vector<double>	MET_pt(200,0);
std::vector<double>	MET_pz(200,0);
std::vector<double>	Muon_energy(200,0);
std::vector<double>	Muon_et(200,0);
std::vector<double>	Muon_eta(200,0);
std::vector<double>	Muon_phi(200,0);
std::vector<double>	Muon_pt(200,0);
std::vector<double>	Muon_pz(200,0);
std::vector<int>	PileupSummaryInfo_getBunchCrossing(10,0);
std::vector<int>	PileupSummaryInfo_getPU_NumInteractions(10,0);
std::vector<float>	PileupSummaryInfo_getTrueNumInteractions(10,0);
std::vector<float>	Tau_againstElectronLoose(200,0);
std::vector<float>	Tau_againstMuonTight(200,0);
std::vector<float>	Tau_byLooseCombinedIsolationDeltaBetaCorr(200,0);
std::vector<float>	Tau_decayModeFinding(200,0);
std::vector<double>	Tau_energy(200,0);
std::vector<double>	Tau_et(200,0);
std::vector<double>	Tau_eta(200,0);
std::vector<double>	Tau_phi(200,0);
std::vector<double>	Tau_pt(200,0);
std::vector<double>	Tau_pz(200,0);
std::vector<double>	Track_caloEMDeltaRp3(2000,0);
std::vector<double>	Track_caloEMDeltaRp4(2000,0);
std::vector<double>	Track_caloEMDeltaRp5(2000,0);
std::vector<double>	Track_caloHadDeltaRp3(2000,0);
std::vector<double>	Track_caloHadDeltaRp4(2000,0);
std::vector<double>	Track_caloHadDeltaRp5(2000,0);
std::vector<double>	Track_eta(2000,0);
std::vector<unsigned short>	Track_hitPattern_trackerLayersWithoutMeasurement(2000,0);
std::vector<unsigned short>	Track_numberOfValidHits(2000,0);
std::vector<double>	Track_phi(2000,0);
std::vector<double>	Track_pt(2000,0);
std::vector<double>	Track_px(2000,0);
std::vector<double>	Track_py(2000,0);
std::vector<double>	Track_pz(2000,0);
std::vector<int>	Track_trackHighPurity(2000,0);
std::vector<double>	Track_trackRelIso03(2000,0);
std::vector<unsigned short>	Track_trackerExpectedHitsInner_numberOfLostHits(2000,0);
std::vector<unsigned short>	Track_trackerExpectedHitsOuter_numberOfHits(2000,0);
std::vector<double>	Track_vx(2000,0);
std::vector<double>	Track_vy(2000,0);
std::vector<double>	Track_vz(2000,0);
std::vector<double>	Vertex_ndof(200,0);
std::vector<double>	Vertex_position_rho(200,0);
std::vector<double>	Vertex_x(200,0);
std::vector<double>	Vertex_y(200,0);
std::vector<double>	Vertex_z(200,0);
int	edmEventHelper_bunchCrossing;
int	edmEventHelper_event;
int	edmEventHelper_isRealData;
int	edmEventHelper_luminosityBlock;
int	edmEventHelper_orbitNumber;
int	edmEventHelper_run;
int	edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v3;
int	edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5;
int	edmTriggerResultsHelper_prescaleHLT_MET120_HBHENoiseCleaned_v3;
int	edmTriggerResultsHelper_prescaleHLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5;
int	nElectron;
int	nGenParticle;
int	nJet;
int	nMET;
int	nMuon;
int	nPileupSummaryInfo;
int	nTau;
int	nTrack;
int	nVertex;
double	sdouble_value;
double weight;

//-----------------------------------------------------------------------------
// --- indexmap keeps track of which objects have been flagged for selection
// --- IMPORTANT: initialize must be called every event to clear selection
std::map<std::string, std::vector<int> > indexmap;
void initialize()
{
  for(std::map<std::string, std::vector<int> >::iterator
    item=indexmap.begin(); 
    item != indexmap.end();
	++item)
	item->second.clear();
}

void select(std::string objname)
{
  indexmap[objname] = std::vector<int>();
}

void select(std::string objname, int index)
{
  try
    {
      indexmap[objname].push_back(index);
    }
  catch (...)
    {
      std::cout << "*** perhaps you failed to call select for " 
                << objname << std::endl;
      assert(0);
    }
}

//-----------------------------------------------------------------------------
// --- Structs can be filled by calling fillObjects()
// --- after the call to stream.read(...)
//-----------------------------------------------------------------------------
struct Electron_s
{
  double	energy;
  double	et;
  double	pz;
  double	pt;
  double	phi;
  double	eta;
  float	mvaNonTrigV0;
};
std::vector<Electron_s> Electron(200);

std::ostream& operator<<(std::ostream& os, const Electron_s& o)
{
  char r[1024];
  os << "Electron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "mvaNonTrigV0", (double)o.mvaNonTrigV0); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct GenParticle_s
{
  int	charge;
  double	p;
  double	energy;
  double	et;
  double	pz;
  double	pt;
  double	phi;
  double	eta;
  double	mass;
  int	pdgId;
};
std::vector<GenParticle_s> GenParticle(2000);

std::ostream& operator<<(std::ostream& os, const GenParticle_s& o)
{
  char r[1024];
  os << "GenParticle" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Jet_s
{
  double	energy;
  double	et;
  double	pt;
  double	pz;
  double	phi;
  double	eta;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
};
std::vector<Jet_s> Jet(200);

std::ostream& operator<<(std::ostream& os, const Jet_s& o)
{
  char r[1024];
  os << "Jet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct MET_s
{
  double	energy;
  double	et;
  double	pz;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<MET_s> MET(200);

std::ostream& operator<<(std::ostream& os, const MET_s& o)
{
  char r[1024];
  os << "MET" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Muon_s
{
  double	energy;
  double	et;
  double	pz;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<Muon_s> Muon(200);

std::ostream& operator<<(std::ostream& os, const Muon_s& o)
{
  char r[1024];
  os << "Muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PileupSummaryInfo_s
{
  int	getBunchCrossing;
  int	getPU_NumInteractions;
  float	getTrueNumInteractions;
};
std::vector<PileupSummaryInfo_s> PileupSummaryInfo(10);

std::ostream& operator<<(std::ostream& os, const PileupSummaryInfo_s& o)
{
  char r[1024];
  os << "PileupSummaryInfo" << std::endl;
  sprintf(r, "  %-32s: %f\n", "getBunchCrossing", (double)o.getBunchCrossing); os << r;
  sprintf(r, "  %-32s: %f\n", "getPU_NumInteractions", (double)o.getPU_NumInteractions); os << r;
  sprintf(r, "  %-32s: %f\n", "getTrueNumInteractions", (double)o.getTrueNumInteractions); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Tau_s
{
  double	energy;
  double	et;
  double	pz;
  double	pt;
  double	phi;
  double	eta;
  float	byLooseCombinedIsolationDeltaBetaCorr;
  float	decayModeFinding;
  float	againstElectronLoose;
  float	againstMuonTight;
};
std::vector<Tau_s> Tau(200);

std::ostream& operator<<(std::ostream& os, const Tau_s& o)
{
  char r[1024];
  os << "Tau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "byLooseCombinedIsolationDeltaBetaCorr", (double)o.byLooseCombinedIsolationDeltaBetaCorr); os << r;
  sprintf(r, "  %-32s: %f\n", "decayModeFinding", (double)o.decayModeFinding); os << r;
  sprintf(r, "  %-32s: %f\n", "againstElectronLoose", (double)o.againstElectronLoose); os << r;
  sprintf(r, "  %-32s: %f\n", "againstMuonTight", (double)o.againstMuonTight); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Track_s
{
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	phi;
  double	eta;
  double	vx;
  double	vy;
  double	vz;
  unsigned short	numberOfValidHits;
  unsigned short	hitPattern_trackerLayersWithoutMeasurement;
  unsigned short	trackerExpectedHitsInner_numberOfLostHits;
  unsigned short	trackerExpectedHitsOuter_numberOfHits;
  int	trackHighPurity;
  double	trackRelIso03;
  double	caloEMDeltaRp3;
  double	caloHadDeltaRp3;
  double	caloEMDeltaRp4;
  double	caloHadDeltaRp4;
  double	caloEMDeltaRp5;
  double	caloHadDeltaRp5;
  int pdgId;
  double beta;
};
std::vector<Track_s> Track(2000);

std::ostream& operator<<(std::ostream& os, const Track_s& o)
{
  char r[1024];
  os << "Track" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "vx", (double)o.vx); os << r;
  sprintf(r, "  %-32s: %f\n", "vy", (double)o.vy); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfValidHits", (double)o.numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "hitPattern_trackerLayersWithoutMeasurement", (double)o.hitPattern_trackerLayersWithoutMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "trackerExpectedHitsInner_numberOfLostHits", (double)o.trackerExpectedHitsInner_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "trackerExpectedHitsOuter_numberOfHits", (double)o.trackerExpectedHitsOuter_numberOfHits); os << r;
  sprintf(r, "  %-32s: %f\n", "trackHighPurity", (double)o.trackHighPurity); os << r;
  sprintf(r, "  %-32s: %f\n", "trackRelIso03", (double)o.trackRelIso03); os << r;
  sprintf(r, "  %-32s: %f\n", "caloEMDeltaRp3", (double)o.caloEMDeltaRp3); os << r;
  sprintf(r, "  %-32s: %f\n", "caloHadDeltaRp3", (double)o.caloHadDeltaRp3); os << r;
  sprintf(r, "  %-32s: %f\n", "caloEMDeltaRp4", (double)o.caloEMDeltaRp4); os << r;
  sprintf(r, "  %-32s: %f\n", "caloHadDeltaRp4", (double)o.caloHadDeltaRp4); os << r;
  sprintf(r, "  %-32s: %f\n", "caloEMDeltaRp5", (double)o.caloEMDeltaRp5); os << r;
  sprintf(r, "  %-32s: %f\n", "caloHadDeltaRp5", (double)o.caloHadDeltaRp5); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Vertex_s
{
  double	ndof;
  double	x;
  double	y;
  double	z;
  double	position_rho;
};
std::vector<Vertex_s> Vertex(200);

std::ostream& operator<<(std::ostream& os, const Vertex_s& o)
{
  char r[1024];
  os << "Vertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  sprintf(r, "  %-32s: %f\n", "position_rho", (double)o.position_rho); os << r;
  return os;
}
//-----------------------------------------------------------------------------

inline void fillElectron()
{
  Electron.resize(Electron_energy.size());
  for(unsigned int i=0; i < Electron.size(); ++i)
    {
      Electron[i].energy	= Electron_energy[i];
      Electron[i].et	= Electron_et[i];
      Electron[i].pz	= Electron_pz[i];
      Electron[i].pt	= Electron_pt[i];
      Electron[i].phi	= Electron_phi[i];
      Electron[i].eta	= Electron_eta[i];
      Electron[i].mvaNonTrigV0	= Electron_mvaNonTrigV0[i];
    }
}

inline void fillGenParticle()
{
  GenParticle.resize(GenParticle_charge.size());
  for(unsigned int i=0; i < GenParticle.size(); ++i)
    {
      GenParticle[i].charge	= GenParticle_charge[i];
      GenParticle[i].p	= GenParticle_p[i];
      GenParticle[i].energy	= GenParticle_energy[i];
      GenParticle[i].et	= GenParticle_et[i];
      GenParticle[i].pz	= GenParticle_pz[i];
      GenParticle[i].pt	= GenParticle_pt[i];
      GenParticle[i].phi	= GenParticle_phi[i];
      GenParticle[i].eta	= GenParticle_eta[i];
      GenParticle[i].mass	= GenParticle_mass[i];
      GenParticle[i].pdgId	= GenParticle_pdgId[i];
    }
}

inline void fillJet()
{
  Jet.resize(Jet_energy.size());
  for(unsigned int i=0; i < Jet.size(); ++i)
    {
      Jet[i].energy	= Jet_energy[i];
      Jet[i].et	= Jet_et[i];
      Jet[i].pt	= Jet_pt[i];
      Jet[i].pz	= Jet_pz[i];
      Jet[i].phi	= Jet_phi[i];
      Jet[i].eta	= Jet_eta[i];
      Jet[i].chargedHadronEnergyFraction	= Jet_chargedHadronEnergyFraction[i];
      Jet[i].neutralHadronEnergyFraction	= Jet_neutralHadronEnergyFraction[i];
      Jet[i].chargedEmEnergyFraction	= Jet_chargedEmEnergyFraction[i];
      Jet[i].neutralEmEnergyFraction	= Jet_neutralEmEnergyFraction[i];
    }
}

inline void fillMET()
{
  MET.resize(MET_energy.size());
  for(unsigned int i=0; i < MET.size(); ++i)
    {
      MET[i].energy	= MET_energy[i];
      MET[i].et	= MET_et[i];
      MET[i].pz	= MET_pz[i];
      MET[i].pt	= MET_pt[i];
      MET[i].phi	= MET_phi[i];
      MET[i].eta	= MET_eta[i];
    }
}

inline void fillMuon()
{
  Muon.resize(Muon_energy.size());
  for(unsigned int i=0; i < Muon.size(); ++i)
    {
      Muon[i].energy	= Muon_energy[i];
      Muon[i].et	= Muon_et[i];
      Muon[i].pz	= Muon_pz[i];
      Muon[i].pt	= Muon_pt[i];
      Muon[i].phi	= Muon_phi[i];
      Muon[i].eta	= Muon_eta[i];
    }
}

inline void fillPileupSummaryInfo()
{
  PileupSummaryInfo.resize(PileupSummaryInfo_getBunchCrossing.size());
  for(unsigned int i=0; i < PileupSummaryInfo.size(); ++i)
    {
      PileupSummaryInfo[i].getBunchCrossing	= PileupSummaryInfo_getBunchCrossing[i];
      PileupSummaryInfo[i].getPU_NumInteractions	= PileupSummaryInfo_getPU_NumInteractions[i];
      PileupSummaryInfo[i].getTrueNumInteractions	= PileupSummaryInfo_getTrueNumInteractions[i];
    }
}

inline void fillTau()
{
  Tau.resize(Tau_energy.size());
  for(unsigned int i=0; i < Tau.size(); ++i)
    {
      Tau[i].energy	= Tau_energy[i];
      Tau[i].et	= Tau_et[i];
      Tau[i].pz	= Tau_pz[i];
      Tau[i].pt	= Tau_pt[i];
      Tau[i].phi	= Tau_phi[i];
      Tau[i].eta	= Tau_eta[i];
      Tau[i].byLooseCombinedIsolationDeltaBetaCorr	= Tau_byLooseCombinedIsolationDeltaBetaCorr[i];
      Tau[i].decayModeFinding	= Tau_decayModeFinding[i];
      Tau[i].againstElectronLoose	= Tau_againstElectronLoose[i];
      Tau[i].againstMuonTight	= Tau_againstMuonTight[i];
    }
}

inline void fillTrack()
{
  Track.resize(Track_pt.size());
  for(unsigned int i=0; i < Track.size(); ++i)
    {
      Track[i].pt	= Track_pt[i];
      Track[i].px	= Track_px[i];
      Track[i].py	= Track_py[i];
      Track[i].pz	= Track_pz[i];
      Track[i].phi	= Track_phi[i];
      Track[i].eta	= Track_eta[i];
      Track[i].vx	= Track_vx[i];
      Track[i].vy	= Track_vy[i];
      Track[i].vz	= Track_vz[i];
      Track[i].numberOfValidHits	= Track_numberOfValidHits[i];
      Track[i].hitPattern_trackerLayersWithoutMeasurement	= Track_hitPattern_trackerLayersWithoutMeasurement[i];
      Track[i].trackerExpectedHitsInner_numberOfLostHits	= Track_trackerExpectedHitsInner_numberOfLostHits[i];
      Track[i].trackerExpectedHitsOuter_numberOfHits	= Track_trackerExpectedHitsOuter_numberOfHits[i];
      Track[i].trackHighPurity	= Track_trackHighPurity[i];
      Track[i].trackRelIso03	= Track_trackRelIso03[i];
      Track[i].caloEMDeltaRp3	= Track_caloEMDeltaRp3[i];
      Track[i].caloHadDeltaRp3	= Track_caloHadDeltaRp3[i];
      Track[i].caloEMDeltaRp4	= Track_caloEMDeltaRp4[i];
      Track[i].caloHadDeltaRp4	= Track_caloHadDeltaRp4[i];
      Track[i].caloEMDeltaRp5	= Track_caloEMDeltaRp5[i];
      Track[i].caloHadDeltaRp5	= Track_caloHadDeltaRp5[i];
    }
}

inline void fillVertex()
{
  Vertex.resize(Vertex_ndof.size());
  for(unsigned int i=0; i < Vertex.size(); ++i)
    {
      Vertex[i].ndof	= Vertex_ndof[i];
      Vertex[i].x	= Vertex_x[i];
      Vertex[i].y	= Vertex_y[i];
      Vertex[i].z	= Vertex_z[i];
      Vertex[i].position_rho	= Vertex_position_rho[i];
    }
}


void fillObjects()
{
  fillElectron();
  fillGenParticle();
  fillJet();
  fillMET();
  fillMuon();
  fillPileupSummaryInfo();
  fillTau();
  fillTrack();
  fillVertex();
}

//-----------------------------------------------------------------------------
// --- Call saveSelectedObjects() just before call to addEvent if
// --- you wish to save only the selected objects
//-----------------------------------------------------------------------------
// Select objects for which the select function was called
void saveSelectedObjects()
{
  int n = 0;

  n = 0;
  try
    {
       n = indexmap["Electron"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Electron"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Electron_energy[i]	= Electron_energy[j];
          Electron_et[i]	= Electron_et[j];
          Electron_pz[i]	= Electron_pz[j];
          Electron_pt[i]	= Electron_pt[j];
          Electron_phi[i]	= Electron_phi[j];
          Electron_eta[i]	= Electron_eta[j];
          Electron_mvaNonTrigV0[i]	= Electron_mvaNonTrigV0[j];
        }
      nElectron = n;
    }

  n = 0;
  try
    {
       n = indexmap["GenParticle"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["GenParticle"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          GenParticle_charge[i]	= GenParticle_charge[j];
          GenParticle_p[i]	= GenParticle_p[j];
          GenParticle_energy[i]	= GenParticle_energy[j];
          GenParticle_et[i]	= GenParticle_et[j];
          GenParticle_pz[i]	= GenParticle_pz[j];
          GenParticle_pt[i]	= GenParticle_pt[j];
          GenParticle_phi[i]	= GenParticle_phi[j];
          GenParticle_eta[i]	= GenParticle_eta[j];
          GenParticle_mass[i]	= GenParticle_mass[j];
          GenParticle_pdgId[i]	= GenParticle_pdgId[j];
        }
      nGenParticle = n;
    }

  n = 0;
  try
    {
       n = indexmap["Jet"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Jet"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Jet_energy[i]	= Jet_energy[j];
          Jet_et[i]	= Jet_et[j];
          Jet_pt[i]	= Jet_pt[j];
          Jet_pz[i]	= Jet_pz[j];
          Jet_phi[i]	= Jet_phi[j];
          Jet_eta[i]	= Jet_eta[j];
          Jet_chargedHadronEnergyFraction[i]	= Jet_chargedHadronEnergyFraction[j];
          Jet_neutralHadronEnergyFraction[i]	= Jet_neutralHadronEnergyFraction[j];
          Jet_chargedEmEnergyFraction[i]	= Jet_chargedEmEnergyFraction[j];
          Jet_neutralEmEnergyFraction[i]	= Jet_neutralEmEnergyFraction[j];
        }
      nJet = n;
    }

  n = 0;
  try
    {
       n = indexmap["MET"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["MET"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          MET_energy[i]	= MET_energy[j];
          MET_et[i]	= MET_et[j];
          MET_pz[i]	= MET_pz[j];
          MET_pt[i]	= MET_pt[j];
          MET_phi[i]	= MET_phi[j];
          MET_eta[i]	= MET_eta[j];
        }
      nMET = n;
    }

  n = 0;
  try
    {
       n = indexmap["Muon"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Muon"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Muon_energy[i]	= Muon_energy[j];
          Muon_et[i]	= Muon_et[j];
          Muon_pz[i]	= Muon_pz[j];
          Muon_pt[i]	= Muon_pt[j];
          Muon_phi[i]	= Muon_phi[j];
          Muon_eta[i]	= Muon_eta[j];
        }
      nMuon = n;
    }

  n = 0;
  try
    {
       n = indexmap["PileupSummaryInfo"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PileupSummaryInfo"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PileupSummaryInfo_getBunchCrossing[i]	= PileupSummaryInfo_getBunchCrossing[j];
          PileupSummaryInfo_getPU_NumInteractions[i]	= PileupSummaryInfo_getPU_NumInteractions[j];
          PileupSummaryInfo_getTrueNumInteractions[i]	= PileupSummaryInfo_getTrueNumInteractions[j];
        }
      nPileupSummaryInfo = n;
    }

  n = 0;
  try
    {
       n = indexmap["Tau"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Tau"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Tau_energy[i]	= Tau_energy[j];
          Tau_et[i]	= Tau_et[j];
          Tau_pz[i]	= Tau_pz[j];
          Tau_pt[i]	= Tau_pt[j];
          Tau_phi[i]	= Tau_phi[j];
          Tau_eta[i]	= Tau_eta[j];
          Tau_byLooseCombinedIsolationDeltaBetaCorr[i]	= Tau_byLooseCombinedIsolationDeltaBetaCorr[j];
          Tau_decayModeFinding[i]	= Tau_decayModeFinding[j];
          Tau_againstElectronLoose[i]	= Tau_againstElectronLoose[j];
          Tau_againstMuonTight[i]	= Tau_againstMuonTight[j];
        }
      nTau = n;
    }

  n = 0;
  try
    {
       n = indexmap["Track"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Track"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Track_pt[i]	= Track_pt[j];
          Track_px[i]	= Track_px[j];
          Track_py[i]	= Track_py[j];
          Track_pz[i]	= Track_pz[j];
          Track_phi[i]	= Track_phi[j];
          Track_eta[i]	= Track_eta[j];
          Track_vx[i]	= Track_vx[j];
          Track_vy[i]	= Track_vy[j];
          Track_vz[i]	= Track_vz[j];
          Track_numberOfValidHits[i]	= Track_numberOfValidHits[j];
          Track_hitPattern_trackerLayersWithoutMeasurement[i]	= Track_hitPattern_trackerLayersWithoutMeasurement[j];
          Track_trackerExpectedHitsInner_numberOfLostHits[i]	= Track_trackerExpectedHitsInner_numberOfLostHits[j];
          Track_trackerExpectedHitsOuter_numberOfHits[i]	= Track_trackerExpectedHitsOuter_numberOfHits[j];
          Track_trackHighPurity[i]	= Track_trackHighPurity[j];
          Track_trackRelIso03[i]	= Track_trackRelIso03[j];
          Track_caloEMDeltaRp3[i]	= Track_caloEMDeltaRp3[j];
          Track_caloHadDeltaRp3[i]	= Track_caloHadDeltaRp3[j];
          Track_caloEMDeltaRp4[i]	= Track_caloEMDeltaRp4[j];
          Track_caloHadDeltaRp4[i]	= Track_caloHadDeltaRp4[j];
          Track_caloEMDeltaRp5[i]	= Track_caloEMDeltaRp5[j];
          Track_caloHadDeltaRp5[i]	= Track_caloHadDeltaRp5[j];
        }
      nTrack = n;
    }

  n = 0;
  try
    {
       n = indexmap["Vertex"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["Vertex"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          Vertex_ndof[i]	= Vertex_ndof[j];
          Vertex_x[i]	= Vertex_x[j];
          Vertex_y[i]	= Vertex_y[j];
          Vertex_z[i]	= Vertex_z[j];
          Vertex_position_rho[i]	= Vertex_position_rho[j];
        }
      nVertex = n;
    }
}

//-----------------------------------------------------------------------------
// --- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("patElectron_patElectronsLoosePFlow.energy", Electron_energy);
  stream.select("patElectron_patElectronsLoosePFlow.et", Electron_et);
  stream.select("patElectron_patElectronsLoosePFlow.eta", Electron_eta);
  stream.select("patElectron_patElectronsLoosePFlow.mvaNonTrigV0", Electron_mvaNonTrigV0);
  stream.select("patElectron_patElectronsLoosePFlow.phi", Electron_phi);
  stream.select("patElectron_patElectronsLoosePFlow.pt", Electron_pt);
  stream.select("patElectron_patElectronsLoosePFlow.pz", Electron_pz);
  stream.select("recoGenParticle_genParticles.charge", GenParticle_charge);
  stream.select("recoGenParticle_genParticles.energy", GenParticle_energy);
  stream.select("recoGenParticle_genParticles.et", GenParticle_et);
  stream.select("recoGenParticle_genParticles.eta", GenParticle_eta);
  stream.select("recoGenParticle_genParticles.mass", GenParticle_mass);
  stream.select("recoGenParticle_genParticles.p", GenParticle_p);
  stream.select("recoGenParticle_genParticles.pdgId", GenParticle_pdgId);
  stream.select("recoGenParticle_genParticles.phi", GenParticle_phi);
  stream.select("recoGenParticle_genParticles.pt", GenParticle_pt);
  stream.select("recoGenParticle_genParticles.pz", GenParticle_pz);
  stream.select("patJet_selectedPatJetsPFlow.chargedEmEnergyFraction", Jet_chargedEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.chargedHadronEnergyFraction", Jet_chargedHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.energy", Jet_energy);
  stream.select("patJet_selectedPatJetsPFlow.et", Jet_et);
  stream.select("patJet_selectedPatJetsPFlow.eta", Jet_eta);
  stream.select("patJet_selectedPatJetsPFlow.neutralEmEnergyFraction", Jet_neutralEmEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.neutralHadronEnergyFraction", Jet_neutralHadronEnergyFraction);
  stream.select("patJet_selectedPatJetsPFlow.phi", Jet_phi);
  stream.select("patJet_selectedPatJetsPFlow.pt", Jet_pt);
  stream.select("patJet_selectedPatJetsPFlow.pz", Jet_pz);
  stream.select("patMET_patMETsPFlow.energy", MET_energy);
  stream.select("patMET_patMETsPFlow.et", MET_et);
  stream.select("patMET_patMETsPFlow.eta", MET_eta);
  stream.select("patMET_patMETsPFlow.phi", MET_phi);
  stream.select("patMET_patMETsPFlow.pt", MET_pt);
  stream.select("patMET_patMETsPFlow.pz", MET_pz);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.energy", Muon_energy);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.et", Muon_et);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.eta", Muon_eta);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.phi", Muon_phi);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.pt", Muon_pt);
  stream.select("patMuon_selectedPatMuonsLoosePFlow.pz", Muon_pz);
  stream.select("PileupSummaryInfo_addPileupInfo.getBunchCrossing", PileupSummaryInfo_getBunchCrossing);
  stream.select("PileupSummaryInfo_addPileupInfo.getPU_NumInteractions", PileupSummaryInfo_getPU_NumInteractions);
  stream.select("PileupSummaryInfo_addPileupInfo.getTrueNumInteractions", PileupSummaryInfo_getTrueNumInteractions);
  stream.select("patTau_selectedPatTaus.againstElectronLoose", Tau_againstElectronLoose);
  stream.select("patTau_selectedPatTaus.againstMuonTight", Tau_againstMuonTight);
  stream.select("patTau_selectedPatTaus.byLooseCombinedIsolationDeltaBetaCorr", Tau_byLooseCombinedIsolationDeltaBetaCorr);
  stream.select("patTau_selectedPatTaus.decayModeFinding", Tau_decayModeFinding);
  stream.select("patTau_selectedPatTaus.energy", Tau_energy);
  stream.select("patTau_selectedPatTaus.et", Tau_et);
  stream.select("patTau_selectedPatTaus.eta", Tau_eta);
  stream.select("patTau_selectedPatTaus.phi", Tau_phi);
  stream.select("patTau_selectedPatTaus.pt", Tau_pt);
  stream.select("patTau_selectedPatTaus.pz", Tau_pz);
  stream.select("recoTrackHelper_generalTracks.caloEMDeltaRp3", Track_caloEMDeltaRp3);
  stream.select("recoTrackHelper_generalTracks.caloEMDeltaRp4", Track_caloEMDeltaRp4);
  stream.select("recoTrackHelper_generalTracks.caloEMDeltaRp5", Track_caloEMDeltaRp5);
  stream.select("recoTrackHelper_generalTracks.caloHadDeltaRp3", Track_caloHadDeltaRp3);
  stream.select("recoTrackHelper_generalTracks.caloHadDeltaRp4", Track_caloHadDeltaRp4);
  stream.select("recoTrackHelper_generalTracks.caloHadDeltaRp5", Track_caloHadDeltaRp5);
  stream.select("recoTrackHelper_generalTracks.eta", Track_eta);
  stream.select("recoTrackHelper_generalTracks.hitPattern_trackerLayersWithoutMeasurement", Track_hitPattern_trackerLayersWithoutMeasurement);
  stream.select("recoTrackHelper_generalTracks.numberOfValidHits", Track_numberOfValidHits);
  stream.select("recoTrackHelper_generalTracks.phi", Track_phi);
  stream.select("recoTrackHelper_generalTracks.pt", Track_pt);
  stream.select("recoTrackHelper_generalTracks.px", Track_px);
  stream.select("recoTrackHelper_generalTracks.py", Track_py);
  stream.select("recoTrackHelper_generalTracks.pz", Track_pz);
  stream.select("recoTrackHelper_generalTracks.trackHighPurity", Track_trackHighPurity);
  stream.select("recoTrackHelper_generalTracks.trackRelIso03", Track_trackRelIso03);
  stream.select("recoTrackHelper_generalTracks.trackerExpectedHitsInner_numberOfLostHits", Track_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("recoTrackHelper_generalTracks.trackerExpectedHitsOuter_numberOfHits", Track_trackerExpectedHitsOuter_numberOfHits);
  stream.select("recoTrackHelper_generalTracks.vx", Track_vx);
  stream.select("recoTrackHelper_generalTracks.vy", Track_vy);
  stream.select("recoTrackHelper_generalTracks.vz", Track_vz);
  stream.select("recoVertex_offlinePrimaryVertices.ndof", Vertex_ndof);
  stream.select("recoVertex_offlinePrimaryVertices.position_rho", Vertex_position_rho);
  stream.select("recoVertex_offlinePrimaryVertices.x", Vertex_x);
  stream.select("recoVertex_offlinePrimaryVertices.y", Vertex_y);
  stream.select("recoVertex_offlinePrimaryVertices.z", Vertex_z);
  stream.select("edmEventHelper_info.bunchCrossing", edmEventHelper_bunchCrossing);
  stream.select("edmEventHelper_info.event", edmEventHelper_event);
  stream.select("edmEventHelper_info.isRealData", edmEventHelper_isRealData);
  stream.select("edmEventHelper_info.luminosityBlock", edmEventHelper_luminosityBlock);
  stream.select("edmEventHelper_info.orbitNumber", edmEventHelper_orbitNumber);
  stream.select("edmEventHelper_info.run", edmEventHelper_run);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_MET120_HBHENoiseCleaned_v3", edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5", edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_MET120_HBHENoiseCleaned_v3", edmTriggerResultsHelper_prescaleHLT_MET120_HBHENoiseCleaned_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5", edmTriggerResultsHelper_prescaleHLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5);
  stream.select("npatElectron_patElectronsLoosePFlow", nElectron);
  stream.select("nrecoGenParticle_genParticles", nGenParticle);
  stream.select("npatJet_selectedPatJetsPFlow", nJet);
  stream.select("npatMET_patMETsPFlow", nMET);
  stream.select("npatMuon_selectedPatMuonsLoosePFlow", nMuon);
  stream.select("nPileupSummaryInfo_addPileupInfo", nPileupSummaryInfo);
  stream.select("npatTau_selectedPatTaus", nTau);
  stream.select("nrecoTrackHelper_generalTracks", nTrack);
  stream.select("nrecoVertex_offlinePrimaryVertices", nVertex);
  stream.select("sdouble_kt6CaloJets_rho.value", sdouble_value);

}
}; // end namespace evt
#endif
