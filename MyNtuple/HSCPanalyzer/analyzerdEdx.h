#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Fri Feb 21 12:20:25 2014 by mkanalyzer.py
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
std::vector<int>	GenParticle_charge(1000,0);
std::vector<double>	GenParticle_energy(1000,0);
std::vector<double>	GenParticle_et(1000,0);
std::vector<double>	GenParticle_eta(1000,0);
std::vector<double>	GenParticle_mass(1000,0);
std::vector<size_t>	GenParticle_numberOfDaughters(1000,0);
std::vector<size_t>	GenParticle_numberOfMothers(1000,0);
std::vector<double>	GenParticle_p(1000,0);
std::vector<int>	GenParticle_pdgId(1000,0);
std::vector<double>	GenParticle_phi(1000,0);
std::vector<double>	GenParticle_pt(1000,0);
std::vector<int>	GenParticle_status(1000,0);
std::vector<double>	GenParticle_vx(1000,0);
std::vector<double>	GenParticle_vy(1000,0);
std::vector<double>	GenParticle_vz(1000,0);
std::vector<float>	PFJet_chargedEmEnergyFraction(200,0);
std::vector<float>	PFJet_chargedHadronEnergyFraction(200,0);
std::vector<double>	PFJet_energy(200,0);
std::vector<double>	PFJet_et(200,0);
std::vector<double>	PFJet_eta(200,0);
std::vector<float>	PFJet_neutralEmEnergyFraction(200,0);
std::vector<float>	PFJet_neutralHadronEnergyFraction(200,0);
std::vector<double>	PFJet_p(200,0);
std::vector<double>	PFJet_phi(200,0);
std::vector<double>	PFJet_pt(200,0);
std::vector<double>	PFMET_energy(200,0);
std::vector<double>	PFMET_et(200,0);
std::vector<double>	PFMET_eta(200,0);
std::vector<double>	PFMET_p(200,0);
std::vector<double>	PFMET_phi(200,0);
std::vector<double>	PFMET_pt(200,0);
std::vector<double>	Track_d0(2000000,0);
std::vector<double>	Track_d0Error(2000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_1(2000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_1000(2000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_2(2000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_3(2000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_5(2000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_7(2000000,0);
std::vector<double>	Track_dEdxHitsNPMedian_1(2000000,0);
std::vector<double>	Track_dEdxHitsNPMedian_1000(2000000,0);
std::vector<double>	Track_dEdxHitsNPMedian_2(2000000,0);
std::vector<double>	Track_dEdxHitsNPMedian_3(2000000,0);
std::vector<double>	Track_dEdxHitsNPMedian_5(2000000,0);
std::vector<double>	Track_dEdxHitsNPMedian_7(2000000,0);
std::vector<double>	Track_dEdxHitsNPTrun40_1(2000000,0);
std::vector<double>	Track_dEdxHitsNPTrun40_1000(2000000,0);
std::vector<double>	Track_dEdxHitsNPTrun40_2(2000000,0);
std::vector<double>	Track_dEdxHitsNPTrun40_3(2000000,0);
std::vector<double>	Track_dEdxHitsNPTrun40_5(2000000,0);
std::vector<double>	Track_dEdxHitsNPTrun40_7(2000000,0);
std::vector<double>	Track_dEdxNPHarm2(2000000,0);
std::vector<unsigned int>	Track_dEdxNPNoM(2000000,0);
std::vector<double>	Track_dEdxNPTru40(2000000,0);
std::vector<double>	Track_dz(2000000,0);
std::vector<double>	Track_dzError(2000000,0);
std::vector<double>	Track_eta(2000000,0);
std::vector<unsigned short>	Track_numberOfLostHits(2000000,0);
std::vector<unsigned short>	Track_numberOfValidHits(2000000,0);
std::vector<double>	Track_phi(2000000,0);
std::vector<double>	Track_pt(2000000,0);
std::vector<double>	Track_ptError(2000000,0);
std::vector<unsigned short>	Track_trackerExpectedHitsInner_numberOfLostHits(2000000,0);
std::vector<unsigned short>	Track_trackerExpectedHitsOuter_numberOfLostHits(2000000,0);
std::vector<double>	Vertex_chi2(200,0);
std::vector<int>	Vertex_isFake(200,0);
std::vector<int>	Vertex_isValid(200,0);
std::vector<double>	Vertex_ndof(200,0);
std::vector<double>	Vertex_normalizedChi2(200,0);
std::vector<size_t>	Vertex_tracksSize(200,0);
std::vector<double>	Vertex_x(200,0);
std::vector<double>	Vertex_xError(200,0);
std::vector<double>	Vertex_y(200,0);
std::vector<double>	Vertex_yError(200,0);
std::vector<double>	Vertex_z(200,0);
std::vector<double>	Vertex_zError(200,0);
int	edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5;
int	edmTriggerResultsHelper_prescaleHLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5;
int	nGenParticle;
int	nPFJet;
int	nPFMET;
int	nTrack;
int	nVertex;

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
struct GenParticle_s
{
  int	charge;
  double	energy;
  double	et;
  double	p;
  double	pt;
  double	phi;
  double	eta;
  size_t	numberOfDaughters;
  size_t	numberOfMothers;
  double	mass;
  double	vx;
  double	vy;
  double	vz;
  int	pdgId;
  int	status;
};
std::vector<GenParticle_s> GenParticle(1000);

std::ostream& operator<<(std::ostream& os, const GenParticle_s& o)
{
  char r[1024];
  os << "GenParticle" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfMothers", (double)o.numberOfMothers); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "vx", (double)o.vx); os << r;
  sprintf(r, "  %-32s: %f\n", "vy", (double)o.vy); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "status", (double)o.status); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFJet_s
{
  double	energy;
  double	et;
  double	p;
  double	pt;
  double	phi;
  double	eta;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	neutralEmEnergyFraction;
  float	chargedEmEnergyFraction;
};
std::vector<PFJet_s> PFJet(200);

std::ostream& operator<<(std::ostream& os, const PFJet_s& o)
{
  char r[1024];
  os << "PFJet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFMET_s
{
  double	energy;
  double	et;
  double	p;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<PFMET_s> PFMET(200);

std::ostream& operator<<(std::ostream& os, const PFMET_s& o)
{
  char r[1024];
  os << "PFMET" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Track_s
{
  double	pt;
  double	eta;
  double	phi;
  unsigned short	numberOfValidHits;
  unsigned short	numberOfLostHits;
  double	d0;
  double	dz;
  double	d0Error;
  double	dzError;
  double	ptError;
  double	dEdxNPHarm2;
  double	dEdxNPTru40;
  unsigned int	dEdxNPNoM;
  double	dEdxHitsNPHarm2_1000;
  double	dEdxHitsNPHarm2_7;
  double	dEdxHitsNPHarm2_5;
  double	dEdxHitsNPHarm2_3;
  double	dEdxHitsNPHarm2_2;
  double	dEdxHitsNPHarm2_1;
  double	dEdxHitsNPTrun40_1000;
  double	dEdxHitsNPTrun40_7;
  double	dEdxHitsNPTrun40_5;
  double	dEdxHitsNPTrun40_3;
  double	dEdxHitsNPTrun40_2;
  double	dEdxHitsNPTrun40_1;
  double	dEdxHitsNPMedian_1000;
  double	dEdxHitsNPMedian_7;
  double	dEdxHitsNPMedian_5;
  double	dEdxHitsNPMedian_3;
  double	dEdxHitsNPMedian_2;
  double	dEdxHitsNPMedian_1;
  unsigned short	trackerExpectedHitsOuter_numberOfLostHits;
  unsigned short	trackerExpectedHitsInner_numberOfLostHits;
  int pdgId;
  double beta;
};
std::vector<Track_s> Track(2000000);

std::ostream& operator<<(std::ostream& os, const Track_s& o)
{
  char r[1024];
  os << "Track" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfValidHits", (double)o.numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfLostHits", (double)o.numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "d0", (double)o.d0); os << r;
  sprintf(r, "  %-32s: %f\n", "dz", (double)o.dz); os << r;
  sprintf(r, "  %-32s: %f\n", "d0Error", (double)o.d0Error); os << r;
  sprintf(r, "  %-32s: %f\n", "dzError", (double)o.dzError); os << r;
  sprintf(r, "  %-32s: %f\n", "ptError", (double)o.ptError); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxNPHarm2", (double)o.dEdxNPHarm2); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxNPTru40", (double)o.dEdxNPTru40); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxNPNoM", (double)o.dEdxNPNoM); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_1000", (double)o.dEdxHitsNPHarm2_1000); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_7", (double)o.dEdxHitsNPHarm2_7); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_5", (double)o.dEdxHitsNPHarm2_5); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_3", (double)o.dEdxHitsNPHarm2_3); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_2", (double)o.dEdxHitsNPHarm2_2); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_1", (double)o.dEdxHitsNPHarm2_1); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPTrun40_1000", (double)o.dEdxHitsNPTrun40_1000); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPTrun40_7", (double)o.dEdxHitsNPTrun40_7); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPTrun40_5", (double)o.dEdxHitsNPTrun40_5); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPTrun40_3", (double)o.dEdxHitsNPTrun40_3); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPTrun40_2", (double)o.dEdxHitsNPTrun40_2); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPTrun40_1", (double)o.dEdxHitsNPTrun40_1); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPMedian_1000", (double)o.dEdxHitsNPMedian_1000); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPMedian_7", (double)o.dEdxHitsNPMedian_7); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPMedian_5", (double)o.dEdxHitsNPMedian_5); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPMedian_3", (double)o.dEdxHitsNPMedian_3); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPMedian_2", (double)o.dEdxHitsNPMedian_2); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPMedian_1", (double)o.dEdxHitsNPMedian_1); os << r;
  sprintf(r, "  %-32s: %f\n", "trackerExpectedHitsOuter_numberOfLostHits", (double)o.trackerExpectedHitsOuter_numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "trackerExpectedHitsInner_numberOfLostHits", (double)o.trackerExpectedHitsInner_numberOfLostHits); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Vertex_s
{
  int	isValid;
  int	isFake;
  size_t	tracksSize;
  double	chi2;
  double	ndof;
  double	normalizedChi2;
  double	x;
  double	y;
  double	z;
  double	xError;
  double	yError;
  double	zError;
};
std::vector<Vertex_s> Vertex(200);

std::ostream& operator<<(std::ostream& os, const Vertex_s& o)
{
  char r[1024];
  os << "Vertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "isValid", (double)o.isValid); os << r;
  sprintf(r, "  %-32s: %f\n", "isFake", (double)o.isFake); os << r;
  sprintf(r, "  %-32s: %f\n", "tracksSize", (double)o.tracksSize); os << r;
  sprintf(r, "  %-32s: %f\n", "chi2", (double)o.chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "normalizedChi2", (double)o.normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  sprintf(r, "  %-32s: %f\n", "xError", (double)o.xError); os << r;
  sprintf(r, "  %-32s: %f\n", "yError", (double)o.yError); os << r;
  sprintf(r, "  %-32s: %f\n", "zError", (double)o.zError); os << r;
  return os;
}
//-----------------------------------------------------------------------------

inline void fillGenParticle()
{
  GenParticle.resize(GenParticle_charge.size());
  for(unsigned int i=0; i < GenParticle.size(); ++i)
    {
      GenParticle[i].charge	= GenParticle_charge[i];
      GenParticle[i].energy	= GenParticle_energy[i];
      GenParticle[i].et	= GenParticle_et[i];
      GenParticle[i].p	= GenParticle_p[i];
      GenParticle[i].pt	= GenParticle_pt[i];
      GenParticle[i].phi	= GenParticle_phi[i];
      GenParticle[i].eta	= GenParticle_eta[i];
      GenParticle[i].numberOfDaughters	= GenParticle_numberOfDaughters[i];
      GenParticle[i].numberOfMothers	= GenParticle_numberOfMothers[i];
      GenParticle[i].mass	= GenParticle_mass[i];
      GenParticle[i].vx	= GenParticle_vx[i];
      GenParticle[i].vy	= GenParticle_vy[i];
      GenParticle[i].vz	= GenParticle_vz[i];
      GenParticle[i].pdgId	= GenParticle_pdgId[i];
      GenParticle[i].status	= GenParticle_status[i];
    }
}

inline void fillPFJet()
{
  PFJet.resize(PFJet_energy.size());
  for(unsigned int i=0; i < PFJet.size(); ++i)
    {
      PFJet[i].energy	= PFJet_energy[i];
      PFJet[i].et	= PFJet_et[i];
      PFJet[i].p	= PFJet_p[i];
      PFJet[i].pt	= PFJet_pt[i];
      PFJet[i].phi	= PFJet_phi[i];
      PFJet[i].eta	= PFJet_eta[i];
      PFJet[i].chargedHadronEnergyFraction	= PFJet_chargedHadronEnergyFraction[i];
      PFJet[i].neutralHadronEnergyFraction	= PFJet_neutralHadronEnergyFraction[i];
      PFJet[i].neutralEmEnergyFraction	= PFJet_neutralEmEnergyFraction[i];
      PFJet[i].chargedEmEnergyFraction	= PFJet_chargedEmEnergyFraction[i];
    }
}

inline void fillPFMET()
{
  PFMET.resize(PFMET_energy.size());
  for(unsigned int i=0; i < PFMET.size(); ++i)
    {
      PFMET[i].energy	= PFMET_energy[i];
      PFMET[i].et	= PFMET_et[i];
      PFMET[i].p	= PFMET_p[i];
      PFMET[i].pt	= PFMET_pt[i];
      PFMET[i].phi	= PFMET_phi[i];
      PFMET[i].eta	= PFMET_eta[i];
    }
}

inline void fillTrack()
{
  Track.resize(Track_pt.size());
  for(unsigned int i=0; i < Track.size(); ++i)
    {
      Track[i].pt	= Track_pt[i];
      Track[i].eta	= Track_eta[i];
      Track[i].phi	= Track_phi[i];
      Track[i].numberOfValidHits	= Track_numberOfValidHits[i];
      Track[i].numberOfLostHits	= Track_numberOfLostHits[i];
      Track[i].d0	= Track_d0[i];
      Track[i].dz	= Track_dz[i];
      Track[i].d0Error	= Track_d0Error[i];
      Track[i].dzError	= Track_dzError[i];
      Track[i].ptError	= Track_ptError[i];
      Track[i].dEdxNPHarm2	= Track_dEdxNPHarm2[i];
      Track[i].dEdxNPTru40	= Track_dEdxNPTru40[i];
      Track[i].dEdxNPNoM	= Track_dEdxNPNoM[i];
      Track[i].dEdxHitsNPHarm2_1000	= Track_dEdxHitsNPHarm2_1000[i];
      Track[i].dEdxHitsNPHarm2_7	= Track_dEdxHitsNPHarm2_7[i];
      Track[i].dEdxHitsNPHarm2_5	= Track_dEdxHitsNPHarm2_5[i];
      Track[i].dEdxHitsNPHarm2_3	= Track_dEdxHitsNPHarm2_3[i];
      Track[i].dEdxHitsNPHarm2_2	= Track_dEdxHitsNPHarm2_2[i];
      Track[i].dEdxHitsNPHarm2_1	= Track_dEdxHitsNPHarm2_1[i];
      Track[i].dEdxHitsNPTrun40_1000	= Track_dEdxHitsNPTrun40_1000[i];
      Track[i].dEdxHitsNPTrun40_7	= Track_dEdxHitsNPTrun40_7[i];
      Track[i].dEdxHitsNPTrun40_5	= Track_dEdxHitsNPTrun40_5[i];
      Track[i].dEdxHitsNPTrun40_3	= Track_dEdxHitsNPTrun40_3[i];
      Track[i].dEdxHitsNPTrun40_2	= Track_dEdxHitsNPTrun40_2[i];
      Track[i].dEdxHitsNPTrun40_1	= Track_dEdxHitsNPTrun40_1[i];
      Track[i].dEdxHitsNPMedian_1000	= Track_dEdxHitsNPMedian_1000[i];
      Track[i].dEdxHitsNPMedian_7	= Track_dEdxHitsNPMedian_7[i];
      Track[i].dEdxHitsNPMedian_5	= Track_dEdxHitsNPMedian_5[i];
      Track[i].dEdxHitsNPMedian_3	= Track_dEdxHitsNPMedian_3[i];
      Track[i].dEdxHitsNPMedian_2	= Track_dEdxHitsNPMedian_2[i];
      Track[i].dEdxHitsNPMedian_1	= Track_dEdxHitsNPMedian_1[i];
      Track[i].trackerExpectedHitsOuter_numberOfLostHits	= Track_trackerExpectedHitsOuter_numberOfLostHits[i];
      Track[i].trackerExpectedHitsInner_numberOfLostHits	= Track_trackerExpectedHitsInner_numberOfLostHits[i];
    }
}

inline void fillVertex()
{
  Vertex.resize(Vertex_isValid.size());
  for(unsigned int i=0; i < Vertex.size(); ++i)
    {
      Vertex[i].isValid	= Vertex_isValid[i];
      Vertex[i].isFake	= Vertex_isFake[i];
      Vertex[i].tracksSize	= Vertex_tracksSize[i];
      Vertex[i].chi2	= Vertex_chi2[i];
      Vertex[i].ndof	= Vertex_ndof[i];
      Vertex[i].normalizedChi2	= Vertex_normalizedChi2[i];
      Vertex[i].x	= Vertex_x[i];
      Vertex[i].y	= Vertex_y[i];
      Vertex[i].z	= Vertex_z[i];
      Vertex[i].xError	= Vertex_xError[i];
      Vertex[i].yError	= Vertex_yError[i];
      Vertex[i].zError	= Vertex_zError[i];
    }
}


void fillObjects()
{
  fillGenParticle();
  fillPFJet();
  fillPFMET();
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
          GenParticle_energy[i]	= GenParticle_energy[j];
          GenParticle_et[i]	= GenParticle_et[j];
          GenParticle_p[i]	= GenParticle_p[j];
          GenParticle_pt[i]	= GenParticle_pt[j];
          GenParticle_phi[i]	= GenParticle_phi[j];
          GenParticle_eta[i]	= GenParticle_eta[j];
          GenParticle_numberOfDaughters[i]	= GenParticle_numberOfDaughters[j];
          GenParticle_numberOfMothers[i]	= GenParticle_numberOfMothers[j];
          GenParticle_mass[i]	= GenParticle_mass[j];
          GenParticle_vx[i]	= GenParticle_vx[j];
          GenParticle_vy[i]	= GenParticle_vy[j];
          GenParticle_vz[i]	= GenParticle_vz[j];
          GenParticle_pdgId[i]	= GenParticle_pdgId[j];
          GenParticle_status[i]	= GenParticle_status[j];
        }
      nGenParticle = n;
    }

  n = 0;
  try
    {
       n = indexmap["PFJet"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PFJet"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PFJet_energy[i]	= PFJet_energy[j];
          PFJet_et[i]	= PFJet_et[j];
          PFJet_p[i]	= PFJet_p[j];
          PFJet_pt[i]	= PFJet_pt[j];
          PFJet_phi[i]	= PFJet_phi[j];
          PFJet_eta[i]	= PFJet_eta[j];
          PFJet_chargedHadronEnergyFraction[i]	= PFJet_chargedHadronEnergyFraction[j];
          PFJet_neutralHadronEnergyFraction[i]	= PFJet_neutralHadronEnergyFraction[j];
          PFJet_neutralEmEnergyFraction[i]	= PFJet_neutralEmEnergyFraction[j];
          PFJet_chargedEmEnergyFraction[i]	= PFJet_chargedEmEnergyFraction[j];
        }
      nPFJet = n;
    }

  n = 0;
  try
    {
       n = indexmap["PFMET"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PFMET"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PFMET_energy[i]	= PFMET_energy[j];
          PFMET_et[i]	= PFMET_et[j];
          PFMET_p[i]	= PFMET_p[j];
          PFMET_pt[i]	= PFMET_pt[j];
          PFMET_phi[i]	= PFMET_phi[j];
          PFMET_eta[i]	= PFMET_eta[j];
        }
      nPFMET = n;
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
          Track_eta[i]	= Track_eta[j];
          Track_phi[i]	= Track_phi[j];
          Track_numberOfValidHits[i]	= Track_numberOfValidHits[j];
          Track_numberOfLostHits[i]	= Track_numberOfLostHits[j];
          Track_d0[i]	= Track_d0[j];
          Track_dz[i]	= Track_dz[j];
          Track_d0Error[i]	= Track_d0Error[j];
          Track_dzError[i]	= Track_dzError[j];
          Track_ptError[i]	= Track_ptError[j];
          Track_dEdxNPHarm2[i]	= Track_dEdxNPHarm2[j];
          Track_dEdxNPTru40[i]	= Track_dEdxNPTru40[j];
          Track_dEdxNPNoM[i]	= Track_dEdxNPNoM[j];
          Track_dEdxHitsNPHarm2_1000[i]	= Track_dEdxHitsNPHarm2_1000[j];
          Track_dEdxHitsNPHarm2_7[i]	= Track_dEdxHitsNPHarm2_7[j];
          Track_dEdxHitsNPHarm2_5[i]	= Track_dEdxHitsNPHarm2_5[j];
          Track_dEdxHitsNPHarm2_3[i]	= Track_dEdxHitsNPHarm2_3[j];
          Track_dEdxHitsNPHarm2_2[i]	= Track_dEdxHitsNPHarm2_2[j];
          Track_dEdxHitsNPHarm2_1[i]	= Track_dEdxHitsNPHarm2_1[j];
          Track_dEdxHitsNPTrun40_1000[i]	= Track_dEdxHitsNPTrun40_1000[j];
          Track_dEdxHitsNPTrun40_7[i]	= Track_dEdxHitsNPTrun40_7[j];
          Track_dEdxHitsNPTrun40_5[i]	= Track_dEdxHitsNPTrun40_5[j];
          Track_dEdxHitsNPTrun40_3[i]	= Track_dEdxHitsNPTrun40_3[j];
          Track_dEdxHitsNPTrun40_2[i]	= Track_dEdxHitsNPTrun40_2[j];
          Track_dEdxHitsNPTrun40_1[i]	= Track_dEdxHitsNPTrun40_1[j];
          Track_dEdxHitsNPMedian_1000[i]	= Track_dEdxHitsNPMedian_1000[j];
          Track_dEdxHitsNPMedian_7[i]	= Track_dEdxHitsNPMedian_7[j];
          Track_dEdxHitsNPMedian_5[i]	= Track_dEdxHitsNPMedian_5[j];
          Track_dEdxHitsNPMedian_3[i]	= Track_dEdxHitsNPMedian_3[j];
          Track_dEdxHitsNPMedian_2[i]	= Track_dEdxHitsNPMedian_2[j];
          Track_dEdxHitsNPMedian_1[i]	= Track_dEdxHitsNPMedian_1[j];
          Track_trackerExpectedHitsOuter_numberOfLostHits[i]	= Track_trackerExpectedHitsOuter_numberOfLostHits[j];
          Track_trackerExpectedHitsInner_numberOfLostHits[i]	= Track_trackerExpectedHitsInner_numberOfLostHits[j];
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
          Vertex_isValid[i]	= Vertex_isValid[j];
          Vertex_isFake[i]	= Vertex_isFake[j];
          Vertex_tracksSize[i]	= Vertex_tracksSize[j];
          Vertex_chi2[i]	= Vertex_chi2[j];
          Vertex_ndof[i]	= Vertex_ndof[j];
          Vertex_normalizedChi2[i]	= Vertex_normalizedChi2[j];
          Vertex_x[i]	= Vertex_x[j];
          Vertex_y[i]	= Vertex_y[j];
          Vertex_z[i]	= Vertex_z[j];
          Vertex_xError[i]	= Vertex_xError[j];
          Vertex_yError[i]	= Vertex_yError[j];
          Vertex_zError[i]	= Vertex_zError[j];
        }
      nVertex = n;
    }
}

//-----------------------------------------------------------------------------
// --- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("recoGenParticle_genParticles.charge", GenParticle_charge);
  stream.select("recoGenParticle_genParticles.energy", GenParticle_energy);
  stream.select("recoGenParticle_genParticles.et", GenParticle_et);
  stream.select("recoGenParticle_genParticles.eta", GenParticle_eta);
  stream.select("recoGenParticle_genParticles.mass", GenParticle_mass);
  stream.select("recoGenParticle_genParticles.numberOfDaughters", GenParticle_numberOfDaughters);
  stream.select("recoGenParticle_genParticles.numberOfMothers", GenParticle_numberOfMothers);
  stream.select("recoGenParticle_genParticles.p", GenParticle_p);
  stream.select("recoGenParticle_genParticles.pdgId", GenParticle_pdgId);
  stream.select("recoGenParticle_genParticles.phi", GenParticle_phi);
  stream.select("recoGenParticle_genParticles.pt", GenParticle_pt);
  stream.select("recoGenParticle_genParticles.status", GenParticle_status);
  stream.select("recoGenParticle_genParticles.vx", GenParticle_vx);
  stream.select("recoGenParticle_genParticles.vy", GenParticle_vy);
  stream.select("recoGenParticle_genParticles.vz", GenParticle_vz);
  stream.select("recoPFJet_ak5PFJetsPt15.chargedEmEnergyFraction", PFJet_chargedEmEnergyFraction);
  stream.select("recoPFJet_ak5PFJetsPt15.chargedHadronEnergyFraction", PFJet_chargedHadronEnergyFraction);
  stream.select("recoPFJet_ak5PFJetsPt15.energy", PFJet_energy);
  stream.select("recoPFJet_ak5PFJetsPt15.et", PFJet_et);
  stream.select("recoPFJet_ak5PFJetsPt15.eta", PFJet_eta);
  stream.select("recoPFJet_ak5PFJetsPt15.neutralEmEnergyFraction", PFJet_neutralEmEnergyFraction);
  stream.select("recoPFJet_ak5PFJetsPt15.neutralHadronEnergyFraction", PFJet_neutralHadronEnergyFraction);
  stream.select("recoPFJet_ak5PFJetsPt15.p", PFJet_p);
  stream.select("recoPFJet_ak5PFJetsPt15.phi", PFJet_phi);
  stream.select("recoPFJet_ak5PFJetsPt15.pt", PFJet_pt);
  stream.select("recoPFMET_pfMet.energy", PFMET_energy);
  stream.select("recoPFMET_pfMet.et", PFMET_et);
  stream.select("recoPFMET_pfMet.eta", PFMET_eta);
  stream.select("recoPFMET_pfMet.p", PFMET_p);
  stream.select("recoPFMET_pfMet.phi", PFMET_phi);
  stream.select("recoPFMET_pfMet.pt", PFMET_pt);
  stream.select("recoTrackHelper_TrackRefitter.d0", Track_d0);
  stream.select("recoTrackHelper_TrackRefitter.d0Error", Track_d0Error);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_1", Track_dEdxHitsNPHarm2_1);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_1000", Track_dEdxHitsNPHarm2_1000);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_2", Track_dEdxHitsNPHarm2_2);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_3", Track_dEdxHitsNPHarm2_3);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_5", Track_dEdxHitsNPHarm2_5);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_7", Track_dEdxHitsNPHarm2_7);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPMedian_1", Track_dEdxHitsNPMedian_1);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPMedian_1000", Track_dEdxHitsNPMedian_1000);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPMedian_2", Track_dEdxHitsNPMedian_2);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPMedian_3", Track_dEdxHitsNPMedian_3);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPMedian_5", Track_dEdxHitsNPMedian_5);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPMedian_7", Track_dEdxHitsNPMedian_7);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPTrun40_1", Track_dEdxHitsNPTrun40_1);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPTrun40_100000", Track_dEdxHitsNPTrun40_1000);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPTrun40_2", Track_dEdxHitsNPTrun40_2);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPTrun40_3", Track_dEdxHitsNPTrun40_3);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPTrun40_5", Track_dEdxHitsNPTrun40_5);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPTrun40_7", Track_dEdxHitsNPTrun40_7);
  stream.select("recoTrackHelper_TrackRefitter.dEdxNPHarm2", Track_dEdxNPHarm2);
  stream.select("recoTrackHelper_TrackRefitter.dEdxNPNoM", Track_dEdxNPNoM);
  stream.select("recoTrackHelper_TrackRefitter.dEdxNPTru40", Track_dEdxNPTru40);
  stream.select("recoTrackHelper_TrackRefitter.dz", Track_dz);
  stream.select("recoTrackHelper_TrackRefitter.dzError", Track_dzError);
  stream.select("recoTrackHelper_TrackRefitter.eta", Track_eta);
  stream.select("recoTrackHelper_TrackRefitter.numberOfLostHits", Track_numberOfLostHits);
  stream.select("recoTrackHelper_TrackRefitter.numberOfValidHits", Track_numberOfValidHits);
  stream.select("recoTrackHelper_TrackRefitter.phi", Track_phi);
  stream.select("recoTrackHelper_TrackRefitter.pt", Track_pt);
  stream.select("recoTrackHelper_TrackRefitter.ptError", Track_ptError);
  stream.select("recoTrackHelper_TrackRefitter.trackerExpectedHitsInner_numberOfLostHits", Track_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("recoTrackHelper_TrackRefitter.trackerExpectedHitsOuter_numberOfLostHits", Track_trackerExpectedHitsOuter_numberOfLostHits);
  stream.select("recoVertex_offlinePrimaryVertices.chi2", Vertex_chi2);
  stream.select("recoVertex_offlinePrimaryVertices.isFake", Vertex_isFake);
  stream.select("recoVertex_offlinePrimaryVertices.isValid", Vertex_isValid);
  stream.select("recoVertex_offlinePrimaryVertices.ndof", Vertex_ndof);
  stream.select("recoVertex_offlinePrimaryVertices.normalizedChi2", Vertex_normalizedChi2);
  stream.select("recoVertex_offlinePrimaryVertices.tracksSize", Vertex_tracksSize);
  stream.select("recoVertex_offlinePrimaryVertices.x", Vertex_x);
  stream.select("recoVertex_offlinePrimaryVertices.xError", Vertex_xError);
  stream.select("recoVertex_offlinePrimaryVertices.y", Vertex_y);
  stream.select("recoVertex_offlinePrimaryVertices.yError", Vertex_yError);
  stream.select("recoVertex_offlinePrimaryVertices.z", Vertex_z);
  stream.select("recoVertex_offlinePrimaryVertices.zError", Vertex_zError);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5", edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5", edmTriggerResultsHelper_prescaleHLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5);
  stream.select("nrecoGenParticle_genParticles", nGenParticle);
  stream.select("nrecoPFJet_ak5PFJetsPt15", nPFJet);
  stream.select("nrecoPFMET_pfMet", nPFMET);
  stream.select("nrecoTrackHelper_TrackRefitter", nTrack);
  stream.select("nrecoVertex_offlinePrimaryVertices", nVertex);

}
}; // end namespace evt
#endif
