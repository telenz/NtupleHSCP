#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Thu Mar 20 12:10:35 2014 by mkanalyzer.py
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
std::vector<double>	CaloCluster1_energy(200,0);
std::vector<double>	CaloCluster1_eta(200,0);
std::vector<double>	CaloCluster1_phi(200,0);
std::vector<double>	CaloCluster1_x(200,0);
std::vector<double>	CaloCluster1_y(200,0);
std::vector<double>	CaloCluster1_z(200,0);
std::vector<double>	CaloCluster2_energy(200,0);
std::vector<double>	CaloCluster2_eta(200,0);
std::vector<double>	CaloCluster2_phi(200,0);
std::vector<double>	CaloCluster2_x(200,0);
std::vector<double>	CaloCluster2_y(200,0);
std::vector<double>	CaloCluster2_z(200,0);
std::vector<double>	CaloCluster3_energy(200,0);
std::vector<double>	CaloCluster3_eta(200,0);
std::vector<double>	CaloCluster3_phi(200,0);
std::vector<double>	CaloCluster3_x(200,0);
std::vector<double>	CaloCluster3_y(200,0);
std::vector<double>	CaloCluster3_z(200,0);
std::vector<double>	CaloCluster_energy(200,0);
std::vector<double>	CaloCluster_eta(200,0);
std::vector<double>	CaloCluster_phi(200,0);
std::vector<double>	CaloCluster_x(200,0);
std::vector<double>	CaloCluster_y(200,0);
std::vector<double>	CaloCluster_z(200,0);
std::vector<int>	GenParticle_charge(200,0);
std::vector<double>	GenParticle_energy(200,0);
std::vector<double>	GenParticle_et(200,0);
std::vector<double>	GenParticle_eta(200,0);
std::vector<double>	GenParticle_mass(200,0);
std::vector<size_t>	GenParticle_numberOfDaughters(200,0);
std::vector<size_t>	GenParticle_numberOfMothers(200,0);
std::vector<double>	GenParticle_p(200,0);
std::vector<int>	GenParticle_pdgId(200,0);
std::vector<double>	GenParticle_phi(200,0);
std::vector<double>	GenParticle_pt(200,0);
std::vector<int>	GenParticle_status(200,0);
std::vector<double>	GenParticle_vx(200,0);
std::vector<double>	GenParticle_vy(200,0);
std::vector<double>	GenParticle_vz(200,0);
std::vector<double>	GsfElectron_energy(200,0);
std::vector<double>	GsfElectron_et(200,0);
std::vector<double>	GsfElectron_eta(200,0);
std::vector<double>	GsfElectron_p(200,0);
std::vector<double>	GsfElectron_phi(200,0);
std::vector<double>	GsfElectron_pt(200,0);
std::vector<double>	Muon_energy(200,0);
std::vector<double>	Muon_et(200,0);
std::vector<double>	Muon_eta(200,0);
std::vector<double>	Muon_p(200,0);
std::vector<double>	Muon_phi(200,0);
std::vector<double>	Muon_pt(200,0);
std::vector<double>	PFCandidate_energy(200,0);
std::vector<double>	PFCandidate_et(200,0);
std::vector<double>	PFCandidate_eta(200,0);
std::vector<double>	PFCandidate_p(200,0);
std::vector<int>	PFCandidate_pdgId(200,0);
std::vector<double>	PFCandidate_phi(200,0);
std::vector<double>	PFCandidate_pt(200,0);
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
std::vector<double>	PFTau1_energy(200,0);
std::vector<double>	PFTau1_et(200,0);
std::vector<double>	PFTau1_eta(200,0);
std::vector<double>	PFTau1_p(200,0);
std::vector<double>	PFTau1_phi(200,0);
std::vector<double>	PFTau1_pt(200,0);
std::vector<double>	PFTau2_energy(200,0);
std::vector<double>	PFTau2_et(200,0);
std::vector<double>	PFTau2_eta(200,0);
std::vector<double>	PFTau2_p(200,0);
std::vector<double>	PFTau2_phi(200,0);
std::vector<double>	PFTau2_pt(200,0);
std::vector<double>	PFTau_energy(200,0);
std::vector<double>	PFTau_et(200,0);
std::vector<double>	PFTau_eta(200,0);
std::vector<double>	PFTau_p(200,0);
std::vector<double>	PFTau_phi(200,0);
std::vector<double>	PFTau_pt(200,0);
std::vector<int>	Track_charge(200,0);
std::vector<double>	Track_d0(200,0);
std::vector<double>	Track_d0Error(200,0);
std::vector<double>	Track_dz(200,0);
std::vector<double>	Track_dzError(200,0);
std::vector<double>	Track_eta(200,0);
std::vector<unsigned short>	Track_numberOfLostHits(200,0);
std::vector<unsigned short>	Track_numberOfValidHits(200,0);
std::vector<double>	Track_p(200,0);
std::vector<double>	Track_phi(200,0);
std::vector<double>	Track_pt(200,0);
std::vector<double>	Track_ptError(200,0);
std::vector<int>	Track_qualityMask(200,0);
std::vector<unsigned short>	Track_trackerExpectedHitsInner_numberOfLostHits(200,0);
std::vector<unsigned short>	Track_trackerExpectedHitsOuter_numberOfLostHits(200,0);
std::vector<double>	Track_vx(200,0);
std::vector<double>	Track_vy(200,0);
std::vector<double>	Track_vz(200,0);
std::vector<double>	Vertex_chi2(200,0);
std::vector<int>	Vertex_hasRefittedTracks(200,0);
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
int	nCaloCluster;
int	nCaloCluster1;
int	nCaloCluster2;
int	nCaloCluster3;
int	nGenParticle;
int	nGsfElectron;
int	nMuon;
int	nPFCandidate;
int	nPFJet;
int	nPFMET;
int	nPFTau;
int	nPFTau1;
int	nPFTau2;
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
struct CaloCluster_s
{
  double	energy;
  double	eta;
  double	phi;
  double	x;
  double	y;
  double	z;
};
std::vector<CaloCluster_s> CaloCluster(200);

std::ostream& operator<<(std::ostream& os, const CaloCluster_s& o)
{
  char r[1024];
  os << "CaloCluster" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct CaloCluster1_s
{
  double	energy;
  double	eta;
  double	phi;
  double	x;
  double	y;
  double	z;
};
std::vector<CaloCluster1_s> CaloCluster1(200);

std::ostream& operator<<(std::ostream& os, const CaloCluster1_s& o)
{
  char r[1024];
  os << "CaloCluster1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct CaloCluster2_s
{
  double	energy;
  double	eta;
  double	phi;
  double	x;
  double	y;
  double	z;
};
std::vector<CaloCluster2_s> CaloCluster2(200);

std::ostream& operator<<(std::ostream& os, const CaloCluster2_s& o)
{
  char r[1024];
  os << "CaloCluster2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct CaloCluster3_s
{
  double	energy;
  double	eta;
  double	phi;
  double	x;
  double	y;
  double	z;
};
std::vector<CaloCluster3_s> CaloCluster3(200);

std::ostream& operator<<(std::ostream& os, const CaloCluster3_s& o)
{
  char r[1024];
  os << "CaloCluster3" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct GenParticle_s
{
  int	charge;
  double	p;
  double	energy;
  double	et;
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
std::vector<GenParticle_s> GenParticle(200);

std::ostream& operator<<(std::ostream& os, const GenParticle_s& o)
{
  char r[1024];
  os << "GenParticle" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
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
struct GsfElectron_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<GsfElectron_s> GsfElectron(200);

std::ostream& operator<<(std::ostream& os, const GsfElectron_s& o)
{
  char r[1024];
  os << "GsfElectron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Muon_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<Muon_s> Muon(200);

std::ostream& operator<<(std::ostream& os, const Muon_s& o)
{
  char r[1024];
  os << "Muon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFCandidate_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  int	pdgId;
};
std::vector<PFCandidate_s> PFCandidate(200);

std::ostream& operator<<(std::ostream& os, const PFCandidate_s& o)
{
  char r[1024];
  os << "PFCandidate" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFJet_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
};
std::vector<PFJet_s> PFJet(200);

std::ostream& operator<<(std::ostream& os, const PFJet_s& o)
{
  char r[1024];
  os << "PFJet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFMET_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<PFMET_s> PFMET(200);

std::ostream& operator<<(std::ostream& os, const PFMET_s& o)
{
  char r[1024];
  os << "PFMET" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFTau_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<PFTau_s> PFTau(200);

std::ostream& operator<<(std::ostream& os, const PFTau_s& o)
{
  char r[1024];
  os << "PFTau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFTau1_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<PFTau1_s> PFTau1(200);

std::ostream& operator<<(std::ostream& os, const PFTau1_s& o)
{
  char r[1024];
  os << "PFTau1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct PFTau2_s
{
  double	p;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<PFTau2_s> PFTau2(200);

std::ostream& operator<<(std::ostream& os, const PFTau2_s& o)
{
  char r[1024];
  os << "PFTau2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct Track_s
{
  int	charge;
  double	p;
  double	pt;
  double	phi;
  double	eta;
  double	d0;
  double	dz;
  double	vx;
  double	vy;
  double	vz;
  double	ptError;
  double	d0Error;
  double	dzError;
  unsigned short	numberOfValidHits;
  unsigned short	numberOfLostHits;
  int	qualityMask;
  unsigned short	trackerExpectedHitsOuter_numberOfLostHits;
  unsigned short	trackerExpectedHitsInner_numberOfLostHits;
  int pdgId;
  double beta;
};
std::vector<Track_s> Track(200);

std::ostream& operator<<(std::ostream& os, const Track_s& o)
{
  char r[1024];
  os << "Track" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "d0", (double)o.d0); os << r;
  sprintf(r, "  %-32s: %f\n", "dz", (double)o.dz); os << r;
  sprintf(r, "  %-32s: %f\n", "vx", (double)o.vx); os << r;
  sprintf(r, "  %-32s: %f\n", "vy", (double)o.vy); os << r;
  sprintf(r, "  %-32s: %f\n", "vz", (double)o.vz); os << r;
  sprintf(r, "  %-32s: %f\n", "ptError", (double)o.ptError); os << r;
  sprintf(r, "  %-32s: %f\n", "d0Error", (double)o.d0Error); os << r;
  sprintf(r, "  %-32s: %f\n", "dzError", (double)o.dzError); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfValidHits", (double)o.numberOfValidHits); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfLostHits", (double)o.numberOfLostHits); os << r;
  sprintf(r, "  %-32s: %f\n", "qualityMask", (double)o.qualityMask); os << r;
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
  int	hasRefittedTracks;
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
  sprintf(r, "  %-32s: %f\n", "hasRefittedTracks", (double)o.hasRefittedTracks); os << r;
  return os;
}
//-----------------------------------------------------------------------------

inline void fillCaloCluster()
{
  CaloCluster.resize(CaloCluster_energy.size());
  for(unsigned int i=0; i < CaloCluster.size(); ++i)
    {
      CaloCluster[i].energy	= CaloCluster_energy[i];
      CaloCluster[i].eta	= CaloCluster_eta[i];
      CaloCluster[i].phi	= CaloCluster_phi[i];
      CaloCluster[i].x	= CaloCluster_x[i];
      CaloCluster[i].y	= CaloCluster_y[i];
      CaloCluster[i].z	= CaloCluster_z[i];
    }
}

inline void fillCaloCluster1()
{
  CaloCluster1.resize(CaloCluster1_energy.size());
  for(unsigned int i=0; i < CaloCluster1.size(); ++i)
    {
      CaloCluster1[i].energy	= CaloCluster1_energy[i];
      CaloCluster1[i].eta	= CaloCluster1_eta[i];
      CaloCluster1[i].phi	= CaloCluster1_phi[i];
      CaloCluster1[i].x	= CaloCluster1_x[i];
      CaloCluster1[i].y	= CaloCluster1_y[i];
      CaloCluster1[i].z	= CaloCluster1_z[i];
    }
}

inline void fillCaloCluster2()
{
  CaloCluster2.resize(CaloCluster2_energy.size());
  for(unsigned int i=0; i < CaloCluster2.size(); ++i)
    {
      CaloCluster2[i].energy	= CaloCluster2_energy[i];
      CaloCluster2[i].eta	= CaloCluster2_eta[i];
      CaloCluster2[i].phi	= CaloCluster2_phi[i];
      CaloCluster2[i].x	= CaloCluster2_x[i];
      CaloCluster2[i].y	= CaloCluster2_y[i];
      CaloCluster2[i].z	= CaloCluster2_z[i];
    }
}

inline void fillCaloCluster3()
{
  CaloCluster3.resize(CaloCluster3_energy.size());
  for(unsigned int i=0; i < CaloCluster3.size(); ++i)
    {
      CaloCluster3[i].energy	= CaloCluster3_energy[i];
      CaloCluster3[i].eta	= CaloCluster3_eta[i];
      CaloCluster3[i].phi	= CaloCluster3_phi[i];
      CaloCluster3[i].x	= CaloCluster3_x[i];
      CaloCluster3[i].y	= CaloCluster3_y[i];
      CaloCluster3[i].z	= CaloCluster3_z[i];
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

inline void fillGsfElectron()
{
  GsfElectron.resize(GsfElectron_p.size());
  for(unsigned int i=0; i < GsfElectron.size(); ++i)
    {
      GsfElectron[i].p	= GsfElectron_p[i];
      GsfElectron[i].energy	= GsfElectron_energy[i];
      GsfElectron[i].et	= GsfElectron_et[i];
      GsfElectron[i].pt	= GsfElectron_pt[i];
      GsfElectron[i].phi	= GsfElectron_phi[i];
      GsfElectron[i].eta	= GsfElectron_eta[i];
    }
}

inline void fillMuon()
{
  Muon.resize(Muon_p.size());
  for(unsigned int i=0; i < Muon.size(); ++i)
    {
      Muon[i].p	= Muon_p[i];
      Muon[i].energy	= Muon_energy[i];
      Muon[i].et	= Muon_et[i];
      Muon[i].pt	= Muon_pt[i];
      Muon[i].phi	= Muon_phi[i];
      Muon[i].eta	= Muon_eta[i];
    }
}

inline void fillPFCandidate()
{
  PFCandidate.resize(PFCandidate_p.size());
  for(unsigned int i=0; i < PFCandidate.size(); ++i)
    {
      PFCandidate[i].p	= PFCandidate_p[i];
      PFCandidate[i].energy	= PFCandidate_energy[i];
      PFCandidate[i].et	= PFCandidate_et[i];
      PFCandidate[i].pt	= PFCandidate_pt[i];
      PFCandidate[i].phi	= PFCandidate_phi[i];
      PFCandidate[i].eta	= PFCandidate_eta[i];
      PFCandidate[i].pdgId	= PFCandidate_pdgId[i];
    }
}

inline void fillPFJet()
{
  PFJet.resize(PFJet_p.size());
  for(unsigned int i=0; i < PFJet.size(); ++i)
    {
      PFJet[i].p	= PFJet_p[i];
      PFJet[i].energy	= PFJet_energy[i];
      PFJet[i].et	= PFJet_et[i];
      PFJet[i].pt	= PFJet_pt[i];
      PFJet[i].phi	= PFJet_phi[i];
      PFJet[i].eta	= PFJet_eta[i];
      PFJet[i].chargedHadronEnergyFraction	= PFJet_chargedHadronEnergyFraction[i];
      PFJet[i].neutralHadronEnergyFraction	= PFJet_neutralHadronEnergyFraction[i];
      PFJet[i].chargedEmEnergyFraction	= PFJet_chargedEmEnergyFraction[i];
      PFJet[i].neutralEmEnergyFraction	= PFJet_neutralEmEnergyFraction[i];
    }
}

inline void fillPFMET()
{
  PFMET.resize(PFMET_p.size());
  for(unsigned int i=0; i < PFMET.size(); ++i)
    {
      PFMET[i].p	= PFMET_p[i];
      PFMET[i].energy	= PFMET_energy[i];
      PFMET[i].et	= PFMET_et[i];
      PFMET[i].pt	= PFMET_pt[i];
      PFMET[i].phi	= PFMET_phi[i];
      PFMET[i].eta	= PFMET_eta[i];
    }
}

inline void fillPFTau()
{
  PFTau.resize(PFTau_p.size());
  for(unsigned int i=0; i < PFTau.size(); ++i)
    {
      PFTau[i].p	= PFTau_p[i];
      PFTau[i].energy	= PFTau_energy[i];
      PFTau[i].et	= PFTau_et[i];
      PFTau[i].pt	= PFTau_pt[i];
      PFTau[i].phi	= PFTau_phi[i];
      PFTau[i].eta	= PFTau_eta[i];
    }
}

inline void fillPFTau1()
{
  PFTau1.resize(PFTau1_p.size());
  for(unsigned int i=0; i < PFTau1.size(); ++i)
    {
      PFTau1[i].p	= PFTau1_p[i];
      PFTau1[i].energy	= PFTau1_energy[i];
      PFTau1[i].et	= PFTau1_et[i];
      PFTau1[i].pt	= PFTau1_pt[i];
      PFTau1[i].phi	= PFTau1_phi[i];
      PFTau1[i].eta	= PFTau1_eta[i];
    }
}

inline void fillPFTau2()
{
  PFTau2.resize(PFTau2_p.size());
  for(unsigned int i=0; i < PFTau2.size(); ++i)
    {
      PFTau2[i].p	= PFTau2_p[i];
      PFTau2[i].energy	= PFTau2_energy[i];
      PFTau2[i].et	= PFTau2_et[i];
      PFTau2[i].pt	= PFTau2_pt[i];
      PFTau2[i].phi	= PFTau2_phi[i];
      PFTau2[i].eta	= PFTau2_eta[i];
    }
}

inline void fillTrack()
{
  Track.resize(Track_charge.size());
  for(unsigned int i=0; i < Track.size(); ++i)
    {
      Track[i].charge	= Track_charge[i];
      Track[i].p	= Track_p[i];
      Track[i].pt	= Track_pt[i];
      Track[i].phi	= Track_phi[i];
      Track[i].eta	= Track_eta[i];
      Track[i].d0	= Track_d0[i];
      Track[i].dz	= Track_dz[i];
      Track[i].vx	= Track_vx[i];
      Track[i].vy	= Track_vy[i];
      Track[i].vz	= Track_vz[i];
      Track[i].ptError	= Track_ptError[i];
      Track[i].d0Error	= Track_d0Error[i];
      Track[i].dzError	= Track_dzError[i];
      Track[i].numberOfValidHits	= Track_numberOfValidHits[i];
      Track[i].numberOfLostHits	= Track_numberOfLostHits[i];
      Track[i].qualityMask	= Track_qualityMask[i];
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
      Vertex[i].hasRefittedTracks	= Vertex_hasRefittedTracks[i];
    }
}


void fillObjects()
{
  fillCaloCluster();
  fillCaloCluster1();
  fillCaloCluster2();
  fillCaloCluster3();
  fillGenParticle();
  fillGsfElectron();
  fillMuon();
  fillPFCandidate();
  fillPFJet();
  fillPFMET();
  fillPFTau();
  fillPFTau1();
  fillPFTau2();
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
       n = indexmap["CaloCluster"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["CaloCluster"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          CaloCluster_energy[i]	= CaloCluster_energy[j];
          CaloCluster_eta[i]	= CaloCluster_eta[j];
          CaloCluster_phi[i]	= CaloCluster_phi[j];
          CaloCluster_x[i]	= CaloCluster_x[j];
          CaloCluster_y[i]	= CaloCluster_y[j];
          CaloCluster_z[i]	= CaloCluster_z[j];
        }
      nCaloCluster = n;
    }

  n = 0;
  try
    {
       n = indexmap["CaloCluster1"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["CaloCluster1"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          CaloCluster1_energy[i]	= CaloCluster1_energy[j];
          CaloCluster1_eta[i]	= CaloCluster1_eta[j];
          CaloCluster1_phi[i]	= CaloCluster1_phi[j];
          CaloCluster1_x[i]	= CaloCluster1_x[j];
          CaloCluster1_y[i]	= CaloCluster1_y[j];
          CaloCluster1_z[i]	= CaloCluster1_z[j];
        }
      nCaloCluster1 = n;
    }

  n = 0;
  try
    {
       n = indexmap["CaloCluster2"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["CaloCluster2"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          CaloCluster2_energy[i]	= CaloCluster2_energy[j];
          CaloCluster2_eta[i]	= CaloCluster2_eta[j];
          CaloCluster2_phi[i]	= CaloCluster2_phi[j];
          CaloCluster2_x[i]	= CaloCluster2_x[j];
          CaloCluster2_y[i]	= CaloCluster2_y[j];
          CaloCluster2_z[i]	= CaloCluster2_z[j];
        }
      nCaloCluster2 = n;
    }

  n = 0;
  try
    {
       n = indexmap["CaloCluster3"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["CaloCluster3"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          CaloCluster3_energy[i]	= CaloCluster3_energy[j];
          CaloCluster3_eta[i]	= CaloCluster3_eta[j];
          CaloCluster3_phi[i]	= CaloCluster3_phi[j];
          CaloCluster3_x[i]	= CaloCluster3_x[j];
          CaloCluster3_y[i]	= CaloCluster3_y[j];
          CaloCluster3_z[i]	= CaloCluster3_z[j];
        }
      nCaloCluster3 = n;
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
       n = indexmap["GsfElectron"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["GsfElectron"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          GsfElectron_p[i]	= GsfElectron_p[j];
          GsfElectron_energy[i]	= GsfElectron_energy[j];
          GsfElectron_et[i]	= GsfElectron_et[j];
          GsfElectron_pt[i]	= GsfElectron_pt[j];
          GsfElectron_phi[i]	= GsfElectron_phi[j];
          GsfElectron_eta[i]	= GsfElectron_eta[j];
        }
      nGsfElectron = n;
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
          Muon_p[i]	= Muon_p[j];
          Muon_energy[i]	= Muon_energy[j];
          Muon_et[i]	= Muon_et[j];
          Muon_pt[i]	= Muon_pt[j];
          Muon_phi[i]	= Muon_phi[j];
          Muon_eta[i]	= Muon_eta[j];
        }
      nMuon = n;
    }

  n = 0;
  try
    {
       n = indexmap["PFCandidate"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PFCandidate"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PFCandidate_p[i]	= PFCandidate_p[j];
          PFCandidate_energy[i]	= PFCandidate_energy[j];
          PFCandidate_et[i]	= PFCandidate_et[j];
          PFCandidate_pt[i]	= PFCandidate_pt[j];
          PFCandidate_phi[i]	= PFCandidate_phi[j];
          PFCandidate_eta[i]	= PFCandidate_eta[j];
          PFCandidate_pdgId[i]	= PFCandidate_pdgId[j];
        }
      nPFCandidate = n;
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
          PFJet_p[i]	= PFJet_p[j];
          PFJet_energy[i]	= PFJet_energy[j];
          PFJet_et[i]	= PFJet_et[j];
          PFJet_pt[i]	= PFJet_pt[j];
          PFJet_phi[i]	= PFJet_phi[j];
          PFJet_eta[i]	= PFJet_eta[j];
          PFJet_chargedHadronEnergyFraction[i]	= PFJet_chargedHadronEnergyFraction[j];
          PFJet_neutralHadronEnergyFraction[i]	= PFJet_neutralHadronEnergyFraction[j];
          PFJet_chargedEmEnergyFraction[i]	= PFJet_chargedEmEnergyFraction[j];
          PFJet_neutralEmEnergyFraction[i]	= PFJet_neutralEmEnergyFraction[j];
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
          PFMET_p[i]	= PFMET_p[j];
          PFMET_energy[i]	= PFMET_energy[j];
          PFMET_et[i]	= PFMET_et[j];
          PFMET_pt[i]	= PFMET_pt[j];
          PFMET_phi[i]	= PFMET_phi[j];
          PFMET_eta[i]	= PFMET_eta[j];
        }
      nPFMET = n;
    }

  n = 0;
  try
    {
       n = indexmap["PFTau"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PFTau"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PFTau_p[i]	= PFTau_p[j];
          PFTau_energy[i]	= PFTau_energy[j];
          PFTau_et[i]	= PFTau_et[j];
          PFTau_pt[i]	= PFTau_pt[j];
          PFTau_phi[i]	= PFTau_phi[j];
          PFTau_eta[i]	= PFTau_eta[j];
        }
      nPFTau = n;
    }

  n = 0;
  try
    {
       n = indexmap["PFTau1"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PFTau1"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PFTau1_p[i]	= PFTau1_p[j];
          PFTau1_energy[i]	= PFTau1_energy[j];
          PFTau1_et[i]	= PFTau1_et[j];
          PFTau1_pt[i]	= PFTau1_pt[j];
          PFTau1_phi[i]	= PFTau1_phi[j];
          PFTau1_eta[i]	= PFTau1_eta[j];
        }
      nPFTau1 = n;
    }

  n = 0;
  try
    {
       n = indexmap["PFTau2"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["PFTau2"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          PFTau2_p[i]	= PFTau2_p[j];
          PFTau2_energy[i]	= PFTau2_energy[j];
          PFTau2_et[i]	= PFTau2_et[j];
          PFTau2_pt[i]	= PFTau2_pt[j];
          PFTau2_phi[i]	= PFTau2_phi[j];
          PFTau2_eta[i]	= PFTau2_eta[j];
        }
      nPFTau2 = n;
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
          Track_charge[i]	= Track_charge[j];
          Track_p[i]	= Track_p[j];
          Track_pt[i]	= Track_pt[j];
          Track_phi[i]	= Track_phi[j];
          Track_eta[i]	= Track_eta[j];
          Track_d0[i]	= Track_d0[j];
          Track_dz[i]	= Track_dz[j];
          Track_vx[i]	= Track_vx[j];
          Track_vy[i]	= Track_vy[j];
          Track_vz[i]	= Track_vz[j];
          Track_ptError[i]	= Track_ptError[j];
          Track_d0Error[i]	= Track_d0Error[j];
          Track_dzError[i]	= Track_dzError[j];
          Track_numberOfValidHits[i]	= Track_numberOfValidHits[j];
          Track_numberOfLostHits[i]	= Track_numberOfLostHits[j];
          Track_qualityMask[i]	= Track_qualityMask[j];
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
          Vertex_hasRefittedTracks[i]	= Vertex_hasRefittedTracks[j];
        }
      nVertex = n;
    }
}

//-----------------------------------------------------------------------------
// --- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("recoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters.energy", CaloCluster1_energy);
  stream.select("recoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters.eta", CaloCluster1_eta);
  stream.select("recoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters.phi", CaloCluster1_phi);
  stream.select("recoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters.x", CaloCluster1_x);
  stream.select("recoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters.y", CaloCluster1_y);
  stream.select("recoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters.z", CaloCluster1_z);
  stream.select("recoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters.energy", CaloCluster2_energy);
  stream.select("recoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters.eta", CaloCluster2_eta);
  stream.select("recoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters.phi", CaloCluster2_phi);
  stream.select("recoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters.x", CaloCluster2_x);
  stream.select("recoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters.y", CaloCluster2_y);
  stream.select("recoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters.z", CaloCluster2_z);
  stream.select("recoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters.energy", CaloCluster3_energy);
  stream.select("recoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters.eta", CaloCluster3_eta);
  stream.select("recoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters.phi", CaloCluster3_phi);
  stream.select("recoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters.x", CaloCluster3_x);
  stream.select("recoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters.y", CaloCluster3_y);
  stream.select("recoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters.z", CaloCluster3_z);
  stream.select("recoCaloCluster_hfEMClusters.energy", CaloCluster_energy);
  stream.select("recoCaloCluster_hfEMClusters.eta", CaloCluster_eta);
  stream.select("recoCaloCluster_hfEMClusters.phi", CaloCluster_phi);
  stream.select("recoCaloCluster_hfEMClusters.x", CaloCluster_x);
  stream.select("recoCaloCluster_hfEMClusters.y", CaloCluster_y);
  stream.select("recoCaloCluster_hfEMClusters.z", CaloCluster_z);
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
  stream.select("recoGsfElectron_gsfElectrons.energy", GsfElectron_energy);
  stream.select("recoGsfElectron_gsfElectrons.et", GsfElectron_et);
  stream.select("recoGsfElectron_gsfElectrons.eta", GsfElectron_eta);
  stream.select("recoGsfElectron_gsfElectrons.p", GsfElectron_p);
  stream.select("recoGsfElectron_gsfElectrons.phi", GsfElectron_phi);
  stream.select("recoGsfElectron_gsfElectrons.pt", GsfElectron_pt);
  stream.select("recoMuon_muons.energy", Muon_energy);
  stream.select("recoMuon_muons.et", Muon_et);
  stream.select("recoMuon_muons.eta", Muon_eta);
  stream.select("recoMuon_muons.p", Muon_p);
  stream.select("recoMuon_muons.phi", Muon_phi);
  stream.select("recoMuon_muons.pt", Muon_pt);
  stream.select("recoPFCandidate_particleFlow.energy", PFCandidate_energy);
  stream.select("recoPFCandidate_particleFlow.et", PFCandidate_et);
  stream.select("recoPFCandidate_particleFlow.eta", PFCandidate_eta);
  stream.select("recoPFCandidate_particleFlow.p", PFCandidate_p);
  stream.select("recoPFCandidate_particleFlow.pdgId", PFCandidate_pdgId);
  stream.select("recoPFCandidate_particleFlow.phi", PFCandidate_phi);
  stream.select("recoPFCandidate_particleFlow.pt", PFCandidate_pt);
  stream.select("recoPFJet_ak5PFJets.chargedEmEnergyFraction", PFJet_chargedEmEnergyFraction);
  stream.select("recoPFJet_ak5PFJets.chargedHadronEnergyFraction", PFJet_chargedHadronEnergyFraction);
  stream.select("recoPFJet_ak5PFJets.energy", PFJet_energy);
  stream.select("recoPFJet_ak5PFJets.et", PFJet_et);
  stream.select("recoPFJet_ak5PFJets.eta", PFJet_eta);
  stream.select("recoPFJet_ak5PFJets.neutralEmEnergyFraction", PFJet_neutralEmEnergyFraction);
  stream.select("recoPFJet_ak5PFJets.neutralHadronEnergyFraction", PFJet_neutralHadronEnergyFraction);
  stream.select("recoPFJet_ak5PFJets.p", PFJet_p);
  stream.select("recoPFJet_ak5PFJets.phi", PFJet_phi);
  stream.select("recoPFJet_ak5PFJets.pt", PFJet_pt);
  stream.select("recoPFMET_pfMet.energy", PFMET_energy);
  stream.select("recoPFMET_pfMet.et", PFMET_et);
  stream.select("recoPFMET_pfMet.eta", PFMET_eta);
  stream.select("recoPFMET_pfMet.p", PFMET_p);
  stream.select("recoPFMET_pfMet.phi", PFMET_phi);
  stream.select("recoPFMET_pfMet.pt", PFMET_pt);
  stream.select("recoPFTau_hpsTancTaus.energy", PFTau1_energy);
  stream.select("recoPFTau_hpsTancTaus.et", PFTau1_et);
  stream.select("recoPFTau_hpsTancTaus.eta", PFTau1_eta);
  stream.select("recoPFTau_hpsTancTaus.p", PFTau1_p);
  stream.select("recoPFTau_hpsTancTaus.phi", PFTau1_phi);
  stream.select("recoPFTau_hpsTancTaus.pt", PFTau1_pt);
  stream.select("recoPFTau_shrinkingConePFTauProducer.energy", PFTau2_energy);
  stream.select("recoPFTau_shrinkingConePFTauProducer.et", PFTau2_et);
  stream.select("recoPFTau_shrinkingConePFTauProducer.eta", PFTau2_eta);
  stream.select("recoPFTau_shrinkingConePFTauProducer.p", PFTau2_p);
  stream.select("recoPFTau_shrinkingConePFTauProducer.phi", PFTau2_phi);
  stream.select("recoPFTau_shrinkingConePFTauProducer.pt", PFTau2_pt);
  stream.select("recoPFTau_hpsPFTauProducer.energy", PFTau_energy);
  stream.select("recoPFTau_hpsPFTauProducer.et", PFTau_et);
  stream.select("recoPFTau_hpsPFTauProducer.eta", PFTau_eta);
  stream.select("recoPFTau_hpsPFTauProducer.p", PFTau_p);
  stream.select("recoPFTau_hpsPFTauProducer.phi", PFTau_phi);
  stream.select("recoPFTau_hpsPFTauProducer.pt", PFTau_pt);
  stream.select("recoTrack_generalTracks.charge", Track_charge);
  stream.select("recoTrack_generalTracks.d0", Track_d0);
  stream.select("recoTrack_generalTracks.d0Error", Track_d0Error);
  stream.select("recoTrack_generalTracks.dz", Track_dz);
  stream.select("recoTrack_generalTracks.dzError", Track_dzError);
  stream.select("recoTrack_generalTracks.eta", Track_eta);
  stream.select("recoTrack_generalTracks.numberOfLostHits", Track_numberOfLostHits);
  stream.select("recoTrack_generalTracks.numberOfValidHits", Track_numberOfValidHits);
  stream.select("recoTrack_generalTracks.p", Track_p);
  stream.select("recoTrack_generalTracks.phi", Track_phi);
  stream.select("recoTrack_generalTracks.pt", Track_pt);
  stream.select("recoTrack_generalTracks.ptError", Track_ptError);
  stream.select("recoTrack_generalTracks.qualityMask", Track_qualityMask);
  stream.select("recoTrack_generalTracks.trackerExpectedHitsInner_numberOfLostHits", Track_trackerExpectedHitsInner_numberOfLostHits);
  stream.select("recoTrack_generalTracks.trackerExpectedHitsOuter_numberOfLostHits", Track_trackerExpectedHitsOuter_numberOfLostHits);
  stream.select("recoTrack_generalTracks.vx", Track_vx);
  stream.select("recoTrack_generalTracks.vy", Track_vy);
  stream.select("recoTrack_generalTracks.vz", Track_vz);
  stream.select("recoVertex_offlinePrimaryVertices.chi2", Vertex_chi2);
  stream.select("recoVertex_offlinePrimaryVertices.hasRefittedTracks", Vertex_hasRefittedTracks);
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
  stream.select("nrecoCaloCluster_hfEMClusters", nCaloCluster);
  stream.select("nrecoCaloCluster_hybridSuperClusters_hybridBarrelBasicClusters", nCaloCluster1);
  stream.select("nrecoCaloCluster_hybridSuperClusters_uncleanOnlyHybridBarrelBasicClusters", nCaloCluster2);
  stream.select("nrecoCaloCluster_multi5x5SuperClusters_multi5x5EndcapBasicClusters", nCaloCluster3);
  stream.select("nrecoGenParticle_genParticles", nGenParticle);
  stream.select("nrecoGsfElectron_gsfElectrons", nGsfElectron);
  stream.select("nrecoMuon_muons", nMuon);
  stream.select("nrecoPFCandidate_particleFlow", nPFCandidate);
  stream.select("nrecoPFJet_ak5PFJets", nPFJet);
  stream.select("nrecoPFMET_pfMet", nPFMET);
  stream.select("nrecoPFTau_hpsPFTauProducer", nPFTau);
  stream.select("nrecoPFTau_hpsTancTaus", nPFTau1);
  stream.select("nrecoPFTau_shrinkingConePFTauProducer", nPFTau2);
  stream.select("nrecoTrack_generalTracks", nTrack);
  stream.select("nrecoVertex_offlinePrimaryVertices", nVertex);

}
}; // end namespace evt
#endif
