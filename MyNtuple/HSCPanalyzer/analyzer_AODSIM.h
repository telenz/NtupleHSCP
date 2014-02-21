#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Fri Feb 21 12:39:41 2014 by mkanalyzer.py
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
double	GenRunInfoProduct_crossSection;
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
  stream.select("GenRunInfoProduct_generator.crossSection", GenRunInfoProduct_crossSection);
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
  stream.select("nrecoPFJet_ak5PFJets", nPFJet);
  stream.select("nrecoPFMET_pfMet", nPFMET);
  stream.select("nrecoVertex_offlinePrimaryVertices", nVertex);

}
}; // end namespace evt
#endif
