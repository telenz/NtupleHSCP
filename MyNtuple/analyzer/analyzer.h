#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Fri Jan 24 05:53:13 2014 by mkanalyzer.py
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
std::vector<int>	GenParticle_charge(100000,0);
std::vector<double>	GenParticle_energy(100000,0);
std::vector<double>	GenParticle_et(100000,0);
std::vector<double>	GenParticle_eta(100000,0);
std::vector<double>	GenParticle_mass(100000,0);
std::vector<size_t>	GenParticle_numberOfDaughters(100000,0);
std::vector<size_t>	GenParticle_numberOfMothers(100000,0);
std::vector<double>	GenParticle_p(100000,0);
std::vector<int>	GenParticle_pdgId(100000,0);
std::vector<double>	GenParticle_phi(100000,0);
std::vector<double>	GenParticle_pt(100000,0);
std::vector<double>	GenParticle_px(100000,0);
std::vector<double>	GenParticle_py(100000,0);
std::vector<double>	GenParticle_pz(100000,0);
std::vector<int>	GenParticle_status(100000,0);
std::vector<double>	GenParticle_vx(100000,0);
std::vector<double>	GenParticle_vy(100000,0);
std::vector<double>	GenParticle_vz(100000,0);
double	GenRunInfoProduct_crossSection;
std::vector<float>	SimTrack_charge(1000000,0);
std::vector<int>	SimTrack_genpartIndex(1000000,0);
std::vector<double>	SimTrack_momentum_energy(1000000,0);
std::vector<double>	SimTrack_momentum_eta(1000000,0);
std::vector<double>	SimTrack_momentum_phi(1000000,0);
std::vector<double>	SimTrack_momentum_pt(1000000,0);
std::vector<int>	SimTrack_noGenpart(1000000,0);
std::vector<int>	SimTrack_noVertex(1000000,0);
std::vector<unsigned int>	SimTrack_trackId(1000000,0);
std::vector<int>	SimTrack_type(1000000,0);
std::vector<int>	SimTrack_vertIndex(1000000,0);
std::vector<int>	SimVertex_noParent(1000000,0);
std::vector<int>	SimVertex_parentIndex(1000000,0);
std::vector<double>	SimVertex_position_t(1000000,0);
std::vector<double>	SimVertex_position_x(1000000,0);
std::vector<double>	SimVertex_position_y(1000000,0);
std::vector<double>	SimVertex_position_z(1000000,0);
std::vector<unsigned int>	SimVertex_vertexId(1000000,0);
int	nGenParticle;
int	nSimTrack;
int	nSimVertex;

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
  double	p;
  double	energy;
  double	et;
  double	px;
  double	py;
  double	pz;
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
std::vector<GenParticle_s> GenParticle(100000);

std::ostream& operator<<(std::ostream& os, const GenParticle_s& o)
{
  char r[1024];
  os << "GenParticle" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "p", (double)o.p); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
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
struct SimTrack_s
{
  float	charge;
  int	vertIndex;
  int	noVertex;
  int	genpartIndex;
  int	noGenpart;
  int	type;
  unsigned int	trackId;
  double	momentum_pt;
  double	momentum_phi;
  double	momentum_eta;
  double	momentum_energy;
};
std::vector<SimTrack_s> SimTrack(1000000);

std::ostream& operator<<(std::ostream& os, const SimTrack_s& o)
{
  char r[1024];
  os << "SimTrack" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "vertIndex", (double)o.vertIndex); os << r;
  sprintf(r, "  %-32s: %f\n", "noVertex", (double)o.noVertex); os << r;
  sprintf(r, "  %-32s: %f\n", "genpartIndex", (double)o.genpartIndex); os << r;
  sprintf(r, "  %-32s: %f\n", "noGenpart", (double)o.noGenpart); os << r;
  sprintf(r, "  %-32s: %f\n", "type", (double)o.type); os << r;
  sprintf(r, "  %-32s: %f\n", "trackId", (double)o.trackId); os << r;
  sprintf(r, "  %-32s: %f\n", "momentum_pt", (double)o.momentum_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "momentum_phi", (double)o.momentum_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "momentum_eta", (double)o.momentum_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "momentum_energy", (double)o.momentum_energy); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct SimVertex_s
{
  int	parentIndex;
  int	noParent;
  unsigned int	vertexId;
  double	position_x;
  double	position_y;
  double	position_z;
  double	position_t;
};
std::vector<SimVertex_s> SimVertex(1000000);

std::ostream& operator<<(std::ostream& os, const SimVertex_s& o)
{
  char r[1024];
  os << "SimVertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "parentIndex", (double)o.parentIndex); os << r;
  sprintf(r, "  %-32s: %f\n", "noParent", (double)o.noParent); os << r;
  sprintf(r, "  %-32s: %f\n", "vertexId", (double)o.vertexId); os << r;
  sprintf(r, "  %-32s: %f\n", "position_x", (double)o.position_x); os << r;
  sprintf(r, "  %-32s: %f\n", "position_y", (double)o.position_y); os << r;
  sprintf(r, "  %-32s: %f\n", "position_z", (double)o.position_z); os << r;
  sprintf(r, "  %-32s: %f\n", "position_t", (double)o.position_t); os << r;
  return os;
}
//-----------------------------------------------------------------------------

inline void fillGenParticle()
{
  GenParticle.resize(GenParticle_charge.size());
  for(unsigned int i=0; i < GenParticle.size(); ++i)
    {
      GenParticle[i].charge	= GenParticle_charge[i];
      GenParticle[i].p	= GenParticle_p[i];
      GenParticle[i].energy	= GenParticle_energy[i];
      GenParticle[i].et	= GenParticle_et[i];
      GenParticle[i].px	= GenParticle_px[i];
      GenParticle[i].py	= GenParticle_py[i];
      GenParticle[i].pz	= GenParticle_pz[i];
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

inline void fillSimTrack()
{
  SimTrack.resize(SimTrack_charge.size());
  for(unsigned int i=0; i < SimTrack.size(); ++i)
    {
      SimTrack[i].charge	= SimTrack_charge[i];
      SimTrack[i].vertIndex	= SimTrack_vertIndex[i];
      SimTrack[i].noVertex	= SimTrack_noVertex[i];
      SimTrack[i].genpartIndex	= SimTrack_genpartIndex[i];
      SimTrack[i].noGenpart	= SimTrack_noGenpart[i];
      SimTrack[i].type	= SimTrack_type[i];
      SimTrack[i].trackId	= SimTrack_trackId[i];
      SimTrack[i].momentum_pt	= SimTrack_momentum_pt[i];
      SimTrack[i].momentum_phi	= SimTrack_momentum_phi[i];
      SimTrack[i].momentum_eta	= SimTrack_momentum_eta[i];
      SimTrack[i].momentum_energy	= SimTrack_momentum_energy[i];
    }
}

inline void fillSimVertex()
{
  SimVertex.resize(SimVertex_parentIndex.size());
  for(unsigned int i=0; i < SimVertex.size(); ++i)
    {
      SimVertex[i].parentIndex	= SimVertex_parentIndex[i];
      SimVertex[i].noParent	= SimVertex_noParent[i];
      SimVertex[i].vertexId	= SimVertex_vertexId[i];
      SimVertex[i].position_x	= SimVertex_position_x[i];
      SimVertex[i].position_y	= SimVertex_position_y[i];
      SimVertex[i].position_z	= SimVertex_position_z[i];
      SimVertex[i].position_t	= SimVertex_position_t[i];
    }
}


void fillObjects()
{
  fillGenParticle();
  fillSimTrack();
  fillSimVertex();
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
          GenParticle_p[i]	= GenParticle_p[j];
          GenParticle_energy[i]	= GenParticle_energy[j];
          GenParticle_et[i]	= GenParticle_et[j];
          GenParticle_px[i]	= GenParticle_px[j];
          GenParticle_py[i]	= GenParticle_py[j];
          GenParticle_pz[i]	= GenParticle_pz[j];
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
       n = indexmap["SimTrack"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["SimTrack"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          SimTrack_charge[i]	= SimTrack_charge[j];
          SimTrack_vertIndex[i]	= SimTrack_vertIndex[j];
          SimTrack_noVertex[i]	= SimTrack_noVertex[j];
          SimTrack_genpartIndex[i]	= SimTrack_genpartIndex[j];
          SimTrack_noGenpart[i]	= SimTrack_noGenpart[j];
          SimTrack_type[i]	= SimTrack_type[j];
          SimTrack_trackId[i]	= SimTrack_trackId[j];
          SimTrack_momentum_pt[i]	= SimTrack_momentum_pt[j];
          SimTrack_momentum_phi[i]	= SimTrack_momentum_phi[j];
          SimTrack_momentum_eta[i]	= SimTrack_momentum_eta[j];
          SimTrack_momentum_energy[i]	= SimTrack_momentum_energy[j];
        }
      nSimTrack = n;
    }

  n = 0;
  try
    {
       n = indexmap["SimVertex"].size();
    }
  catch (...)
    {}
  if ( n > 0 )
    {
      std::vector<int>& index = indexmap["SimVertex"];
      for(int i=0; i < n; ++i)
        {
          int j = index[i];
          SimVertex_parentIndex[i]	= SimVertex_parentIndex[j];
          SimVertex_noParent[i]	= SimVertex_noParent[j];
          SimVertex_vertexId[i]	= SimVertex_vertexId[j];
          SimVertex_position_x[i]	= SimVertex_position_x[j];
          SimVertex_position_y[i]	= SimVertex_position_y[j];
          SimVertex_position_z[i]	= SimVertex_position_z[j];
          SimVertex_position_t[i]	= SimVertex_position_t[j];
        }
      nSimVertex = n;
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
  stream.select("recoGenParticle_genParticles.px", GenParticle_px);
  stream.select("recoGenParticle_genParticles.py", GenParticle_py);
  stream.select("recoGenParticle_genParticles.pz", GenParticle_pz);
  stream.select("recoGenParticle_genParticles.status", GenParticle_status);
  stream.select("recoGenParticle_genParticles.vx", GenParticle_vx);
  stream.select("recoGenParticle_genParticles.vy", GenParticle_vy);
  stream.select("recoGenParticle_genParticles.vz", GenParticle_vz);
  stream.select("GenRunInfoProduct_generator.crossSection", GenRunInfoProduct_crossSection);
  stream.select("SimTrack_g4SimHits.charge", SimTrack_charge);
  stream.select("SimTrack_g4SimHits.genpartIndex", SimTrack_genpartIndex);
  stream.select("SimTrack_g4SimHits.momentum_energy", SimTrack_momentum_energy);
  stream.select("SimTrack_g4SimHits.momentum_eta", SimTrack_momentum_eta);
  stream.select("SimTrack_g4SimHits.momentum_phi", SimTrack_momentum_phi);
  stream.select("SimTrack_g4SimHits.momentum_pt", SimTrack_momentum_pt);
  stream.select("SimTrack_g4SimHits.noGenpart", SimTrack_noGenpart);
  stream.select("SimTrack_g4SimHits.noVertex", SimTrack_noVertex);
  stream.select("SimTrack_g4SimHits.trackId", SimTrack_trackId);
  stream.select("SimTrack_g4SimHits.type", SimTrack_type);
  stream.select("SimTrack_g4SimHits.vertIndex", SimTrack_vertIndex);
  stream.select("SimVertex_g4SimHits.noParent", SimVertex_noParent);
  stream.select("SimVertex_g4SimHits.parentIndex", SimVertex_parentIndex);
  stream.select("SimVertex_g4SimHits.position_t", SimVertex_position_t);
  stream.select("SimVertex_g4SimHits.position_x", SimVertex_position_x);
  stream.select("SimVertex_g4SimHits.position_y", SimVertex_position_y);
  stream.select("SimVertex_g4SimHits.position_z", SimVertex_position_z);
  stream.select("SimVertex_g4SimHits.vertexId", SimVertex_vertexId);
  stream.select("nrecoGenParticle_genParticles", nGenParticle);
  stream.select("nSimTrack_g4SimHits", nSimTrack);
  stream.select("nSimVertex_g4SimHits", nSimVertex);

}
}; // end namespace evt
#endif
