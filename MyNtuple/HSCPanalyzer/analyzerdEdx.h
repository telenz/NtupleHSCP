#ifndef ANALYZER_H
#define ANALYZER_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Tue Jan 21 11:58:54 2014 by mkanalyzer.py
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
std::vector<double>	GenParticle_px(1000,0);
std::vector<double>	GenParticle_py(1000,0);
std::vector<double>	GenParticle_pz(1000,0);
std::vector<int>	GenParticle_status(1000,0);
std::vector<double>	GenParticle_vx(1000,0);
std::vector<double>	GenParticle_vy(1000,0);
std::vector<double>	GenParticle_vz(1000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_1000(20000000,0);
std::vector<double>	Track_dEdxHitsNPHarm2_7(20000000,0);
std::vector<double>	Track_dEdxNPHarm2(20000000,0);
std::vector<unsigned int>	Track_dEdxNPNoM(20000000,0);
std::vector<double>	Track_eta(20000000,0);
std::vector<double>	Track_phi(20000000,0);
std::vector<double>	Track_pt(20000000,0);
std::vector<double>	Track_px(20000000,0);
std::vector<double>	Track_py(20000000,0);
std::vector<double>	Track_pz(20000000,0);
int	nGenParticle;
int	nTrack;

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
std::vector<GenParticle_s> GenParticle(1000);

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
struct Track_s
{
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	eta;
  double	phi;
  double	dEdxNPHarm2;
  unsigned int	dEdxNPNoM;
  double	dEdxHitsNPHarm2_1000;
  double	dEdxHitsNPHarm2_7;
};
std::vector<Track_s> Track(20000000);

std::ostream& operator<<(std::ostream& os, const Track_s& o)
{
  char r[1024];
  os << "Track" << std::endl;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "px", (double)o.px); os << r;
  sprintf(r, "  %-32s: %f\n", "py", (double)o.py); os << r;
  sprintf(r, "  %-32s: %f\n", "pz", (double)o.pz); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxNPHarm2", (double)o.dEdxNPHarm2); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxNPNoM", (double)o.dEdxNPNoM); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_1000", (double)o.dEdxHitsNPHarm2_1000); os << r;
  sprintf(r, "  %-32s: %f\n", "dEdxHitsNPHarm2_7", (double)o.dEdxHitsNPHarm2_7); os << r;
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

inline void fillTrack()
{
  Track.resize(Track_pt.size());
  for(unsigned int i=0; i < Track.size(); ++i)
    {
      Track[i].pt	= Track_pt[i];
      Track[i].px	= Track_px[i];
      Track[i].py	= Track_py[i];
      Track[i].pz	= Track_pz[i];
      Track[i].eta	= Track_eta[i];
      Track[i].phi	= Track_phi[i];
      Track[i].dEdxNPHarm2	= Track_dEdxNPHarm2[i];
      Track[i].dEdxNPNoM	= Track_dEdxNPNoM[i];
      Track[i].dEdxHitsNPHarm2_1000	= Track_dEdxHitsNPHarm2_1000[i];
      Track[i].dEdxHitsNPHarm2_7	= Track_dEdxHitsNPHarm2_7[i];
    }
}


void fillObjects()
{
  fillGenParticle();
  fillTrack();
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
          Track_eta[i]	= Track_eta[j];
          Track_phi[i]	= Track_phi[j];
          Track_dEdxNPHarm2[i]	= Track_dEdxNPHarm2[j];
          Track_dEdxNPNoM[i]	= Track_dEdxNPNoM[j];
          Track_dEdxHitsNPHarm2_1000[i]	= Track_dEdxHitsNPHarm2_1000[j];
          Track_dEdxHitsNPHarm2_7[i]	= Track_dEdxHitsNPHarm2_7[j];
        }
      nTrack = n;
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
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_1000", Track_dEdxHitsNPHarm2_1000);
  stream.select("recoTrackHelper_TrackRefitter.dEdxHitsNPHarm2_7", Track_dEdxHitsNPHarm2_7);
  stream.select("recoTrackHelper_TrackRefitter.dEdxNPHarm2", Track_dEdxNPHarm2);
  stream.select("recoTrackHelper_TrackRefitter.dEdxNPNoM", Track_dEdxNPNoM);
  stream.select("recoTrackHelper_TrackRefitter.eta", Track_eta);
  stream.select("recoTrackHelper_TrackRefitter.phi", Track_phi);
  stream.select("recoTrackHelper_TrackRefitter.pt", Track_pt);
  stream.select("recoTrackHelper_TrackRefitter.px", Track_px);
  stream.select("recoTrackHelper_TrackRefitter.py", Track_py);
  stream.select("recoTrackHelper_TrackRefitter.pz", Track_pz);
  stream.select("nrecoGenParticle_genParticles", nGenParticle);
  stream.select("nrecoTrackHelper_TrackRefitter", nTrack);

}
}; // end namespace evt
#endif
