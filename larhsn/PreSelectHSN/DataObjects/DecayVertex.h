/******************************************************************************
 * @file DecayVertex.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.cxx
 * ****************************************************************************/

#ifndef DECAYVERTEX_H
#define DECAYVERTEX_H

#include "TMatrixD.h"
#include "TVector3.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


namespace AuxVertex
{

  // Decay vertex class and functions
  class DecayVertex
  {
  public:
    // Constructor and destructor
    DecayVertex();
    virtual ~DecayVertex();

    DecayVertex(
            const art::Ptr<recob::Vertex> &nuVertex,
            const art::Ptr<recob::Vertex> &t1Vertex,
            const art::Ptr<recob::Vertex> &t2Vertex,
            const art::Ptr<recob::Track> &t1Track,
            const art::Ptr<recob::Track> &t2Track,
            const std::vector<art::Ptr<recob::Hit>> &t1Hits,
            const std::vector<art::Ptr<recob::Hit>> &t2Hits);

    // New Getters
    art::Ptr<recob::Vertex> GetNuVertex() const;
    art::Ptr<recob::Vertex> GetProngVertex(int prong) const;
    art::Ptr<recob::Track> GetProngTrack(int prong) const;
    std::vector<art::Ptr<recob::Hit>> GetProngHits(int prong) const;
    std::vector<art::Ptr<recob::Hit>> GetTotHits() const;
    float GetX() const;
    float GetY() const;
    float GetZ() const;
    float GetProngX(int par) const;
    float GetProngY(int par) const;
    float GetProngZ(int par) const;
    float GetProngDirPx(int prong) const;
    float GetProngDirPy(int prong) const;
    float GetProngDirPz(int prong) const;
    float GetProngMagP(int prong) const;
    float GetProngLength(int prong) const;
    float GetProngTheta(int prong) const;
    float GetProngPhi(int prong) const;
    int GetProngNumHits(int prong) const;
    float GetProngStartToNeutrinoDistance(int prong) const;
    float GetOpeningAngle() const;
    int GetChannelLoc(int plane) const;
    float GetTickLoc(int plane) const;
    int GetProngChannelLoc(int par,int plane) const;
    float GetProngTickLoc(int par,int plane) const;
    bool IsInsideTPC() const;
    bool IsDetLocAssigned() const;

    // Setters
    void SetDetectorCoordinates(
      const std::vector<double>& minTpcBound,
      const std::vector<double>& maxTpcBound,
      geo::GeometryCore const* geometry,
      detinfo::DetectorProperties const* detectorProperties);
    void SetChannelLoc(int channel0, int channel1, int channel2);
    void SetTickLoc(float tick0, float tick1, float tick2);
    void SetProngChannelLoc(int par, int channel0, int channel1, int channel2);
    void SetProngTickLoc(int par, float tick0, float tick1, float tick2);
    void SetProngXYZ(int par, float x, float y, float z);
    void SetIsInsideTPC(bool val);
    void SetIsDetLocAssigned(bool val);
    void SetTotHits(std::vector<art::Ptr<recob::Hit>> totHitsInMaxRadius);
    void DetermineMomenta();

    // Printers
    void PrintInformation() const;

    private:
      // Data products pointers
      art::Ptr<recob::Vertex> fNuVertex;
      std::vector<art::Ptr<recob::Vertex>> fProngVertex;
      std::vector<art::Ptr<recob::Track>> fProngTrack;
      std::vector<std::vector<art::Ptr<recob::Hit>>> fProngHits;
      std::vector<art::Ptr<recob::Hit>> fTotHitsInMaxRadius;

      // Coordinates
      float fX, fY, fZ; // Spatial coordinates of the vertex inside the detector.
      std::vector<float> fProngX, fProngY, fProngZ; // Spatial coordinates of the parent of the vertex inside the detector.
      std::vector<int> fChannelLoc; // Nearest channel in each plane.
      std::vector<float> fTickLoc; // Nearest time tick in each plane.
      std::vector<std::vector<int>> fProngChannelLoc; // Nearest channel in each plane for the vertex parent.
      std::vector<std::vector<float>> fProngTickLoc; // Nearest time tick in each plane for the vertex parent.

      // Physical variables
      float fOpeningAngle; // Opening angle between the two prongs

      std::vector<std::vector<float>> fProngMomentumDir; // Momentum direction of each prong
      std::vector<float> fProngMomentumMag; // Momentum magnitude of each prong
      std::vector<float> fProngLength; // Length of each prong
      std::vector<float> fProngTheta; // Theta angle of each prong
      std::vector<float> fProngPhi; // Phi angle of each prong
      std::vector<float> fProngStartToNeutrinoDistance; // Distance from start point to neutrino vertex for each prong
      std::vector<int> fProngNumHits; // Number of hits associated with each prong

      // Status
      bool fIsInsideTPC; // Whetehr the vertex is inside the TPC.
      bool fIsDetLocAssigned; // Whether channel/tick coordinates have been determined.
  };
  
} //END namespace AuxVertex

#endif