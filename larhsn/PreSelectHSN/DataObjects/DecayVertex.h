/******************************************************************************
 * @file DecayVertex.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.cxx
 * ****************************************************************************/

#ifndef DECAYVERTEX_H
#define DECAYVERTEX_H

// C++ standard libraries
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
            recob::Vertex const * nuVertex,
            recob::Vertex const * t1Vertex,
            recob::Vertex const * t2Vertex,
            recob::Track const * t1Track,
            recob::Track const * t2Track,
            std::vector<recob::Hit const *> t1Hits,
            std::vector<recob::Hit const *> t2Hits);

    // New Getters
    recob::Vertex const * GetNuVertex() const;
    recob::Vertex const * GetProngVertex(int prong) const;
    recob::Track const * GetProngTrack(int prong) const;
    std::vector<recob::Hit const *> GetProngHits(int prong) const;
    std::vector<recob::Hit const *> GetTotHits() const;
    float GetX() const;
    float GetY() const;
    float GetZ() const;
    float GetProngX(int par) const;
    float GetProngY(int par) const;
    float GetProngZ(int par) const;
    float GetProngLength(int prong) const;
    float GetProngTheta(int prong) const;
    float GetProngPhi(int prong) const;
    int GetProngNumHits(int prong) const;
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
    void SetTotHits(std::vector<recob::Hit const*> totHitsInMaxRadius);

    // Printers
    void PrintInformation() const;

    private:
      // Data products pointers
      recob::Vertex const * fNuVertex;
      std::vector<recob::Vertex const *> fProngVertex;
      std::vector<recob::Track const *> fProngTrack;
      std::vector<std::vector<recob::Hit const *>> fProngHits;
      std::vector<recob::Hit const *> fTotHitsInMaxRadius;

      // Coordinates
      float fX, fY, fZ; // Spatial coordinates of the vertex inside the detector.
      std::vector<float> fProngX, fProngY, fProngZ; // Spatial coordinates of the parent of the vertex inside the detector.
      std::vector<int> fChannelLoc; // Nearest channel in each plane.
      std::vector<float> fTickLoc; // Nearest time tick in each plane.
      std::vector<std::vector<int>> fProngChannelLoc; // Nearest channel in each plane for the vertex parent.
      std::vector<std::vector<float>> fProngTickLoc; // Nearest time tick in each plane for the vertex parent.

      // Physical variables
      std::vector<float> fProngLength;
      std::vector<float> fProngTheta;
      std::vector<float> fProngPhi;
      std::vector<int> fProngNumHits;

      // Status
      bool fIsInsideTPC; // Whetehr the vertex is inside the TPC.
      bool fIsDetLocAssigned; // Whether channel/tick coordinates have been determined.
  };
  
} //END namespace AuxVertex

#endif