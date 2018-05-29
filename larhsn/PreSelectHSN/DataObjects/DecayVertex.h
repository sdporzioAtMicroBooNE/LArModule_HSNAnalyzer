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
// #include "larreco/RecoAlg/TrackMomentumCalculator.h"

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
    // Coordinates
    float GetX() const;
    float GetY() const;
    float GetZ() const;
    float GetProngX(int par) const;
    float GetProngY(int par) const;
    float GetProngZ(int par) const;
    float GetProngStartX(int par) const;
    float GetProngStartY(int par) const;
    float GetProngStartZ(int par) const;
    float GetProngEndX(int par) const;
    float GetProngEndY(int par) const;
    float GetProngEndZ(int par) const;
    // Wire-Tick coordinates
    int GetChannelLoc(int plane) const;
    float GetTickLoc(int plane) const;
    int GetProngChannelLoc(int par,int plane) const;
    float GetProngTickLoc(int par,int plane) const;
    // Direction
    float GetProngDirX(int prong) const;
    float GetProngDirY(int prong) const;
    float GetProngDirZ(int prong) const;
    float GetProngTheta(int prong) const;
    float GetProngPhi(int prong) const;
    // Prong momentum (by range, assuming both muons)
    float GetProngEnergy_ByRange_AssMuon(int prong) const;
    float GetProngMomMag_ByRange_AssMuon(int prong) const;
    float GetProngMom_ByRange_AssMuonX(int prong) const;
    float GetProngMom_ByRange_AssMuonY(int prong) const;
    float GetProngMom_ByRange_AssMuonZ(int prong) const;
    // Tot momentum (by range, assuming both muons)
    float GetTotEnergy_ByRange_AssMuon() const;
    float GetInvMass_ByRange_AssMuon() const;
    float GetTotMomMag_ByRange_AssMuon() const;
    float GetTotMom_ByRange_AssMuonX() const;
    float GetTotMom_ByRange_AssMuonY() const;
    float GetTotMom_ByRange_AssMuonZ() const;
    // Tot momentum direction (by range, assuming both muons)
    float GetTotDir_ByRange_AssMuonX() const;
    float GetTotDir_ByRange_AssMuonY() const;
    float GetTotDir_ByRange_AssMuonZ() const;
    float GetTotTheta_ByRange_AssMuon() const;
    float GetTotPhi_ByRange_AssMuon() const;
    // Others
    float GetProngLength(int prong) const;
    float GetOpeningAngle() const;
    int GetProngNumHits(int prong) const;
    float GetProngStartToNeutrinoDistance(int prong) const;
    // Status
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

    // Printers
    void PrintInformation() const;

    private:
      // Data products pointers
      art::Ptr<recob::Vertex> fNuVertex;
      std::vector<art::Ptr<recob::Vertex>> fProngVertex;
      std::vector<art::Ptr<recob::Track>> fProngTrack;
      std::vector<std::vector<art::Ptr<recob::Hit>>> fProngHits;
      std::vector<art::Ptr<recob::Hit>> fTotHitsInMaxRadius;

      // Coordinates of the pandora neutrino recob::Vertex object
      float fX, fY, fZ; // Spatial coordinates of the vertex inside the detector.
      std::vector<int> fChannelLoc; // Nearest channel in each plane.
      std::vector<float> fTickLoc; // Nearest time tick in each plane.
      /**/
      // Coordinates of the two prongs recob::Vertex objects
      std::vector<float> fProngX, fProngY, fProngZ; // Spatial coordinates of the vertex of the track inside the detector.
      std::vector<std::vector<int>> fProngChannelLoc; // Nearest channel in each plane for the vertex parent.
      std::vector<std::vector<float>> fProngTickLoc; // Nearest time tick in each plane for the vertex parent.
      /**/
      // Coordinates of the two prongs start and end points for the recob::Track objects
      std::vector<float> fProngStartX, fProngStartY, fProngStartZ; // Spatial coordinates for the start of the track.
      std::vector<float> fProngEndX, fProngEndY, fProngEndZ; // Spatial coordinates for the end of the track.

      // Track direction (no calorimetry data)
      std::vector<float> fProngDirX, fProngDirY, fProngDirZ; // Direction of momentum for each track.
      std::vector<float> fProngTheta, fProngPhi; // Direction angles of each prong.

      // Momentum information (measured by range and assuming both particles are muons)
      std::vector<float> fProngMom_ByRange_AssMuonX, fProngMom_ByRange_AssMuonY, fProngMom_ByRange_AssMuonZ; // Components of momentum.
      std::vector<float> fProngMomMag_ByRange_AssMuon, fProngEnergy_ByRange_AssMuon; // Momentum and energy of each prong.
      /**/
      float fTotMom_ByRange_AssMuonX, fTotMom_ByRange_AssMuonY, fTotMom_ByRange_AssMuonZ; // Momentum component for neutrino.
      float fTotDir_ByRange_AssMuonX, fTotDir_ByRange_AssMuonY, fTotDir_ByRange_AssMuonZ; // Direction components of neutrino
      float fTotTheta_ByRange_AssMuon, fTotPhi_ByRange_AssMuon; // Direction angles of neutrino
      float fTotMomMag_ByRange_AssMuon, fTotEnergy_ByRange_AssMuon, fInvMass_ByRange_AssMuon; // Total momentum, total energy and invariant mass.

      // Other variables
      std::vector<float> fProngLength; // Length of each prong
      float fOpeningAngle; // Opening angle between the two prongs
      std::vector<float> fProngStartToNeutrinoDistance; // Distance from start point to neutrino vertex for each prong
      std::vector<int> fProngNumHits; // Number of hits associated with each prong

      // Status
      bool fIsInsideTPC; // Whetehr the vertex is inside the TPC.
      bool fIsDetLocAssigned; // Whether channel/tick coordinates have been determined.
  };
  
} //END namespace AuxVertex

#endif