/******************************************************************************
 * @file DecayVertex.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.h
 * ****************************************************************************/

// Decay vertex header
#include "DecayVertex.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

namespace AuxVertex
{
  DecayVertex::DecayVertex()
  {}
  DecayVertex::~DecayVertex()
  {}

  DecayVertex::DecayVertex(
            const art::Ptr<recob::Vertex> &nuVertex,
            const art::Ptr<recob::Vertex> &t1Vertex,
            const art::Ptr<recob::Vertex> &t2Vertex,
            const art::Ptr<recob::Track> &t1Track,
            const art::Ptr<recob::Track> &t2Track,
            const std::vector<art::Ptr<recob::Hit>> &t1Hits,
            const std::vector<art::Ptr<recob::Hit>> &t2Hits)
  {
    /*
    This function creates a HSN decay vertex candidate object.
    It is assigned when a neutrino with two (and only two) track daughters has been found.
    This object is build using the vertex pointers to the neutrino and the two start points of the tracks, the actual track objects and all the hits associated with the two track objects.
    */

    // Set default mock attributes for the Decay Vertex candidate
    // The real values for these attributes will be assigned later
    fChannelLoc = {-1,-1,-1};
    fTickLoc = {-1.,-1.,-1.};
    fProngChannelLoc = {{-1,-1,-1}, {-1,-1,-1}};
    fProngTickLoc = {{-1.,-1.,-1.}, {-1.,-1.,-1.}};
    fIsDetLocAssigned = false;
    fIsInsideTPC = false;

    // Store pointers to reconstructed objects provided as input in internal attributes
    fNuVertex = nuVertex;
    fProngVertex = {t1Vertex, t2Vertex};
    fProngTrack = {t1Track, t2Track};
    fProngHits = {t1Hits, t2Hits};

    // Use pointers to reconstructed objects to obtain start coordinates for vertices.
    double nuVertexPosition[3], t1VertexPosition[3], t2VertexPosition[3];
    nuVertex->XYZ(nuVertexPosition);
    t1Vertex->XYZ(t1VertexPosition);
    t2Vertex->XYZ(t2VertexPosition);
    fX = nuVertexPosition[0];
    fY = nuVertexPosition[1];
    fZ = nuVertexPosition[2];
    fProngX = {(float) t1VertexPosition[0], (float) t2VertexPosition[0]};
    fProngY = {(float) t1VertexPosition[1], (float) t2VertexPosition[1]};
    fProngZ = {(float) t1VertexPosition[2], (float) t2VertexPosition[2]};
    // Do the same for tracks.
    fProngStartX = {(float) t1Track->Start().X(), (float) t2Track->Start().X()};
    fProngStartY = {(float) t1Track->Start().Y(), (float) t2Track->Start().Y()};
    fProngStartZ = {(float) t1Track->Start().Z(), (float) t2Track->Start().Z()};
    fProngEndX = {(float) t1Track->End()[0], (float) t2Track->End()[0]};
    fProngEndY = {(float) t1Track->End()[1], (float) t2Track->End()[1]};
    fProngEndZ = {(float) t1Track->End()[2], (float) t2Track->End()[2]};

    // Calculate length, theta, phi and number of hits
    fProngLength = {(float) t1Track->Length(), (float) t2Track->Length()};
    fProngTheta = {(float) t1Track->Theta(), (float) t2Track->Theta()};
    fProngPhi = {(float) t1Track->Phi(), (float) t2Track->Phi()};
    fProngNumHits = {(int) t1Hits.size(), (int) t2Hits.size()};

    // Calculate track directions
    fProngDirX = {(float) t1Track->VertexMomentumVector().X(), (float) t2Track->VertexMomentumVector().X()};
    fProngDirY = {(float) t1Track->VertexMomentumVector().Y(), (float) t2Track->VertexMomentumVector().Y()};
    fProngDirZ = {(float) t1Track->VertexMomentumVector().Z(), (float) t2Track->VertexMomentumVector().Z()};

    // Calculate kinematic quantities for each prong
    trkf::TrackMomentumCalculator tmc;
    // Prong momentum magnitude
    fProngMomMag_ByRange_AssMuon = {(float) tmc.GetTrackMomentum(fProngLength[0],13), (float) tmc.GetTrackMomentum(fProngLength[1],13)};
    // Prong momentum components
    fProngMom_ByRange_AssMuonX = {(float) fProngDirX[0]*fProngMomMag_ByRange_AssMuon[0], (float) fProngDirX[1]*fProngMomMag_ByRange_AssMuon[1]};
    fProngMom_ByRange_AssMuonY = {(float) fProngDirY[0]*fProngMomMag_ByRange_AssMuon[0], (float) fProngDirY[1]*fProngMomMag_ByRange_AssMuon[1]};
    fProngMom_ByRange_AssMuonZ = {(float) fProngDirZ[0]*fProngMomMag_ByRange_AssMuon[0], (float) fProngDirZ[1]*fProngMomMag_ByRange_AssMuon[1]};
    // Prong energy
    float muonMass = 0.10566;
    float e1 = sqrt(pow(muonMass,2.) + pow(fProngMomMag_ByRange_AssMuon[0],2.));
    float e2 = sqrt(pow(muonMass,2.) + pow(fProngMomMag_ByRange_AssMuon[1],2.));
    fProngEnergy_ByRange_AssMuon = {(float) e1, (float) e2};
    /**/
    // Calculate kinematic quantities for parent neutrino
    // Total momentum components
    fTotMom_ByRange_AssMuonX = fProngMom_ByRange_AssMuonX[0] + fProngMom_ByRange_AssMuonX[1];
    fTotMom_ByRange_AssMuonY = fProngMom_ByRange_AssMuonY[0] + fProngMom_ByRange_AssMuonY[1];
    fTotMom_ByRange_AssMuonZ = fProngMom_ByRange_AssMuonZ[0] + fProngMom_ByRange_AssMuonZ[1];
    // Total momentum magnitude
    fTotMomMag_ByRange_AssMuon = sqrt(pow(fTotMom_ByRange_AssMuonX,2.) + pow(fTotMom_ByRange_AssMuonY,2.) + pow(fTotMom_ByRange_AssMuonZ,2.));
    // Total direction components
    fTotDir_ByRange_AssMuonX = fTotMom_ByRange_AssMuonX/fTotMomMag_ByRange_AssMuon;
    fTotDir_ByRange_AssMuonY = fTotMom_ByRange_AssMuonY/fTotMomMag_ByRange_AssMuon;
    fTotDir_ByRange_AssMuonZ = fTotMom_ByRange_AssMuonZ/fTotMomMag_ByRange_AssMuon;
    // Total direction angles
    fTotTheta_ByRange_AssMuon = acos(fTotDir_ByRange_AssMuonZ);
    fTotPhi_ByRange_AssMuon = atan2(fTotDir_ByRange_AssMuonY,fTotDir_ByRange_AssMuonX);
    // Total energy
    fTotEnergy_ByRange_AssMuon = fProngEnergy_ByRange_AssMuon[0] + fProngEnergy_ByRange_AssMuon[1];
    // Invariant mass
    fInvMass_ByRange_AssMuon = sqrt(pow(fTotEnergy_ByRange_AssMuon,2.) - pow(fTotMomMag_ByRange_AssMuon,2.));

    // Calculate opening angle
    std::vector<double> startDirection1 = {
      t1Track->StartDirection().X(),
      t1Track->StartDirection().Y(),
      t1Track->StartDirection().Z()
    };
    std::vector<double> startDirection2 = {
      t2Track->StartDirection().X(),
      t2Track->StartDirection().Y(),
      t2Track->StartDirection().Z()
    };
    float magnitude1 = sqrt(startDirection1[0]*startDirection1[0] + startDirection1[1]*startDirection1[1] + startDirection1[2]*startDirection1[2]);
    float magnitude2 = sqrt(startDirection2[0]*startDirection2[0] + startDirection2[1]*startDirection2[1] + startDirection2[2]*startDirection2[2]);
    float dotProduct = startDirection1[0]*startDirection2[0] + startDirection1[1]*startDirection2[1] + startDirection1[2]*startDirection2[2];
    fOpeningAngle = acos(dotProduct / (magnitude1*magnitude2));
    // Calculate start point to neutrino vertex distance
    float prong1_distance = sqrt(pow(fX - fProngX[0],2.) + pow(fY - fProngY[0],2.) + pow(fZ - fProngZ[0],2.));
    float prong2_distance = sqrt(pow(fX - fProngX[1],2.) + pow(fY - fProngY[1],2.) + pow(fZ - fProngZ[1],2.));
    fProngStartToNeutrinoDistance = {prong1_distance,prong2_distance};
  } //  END constructor DecayVertex

  // Getters
  // Pointers
  art::Ptr<recob::Vertex> DecayVertex::GetNuVertex() const {return fNuVertex;}
  art::Ptr<recob::Vertex> DecayVertex::GetProngVertex(int prong) const {return fProngVertex[prong];}
  art::Ptr<recob::Track> DecayVertex::GetProngTrack(int prong) const {return fProngTrack[prong];}
  std::vector<art::Ptr<recob::Hit>> DecayVertex::GetProngHits(int prong) const {return fProngHits[prong];}
  std::vector<art::Ptr<recob::Hit>> DecayVertex::GetTotHits() const {return fTotHitsInMaxRadius;}
  // Coordinates
  float DecayVertex::GetX() const {return fX;}
  float DecayVertex::GetY() const {return fY;}
  float DecayVertex::GetZ() const {return fZ;}
  float DecayVertex::GetProngX(int prong) const {return fProngX[prong];}
  float DecayVertex::GetProngY(int prong) const {return fProngY[prong];}
  float DecayVertex::GetProngZ(int prong) const {return fProngZ[prong];}
  float DecayVertex::GetProngStartX(int prong) const {return fProngStartX[prong];}
  float DecayVertex::GetProngStartY(int prong) const {return fProngStartY[prong];}
  float DecayVertex::GetProngStartZ(int prong) const {return fProngStartZ[prong];}
  float DecayVertex::GetProngEndX(int prong) const {return fProngEndX[prong];}
  float DecayVertex::GetProngEndY(int prong) const {return fProngEndY[prong];}
  float DecayVertex::GetProngEndZ(int prong) const {return fProngEndZ[prong];}
  // Wire-Tick coordinates
  int DecayVertex::GetChannelLoc(int plane) const {return fChannelLoc[plane];}
  float DecayVertex::GetTickLoc(int plane) const {return fTickLoc[plane];}
  int DecayVertex::GetProngChannelLoc(int prong,int plane) const {return fProngChannelLoc[prong][plane];}
  float DecayVertex::GetProngTickLoc(int prong,int plane) const {return fProngTickLoc[prong][plane];}
  // Direction
  float DecayVertex::GetProngDirX(int prong) const {return fProngDirX[prong];}
  float DecayVertex::GetProngDirY(int prong) const {return fProngDirY[prong];}
  float DecayVertex::GetProngDirZ(int prong) const {return fProngDirZ[prong];}
  float DecayVertex::GetProngTheta(int prong) const {return fProngTheta[prong];}
  float DecayVertex::GetProngPhi(int prong) const {return fProngPhi[prong];}
  // Prong momentum (by range, assuming both muons)
  float DecayVertex::GetProngMomMag_ByRange_AssMuon(int prong) const {return fProngMomMag_ByRange_AssMuon[prong];}
  float DecayVertex::GetProngMom_ByRange_AssMuonX(int prong) const {return fProngMom_ByRange_AssMuonX[prong];}
  float DecayVertex::GetProngMom_ByRange_AssMuonY(int prong) const {return fProngMom_ByRange_AssMuonY[prong];}
  float DecayVertex::GetProngMom_ByRange_AssMuonZ(int prong) const {return fProngMom_ByRange_AssMuonZ[prong];}
  float DecayVertex::GetProngEnergy_ByRange_AssMuon(int prong) const {return fProngEnergy_ByRange_AssMuon[prong];}
  // Tot momentum (by range, assuming both muons)
  float DecayVertex::GetTotMomMag_ByRange_AssMuon() const {return fTotMomMag_ByRange_AssMuon;}
  float DecayVertex::GetTotMom_ByRange_AssMuonX() const {return fTotMom_ByRange_AssMuonX;}
  float DecayVertex::GetTotMom_ByRange_AssMuonY() const {return fTotMom_ByRange_AssMuonY;}
  float DecayVertex::GetTotMom_ByRange_AssMuonZ() const {return fTotMom_ByRange_AssMuonZ;}
  float DecayVertex::GetTotEnergy_ByRange_AssMuon() const {return fTotEnergy_ByRange_AssMuon;}
  float DecayVertex::GetInvMass_ByRange_AssMuon() const {return fInvMass_ByRange_AssMuon;}
  // Tot momentum direction (by range, assuming both muons)
  float DecayVertex::GetTotDir_ByRange_AssMuonX() const {return fTotDir_ByRange_AssMuonX;}
  float DecayVertex::GetTotDir_ByRange_AssMuonY() const {return fTotDir_ByRange_AssMuonY;}
  float DecayVertex::GetTotDir_ByRange_AssMuonZ() const {return fTotDir_ByRange_AssMuonZ;}
  float DecayVertex::GetTotTheta_ByRange_AssMuon() const {return fTotTheta_ByRange_AssMuon;}
  float DecayVertex::GetTotPhi_ByRange_AssMuon() const {return fTotPhi_ByRange_AssMuon;}
  // Others
  float DecayVertex::GetProngLength(int prong) const {return fProngLength[prong];}
  int DecayVertex::GetProngNumHits(int prong) const {return fProngNumHits[prong];}
  float DecayVertex::GetProngStartToNeutrinoDistance(int prong) const {return fProngStartToNeutrinoDistance[prong];}
  float DecayVertex::GetOpeningAngle() const {return fOpeningAngle;}
  bool DecayVertex::IsInsideTPC() const {return fIsInsideTPC;}
  bool DecayVertex::IsDetLocAssigned() const {return fIsDetLocAssigned;}


  // Setters
  void DecayVertex::SetChannelLoc(int channel0, int channel1, int channel2) {fChannelLoc = {channel0,channel1,channel2}; return;}
  void DecayVertex::SetTickLoc(float tick0, float tick1, float tick2) {fTickLoc = {tick0, tick1, tick2}; return;}
  void DecayVertex::SetProngChannelLoc(int prong, int channel0, int channel1, int channel2) {fProngChannelLoc[prong] =  {channel0,channel1,channel2}; return;}
  void DecayVertex::SetProngTickLoc(int prong, float tick0, float tick1, float tick2) {fProngTickLoc[prong] = {tick0, tick1, tick2}; return;}
  void DecayVertex::SetProngXYZ(int prong, float x, float y, float z) {fProngX[prong] = x; fProngY[prong] = y; fProngZ[prong] = z; return;}
  void DecayVertex::SetIsInsideTPC(bool val) {fIsInsideTPC = val; return;}
  void DecayVertex::SetIsDetLocAssigned(bool val) {fIsDetLocAssigned = val; return;}
  void DecayVertex::SetTotHits(std::vector<art::Ptr<recob::Hit>> totHitsInMaxRadius) {fTotHitsInMaxRadius = totHitsInMaxRadius; return;}

  void DecayVertex::SetDetectorCoordinates(
    const std::vector<double>& minTpcBound,
    const std::vector<double>& maxTpcBound,
    geo::GeometryCore const* geometry,
    detinfo::DetectorProperties const* detectorProperties)
  {
    /* This functions call special geometry and detectorProperties objects
    in order to translate x,y,z coordinates in wire,tick coordinates.
    This allows the determination of the vertices location in the event display
    */
    // Get spatial coordinates and mark vertex as assigned
    float xyz[3] = {fX,fY,fZ};
    float prong1_xyz[3] = {fProngX[0],fProngY[0],fProngZ[0]};
    float prong2_xyz[3] = {fProngX[1],fProngY[1],fProngZ[1]};

    fIsDetLocAssigned = true;

    // Check whether coordinates are inside TPC
    double extraEdge = 0;
    bool isInsideX = (xyz[0]>minTpcBound[0]+extraEdge &&
      xyz[0]<maxTpcBound[0]-extraEdge);
    bool isInsideY = (xyz[1]>minTpcBound[1]+extraEdge &&
      xyz[1]<maxTpcBound[1]-extraEdge);
    bool isInsideZ = (xyz[2]>minTpcBound[2]+extraEdge &&
      xyz[2]<maxTpcBound[2]-extraEdge);

    // If vertex is inside TPC, determine channel/tick coordinates and assign them
    if (isInsideX && isInsideY && isInsideZ)
    {
      fIsInsideTPC = true;
      raw::ChannelID_t channel0 = geometry->NearestChannel(xyz,0);
      raw::ChannelID_t channel1 = geometry->NearestChannel(xyz,1);
      raw::ChannelID_t channel2 = geometry->NearestChannel(xyz,2);
      double tick0 = detectorProperties->ConvertXToTicks(xyz[0], 0, 0, 0);
      double tick1 = detectorProperties->ConvertXToTicks(xyz[0], 1, 0, 0);
      double tick2 = detectorProperties->ConvertXToTicks(xyz[0], 2, 0, 0);
      fChannelLoc = {(int) channel0,(int) channel1,(int) channel2};
      fTickLoc = { (float) tick0, (float) tick1, (float) tick2};

      raw::ChannelID_t prong1_channel0 = geometry->NearestChannel(prong1_xyz,0);
      raw::ChannelID_t prong1_channel1 = geometry->NearestChannel(prong1_xyz,1);
      raw::ChannelID_t prong1_channel2 = geometry->NearestChannel(prong1_xyz,2);
      double prong1_tick0 = detectorProperties->ConvertXToTicks(prong1_xyz[0], 0, 0, 0);
      double prong1_tick1 = detectorProperties->ConvertXToTicks(prong1_xyz[0], 1, 0, 0);
      double prong1_tick2 = detectorProperties->ConvertXToTicks(prong1_xyz[0], 2, 0, 0);
      raw::ChannelID_t prong2_channel0 = geometry->NearestChannel(prong2_xyz,0);
      raw::ChannelID_t prong2_channel1 = geometry->NearestChannel(prong2_xyz,1);
      raw::ChannelID_t prong2_channel2 = geometry->NearestChannel(prong2_xyz,2);
      double prong2_tick0 = detectorProperties->ConvertXToTicks(prong2_xyz[0], 0, 0, 0);
      double prong2_tick1 = detectorProperties->ConvertXToTicks(prong2_xyz[0], 1, 0, 0);
      double prong2_tick2 = detectorProperties->ConvertXToTicks(prong2_xyz[0], 2, 0, 0);

      fProngChannelLoc = {{ (int) prong1_channel0, (int) prong1_channel1, (int) prong1_channel2}, { (int) prong2_channel0, (int) prong2_channel1, (int) prong2_channel2}};
      fProngTickLoc = {{ (float) prong1_tick0, (float) prong1_tick1, (float) prong1_tick2}, { (float) prong2_tick0, (float) prong2_tick1, (float) prong2_tick2}};
      return;
    }

    // Else flag it as outside the TPC and exit function
    else
      {
        fIsInsideTPC = false;
        return;
      }
  } // END function SetDetectorCoordinates


  // Printers
  void DecayVertex::PrintInformation() const{
    int fStartWire[3] = {0,2399,4798};
    printf("\n-|Vertex information|\n");
    printf("|_Vertex inside TPC: %d\n", fIsInsideTPC);
    printf("|_Spatial Coordinates: [%.1f, %.1f, %.1f]\n",fX,fY,fZ);
    printf("|_Prongs lengths: [%.1f,%.1f]\n", fProngLength[0], fProngLength[1]);
    printf("|_Prongs theta: [%.1f,%.1f]\n", fProngTheta[0], fProngTheta[1]);
    printf("|_Prongs phi: [%.1f,%.1f]\n", fProngPhi[0], fProngPhi[1]);
    printf("|_Prongs hit number: [%i,%i]\n", fProngNumHits[0], fProngNumHits[1]);
    printf("|_Detector location assigned: %d\n", fIsDetLocAssigned);
    if (fIsDetLocAssigned)
    {
      printf("|_Channel coordinates: [%i, %i, %i]\n",fChannelLoc[0]-fStartWire[0],fChannelLoc[1]-fStartWire[1],fChannelLoc[2]-fStartWire[2]);
      printf("|_Ticks coordinates: [%.1f, %.1f, %.1f]\n",fTickLoc[0],fTickLoc[1],fTickLoc[2]);
    }
    return;
  }
} // END namespace AuxVertex
