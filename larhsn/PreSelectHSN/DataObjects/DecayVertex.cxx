/******************************************************************************
 * @file DecayVertex.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DecayVertex.h
 * ****************************************************************************/

// Decay vertex header
#include "DecayVertex.h"

namespace AuxVertex
{
  DecayVertex::DecayVertex()
  {}
  DecayVertex::~DecayVertex()
  {}

  DecayVertex::DecayVertex(
            recob::Vertex const * nuVertex,
            recob::Vertex const * t1Vertex,
            recob::Vertex const * t2Vertex,
            recob::Track const * t1Track,
            recob::Track const * t2Track,
            std::vector<recob::Hit const *> t1Hits,
            std::vector<recob::Hit const *> t2Hits)
  {
    // Set default attributes
    fChannelLoc = {-1,-1,-1};
    fTickLoc = {-1.,-1.,-1.};
    fProngChannelLoc = {{-1,-1,-1}, {-1,-1,-1}};
    fProngTickLoc = {{-1.,-1.,-1.}, {-1.,-1.,-1.}};
    fIsDetLocAssigned = false;
    fIsInsideTPC = false;

    // Assign arguments to attributes
    fNuVertex = nuVertex;
    fProngVertex = {t1Vertex, t2Vertex};
    fProngTrack = {t1Track, t2Track};
    fProngHits = {t1Hits, t2Hits};

    //std::cout << "\tThere are " << t1Hits.size() << " associated hits." << std::endl;
    //std::cout << "\tThere are " << t2Hits.size() << " associated hits." << std::endl;
    // Get vertex and daughters coordinates and assign them to attributes
    double nuVertexPosition[3], t1VertexPosition[3], t2VertexPosition[3];
    fNuVertex->XYZ(nuVertexPosition);
    fProngVertex[0]->XYZ(t1VertexPosition);
    fProngVertex[1]->XYZ(t2VertexPosition);
    fX = nuVertexPosition[0];
    fY = nuVertexPosition[1];
    fZ = nuVertexPosition[2];
    fProngX = {(float) t1VertexPosition[0], (float) t2VertexPosition[0]};
    fProngY = {(float) t1VertexPosition[1], (float) t2VertexPosition[1]};
    fProngZ = {(float) t1VertexPosition[2], (float) t2VertexPosition[2]};

    // Get physics variables
    fProngLength = {(float) t1Track->Length(), (float) t2Track->Length()};
    fProngTheta = {(float) t1Track->Theta(), (float) t2Track->Theta()};
    fProngPhi = {(float) t1Track->Phi(), (float) t2Track->Phi()};
    fProngNumHits = {(int) t1Hits.size(), (int) t2Hits.size()};


  }

  // Getters
  recob::Vertex const * DecayVertex::GetNuVertex() const {return fNuVertex;}
  recob::Vertex const * DecayVertex::GetProngVertex(int prong) const {return fProngVertex[prong];}
  recob::Track const * DecayVertex::GetProngTrack(int prong) const {return fProngTrack[prong];}
  std::vector<recob::Hit const *> DecayVertex::GetProngHits(int prong) const {return fProngHits[prong];}
  std::vector<recob::Hit const *> DecayVertex::GetTotHits() const {return fTotHitsInMaxRadius;}
  float DecayVertex::GetX() const {return fX;}
  float DecayVertex::GetY() const {return fY;}
  float DecayVertex::GetZ() const {return fZ;}
  float DecayVertex::GetProngX(int prong) const {return fProngX[prong];}
  float DecayVertex::GetProngY(int prong) const {return fProngY[prong];}
  float DecayVertex::GetProngZ(int prong) const {return fProngZ[prong];}
  int DecayVertex::GetChannelLoc(int plane) const {return fChannelLoc[plane];}
  float DecayVertex::GetTickLoc(int plane) const {return fTickLoc[plane];}
  int DecayVertex::GetProngChannelLoc(int prong,int plane) const {return fProngChannelLoc[prong][plane];}
  float DecayVertex::GetProngTickLoc(int prong,int plane) const {return fProngTickLoc[prong][plane];}
  float DecayVertex::GetProngLength(int prong) const {return fProngLength[prong];}
  float DecayVertex::GetProngTheta(int prong) const {return fProngTheta[prong];}
  float DecayVertex::GetProngPhi(int prong) const {return fProngPhi[prong];}
  int DecayVertex::GetProngNumHits(int prong) const {return fProngNumHits[prong];}
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
  void DecayVertex::SetTotHits(std::vector<recob::Hit const*> totHitsInMaxRadius) {fTotHitsInMaxRadius = totHitsInMaxRadius; return;}

  void DecayVertex::SetDetectorCoordinates(
    const std::vector<double>& minTpcBound,
    const std::vector<double>& maxTpcBound,
    geo::GeometryCore const* geometry,
    detinfo::DetectorProperties const* detectorProperties)
  {
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
