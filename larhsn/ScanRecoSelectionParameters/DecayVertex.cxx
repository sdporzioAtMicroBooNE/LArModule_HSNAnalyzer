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
  DecayVertex::DecayVertex(double x,
                           double y,
                           double z,
                           int parIdx1,
                           int parIdx2,
                           std::string parType1,
                           std::string parType2)
  {
    fX = x;
    fY = y;
    fZ = z;
    fParIdx1 = parIdx1;
    fParIdx2 = parIdx2;
    fParType1 = parType1;
    fParType2 = parType2;
  }

  DecayVertex::~DecayVertex()
  {}

  double DecayVertex::GetX() {return fX;}
  double DecayVertex::GetY() {return fY;}
  double DecayVertex::GetZ() {return fZ;}
  int DecayVertex::GetParIdx1() {return fParIdx1;}
  int DecayVertex::GetParIdx2() {return fParIdx2;}
  std::string DecayVertex::GetParType1() {return fParType1;}
  std::string DecayVertex::GetParType2() {return fParType2;}
  bool DecayVertex::IsInsideTPC() {return fIsInsideTPC;}
  int DecayVertex::GetWireLoc(int plane) {return fWireLoc[plane];}
  double DecayVertex::GetTickLoc(int plane) {return fTickLoc[plane];}

  double Distance(DecayVertex v1, DecayVertex v2)
  {
    // Calculate distance between vertices
    double dist = sqrt(
      pow((v1.GetX() - v2.GetX()),2.) +
      pow((v1.GetY() - v2.GetY()),2.) +
      pow((v1.GetZ() - v2.GetZ()),2.)
    );
    return dist;
  } // END function Distance

  DecayVertex MeanVertex(DecayVertex v1, DecayVertex v2)
  {
    // Find vertex half-way between two ORIGIN vertices
    if (v1.GetParIdx1()!=v1.GetParIdx2() || v2.GetParIdx1()!=v2.GetParIdx2()){
      throw std::runtime_error("Mean vertex is supposed to be calculated only from two origin vertices (vertices which correspond to the actual start or end of a track or shower). Calculating the mean vertex between two mean vertices is not implemented yet.");
    }
    double x = (v1.GetX() + v2.GetX())/2.;
    double y = (v1.GetY() + v2.GetY())/2.;
    double z = (v1.GetZ() + v2.GetZ())/2.;
    int pidx1 = v1.GetParIdx1();
    int pidx2 = v2.GetParIdx1();
    std::string ptype1 = v1.GetParType1();
    std::string ptype2 = v2.GetParType1();
    DecayVertex meanVertex(x,y,z,pidx1,pidx2,ptype1,ptype2);
    return meanVertex;
  } // END function MeanVertex

} // END namespace AuxVertex




// // Geometry functions
// void CoordinatesToChannels(geo::GeometryCore const* fGeometry, double xyz[], int channels[])
// {
//   // Find nearest wires to provided coordinates
//   raw::ChannelID_t channels_info_0 = fGeometry->NearestChannel(xyz,0);
//   raw::ChannelID_t channels_info_1 = fGeometry->NearestChannel(xyz,1);
//   raw::ChannelID_t channels_info_2 = fGeometry->NearestChannel(xyz,2);
//   channels[0] = channels_info_0;
//   channels[1] = channels_info_1;
//   channels[2] = channels_info_2;
//   return;
// } // END function CoordinatesToChannels

// void XToTicks(detinfo::DetectorProperties const* fDetectorProperties, double x, double ticks[])
// {
//   // Find nearest ticks to provided coordinates
//   ticks[0] = fDetectorProperties->ConvertXToTicks(x, 0, 0, 0);
//   ticks[1] = fDetectorProperties->ConvertXToTicks(x, 1, 0, 0);
//   ticks[2] = fDetectorProperties->ConvertXToTicks(x, 2, 0, 0);
//   return;
// } // END function XToTicks

// bool IsInsideTPC(double xyz[], double fDistanceCut, std::vector<double> fMinTpcBound, std::vector<double> fMaxTpcBound)
// {
//   // Check whether coordinates are inside TPC 
//   bool isInsideX = (xyz[0]>fMinTpcBound[0]+fDistanceCut &&
//     xyz[0]<fMaxTpcBound[0]-fDistanceCut);
//   bool isInsideY = (xyz[1]>fMinTpcBound[1]+fDistanceCut &&
//     xyz[1]<fMaxTpcBound[1]-fDistanceCut);
//   bool isInsideZ = (xyz[2]>fMinTpcBound[2]+fDistanceCut &&
//     xyz[2]<fMaxTpcBound[2]-fDistanceCut);
//   if (isInsideX && isInsideY && isInsideZ) return true;
//   else return false;
// } // END function IsInsideTPC