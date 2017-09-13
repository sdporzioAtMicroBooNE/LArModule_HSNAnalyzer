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

  DecayVertex::DecayVertex(double x,
                           double y,
                           double z,
                           int parIdx1,
                           int parIdx2,
                           std::string parType1,
                           std::string parType2,
                           std::string direction1,
                           std::string direction2)
  {
    fIsDetLocAssigned = false;
    fX = x;
    fY = y;
    fZ = z;
    fParIdx1 = parIdx1;
    fParIdx2 = parIdx2;
    fParType1 = parType1;
    fParType2 = parType2;
    fDirection1 = direction1;
    fDirection2 = direction2;
    fChannelLoc = {-1,-1,-1};
    fTickLoc = {-1.,-1.,-1.};
  }

  // Getters
  double DecayVertex::GetX() {return fX;}
  double DecayVertex::GetY() {return fY;}
  double DecayVertex::GetZ() {return fZ;}
  int DecayVertex::GetParIdx1() {return fParIdx1;}
  int DecayVertex::GetParIdx2() {return fParIdx2;}
  std::string DecayVertex::GetParType1() {return fParType1;}
  std::string DecayVertex::GetParType2() {return fParType2;}
  std::string DecayVertex::GetDirection1() {return fDirection1;}
  std::string DecayVertex::GetDirection2() {return fDirection2;}
  bool DecayVertex::IsInsideTPC() {return fIsInsideTPC;}
  bool DecayVertex::IsDetLocAssigned() {return fIsDetLocAssigned;}
  int DecayVertex::GetChannelLoc(int plane) {return fChannelLoc[plane];}
  double DecayVertex::GetTickLoc(int plane) {return fTickLoc[plane];}

  // Setters
  void DecayVertex::SetChannelLoc(int channel0, int channel1, int channel2) {fChannelLoc = {channel0,channel1,channel2};return;}
  void DecayVertex::SetTickLoc(double tick0, double tick1, double tick2) {fTickLoc = {tick0, tick1, tick2};return;}
  void DecayVertex::SetIsInsideTPC(bool val) {fIsInsideTPC = val;return;}
  void DecayVertex::SetIsDetLocAssigned(bool val) {fIsDetLocAssigned = val;return;}

  // Printers
  void DecayVertex::PrintInformation(){
    printf("|Vertex information|\n");
    printf("Parents type: [%s,%s]\n", fParType1.c_str(), fParType2.c_str());
    printf("Parents direction: [%s,%s]\n", fDirection1.c_str(), fDirection2.c_str());
    printf("Parents index: [%i,%i]\n", fParIdx1, fParIdx2);
    printf("Vertex inside TPC: %d\n", fIsInsideTPC);
    printf("Spatial Coordinates: [%.1f, %.1f, %.1f]\n",fX,fY,fZ);
    printf("Detector location assigned: %d\n", fIsDetLocAssigned);
    if (fIsDetLocAssigned)
    {
      printf("Channel coordinates: [%i, %i, %i]\n",fChannelLoc[0],fChannelLoc[1],fChannelLoc[2]);
      printf("Ticks coordinates: [%.1f, %.1f, %.1f]\n",fTickLoc[0],fTickLoc[1],fTickLoc[2]);
    }
    printf("\n\n");
    return;
  }


  // Aux Functions
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
    std::string direction1 = v1.GetDirection1();
    std::string direction2 = v2.GetDirection1();
    DecayVertex meanVertex(x,y,z,pidx1,pidx2,ptype1,ptype2,direction1,direction2);
    return meanVertex;
  } // END function MeanVertex

} // END namespace AuxVertex
