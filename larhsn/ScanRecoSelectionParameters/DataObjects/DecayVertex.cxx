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

  DecayVertex::DecayVertex(float x,
                           float y,
                           float z,
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
    fParX = {x,x};
    fParY = {y,y};
    fParZ = {z,z};
    fParIdx1 = parIdx1;
    fParIdx2 = parIdx2;
    fParType1 = parType1;
    fParType2 = parType2;
    fDirection1 = direction1;
    fDirection2 = direction2;
    fChannelLoc = {-1,-1,-1};
    fTickLoc = {-1.,-1.,-1.};
    fParChannelLoc = {{-1,-1,-1}, {-1,-1,-1}};
    fParTickLoc = {{-1.,-1.,-1.}, {-1.,-1.,-1.}};
  }

  // Getters
  float DecayVertex::GetX() const {return fX;}
  float DecayVertex::GetY() const {return fY;}
  float DecayVertex::GetZ() const {return fZ;}
  float DecayVertex::GetParX(int par) const {return fParX[par];}
  float DecayVertex::GetParY(int par) const {return fParY[par];}
  float DecayVertex::GetParZ(int par) const {return fParZ[par];}
  int DecayVertex::GetParIdx1() const {return fParIdx1;}
  int DecayVertex::GetParIdx2() const {return fParIdx2;}
  std::string DecayVertex::GetParType1() const {return fParType1;}
  std::string DecayVertex::GetParType2() const {return fParType2;}
  std::string DecayVertex::GetDirection1() const {return fDirection1;}
  std::string DecayVertex::GetDirection2() const {return fDirection2;}
  bool DecayVertex::IsInsideTPC() const {return fIsInsideTPC;}
  bool DecayVertex::IsDetLocAssigned() const {return fIsDetLocAssigned;}
  int DecayVertex::GetChannelLoc(int plane) const {return fChannelLoc[plane];}
  float DecayVertex::GetTickLoc(int plane) const {return fTickLoc[plane];}
  int DecayVertex::GetParChannelLoc(int par,int plane) const {return fParChannelLoc[par][plane];}
  float DecayVertex::GetParTickLoc(int par,int plane) const {return fParTickLoc[par][plane];}


  // Setters
  void DecayVertex::SetChannelLoc(int channel0, int channel1, int channel2) {fChannelLoc = {channel0,channel1,channel2}; return;}
  void DecayVertex::SetTickLoc(float tick0, float tick1, float tick2) {fTickLoc = {tick0, tick1, tick2}; return;}
  void DecayVertex::SetParChannelLoc(int par, int channel0, int channel1, int channel2) {fParChannelLoc[par] =  {channel0,channel1,channel2}; return;}
  void DecayVertex::SetParTickLoc(int par, float tick0, float tick1, float tick2) {fParTickLoc[par] = {tick0, tick1, tick2}; return;}
  void DecayVertex::SetParXYZ(int par, float x, float y, float z) {fParX[par] = x; fParY[par] = y; fParZ[par] = z; return;}
  void DecayVertex::SetIsInsideTPC(bool val) {fIsInsideTPC = val; return;}
  void DecayVertex::SetIsDetLocAssigned(bool val) {fIsDetLocAssigned = val; return;}

  // Printers
  void DecayVertex::PrintInformation() const{
    int fStartWire[3] = {0,2399,4798};
    printf("\n-|Vertex information|\n");
    printf("|_Parents type: [%s,%s]\n", fParType1.c_str(), fParType2.c_str());
    printf("|_Parents direction: [%s,%s]\n", fDirection1.c_str(), fDirection2.c_str());
    printf("|_Parents index: [%i,%i]\n", fParIdx1, fParIdx2);
    printf("|_Vertex inside TPC: %d\n", fIsInsideTPC);
    printf("|_Spatial Coordinates: [%.1f, %.1f, %.1f]\n",fX,fY,fZ);
    printf("|_Detector location assigned: %d\n", fIsDetLocAssigned);
    if (fIsDetLocAssigned)
    {
      printf("|_Channel coordinates: [%i, %i, %i]\n",fChannelLoc[0]-fStartWire[0],fChannelLoc[1]-fStartWire[1],fChannelLoc[2]-fStartWire[2]);
      printf("|_Ticks coordinates: [%.1f, %.1f, %.1f]\n",fTickLoc[0],fTickLoc[1],fTickLoc[2]);
    }
    return;
  }


  // Aux Functions
  float Distance(DecayVertex v1, DecayVertex v2)
  {
    // Calculate distance between vertices
    float dist = sqrt(
      pow((v1.GetX() - v2.GetX()),2.) +
      pow((v1.GetY() - v2.GetY()),2.) +
      pow((v1.GetZ() - v2.GetZ()),2.)
    );
    return dist;
  } // END function Distance

  DecayVertex MeanVertex(DecayVertex v1, DecayVertex v2)
  {
    // Find vertex half-way between two ORIGIN vertices
    // if (v1.GetParIdx1()!=v1.GetParIdx2() || v2.GetParIdx1()!=v2.GetParIdx2())
    // {
    //   throw std::runtime_error("Mean vertex is supposed to be calculated only from two origin vertices (vertices which correspond to the actual start or end of a track or shower). Calculating the mean vertex between two mean vertices is not implemented yet.");
    // }
    float x = (v1.GetX() + v2.GetX())/2.;
    float y = (v1.GetY() + v2.GetY())/2.;
    float z = (v1.GetZ() + v2.GetZ())/2.;
    int pidx1 = v1.GetParIdx1();
    int pidx2 = v2.GetParIdx1();
    std::string ptype1 = v1.GetParType1();
    std::string ptype2 = v2.GetParType1();
    std::string direction1 = v1.GetDirection1();
    std::string direction2 = v2.GetDirection1();
    DecayVertex meanVertex(x,y,z,pidx1,pidx2,ptype1,ptype2,direction1,direction2);
    meanVertex.SetParXYZ(0,v1.GetX(),v1.GetY(),v1.GetZ());
    meanVertex.SetParXYZ(1,v2.GetX(),v2.GetY(),v2.GetZ());
    return meanVertex;
  } // END function MeanVertex

} // END namespace AuxVertex
