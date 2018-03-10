/******************************************************************************
 * @file DrawTreeDescriptor.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DrawTreeDescriptor.h
 * ****************************************************************************/

// Decay vertex header
#include "DrawTreeDescriptor.h"

namespace AuxEvent
{
  DrawTreeDescriptor::DrawTreeDescriptor()
  {}
  DrawTreeDescriptor::~DrawTreeDescriptor()
  {}

  void DrawTreeDescriptor::Initialize()
  {
    dv_xyzCoordinates.clear();
    dv_wireCoordinates.clear();
    dv_tickCoordinates.clear();
    par1_xyzCoordinates.clear();
    par1_wireCoordinates.clear();
    par1_tickCoordinates.clear();
    par2_xyzCoordinates.clear();
    par2_wireCoordinates.clear();
    par2_tickCoordinates.clear();
    par1_hits_p0_wireCoordinates.clear();
    par1_hits_p1_wireCoordinates.clear();
    par1_hits_p2_wireCoordinates.clear();
    par1_hits_p0_tickCoordinates.clear();
    par1_hits_p1_tickCoordinates.clear();
    par1_hits_p2_tickCoordinates.clear();
    par2_hits_p0_wireCoordinates.clear();
    par2_hits_p1_wireCoordinates.clear();
    par2_hits_p2_wireCoordinates.clear();
    par2_hits_p0_tickCoordinates.clear();
    par2_hits_p1_tickCoordinates.clear();
    par2_hits_p2_tickCoordinates.clear();
    tot_hits_p0_wireCoordinates.clear();
    tot_hits_p1_wireCoordinates.clear();
    tot_hits_p2_wireCoordinates.clear();
    tot_hits_p0_tickCoordinates.clear();
    tot_hits_p1_tickCoordinates.clear();
    tot_hits_p2_tickCoordinates.clear();
  }

  void DrawTreeDescriptor::FillDrawTreeVariables(const std::vector<AuxVertex::DecayVertex>& cleanVertices, const std::vector<std::vector<recob::Hit const*>>& totHitsInMaxRadius, const std::vector<std::vector<recob::Hit const*>>& trackHits, const std::vector<std::vector<recob::Hit const*>>& showerHits)
  {
    for (std::vector<int>::size_type i=0; i!=cleanVertices.size(); i++)
    {
      // Get clean vertex
      auto dv = cleanVertices[i];
      int parIdx1 = dv.GetParIdx1();
      int parIdx2 = dv.GetParIdx2();
      std::vector<recob::Hit const*> par1_hits = trackHits[parIdx1];
      std::vector<recob::Hit const*> par2_hits = trackHits[parIdx2];
      std::vector<recob::Hit const*> thisTot_hits = totHitsInMaxRadius[i];
      std::vector<float> thisPar1_hits_p0_tickCoordinates,
        thisPar1_hits_p1_tickCoordinates,
        thisPar1_hits_p2_tickCoordinates,
        thisPar2_hits_p0_tickCoordinates,
        thisPar2_hits_p1_tickCoordinates,
        thisPar2_hits_p2_tickCoordinates,
        thisTot_hits_p0_tickCoordinates,
        thisTot_hits_p1_tickCoordinates,
        thisTot_hits_p2_tickCoordinates;
      std::vector<int> thisPar1_hits_p0_wireCoordinates,
        thisPar1_hits_p1_wireCoordinates,
        thisPar1_hits_p2_wireCoordinates,
        thisPar2_hits_p0_wireCoordinates,
        thisPar2_hits_p1_wireCoordinates,
        thisPar2_hits_p2_wireCoordinates,
        thisTot_hits_p0_wireCoordinates,
        thisTot_hits_p1_wireCoordinates,
        thisTot_hits_p2_wireCoordinates;

      for (auto hit : par1_hits)
      {
        if (hit->View() == 0) {
          thisPar1_hits_p0_wireCoordinates.push_back(hit->Channel());
          thisPar1_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 1) {
          thisPar1_hits_p1_wireCoordinates.push_back(hit->Channel());
          thisPar1_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 2) {
          thisPar1_hits_p2_wireCoordinates.push_back(hit->Channel());
          thisPar1_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
      }

      for (auto hit : par2_hits)
      {
        if (hit->View() == 0) {
          thisPar2_hits_p0_wireCoordinates.push_back(hit->Channel());
          thisPar2_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 1) {
          thisPar2_hits_p1_wireCoordinates.push_back(hit->Channel());
          thisPar2_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 2) {
          thisPar2_hits_p2_wireCoordinates.push_back(hit->Channel());
          thisPar2_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
      }

      for (auto hit : thisTot_hits)
      {
        if (hit->View() == 0) {
          thisTot_hits_p0_wireCoordinates.push_back(hit->Channel());
          thisTot_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 1) {
          thisTot_hits_p1_wireCoordinates.push_back(hit->Channel());
          thisTot_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 2) {
          thisTot_hits_p2_wireCoordinates.push_back(hit->Channel());
          thisTot_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
      }

      // Get coordinates
      dv_xyzCoordinates.push_back({dv.GetX(),dv.GetY(), (float)dv.GetZ()});
      dv_wireCoordinates.push_back({dv.GetChannelLoc(0),dv.GetChannelLoc(1),dv.GetChannelLoc(2)});
      dv_tickCoordinates.push_back({dv.GetTickLoc(0),dv.GetTickLoc(1),dv.GetTickLoc(2)});
      par1_xyzCoordinates.push_back({dv.GetParX(0),dv.GetParY(0),dv.GetParZ(0)});
      par1_wireCoordinates.push_back({dv.GetParChannelLoc(0,0),dv.GetParChannelLoc(0,1),dv.GetParChannelLoc(0,2)});
      par1_tickCoordinates.push_back({dv.GetParTickLoc(0,0),dv.GetParTickLoc(0,1),dv.GetParTickLoc(0,2)});
      par2_xyzCoordinates.push_back({dv.GetParX(1),dv.GetParY(1),dv.GetParZ(1)});
      par2_wireCoordinates.push_back({dv.GetParChannelLoc(1,0),dv.GetParChannelLoc(1,1),dv.GetParChannelLoc(1,2)});
      par2_tickCoordinates.push_back({dv.GetParTickLoc(1,0),dv.GetParTickLoc(1,1),dv.GetParTickLoc(1,2)});

      par1_hits_p0_wireCoordinates.push_back(thisPar1_hits_p0_wireCoordinates);
      par1_hits_p1_wireCoordinates.push_back(thisPar1_hits_p1_wireCoordinates);
      par1_hits_p2_wireCoordinates.push_back(thisPar1_hits_p2_wireCoordinates);
      par1_hits_p0_tickCoordinates.push_back(thisPar1_hits_p0_tickCoordinates);
      par1_hits_p1_tickCoordinates.push_back(thisPar1_hits_p1_tickCoordinates);
      par1_hits_p2_tickCoordinates.push_back(thisPar1_hits_p2_tickCoordinates);

      par2_hits_p0_wireCoordinates.push_back(thisPar2_hits_p0_wireCoordinates);
      par2_hits_p1_wireCoordinates.push_back(thisPar2_hits_p1_wireCoordinates);
      par2_hits_p2_wireCoordinates.push_back(thisPar2_hits_p2_wireCoordinates);
      par2_hits_p0_tickCoordinates.push_back(thisPar2_hits_p0_tickCoordinates);
      par2_hits_p1_tickCoordinates.push_back(thisPar2_hits_p1_tickCoordinates);
      par2_hits_p2_tickCoordinates.push_back(thisPar2_hits_p2_tickCoordinates);

      tot_hits_p0_wireCoordinates.push_back(thisTot_hits_p0_wireCoordinates);
      tot_hits_p1_wireCoordinates.push_back(thisTot_hits_p1_wireCoordinates);
      tot_hits_p2_wireCoordinates.push_back(thisTot_hits_p2_wireCoordinates);
      tot_hits_p0_tickCoordinates.push_back(thisTot_hits_p0_tickCoordinates);
      tot_hits_p1_tickCoordinates.push_back(thisTot_hits_p1_tickCoordinates);
      tot_hits_p2_tickCoordinates.push_back(thisTot_hits_p2_tickCoordinates);
    }
    return;
  } // END function FillDrawTree
} // END namespace DrawTreeDescriptor 
