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

  void DrawTreeDescriptor::Initialize(int i_run, int i_subrun, int i_event)
  {
    run = i_run;
    subrun = i_subrun;
    event = i_event;
    
    dv_xyzCoordinates.clear();
    dv_wireCoordinates.clear();
    dv_tickCoordinates.clear();
    prong1_xyzCoordinates.clear();
    prong1_wireCoordinates.clear();
    prong1_tickCoordinates.clear();
    prong2_xyzCoordinates.clear();
    prong2_wireCoordinates.clear();
    prong2_tickCoordinates.clear();
    prong1_hits_p0_wireCoordinates.clear();
    prong1_hits_p1_wireCoordinates.clear();
    prong1_hits_p2_wireCoordinates.clear();
    prong1_hits_p0_tickCoordinates.clear();
    prong1_hits_p1_tickCoordinates.clear();
    prong1_hits_p2_tickCoordinates.clear();
    prong2_hits_p0_wireCoordinates.clear();
    prong2_hits_p1_wireCoordinates.clear();
    prong2_hits_p2_wireCoordinates.clear();
    prong2_hits_p0_tickCoordinates.clear();
    prong2_hits_p1_tickCoordinates.clear();
    prong2_hits_p2_tickCoordinates.clear();
    tot_hits_p0_wireCoordinates.clear();
    tot_hits_p1_wireCoordinates.clear();
    tot_hits_p2_wireCoordinates.clear();
    tot_hits_p0_tickCoordinates.clear();
    tot_hits_p1_tickCoordinates.clear();
    tot_hits_p2_tickCoordinates.clear();
  }

  void DrawTreeDescriptor::FillDrawTreeVariables(
          const std::vector<AuxVertex::DecayVertex>& decayVertices)
  {
    for (std::vector<int>::size_type i=0; i!=decayVertices.size(); i++)
    {
      // Get decay vertex
      AuxVertex::DecayVertex currentVertex = decayVertices[i];

      std::vector<art::Ptr<recob::Hit>> prong1_hits = currentVertex.GetProngHits(0);
      std::vector<art::Ptr<recob::Hit>> prong2_hits = currentVertex.GetProngHits(1);
      std::vector<art::Ptr<recob::Hit>> thisTot_hits = currentVertex.GetTotHits();

      std::vector<float> thisProng1_hits_p0_tickCoordinates,
        thisProng1_hits_p1_tickCoordinates,
        thisProng1_hits_p2_tickCoordinates,
        thisProng2_hits_p0_tickCoordinates,
        thisProng2_hits_p1_tickCoordinates,
        thisProng2_hits_p2_tickCoordinates,
        thisTot_hits_p0_tickCoordinates,
        thisTot_hits_p1_tickCoordinates,
        thisTot_hits_p2_tickCoordinates;
      std::vector<int> thisProng1_hits_p0_wireCoordinates,
        thisProng1_hits_p1_wireCoordinates,
        thisProng1_hits_p2_wireCoordinates,
        thisProng2_hits_p0_wireCoordinates,
        thisProng2_hits_p1_wireCoordinates,
        thisProng2_hits_p2_wireCoordinates,
        thisTot_hits_p0_wireCoordinates,
        thisTot_hits_p1_wireCoordinates,
        thisTot_hits_p2_wireCoordinates;

      for (auto hit : prong1_hits)
      {
        if (hit->View() == 0) {
          thisProng1_hits_p0_wireCoordinates.push_back(hit->Channel());
          thisProng1_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 1) {
          thisProng1_hits_p1_wireCoordinates.push_back(hit->Channel());
          thisProng1_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 2) {
          thisProng1_hits_p2_wireCoordinates.push_back(hit->Channel());
          thisProng1_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
      }

      for (auto hit : prong2_hits)
      {
        if (hit->View() == 0) {
          thisProng2_hits_p0_wireCoordinates.push_back(hit->Channel());
          thisProng2_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 1) {
          thisProng2_hits_p1_wireCoordinates.push_back(hit->Channel());
          thisProng2_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
        }
        if (hit->View() == 2) {
          thisProng2_hits_p2_wireCoordinates.push_back(hit->Channel());
          thisProng2_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
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
      dv_xyzCoordinates.push_back({currentVertex.GetX(),currentVertex.GetY(),currentVertex.GetZ()});
      dv_wireCoordinates.push_back({currentVertex.GetChannelLoc(0),currentVertex.GetChannelLoc(1),currentVertex.GetChannelLoc(2)});
      dv_tickCoordinates.push_back({currentVertex.GetTickLoc(0),currentVertex.GetTickLoc(1),currentVertex.GetTickLoc(2)});
      prong1_xyzCoordinates.push_back({currentVertex.GetProngX(0),currentVertex.GetProngY(0),currentVertex.GetProngZ(0)});
      prong1_wireCoordinates.push_back({currentVertex.GetProngChannelLoc(0,0),currentVertex.GetProngChannelLoc(0,1),currentVertex.GetProngChannelLoc(0,2)});
      prong1_tickCoordinates.push_back({currentVertex.GetProngTickLoc(0,0),currentVertex.GetProngTickLoc(0,1),currentVertex.GetProngTickLoc(0,2)});
      prong2_xyzCoordinates.push_back({currentVertex.GetProngX(1),currentVertex.GetProngY(1),currentVertex.GetProngZ(1)});
      prong2_wireCoordinates.push_back({currentVertex.GetProngChannelLoc(1,0),currentVertex.GetProngChannelLoc(1,1),currentVertex.GetProngChannelLoc(1,2)});
      prong2_tickCoordinates.push_back({currentVertex.GetProngTickLoc(1,0),currentVertex.GetProngTickLoc(1,1),currentVertex.GetProngTickLoc(1,2)});

      prong1_hits_p0_wireCoordinates.push_back(thisProng1_hits_p0_wireCoordinates);
      prong1_hits_p1_wireCoordinates.push_back(thisProng1_hits_p1_wireCoordinates);
      prong1_hits_p2_wireCoordinates.push_back(thisProng1_hits_p2_wireCoordinates);
      prong1_hits_p0_tickCoordinates.push_back(thisProng1_hits_p0_tickCoordinates);
      prong1_hits_p1_tickCoordinates.push_back(thisProng1_hits_p1_tickCoordinates);
      prong1_hits_p2_tickCoordinates.push_back(thisProng1_hits_p2_tickCoordinates);

      prong2_hits_p0_wireCoordinates.push_back(thisProng2_hits_p0_wireCoordinates);
      prong2_hits_p1_wireCoordinates.push_back(thisProng2_hits_p1_wireCoordinates);
      prong2_hits_p2_wireCoordinates.push_back(thisProng2_hits_p2_wireCoordinates);
      prong2_hits_p0_tickCoordinates.push_back(thisProng2_hits_p0_tickCoordinates);
      prong2_hits_p1_tickCoordinates.push_back(thisProng2_hits_p1_tickCoordinates);
      prong2_hits_p2_tickCoordinates.push_back(thisProng2_hits_p2_tickCoordinates);

      tot_hits_p0_wireCoordinates.push_back(thisTot_hits_p0_wireCoordinates);
      tot_hits_p1_wireCoordinates.push_back(thisTot_hits_p1_wireCoordinates);
      tot_hits_p2_wireCoordinates.push_back(thisTot_hits_p2_wireCoordinates);
      tot_hits_p0_tickCoordinates.push_back(thisTot_hits_p0_tickCoordinates);
      tot_hits_p1_tickCoordinates.push_back(thisTot_hits_p1_tickCoordinates);
      tot_hits_p2_tickCoordinates.push_back(thisTot_hits_p2_tickCoordinates);
    } // END loop for each decay vertex
    return;
  } // END function FillDrawTree
} // END namespace DrawTreeDescriptor 
