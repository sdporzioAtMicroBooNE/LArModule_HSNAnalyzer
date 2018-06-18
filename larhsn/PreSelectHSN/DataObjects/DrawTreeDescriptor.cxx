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
      dv_xyzCoordinates.push_back({currentVertex.fX,currentVertex.fY,currentVertex.fZ});
      dv_wireCoordinates.push_back({currentVertex.fChannelLoc[0],currentVertex.fChannelLoc[1],currentVertex.fChannelLoc[2]});
      dv_tickCoordinates.push_back({currentVertex.fTickLoc[0],currentVertex.fTickLoc[1],currentVertex.fTickLoc[2]});
      prong1_xyzCoordinates.push_back({currentVertex.fProngX[0],currentVertex.fProngY[0],currentVertex.fProngZ[0]});
      prong1_wireCoordinates.push_back({currentVertex.fProngChannelLoc[0][0],currentVertex.fProngChannelLoc[0][1],currentVertex.fProngChannelLoc[0][2]});
      prong1_tickCoordinates.push_back({currentVertex.fProngTickLoc[0][0],currentVertex.fProngTickLoc[0][1],currentVertex.fProngTickLoc[0][2]});
      prong2_xyzCoordinates.push_back({currentVertex.fProngX[1],currentVertex.fProngY[1],currentVertex.fProngZ[1]});
      prong2_wireCoordinates.push_back({currentVertex.fProngChannelLoc[1][0],currentVertex.fProngChannelLoc[1][1],currentVertex.fProngChannelLoc[1][2]});
      prong2_tickCoordinates.push_back({currentVertex.fProngTickLoc[1][0],currentVertex.fProngTickLoc[1][1],currentVertex.fProngTickLoc[1][2]});

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
