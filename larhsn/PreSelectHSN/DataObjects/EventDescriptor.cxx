/******************************************************************************
 * @file EventDescriptor.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  EventDescriptor.h
 * ****************************************************************************/

// Decay vertex header
#include "EventDescriptor.h"

namespace AuxEvent
{
  EventDescriptor::EventDescriptor()
  {}
  EventDescriptor::~EventDescriptor()
  {}

  void EventDescriptor::Initialize(int i_run, int i_subrun, int i_event)
  {
    run = i_run;
    subrun = i_subrun;
    event = i_event;

    // Pandora search
    nNeutrinos = -999;
    nTwoProngedNeutrinos = -999;
    nContainedTwoProngedNeutrinos = -999;
    neutrinoPdgCode.clear();
    neutrinoNumDaughters.clear();
    neutrinoNumTracks.clear();
    neutrinoNumShowers.clear();
    // Pandora calo
    calo_totChargeInRadius.clear();
    calo_prong1ChargeInRadius.clear();
    calo_prong2ChargeInRadius.clear();
    calo_caloRatio.clear();
    // Pandora physics
    phys_prongLength.clear();
    phys_prongTheta.clear();
    phys_prongPhi.clear();
    phys_prongStartToNeutrinoDistance.clear();
    phys_prongNumHits.clear();
    phys_openingAngle.clear();
    // Pandora diagnostic
    diag_nuWithMissingAssociatedVertex = -999;
    diag_nuWithMissingAssociatedTrack = -999;
    diag_nuProngWithMissingAssociatedHits = -999;
    // Truth information
    prong_matchedPDG.clear();
  } // END function Initialize

  void EventDescriptor::ExtractVertexPhysics(const std::vector<AuxVertex::DecayVertex> & decayVertices)
  {
    // Clear vectors that will be filled
    phys_prongLength.clear();
    phys_prongTheta.clear();
    phys_prongPhi.clear();
    phys_prongNumHits.clear();

    for (std::vector<int>::size_type i=0; i!=decayVertices.size(); i++)
    {
      AuxVertex::DecayVertex currentVertex = decayVertices[i];
      phys_prongLength.push_back({currentVertex.GetProngLength(0),currentVertex.GetProngLength(1)});
      phys_prongTheta.push_back({currentVertex.GetProngTheta(0),currentVertex.GetProngTheta(1)});
      phys_prongPhi.push_back({currentVertex.GetProngPhi(0),currentVertex.GetProngPhi(1)});
      phys_prongNumHits.push_back({currentVertex.GetProngNumHits(0),currentVertex.GetProngNumHits(1)});
      phys_openingAngle.push_back(currentVertex.GetOpeningAngle());
      phys_prongStartToNeutrinoDistance.push_back({currentVertex.GetProngStartToNeutrinoDistance(0),currentVertex.GetProngStartToNeutrinoDistance(1)});
    }
    return;
  } // END function ExtractVertexPhysics

} // END namespace EventDescriptor 
