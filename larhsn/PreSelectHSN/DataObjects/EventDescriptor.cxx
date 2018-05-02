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
    phys_nuStartPosition.clear();
    phys_prong1StartPosition.clear();
    phys_prong2StartPosition.clear();
    phys_prong1MomentumDir.clear();
    phys_prong2MomentumDir.clear();
    phys_prong1MomentumMag.clear();
    phys_prong2MomentumMag.clear();
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
    match_pdgCode.clear();
    match_mass.clear();
    match_energy.clear();
    match_prong1StartPosition.clear();
    match_prong2StartPosition.clear();
    match_prong1Momentum.clear();
    match_prong2Momentum.clear();

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
      phys_nuStartPosition.push_back(
      {
        currentVertex.GetX(),
        currentVertex.GetY(),
        currentVertex.GetZ()
      });
      phys_prong1StartPosition.push_back(
      {
        currentVertex.GetProngX(0),
        currentVertex.GetProngY(0), 
        currentVertex.GetProngZ(0)
      });
      phys_prong2StartPosition.push_back(
      {
        currentVertex.GetProngX(1),
        currentVertex.GetProngY(1), 
        currentVertex.GetProngZ(1)
      });
      phys_prong1MomentumDir.push_back(
      {
        currentVertex.GetProngDirPx(0),
        currentVertex.GetProngDirPy(0), 
        currentVertex.GetProngDirPz(0)
      });
      phys_prong2MomentumDir.push_back(
      {
        currentVertex.GetProngDirPx(1),
        currentVertex.GetProngDirPy(1), 
        currentVertex.GetProngDirPz(1)
      });
      phys_prong1MomentumMag.push_back(currentVertex.GetProngMagP(0));
      phys_prong2MomentumMag.push_back(currentVertex.GetProngMagP(1));
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
