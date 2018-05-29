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
    // Coordinates
    phys_nuPosX.clear();
    phys_nuPosY.clear();
    phys_nuPosZ.clear();
    phys_prongPosX.clear();
    phys_prongPosY.clear();
    phys_prongPosZ.clear();
    phys_prongStartPosX.clear();
    phys_prongStartPosY.clear();
    phys_prongStartPosZ.clear();
    phys_prongEndPosX.clear();
    phys_prongEndPosY.clear();
    phys_prongEndPosZ.clear();
    // Direction
    phys_prongDirX.clear();
    phys_prongDirY.clear();
    phys_prongDirZ.clear();
    phys_prongTheta.clear();
    phys_prongPhi.clear();
    // Prong momentum (by range, assuming both muons)
    phys_prongMomMag_ByRange_AssMuon.clear();
    phys_prongEnergy_ByRange_AssMuon.clear();
    phys_prongMom_ByRange_AssMuonX.clear();
    phys_prongMom_ByRange_AssMuonY.clear();
    phys_prongMom_ByRange_AssMuonZ.clear();
    // Tot momentum (by range, assuming both muons)
    phys_totMomMag_ByRange_AssMuon.clear();
    phys_totEnergy_ByRange_AssMuon.clear();
    phys_invariantMass_ByRange_AssMuon.clear();
    phys_totMom_ByRange_AssMuonX.clear();
    phys_totMom_ByRange_AssMuonY.clear();
    phys_totMom_ByRange_AssMuonZ.clear();
    // Tot momentum direction (by range, assuming both muons)
    phys_totDir_ByRange_AssMuonX.clear();
    phys_totDir_ByRange_AssMuonY.clear();
    phys_totDir_ByRange_AssMuonZ.clear();
    phys_totTheta_ByRange_AssMuon.clear();
    phys_totPhi_ByRange_AssMuon.clear();
    // Others
    phys_prongLength.clear();
    phys_prongStartToNeutrinoDistance.clear();
    phys_prongNumHits.clear();
    phys_openingAngle.clear();
    // Pandora statusnostic
    status_nuWithMissingAssociatedVertex = -999;
    status_nuWithMissingAssociatedTrack = -999;
    status_nuProngWithMissingAssociatedHits = -999;
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
      AuxVertex::DecayVertex cv = decayVertices[i];
      // Coordinates
      phys_nuPosX.push_back(cv.GetX());
      phys_nuPosY.push_back(cv.GetY());
      phys_nuPosZ.push_back(cv.GetZ());
      phys_prongPosX.push_back({cv.GetProngX(0),cv.GetProngX(1)});
      phys_prongPosY.push_back({cv.GetProngY(0),cv.GetProngY(1)});
      phys_prongPosZ.push_back({cv.GetProngZ(0),cv.GetProngZ(1)});
      phys_prongStartPosX.push_back({cv.GetProngStartX(0),cv.GetProngStartX(1)});
      phys_prongStartPosY.push_back({cv.GetProngStartY(0),cv.GetProngStartY(1)});
      phys_prongStartPosZ.push_back({cv.GetProngStartZ(0),cv.GetProngStartZ(1)});
      phys_prongEndPosX.push_back({cv.GetProngEndX(0),cv.GetProngEndX(1)});
      phys_prongEndPosY.push_back({cv.GetProngEndY(0),cv.GetProngEndY(1)});
      phys_prongEndPosZ.push_back({cv.GetProngEndZ(0),cv.GetProngEndZ(1)});
      // Direction
      phys_prongDirX.push_back({cv.GetProngDirX(0),cv.GetProngDirX(1)});
      phys_prongDirY.push_back({cv.GetProngDirY(0),cv.GetProngDirY(1)});
      phys_prongDirZ.push_back({cv.GetProngDirZ(0),cv.GetProngDirZ(1)});
      phys_prongTheta.push_back({cv.GetProngTheta(0),cv.GetProngTheta(1)});
      phys_prongPhi.push_back({cv.GetProngPhi(0),cv.GetProngPhi(1)});
      // Prong momentum (by range, assuming both muons)
      phys_prongMomMag_ByRange_AssMuon.push_back({cv.GetProngMomMag_ByRange_AssMuon(0),cv.GetProngMomMag_ByRange_AssMuon(1)});
      phys_prongEnergy_ByRange_AssMuon.push_back({cv.GetProngEnergy_ByRange_AssMuon(0),cv.GetProngEnergy_ByRange_AssMuon(1)});
      phys_prongMom_ByRange_AssMuonX.push_back({cv.GetProngMom_ByRange_AssMuonX(0),cv.GetProngMom_ByRange_AssMuonX(1)});
      phys_prongMom_ByRange_AssMuonY.push_back({cv.GetProngMom_ByRange_AssMuonY(0),cv.GetProngMom_ByRange_AssMuonY(1)});
      phys_prongMom_ByRange_AssMuonZ.push_back({cv.GetProngMom_ByRange_AssMuonZ(0),cv.GetProngMom_ByRange_AssMuonZ(1)});
      // Tot momentum (by direction, assuming both muons)
      phys_totMomMag_ByRange_AssMuon.push_back(cv.GetTotMomMag_ByRange_AssMuon());
      phys_totEnergy_ByRange_AssMuon.push_back(cv.GetTotEnergy_ByRange_AssMuon());
      phys_invariantMass_ByRange_AssMuon.push_back(cv.GetInvMass_ByRange_AssMuon());
      phys_totMom_ByRange_AssMuonX.push_back(cv.GetTotMom_ByRange_AssMuonX());
      phys_totMom_ByRange_AssMuonY.push_back(cv.GetTotMom_ByRange_AssMuonY());
      phys_totMom_ByRange_AssMuonZ.push_back(cv.GetTotMom_ByRange_AssMuonZ());
      // Tot momentum direction (by direction, assuming both muons)
      phys_totDir_ByRange_AssMuonX.push_back(cv.GetTotDir_ByRange_AssMuonX());
      phys_totDir_ByRange_AssMuonY.push_back(cv.GetTotDir_ByRange_AssMuonY());
      phys_totDir_ByRange_AssMuonZ.push_back(cv.GetTotDir_ByRange_AssMuonZ());
      phys_totTheta_ByRange_AssMuon.push_back(cv.GetTotTheta_ByRange_AssMuon());
      phys_totPhi_ByRange_AssMuon.push_back(cv.GetTotPhi_ByRange_AssMuon());


      // Others
      phys_prongLength.push_back({cv.GetProngLength(0),cv.GetProngLength(1)});
      phys_prongNumHits.push_back({cv.GetProngNumHits(0),cv.GetProngNumHits(1)});
      phys_openingAngle.push_back(cv.GetOpeningAngle());
      phys_prongStartToNeutrinoDistance.push_back({cv.GetProngStartToNeutrinoDistance(0),cv.GetProngStartToNeutrinoDistance(1)});
    }
    return;
  } // END function ExtractVertexPhysics

} // END namespace EventDescriptor 
