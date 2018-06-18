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
    // Hypothesis info
    phys_prongPdgCode_h1.clear();
    phys_prongMass_h1.clear();
    phys_prongPdgCode_h2.clear();
    phys_prongMass_h2.clear();
    // Prong momentum (by range, assuming both muons)
    phys_prongMomMag_ByRange_h1.clear();
    phys_prongEnergy_ByRange_h1.clear();
    phys_prongMom_ByRange_h1_X.clear();
    phys_prongMom_ByRange_h1_Y.clear();
    phys_prongMom_ByRange_h1_Z.clear();
    //
    phys_prongMomMag_ByRange_h2.clear();
    phys_prongEnergy_ByRange_h2.clear();
    phys_prongMom_ByRange_h2_X.clear();
    phys_prongMom_ByRange_h2_Y.clear();
    phys_prongMom_ByRange_h2_Z.clear();
    // Tot momentum (by range, assuming both muons)
    phys_totMomMag_ByRange_h1.clear();
    phys_totEnergy_ByRange_h1.clear();
    phys_invariantMass_ByRange_h1.clear();
    phys_totMom_ByRange_h1_X.clear();
    phys_totMom_ByRange_h1_Y.clear();
    phys_totMom_ByRange_h1_Z.clear();
    //
    phys_totMomMag_ByRange_h2.clear();
    phys_totEnergy_ByRange_h2.clear();
    phys_invariantMass_ByRange_h2.clear();
    phys_totMom_ByRange_h2_X.clear();
    phys_totMom_ByRange_h2_Y.clear();
    phys_totMom_ByRange_h2_Z.clear();
    // Tot momentum direction (by range, assuming both muons)
    phys_totDir_ByRange_h1_X.clear();
    phys_totDir_ByRange_h1_Y.clear();
    phys_totDir_ByRange_h1_Z.clear();
    phys_totTheta_ByRange_h1.clear();
    phys_totPhi_ByRange_h1.clear();
    //
    phys_totDir_ByRange_h2_X.clear();
    phys_totDir_ByRange_h2_Y.clear();
    phys_totDir_ByRange_h2_Z.clear();
    phys_totTheta_ByRange_h2.clear();
    phys_totPhi_ByRange_h2.clear();
    // Prong momentum (by MCS)
    phys_prongPdgCodeHypothesis_ByMcs.clear();
    phys_prongIsBestFwd_ByMcs.clear();
    // Prong Momentum (By Mcs, forward)
    phys_prongMomMag_ByMcs_fwd_h1.clear();
    phys_prongEnergy_ByMcs_fwd_h1.clear();
    phys_prongMom_ByMcs_fwd_h1_X.clear();
    phys_prongMom_ByMcs_fwd_h1_Y.clear();
    phys_prongMom_ByMcs_fwd_h1_Z.clear();
    phys_prongMomMag_ByMcs_fwd_h2.clear();
    phys_prongEnergy_ByMcs_fwd_h2.clear();
    phys_prongMom_ByMcs_fwd_h2_X.clear();
    phys_prongMom_ByMcs_fwd_h2_Y.clear();
    phys_prongMom_ByMcs_fwd_h2_Z.clear();
    // Tot momentum (by range, assuming both muons, forward)
    phys_totMomMag_ByMcs_fwd_h1.clear();
    phys_totEnergy_ByMcs_fwd_h1.clear();
    phys_invariantMass_ByMcs_fwd_h1.clear();
    phys_totMom_ByMcs_fwd_h1_X.clear();
    phys_totMom_ByMcs_fwd_h1_Y.clear();
    phys_totMom_ByMcs_fwd_h1_Z.clear();
    phys_totMomMag_ByMcs_fwd_h2.clear();
    phys_totEnergy_ByMcs_fwd_h2.clear();
    phys_invariantMass_ByMcs_fwd_h2.clear();
    phys_totMom_ByMcs_fwd_h2_X.clear();
    phys_totMom_ByMcs_fwd_h2_Y.clear();
    phys_totMom_ByMcs_fwd_h2_Z.clear();
    // Tot momentum direction (by range, assuming both muons, forward)
    phys_totTheta_ByMcs_fwd_h1.clear();
    phys_totPhi_ByMcs_fwd_h1.clear();
    phys_totDir_ByMcs_fwd_h1_X.clear();
    phys_totDir_ByMcs_fwd_h1_Y.clear();
    phys_totDir_ByMcs_fwd_h1_Z.clear();
    phys_totTheta_ByMcs_fwd_h2.clear();
    phys_totPhi_ByMcs_fwd_h2.clear();
    phys_totDir_ByMcs_fwd_h2_X.clear();
    phys_totDir_ByMcs_fwd_h2_Y.clear();
    phys_totDir_ByMcs_fwd_h2_Z.clear();
    // Prong Momentum (By Mcs, best)
    phys_prongMomMag_ByMcs_best_h1.clear();
    phys_prongEnergy_ByMcs_best_h1.clear();
    phys_prongMom_ByMcs_best_h1_X.clear();
    phys_prongMom_ByMcs_best_h1_Y.clear();
    phys_prongMom_ByMcs_best_h1_Z.clear();
    phys_prongMomMag_ByMcs_best_h2.clear();
    phys_prongEnergy_ByMcs_best_h2.clear();
    phys_prongMom_ByMcs_best_h2_X.clear();
    phys_prongMom_ByMcs_best_h2_Y.clear();
    phys_prongMom_ByMcs_best_h2_Z.clear();
    // Tot momentum (by range, assuming both muons, best)
    phys_totMomMag_ByMcs_best_h1.clear();
    phys_totEnergy_ByMcs_best_h1.clear();
    phys_invariantMass_ByMcs_best_h1.clear();
    phys_totMom_ByMcs_best_h1_X.clear();
    phys_totMom_ByMcs_best_h1_Y.clear();
    phys_totMom_ByMcs_best_h1_Z.clear();
    phys_totMomMag_ByMcs_best_h2.clear();
    phys_totEnergy_ByMcs_best_h2.clear();
    phys_invariantMass_ByMcs_best_h2.clear();
    phys_totMom_ByMcs_best_h2_X.clear();
    phys_totMom_ByMcs_best_h2_Y.clear();
    phys_totMom_ByMcs_best_h2_Z.clear();
    // Tot momentum direction (by range, assuming both muons, best)
    phys_totTheta_ByMcs_best_h1.clear();
    phys_totPhi_ByMcs_best_h1.clear();
    phys_totDir_ByMcs_best_h1_X.clear();
    phys_totDir_ByMcs_best_h1_Y.clear();
    phys_totDir_ByMcs_best_h1_Z.clear();
    phys_totTheta_ByMcs_best_h2.clear();
    phys_totPhi_ByMcs_best_h2.clear();
    phys_totDir_ByMcs_best_h2_X.clear();
    phys_totDir_ByMcs_best_h2_Y.clear();
    phys_totDir_ByMcs_best_h2_Z.clear();
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
      phys_nuPosX.push_back(cv.fX);
      phys_nuPosY.push_back(cv.fY);
      phys_nuPosZ.push_back(cv.fZ);
      phys_prongPosX.push_back({cv.fProngX[0],cv.fProngX[1]});
      phys_prongPosY.push_back({cv.fProngY[0],cv.fProngY[1]});
      phys_prongPosZ.push_back({cv.fProngZ[0],cv.fProngZ[1]});
      phys_prongStartPosX.push_back({cv.fProngStartX[0],cv.fProngStartX[1]});
      phys_prongStartPosY.push_back({cv.fProngStartY[0],cv.fProngStartY[1]});
      phys_prongStartPosZ.push_back({cv.fProngStartZ[0],cv.fProngStartZ[1]});
      phys_prongEndPosX.push_back({cv.fProngEndX[0],cv.fProngEndX[1]});
      phys_prongEndPosY.push_back({cv.fProngEndY[0],cv.fProngEndY[1]});
      phys_prongEndPosZ.push_back({cv.fProngEndZ[0],cv.fProngEndZ[1]});
      // Direction
      phys_prongDirX.push_back({cv.fProngDirX[0],cv.fProngDirX[1]});
      phys_prongDirY.push_back({cv.fProngDirY[0],cv.fProngDirY[1]});
      phys_prongDirZ.push_back({cv.fProngDirZ[0],cv.fProngDirZ[1]});
      phys_prongTheta.push_back({cv.fProngTheta[0],cv.fProngTheta[1]});
      phys_prongPhi.push_back({cv.fProngPhi[0],cv.fProngPhi[1]});
      // Prong momentum (by range, assuming both muons)
      phys_prongPdgCode_h1.push_back({cv.fProngPdgCode_h1[0],cv.fProngPdgCode_h1[1]});
      phys_prongMass_h1.push_back({cv.fProngMass_h1[0],cv.fProngMass_h1[1]});
      phys_prongMomMag_ByRange_h1.push_back({cv.fProngMomMag_ByRange_h1[0],cv.fProngMomMag_ByRange_h1[1]});
      phys_prongEnergy_ByRange_h1.push_back({cv.fProngEnergy_ByRange_h1[0],cv.fProngEnergy_ByRange_h1[1]});
      phys_prongMom_ByRange_h1_X.push_back({cv.fProngMom_ByRange_h1_X[0],cv.fProngMom_ByRange_h1_X[1]});
      phys_prongMom_ByRange_h1_Y.push_back({cv.fProngMom_ByRange_h1_Y[0],cv.fProngMom_ByRange_h1_Y[1]});
      phys_prongMom_ByRange_h1_Z.push_back({cv.fProngMom_ByRange_h1_Z[0],cv.fProngMom_ByRange_h1_Z[1]});
      //
      phys_prongPdgCode_h2.push_back({cv.fProngPdgCode_h2[0],cv.fProngPdgCode_h2[1]});
      phys_prongMass_h2.push_back({cv.fProngMass_h2[0],cv.fProngMass_h2[1]});
      phys_prongMomMag_ByRange_h2.push_back({cv.fProngMomMag_ByRange_h2[0],cv.fProngMomMag_ByRange_h2[1]});
      phys_prongEnergy_ByRange_h2.push_back({cv.fProngEnergy_ByRange_h2[0],cv.fProngEnergy_ByRange_h2[1]});
      phys_prongMom_ByRange_h2_X.push_back({cv.fProngMom_ByRange_h2_X[0],cv.fProngMom_ByRange_h2_X[1]});
      phys_prongMom_ByRange_h2_Y.push_back({cv.fProngMom_ByRange_h2_Y[0],cv.fProngMom_ByRange_h2_Y[1]});
      phys_prongMom_ByRange_h2_Z.push_back({cv.fProngMom_ByRange_h2_Z[0],cv.fProngMom_ByRange_h2_Z[1]});
      // Tot momentum (by direction, assuming both muons)
      phys_totMomMag_ByRange_h1.push_back(cv.fTotMomMag_ByRange_h1);
      phys_totEnergy_ByRange_h1.push_back(cv.fTotEnergy_ByRange_h1);
      phys_invariantMass_ByRange_h1.push_back(cv.fInvMass_ByRange_h1);
      phys_totMom_ByRange_h1_X.push_back(cv.fTotMom_ByRange_h1_X);
      phys_totMom_ByRange_h1_Y.push_back(cv.fTotMom_ByRange_h1_Y);
      phys_totMom_ByRange_h1_Z.push_back(cv.fTotMom_ByRange_h1_Z);
      //
      phys_totMomMag_ByRange_h2.push_back(cv.fTotMomMag_ByRange_h2);
      phys_totEnergy_ByRange_h2.push_back(cv.fTotEnergy_ByRange_h2);
      phys_invariantMass_ByRange_h2.push_back(cv.fInvMass_ByRange_h2);
      phys_totMom_ByRange_h2_X.push_back(cv.fTotMom_ByRange_h2_X);
      phys_totMom_ByRange_h2_Y.push_back(cv.fTotMom_ByRange_h2_Y);
      phys_totMom_ByRange_h2_Z.push_back(cv.fTotMom_ByRange_h2_Z);
      // Tot momentum direction (by direction, assuming both muons)
      phys_totTheta_ByRange_h1.push_back(cv.fTotTheta_ByRange_h1);
      phys_totPhi_ByRange_h1.push_back(cv.fTotPhi_ByRange_h1);
      phys_totDir_ByRange_h1_X.push_back(cv.fTotDir_ByRange_h1_X);
      phys_totDir_ByRange_h1_Y.push_back(cv.fTotDir_ByRange_h1_Y);
      phys_totDir_ByRange_h1_Z.push_back(cv.fTotDir_ByRange_h1_Z);
      //
      phys_totTheta_ByRange_h2.push_back(cv.fTotTheta_ByRange_h2);
      phys_totPhi_ByRange_h2.push_back(cv.fTotPhi_ByRange_h2);
      phys_totDir_ByRange_h2_X.push_back(cv.fTotDir_ByRange_h2_X);
      phys_totDir_ByRange_h2_Y.push_back(cv.fTotDir_ByRange_h2_Y);
      phys_totDir_ByRange_h2_Z.push_back(cv.fTotDir_ByRange_h2_Z);
      // Momentum (By Mcs)
      phys_prongPdgCodeHypothesis_ByMcs.push_back(cv.fProngPdgCodeHypothesis_ByMcs);
      phys_prongIsBestFwd_ByMcs.push_back(cv.fProngIsBestFwd_ByMcs);
      // Prong Momentum (By Mcs, best)
      phys_prongMomMag_ByMcs_best_h1.push_back(cv.fProngMomMag_ByMcs_best_h1);
      phys_prongEnergy_ByMcs_best_h1.push_back(cv.fProngEnergy_ByMcs_best_h1);
      phys_prongMom_ByMcs_best_h1_X.push_back(cv.fProngMom_ByMcs_best_h1_X);
      phys_prongMom_ByMcs_best_h1_Y.push_back(cv.fProngMom_ByMcs_best_h1_Y);
      phys_prongMom_ByMcs_best_h1_Z.push_back(cv.fProngMom_ByMcs_best_h1_Z);
      phys_prongMomMag_ByMcs_best_h2.push_back(cv.fProngMomMag_ByMcs_best_h2);
      phys_prongEnergy_ByMcs_best_h2.push_back(cv.fProngEnergy_ByMcs_best_h2);
      phys_prongMom_ByMcs_best_h2_X.push_back(cv.fProngMom_ByMcs_best_h2_X);
      phys_prongMom_ByMcs_best_h2_Y.push_back(cv.fProngMom_ByMcs_best_h2_Y);
      phys_prongMom_ByMcs_best_h2_Z.push_back(cv.fProngMom_ByMcs_best_h2_Z);
      // Tot momentum (by range, assuming both muons, best)
      phys_totMomMag_ByMcs_best_h1.push_back(cv.fTotMomMag_ByMcs_best_h1);
      phys_totEnergy_ByMcs_best_h1.push_back(cv.fTotEnergy_ByMcs_best_h1);
      phys_invariantMass_ByMcs_best_h1.push_back(cv.fInvMass_ByMcs_best_h1);
      phys_totMom_ByMcs_best_h1_X.push_back(cv.fTotMom_ByMcs_best_h1_X);
      phys_totMom_ByMcs_best_h1_Y.push_back(cv.fTotMom_ByMcs_best_h1_Y);
      phys_totMom_ByMcs_best_h1_Z.push_back(cv.fTotMom_ByMcs_best_h1_Z);
      phys_totMomMag_ByMcs_best_h2.push_back(cv.fTotMomMag_ByMcs_best_h2);
      phys_totEnergy_ByMcs_best_h2.push_back(cv.fTotEnergy_ByMcs_best_h2);
      phys_invariantMass_ByMcs_best_h2.push_back(cv.fInvMass_ByMcs_best_h2);
      phys_totMom_ByMcs_best_h2_X.push_back(cv.fTotMom_ByMcs_best_h2_X);
      phys_totMom_ByMcs_best_h2_Y.push_back(cv.fTotMom_ByMcs_best_h2_Y);
      phys_totMom_ByMcs_best_h2_Z.push_back(cv.fTotMom_ByMcs_best_h2_Z);
      // Tot momentum direction (by range, assuming both muons, best)
      phys_totTheta_ByMcs_best_h1.push_back(cv.fTotTheta_ByMcs_best_h1);
      phys_totPhi_ByMcs_best_h1.push_back(cv.fTotPhi_ByMcs_best_h1);
      phys_totDir_ByMcs_best_h1_X.push_back(cv.fTotDir_ByMcs_best_h1_X);
      phys_totDir_ByMcs_best_h1_Y.push_back(cv.fTotDir_ByMcs_best_h1_Y);
      phys_totDir_ByMcs_best_h1_Z.push_back(cv.fTotDir_ByMcs_best_h1_Z);
      phys_totTheta_ByMcs_best_h2.push_back(cv.fTotTheta_ByMcs_best_h2);
      phys_totPhi_ByMcs_best_h2.push_back(cv.fTotPhi_ByMcs_best_h2);
      phys_totDir_ByMcs_best_h2_X.push_back(cv.fTotDir_ByMcs_best_h2_X);
      phys_totDir_ByMcs_best_h2_Y.push_back(cv.fTotDir_ByMcs_best_h2_Y);
      phys_totDir_ByMcs_best_h2_Z.push_back(cv.fTotDir_ByMcs_best_h2_Z);

      // Others
      phys_prongLength.push_back({cv.fProngLength[0],cv.fProngLength[1]});
      phys_prongNumHits.push_back({cv.fProngNumHits[0],cv.fProngNumHits[1]});
      phys_openingAngle.push_back(cv.fOpeningAngle);
      phys_prongStartToNeutrinoDistance.push_back({cv.fProngStartToNeutrinoDistance[0],cv.fProngStartToNeutrinoDistance[1]});
    }
    return;
  } // END function ExtractVertexPhysics

} // END namespace EventDescriptor 
