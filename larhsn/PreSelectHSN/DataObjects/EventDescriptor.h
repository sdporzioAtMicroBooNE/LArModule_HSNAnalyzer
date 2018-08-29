/******************************************************************************
 * @file EventDescriptor.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  EventDescriptor.cxx
 * ****************************************************************************/

#ifndef EVENTDESCRIPTOR_H
#define EVENTDESCRIPTOR_H

// C++ standard libraries
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <stdexcept>
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcorealg/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larhsn/PreSelectHSN/DataObjects/DecayVertex.h"

namespace AuxEvent
{

  // EventDescriptor class and functions
  class EventDescriptor
  {
  public:
    // Constructor and destructor
    EventDescriptor();
    virtual ~EventDescriptor();

    void Initialize(int i_run, int i_subrun, int i_event);
    void ExtractVertexPhysics(const std::vector<AuxVertex::DecayVertex> & decayVertices);

    // General
    int run;
    int subrun;
    int event;

    // Pandora search
    int nNeutrinos;
    int nTwoProngedNeutrinos;
    int nContainedTwoProngedNeutrinos;
    std::vector<int> neutrinoPdgCode;
    std::vector<int> neutrinoNumDaughters;
    std::vector<int> neutrinoNumTracks;
    std::vector<int> neutrinoNumShowers;
    // Pandora calo
    std::vector<std::vector<float>> calo_totChargeInRadius;
    std::vector<std::vector<float>> calo_prong1ChargeInRadius;
    std::vector<std::vector<float>> calo_prong2ChargeInRadius;
    std::vector<std::vector<float>> calo_caloRatio;
    // Pandora physics
    // Coordinates
    std::vector<float> phys_nuPosX, phys_nuPosY, phys_nuPosZ;
    std::vector<std::vector<float>> phys_prongPosX, phys_prongPosY, phys_prongPosZ;
    std::vector<std::vector<float>> phys_prongStartPosX, phys_prongStartPosY, phys_prongStartPosZ;
    std::vector<std::vector<float>> phys_prongEndPosX, phys_prongEndPosY, phys_prongEndPosZ;
    // Direction
    std::vector<std::vector<float>> phys_prongDirX, phys_prongDirY, phys_prongDirZ;
    std::vector<std::vector<float>> phys_prongTheta, phys_prongPhi;
    // Hypothesis information
    std::vector<std::vector<int>> phys_prongPdgCode_h1, phys_prongPdgCode_h2;
    std::vector<std::vector<float>> phys_prongMass_h1, phys_prongMass_h2;
    // Prong momentum (by range, assuming both muons)
    std::vector<std::vector<float>> phys_prongMomMag_ByRange_h1, phys_prongEnergy_ByRange_h1;
    std::vector<std::vector<float>> phys_prongMom_ByRange_h1_X, phys_prongMom_ByRange_h1_Y, phys_prongMom_ByRange_h1_Z;
    std::vector<std::vector<float>> phys_prongMomMag_ByRange_h2, phys_prongEnergy_ByRange_h2;
    std::vector<std::vector<float>> phys_prongMom_ByRange_h2_X, phys_prongMom_ByRange_h2_Y, phys_prongMom_ByRange_h2_Z;
    // Tot momentum (by range, assuming both muons)
    std::vector<float> phys_totMomMag_ByRange_h1, phys_totEnergy_ByRange_h1, phys_invariantMass_ByRange_h1;
    std::vector<float> phys_totMom_ByRange_h1_X, phys_totMom_ByRange_h1_Y, phys_totMom_ByRange_h1_Z;
    std::vector<float> phys_totMomMag_ByRange_h2, phys_totEnergy_ByRange_h2, phys_invariantMass_ByRange_h2;
    std::vector<float> phys_totMom_ByRange_h2_X, phys_totMom_ByRange_h2_Y, phys_totMom_ByRange_h2_Z;
    // Tot momentum direction (by range, assuming both muons)
    std::vector<float> phys_totTheta_ByRange_h1, phys_totPhi_ByRange_h1;
    std::vector<float> phys_totDir_ByRange_h1_X, phys_totDir_ByRange_h1_Y, phys_totDir_ByRange_h1_Z;
    std::vector<float> phys_totTheta_ByRange_h2, phys_totPhi_ByRange_h2;
    std::vector<float> phys_totDir_ByRange_h2_X, phys_totDir_ByRange_h2_Y, phys_totDir_ByRange_h2_Z;
    // MCS variables
    std::vector<std::vector<int>> phys_prongPdgCodeHypothesis_ByMcs;
    std::vector<std::vector<bool>> phys_prongIsBestFwd_ByMcs;
    // Prong Momentum (By Mcs, forward)
    std::vector<std::vector<float>> phys_prongMomMag_ByMcs_fwd_h1, phys_prongEnergy_ByMcs_fwd_h1;
    std::vector<std::vector<float>> phys_prongMom_ByMcs_fwd_h1_X, phys_prongMom_ByMcs_fwd_h1_Y, phys_prongMom_ByMcs_fwd_h1_Z;
    std::vector<std::vector<float>> phys_prongMomMag_ByMcs_fwd_h2, phys_prongEnergy_ByMcs_fwd_h2;
    std::vector<std::vector<float>> phys_prongMom_ByMcs_fwd_h2_X, phys_prongMom_ByMcs_fwd_h2_Y, phys_prongMom_ByMcs_fwd_h2_Z;
    // Tot momentum (by range, assuming both muons, forward)
    std::vector<float> phys_totMomMag_ByMcs_fwd_h1, phys_totEnergy_ByMcs_fwd_h1, phys_invariantMass_ByMcs_fwd_h1;
    std::vector<float> phys_totMom_ByMcs_fwd_h1_X, phys_totMom_ByMcs_fwd_h1_Y, phys_totMom_ByMcs_fwd_h1_Z;
    std::vector<float> phys_totMomMag_ByMcs_fwd_h2, phys_totEnergy_ByMcs_fwd_h2, phys_invariantMass_ByMcs_fwd_h2;
    std::vector<float> phys_totMom_ByMcs_fwd_h2_X, phys_totMom_ByMcs_fwd_h2_Y, phys_totMom_ByMcs_fwd_h2_Z;
    // Tot momentum direction (by range, assuming both muons, forward)
    std::vector<float> phys_totTheta_ByMcs_fwd_h1, phys_totPhi_ByMcs_fwd_h1;
    std::vector<float> phys_totDir_ByMcs_fwd_h1_X, phys_totDir_ByMcs_fwd_h1_Y, phys_totDir_ByMcs_fwd_h1_Z;
    std::vector<float> phys_totTheta_ByMcs_fwd_h2, phys_totPhi_ByMcs_fwd_h2;
    std::vector<float> phys_totDir_ByMcs_fwd_h2_X, phys_totDir_ByMcs_fwd_h2_Y, phys_totDir_ByMcs_fwd_h2_Z;
    // Prong Momentum (By Mcs, best)
    std::vector<std::vector<float>> phys_prongMomMag_ByMcs_best_h1, phys_prongEnergy_ByMcs_best_h1;
    std::vector<std::vector<float>> phys_prongMom_ByMcs_best_h1_X, phys_prongMom_ByMcs_best_h1_Y, phys_prongMom_ByMcs_best_h1_Z;
    std::vector<std::vector<float>> phys_prongMomMag_ByMcs_best_h2, phys_prongEnergy_ByMcs_best_h2;
    std::vector<std::vector<float>> phys_prongMom_ByMcs_best_h2_X, phys_prongMom_ByMcs_best_h2_Y, phys_prongMom_ByMcs_best_h2_Z;
    // Tot momentum (by range, assuming both muons, best)
    std::vector<float> phys_totMomMag_ByMcs_best_h1, phys_totEnergy_ByMcs_best_h1, phys_invariantMass_ByMcs_best_h1;
    std::vector<float> phys_totMom_ByMcs_best_h1_X, phys_totMom_ByMcs_best_h1_Y, phys_totMom_ByMcs_best_h1_Z;
    std::vector<float> phys_totMomMag_ByMcs_best_h2, phys_totEnergy_ByMcs_best_h2, phys_invariantMass_ByMcs_best_h2;
    std::vector<float> phys_totMom_ByMcs_best_h2_X, phys_totMom_ByMcs_best_h2_Y, phys_totMom_ByMcs_best_h2_Z;
    // Tot momentum direction (by range, assuming both muons, best)
    std::vector<float> phys_totTheta_ByMcs_best_h1, phys_totPhi_ByMcs_best_h1;
    std::vector<float> phys_totDir_ByMcs_best_h1_X, phys_totDir_ByMcs_best_h1_Y, phys_totDir_ByMcs_best_h1_Z;
    std::vector<float> phys_totTheta_ByMcs_best_h2, phys_totPhi_ByMcs_best_h2;
    std::vector<float> phys_totDir_ByMcs_best_h2_X, phys_totDir_ByMcs_best_h2_Y, phys_totDir_ByMcs_best_h2_Z;
    // Others
    std::vector<std::vector<float>> phys_prongLength;
    std::vector<std::vector<float>> phys_prongStartToNeutrinoDistance;
    std::vector<std::vector<int>> phys_prongNumHits;
    std::vector<float> phys_openingAngle;
    // Pandora statusnostic
    int status_nuWithMissingAssociatedVertex;
    int status_nuWithMissingAssociatedTrack;
    int status_nuProngWithMissingAssociatedHits;
    // Truth information
    std::vector<std::vector<int>> match_pdgCode;
    std::vector<std::vector<float>> match_mass;
    std::vector<std::vector<float>> match_energy;
    std::vector<std::vector<float>> match_prong1StartPosition;
    std::vector<std::vector<float>> match_prong2StartPosition;
    std::vector<std::vector<float>> match_prong1Momentum;
    std::vector<std::vector<float>> match_prong2Momentum;

  };


} //END namespace AuxEvent

#endif
