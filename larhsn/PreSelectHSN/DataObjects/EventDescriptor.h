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
#include "larcore/Geometry/geo.h"
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

    void Initialize(int run, int subrun, int event);
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
    int calo_NumTotHits;
    std::vector<std::vector<float>> calo_totChargeInRadius;
    std::vector<std::vector<float>> calo_prong1ChargeInRadius;
    std::vector<std::vector<float>> calo_prong2ChargeInRadius;
    std::vector<std::vector<float>> calo_caloRatio;
    // Pandora physics
    std::vector<std::vector<float>> phys_prongLength;
    std::vector<std::vector<float>> phys_prongTheta;
    std::vector<std::vector<float>> phys_prongPhi;
    std::vector<std::vector<int>> phys_prongNumHits;
    // Pandora diagnostic
    int diag_nuWithMissingAssociatedVertex;
    int diag_nuWithMissingAssociatedTrack;
    int diag_nuProngWithMissingAssociatedHits;
  };
  

} //END namespace AuxEvent

#endif