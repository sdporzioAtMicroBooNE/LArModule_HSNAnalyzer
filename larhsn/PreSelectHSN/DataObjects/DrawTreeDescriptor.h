/******************************************************************************
 * @file DrawTreeDescriptor.h
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DrawTreeDescriptor.cxx
 * ****************************************************************************/

#ifndef DRAWTREEDESCRIPTOR_H
#define DRAWTREEDESCRIPTOR_H

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

  // DrawTreeDescriptor class and functions
  class DrawTreeDescriptor
  {
  public:
    // Constructor and destructor
    DrawTreeDescriptor();
    virtual ~DrawTreeDescriptor();
    void Initialize();

    void FillDrawTreeVariables(
          const std::vector<AuxVertex::DecayVertex>& decayVertices);

    // Declare drawTree variables
    std::vector<std::vector<int>> dv_wireCoordinates;
    std::vector<std::vector<int>> prong1_wireCoordinates;
    std::vector<std::vector<int>> prong2_wireCoordinates;

    std::vector<std::vector<int>> prong1_hits_p0_wireCoordinates;
    std::vector<std::vector<int>> prong1_hits_p1_wireCoordinates;
    std::vector<std::vector<int>> prong1_hits_p2_wireCoordinates;
    std::vector<std::vector<int>> prong2_hits_p0_wireCoordinates; 
    std::vector<std::vector<int>> prong2_hits_p1_wireCoordinates; 
    std::vector<std::vector<int>> prong2_hits_p2_wireCoordinates;
    std::vector<std::vector<int>> tot_hits_p0_wireCoordinates; 
    std::vector<std::vector<int>> tot_hits_p1_wireCoordinates; 
    std::vector<std::vector<int>> tot_hits_p2_wireCoordinates;

    std::vector<std::vector<float>> dv_tickCoordinates;
    std::vector<std::vector<float>> prong1_tickCoordinates;
    std::vector<std::vector<float>> prong2_tickCoordinates;

    std::vector<std::vector<float>> prong1_hits_p0_tickCoordinates; 
    std::vector<std::vector<float>> prong1_hits_p1_tickCoordinates; 
    std::vector<std::vector<float>> prong1_hits_p2_tickCoordinates;
    std::vector<std::vector<float>> prong2_hits_p0_tickCoordinates; 
    std::vector<std::vector<float>> prong2_hits_p1_tickCoordinates; 
    std::vector<std::vector<float>> prong2_hits_p2_tickCoordinates;
    std::vector<std::vector<float>> tot_hits_p0_tickCoordinates; 
    std::vector<std::vector<float>> tot_hits_p1_tickCoordinates; 
    std::vector<std::vector<float>> tot_hits_p2_tickCoordinates;

    std::vector<std::vector<float>> dv_xyzCoordinates; 
    std::vector<std::vector<float>> prong1_xyzCoordinates; 
    std::vector<std::vector<float>> prong2_xyzCoordinates; 
  };
  

} //END namespace AuxEvent

#endif