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

    // General
    int run;
    int subrun;
    int event;

    // Manual search
    int manual_NumPrimaries;
    int manual_NumTracks;
    int manual_NumShowers;
    int manual_NumTrackVertices;
    int manual_NumShowerVertices;
    int manual_NumVertexPairs;
    int manual_NumPotVertices;
    int manual_NumCleanVertices;
    int manual_NumUnassociatedTracks;
    int manual_NumUnassociatedShowers;
    std::vector<float> manual_pairDistances;
    std::vector<float> manual_potPairDistances;
  };
  

} //END namespace AuxEvent

#endif