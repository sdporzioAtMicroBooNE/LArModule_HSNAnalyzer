#ifndef FINDDECAYVERTEXALG_H
#define FINDDECAYVERTEXALG_H

// c++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <chrono>
#include <exception>

// root includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2I.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TGraph.h"

// framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"

// art includes
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindOneP.h"

// larsoft object includes
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
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Auxiliary objects includes
#include "larhsn/PreSelectHSN/DataObjects/DecayVertex.h"
#include "larhsn/PreSelectHSN/DataObjects/EventDescriptor.h"

namespace FindDecayVertex
{

  class FindDecayVertexAlg
  {
  public:
    FindDecayVertexAlg(fhicl::ParameterSet const & pset);
    ~FindDecayVertexAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

    // algorithms
    void GetTrackShowerVectors(
            AuxEvent::EventDescriptor & evd,
            art::Event const & evt);

    void GetOriginVertices(
            AuxEvent::EventDescriptor & evd,
            art::Event const & evt,
            const std::vector<recob::PFParticle const*>& tracks,
            const std::vector<recob::PFParticle const*>& showers);

    void GetDecayVertices(
            AuxEvent::EventDescriptor & evd,
            const std::vector<AuxVertex::DecayVertex>& trackVertices,
            const std::vector<AuxVertex::DecayVertex>& showerVertices);

    // analysis variables
    std::vector<recob::PFParticle const*>  ana_manual_primaries, ana_manual_tracks, ana_manual_showers;
    std::vector<AuxVertex::DecayVertex> ana_manual_trackVertices, ana_manual_showerVertices;
    std::vector<AuxVertex::DecayVertex> ana_manual_potVertices, ana_manual_cleanVertices;

  private:
    // fhicl parameters
    std::string fPfpLabel;
    std::string fAnaType;
    std::vector<double> fMinTpcBound;
    std::vector<double> fMaxTpcBound;
    double fDistanceCut;
    bool fPrimaryOnly;
    bool fEndVerticesAlso;
    bool fVerbose;

    // microboone services
    const geo::GeometryCore* fGeometry;
    const detinfo::DetectorProperties* fDetectorProperties;
  };

} // END namespace FindDecayVertex

#endif