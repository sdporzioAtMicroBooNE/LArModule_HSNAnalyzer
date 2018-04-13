#ifndef TRUTHMATCHINGALG_H
#define TRUTHMATCHINGALG_H

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
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"
//s#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"


// Auxiliary objects includes
#include "larhsn/PreSelectHSN/DataObjects/DecayVertex.h"
#include "larhsn/PreSelectHSN/DataObjects/EventDescriptor.h"

namespace TruthMatching
{

  class TruthMatchingAlg
  {
  public:
    TruthMatchingAlg(fhicl::ParameterSet const & pset);
    ~TruthMatchingAlg();
    void reconfigure(fhicl::ParameterSet const & pset);

  // Algorithms
  void PerformTruthMatching(
          art::Event const & evt,
          AuxEvent::EventDescriptor & evd,
          std::vector<AuxVertex::DecayVertex>& decayVertices);

  // PerformCalorimetry returns
  std::vector<std::vector<int>> prong_matchedPDG;

  private:
    // fhicl parameters
    std::string fPfpLabel;
    std::string fHitLabel;
    bool fVerbose;

    // microboone services
    const geo::GeometryCore* fGeometry;
    const detinfo::DetectorProperties* fDetectorProperties;
  };

} // END namespace CalorimetryRadius

#endif
