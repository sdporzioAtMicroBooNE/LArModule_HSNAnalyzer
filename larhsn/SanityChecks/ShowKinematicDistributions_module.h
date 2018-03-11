#ifndef SHOWKINEMATICDISTRIBUTIONS_H
#define SHOWKINEMATICDISTRIBUTIONS_H

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
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/geo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Analyzer class
class ShowKinematicDistributions : public art::EDAnalyzer
{
public:
  explicit ShowKinematicDistributions(fhicl::ParameterSet const & pset);
  virtual ~ShowKinematicDistributions();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
  void GetTruthParticles(art::Event const & evt, std::vector<simb::MCParticle const*>& mcps);
  void ExtractKinematic(const std::vector<simb::MCParticle const*>& mcps);
private:

  // Declare fhiclcpp variables
  std::string fMCLabel;

  // Declare trees and tree variables
  TTree *tDataTree;
  std::vector<int> pdgCode;
  std::vector<double> Vx, Vy, Vz, T, EndX, EndY, EndZ, EndT, Px, Py, Pz, E, P, Pt, Theta, Phi;
  double Nu_E, Nu_Px, Nu_Py, Nu_Pz, Nu_P, Nu_Theta, Nu_Phi;
  double OpeningAngle, InvariantMass;


  // Declare analysis variables
  int run, subrun, event;

  // Declare analysis functions
  void ClearData();
}; // End class ShowKinematicDistributions

#endif // END def ShowKinematicDistributions header