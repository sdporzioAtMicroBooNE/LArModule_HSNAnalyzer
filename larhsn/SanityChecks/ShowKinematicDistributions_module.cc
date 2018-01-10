#ifndef ShowKinematicDistributions_module
#define ShowKinematicDistributions_module

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
// #include "lardataobj/RecoBase/Track.h"
// #include "lardataobj/RecoBase/Shower.h"
// #include "lardataobj/RecoBase/Vertex.h"
// #include "lardataobj/RecoBase/PFParticle.h"
// #include "lardataobj/RecoBase/Wire.h"
// #include "lardataobj/RecoBase/Hit.h"
// #include "lardataobj/RecoBase/TrackingTypes.h"
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

  // Declare services
  geo::GeometryCore const* fGeometry; // Pointer to the Geometry service
  detinfo::DetectorProperties const* fDetectorProperties; // Pointer to the Detector Properties

  // Declare trees
  TTree *tDataTree;
  std::vector<int> pdgCode;
  std::vector<double> Vx, Vy, Vz, T, EndX, EndY, EndZ, EndT, Px, Py, Pz, E, P, Pt;
  double OpeningAngle, InvariantMass;

  // Declare analysis variables
  int run, subrun, event;

  // Declare analysis functions
  void ClearData();
}; // End class ShowKinematicDistributions

ShowKinematicDistributions::ShowKinematicDistributions(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMCLabel(pset.get<std::string>("mcLabel"))  
{} // END constructor ShowKinematicDistributions

ShowKinematicDistributions::~ShowKinematicDistributions()
{} // END destructor ShowKinematicDistributions

void ShowKinematicDistributions::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run,"run/I");
  tDataTree->Branch("subrun",&subrun,"subrun/I");
  tDataTree->Branch("event",&event,"event/I");
  tDataTree->Branch("Vx",&Vx);
  tDataTree->Branch("Vy",&Vy);
  tDataTree->Branch("Vz",&Vz);
  tDataTree->Branch("T",&T);
  tDataTree->Branch("EndX",&EndX);
  tDataTree->Branch("EndY",&EndY);
  tDataTree->Branch("EndZ",&EndZ);
  tDataTree->Branch("EndT",&EndT);
  tDataTree->Branch("Px",&Px);
  tDataTree->Branch("Py",&Py);
  tDataTree->Branch("Pz",&Pz);
  tDataTree->Branch("E",&E);
  tDataTree->Branch("P",&P);
  tDataTree->Branch("Pt",&Pt);
  tDataTree->Branch("OpeningAngle",&OpeningAngle);
  tDataTree->Branch("InvariantMass",&InvariantMass);

  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
} // END function beginJob

void ShowKinematicDistributions::endJob()
{
} // END function endJob

void ShowKinematicDistributions::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  OpeningAngle = -999;
  InvariantMass = -999;
  pdgCode.clear();
  Vx.clear();
  Vy.clear();
  Vz.clear();
  T.clear();
  EndX.clear();
  EndY.clear();
  EndZ.clear();
  EndT.clear();
  Px.clear();
  Py.clear();
  Pz.clear();
  E.clear();
  P.clear();
  Pt.clear();
} // END function ClearData

void ShowKinematicDistributions::GetTruthParticles(art::Event const & evt, std::vector<simb::MCParticle const*>& mcps)
{
  // Prepare handle labels
  art::InputTag mcTag {fMCLabel};

  // Find mcTruth
  const auto& mctHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTag);
  for (auto const& mct : (*mctHandle)){
    int nParticles = mct.NParticles();
    printf("|_Number of particles: %i\n", nParticles);
    for (int i=0; i<nParticles; i++)
    {
      const simb::MCParticle & mcp = mct.GetParticle(i);
      pdgCode.push_back(mcp.PdgCode());
      Vx.push_back(mcp.Vx());
      Vy.push_back(mcp.Vy());
      Vz.push_back(mcp.Vz());
      T.push_back(mcp.T());
      EndX.push_back(mcp.EndX());
      EndY.push_back(mcp.EndY());
      EndZ.push_back(mcp.EndZ());
      EndT.push_back(mcp.EndT());
      Px.push_back(mcp.Px());
      Py.push_back(mcp.Py());
      Pz.push_back(mcp.Pz());
      E.push_back(mcp.E());
      P.push_back(mcp.P());
      Pt.push_back(mcp.Pt());
      mcps.push_back(&mcp);
    }
    if (nParticles==2)
    {
      double dotProduct = Px[0]*Px[1] + Py[0]*Py[1] + Pz[0]*Pz[1];
      OpeningAngle = dotProduct / float(P[0]*P[1]);
      double eTerm = pow((E[0] + E[1]),2.);
      double pTerm = pow(P[0],2.) + pow(P[1],2.) + 2.*dotProduct;
      InvariantMass = eTerm - pTerm;
    }
    else
    {
      OpeningAngle = -999;
      InvariantMass = -999;
    }

  } // End of pfp loop

  return;
} // END function GetTruthParticles


void ShowKinematicDistributions::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  printf("\n-------------------------------------------------------\n");
  
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

  // Get vector of primaries and secondaries pfps
  std::vector<simb::MCParticle const*> mcps;
  GetTruthParticles(evt, mcps);

  // Fill tree and finish event loop
  tDataTree->Fill();
  printf("-------------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ShowKinematicDistributions)

#endif // END def ShowKinematicDistributions_module