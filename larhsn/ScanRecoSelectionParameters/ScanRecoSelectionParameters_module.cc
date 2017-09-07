#ifndef ScanRecoSelectionParameters_module
#define ScanRecoSelectionParameters_module

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

// larsoft object includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Wire.h"

// Auxiliary objects includes
// #include "DecayVertex.h"


// Analyzer class
class ScanRecoSelectionParameters : public art::EDAnalyzer
{
public:
  explicit ScanRecoSelectionParameters(fhicl::ParameterSet const & pset);
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Declare fhiclcpp variables
  std::string fInstanceName;
  int fIteration;
  std::vector<double> fMinTpcBound;
  std::vector<double> fMaxTpcBound;
  std::string fDataType;
  std::string fTrackLabel;
  std::string fPfpLabel;
  int fDecayChannel;
  double fSterileMass;
  double fDistanceCut;
  bool fPrimaryOnly;
  bool fEndVerticesAlso;
  std::string fAnaType;

  // Declare services
  geo::GeometryCore const* fGeometry; // Pointer to the Geometry service
  detinfo::DetectorProperties const* fDetectorProperties; // Pointer to the Detector Properties

  // Declare script variables
  TTree *tTree;
  TTree *metaTree;
  Int_t run, subrun, event, nTracks, nShowers, nPairs, nTrackVertices, nShowerVertices, nPotVertices, nCleanVertices, pandora_nPrimaryVertices, pandora_nCleanVertices, nCleanVerticesOutsideTPC;
  std::vector<float> pairDistance, potPairDistance;
  std::vector<int> pandora_primaryVertexPDG, pandora_nDaughters, pandora_nTracks, pandora_nShowers, pandora_nNearTracks, pandora_nNearShowers;

  void ClearData();
  // void GetTrackShowerVectors(art::Event const & evt, std::vector<recob::PFParticle>& pandora_primaryPFP, std::vector<recob::Track const*>& tracks, std::vector<recob::Shower const*>& showers);
  // void GetOriginVertices(const std::vector<recob::Track const*>& tracks, const std::vector<recob::Shower const*>& showers, std::vector<ScanRecoSelectionParameters::DecayVertex>& trackVertices, std::vector<ScanRecoSelectionParameters::DecayVertex>& showerVertices);
  // void GetDecayVertices(const std::vector<ScanRecoSelectionParameters::DecayVertex>& trackVertices, const std::vector<ScanRecoSelectionParameters::DecayVertex>& showerVertices, std::vector<ScanRecoSelectionParameters::DecayVertex>& potVertices, std::vector<ScanRecoSelectionParameters::DecayVertex>& cleanVertices);
  // void PerformPandoraAnalysis(art::Event const & evt, const std::vector<recob::PFParticle>& pandora_primaryPFP);
}; // End class ScanRecoSelectionParameters

ScanRecoSelectionParameters::ScanRecoSelectionParameters(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fInstanceName(pset.get<std::string>("InstanceName")),
    fIteration(pset.get<int>("Iteration")),
    fMinTpcBound(pset.get<std::vector<double>>("MinTpcBound")),
    fMaxTpcBound(pset.get<std::vector<double>>("MaxTpcBound")),
    fDataType(pset.get<std::string>("DataType")),
    fTrackLabel(pset.get<std::string>("TrackLabel")),
    fPfpLabel(pset.get<std::string>("PfpLabel")),
    fDecayChannel(pset.get<int>("DecayChannel")),
    fSterileMass(pset.get<double>("SterileMass")),
    fDistanceCut(pset.get<double>("DistanceCut")),
    fPrimaryOnly(pset.get<bool>("PrimaryOnly")),
    fEndVerticesAlso(pset.get<bool>("EndVerticesAlso")),
    fAnaType(pset.get<std::string>("AnalysisType"))
{} // END function ScanRecoSelectionParameters

void ScanRecoSelectionParameters::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  metaTree = tfs->make<TTree>("Metadata","");
  metaTree->Branch("instanceName",&fInstanceName);
  metaTree->Branch("iteration",&fIteration,"iteration/I");  
  metaTree->Branch("dataType",&fDataType);
  metaTree->Branch("trackLabel",&fTrackLabel);
  metaTree->Branch("pfpLabel",&fPfpLabel);
  metaTree->Branch("decayChannel",&fDecayChannel,"decayChannel/I");
  metaTree->Branch("distanceCut",&fDistanceCut,"distanceCut/D");
  metaTree->Branch("sterileMass",&fSterileMass,"sterileMass/D");
  metaTree->Branch("primaryOnly",&fPrimaryOnly,"primaryOnly/O");
  metaTree->Branch("endVerticesAlso",&fEndVerticesAlso,"endVerticesAlso/O");
  metaTree->Branch("anaType",&fAnaType);
  metaTree->Fill();

  tTree = tfs->make<TTree>("Data","");
  tTree->Branch("run",&run,"run/I");
  tTree->Branch("subrun",&subrun,"subrun/I");
  tTree->Branch("event",&event,"event/I");
  tTree->Branch("nTracks",&nTracks,"nTracks/I");
  tTree->Branch("nShowers",&nShowers,"nShowers/I");
  tTree->Branch("pairDistance",&pairDistance);
  tTree->Branch("nPairs",&nPairs,"nPairs/I");
  tTree->Branch("nTrackVertices",&nTrackVertices,"nTrackVertices/I");
  tTree->Branch("nShowerVertices",&nShowerVertices,"nShowerVertices/I");
  tTree->Branch("potPairDistance",&potPairDistance);
  tTree->Branch("nPotVertices",&nPotVertices,"nPotVertices/I");
  tTree->Branch("nCleanVertices",&nCleanVertices,"nCleanVertices/I");
  tTree->Branch("pandora_nPrimaryVertices",&pandora_nPrimaryVertices,"pandora_nPrimaryVertices/I");
  tTree->Branch("pandora_primaryVertexPDG",&pandora_primaryVertexPDG);
  tTree->Branch("pandora_nDaughters",&pandora_nDaughters);
  tTree->Branch("pandora_nTracks",&pandora_nTracks);
  tTree->Branch("pandora_nShowers",&pandora_nShowers);
  tTree->Branch("pandora_nNearTracks",&pandora_nNearTracks);
  tTree->Branch("pandora_nNearShowers",&pandora_nNearShowers);
  tTree->Branch("pandora_nCleanVertices",&pandora_nCleanVertices,"pandora_nCleanVertices/I");
} // END function beginJob

void ScanRecoSelectionParameters::endJob()
{
} // END function endJob

void ScanRecoSelectionParameters::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  nTracks = 0;
  nShowers= 0;
  nPairs = 0;
  nTrackVertices = 0;
  nShowerVertices = 0;
  nPotVertices = 0;
  nCleanVertices = 0;
  pandora_nPrimaryVertices = 0;
  pandora_nCleanVertices = 0;
  nCleanVerticesOutsideTPC = 0;
  pairDistance.clear();
  potPairDistance.clear();
  pandora_primaryVertexPDG.clear();
  pandora_nDaughters.clear();
  pandora_nTracks.clear();
  pandora_nShowers.clear();
  pandora_nNearTracks.clear();
  pandora_nNearShowers.clear();
} // END function ClearData

// void ScanRecoSelectionParameters::GetTrackShowerVectors(art::Event const & evt, std::vector<recob::PFParticle>& pandora_primaryPFP, std::vector<recob::Track const*>& tracks, std::vector<recob::Shower const*>& showers)
// {
//   // Take an event (evt), then fill the (empty!) track, shower and primaries pfp vectors (tracks,showers,pandora_primaryPFP) with all the tracks, showers and primaries pfps in the event. These vectors will then be used to determine all the possible vertices.

//   // Prepare handle labels
//   art::InputTag trackTag {fTrackLabel};
//   art::InputTag pfpTag {fPfpLabel};

//   // Find associations from PFP to tracks and showers
//   const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
//   art::FindMany<recob::Track> pfp_to_track_assns(pfpHandle,evt,trackTag);
//   art::FindMany<recob::Shower> pfp_to_shower_assns(pfpHandle,evt,trackTag);

//   // Loop through each PFP
//   int index = 0;
//   for (auto const& pfp : (*pfpHandle)){
//     // Generate temporary vectors to contain objects associated with pfp
//     std::vector<recob::Track const*> tempTrackVector;
//     std::vector<recob::Shower const*> tempShowerVector;

//     // Fill temporary vectors with objects associated with the pfp
//     pfp_to_track_assns.get(index,tempTrackVector);
//     pfp_to_shower_assns.get(index,tempShowerVector);

//     // Group all the primaries to analyze them later
//     if (pfp.IsPrimary()){
//       pandora_primaryPFP.push_back(pfp);
//     }
//     // Ignore (unphysical) reconstructed parents
//     // (e.g. "ghost" neutrinos created by PandoraNu)
//     else {
//       // Find current pfp parent and check if parent is primary (so that pfp are only secondaries, and not tertiaries, etc.), then fill the correspondending vector with the type of particle found
//       // tempTrackVector, tempShowerVector should always contain only ONE element; that exact track or shower (which is why .at(0)). No idea why it's a vector.
//       if (fPrimaryOnly){
//         auto const& parent = (*pfpHandle).at(pfp.Parent());
//         if (parent.IsPrimary()){
//           if (tempTrackVector.size()!=0) tracks.push_back(tempTrackVector.at(0));
//           if (tempShowerVector.size()!=0) showers.push_back(tempShowerVector.at(0));
//         }
//       }
//       // Fill the track and shower vectors with the products associated with the pfp (in any case, regardless of their "genealogy")
//       else {
//         if (tempTrackVector.size()!=0) tracks.push_back(tempTrackVector.at(0));
//         if (tempShowerVector.size()!=0) showers.push_back(tempShowerVector.at(0));
//       }
//       index++; // Index loops through pfps and MUST always keep on (no ifs)
//     }
//   } // End of pfp loop

//   return;
// } // END function GetTrackShowerVectors

// void ScanRecoSelectionParameters::GetOriginVertices(const std::vector<recob::Track const*>& tracks, const std::vector<recob::Shower const*>& showers, std::vector<ScanRecoSelectionParameters::DecayVertex>& trackVertices, std::vector<ScanRecoSelectionParameters::DecayVertex>& showerVertices)
// {
//     // Given track and shower vectors (tracks,showers), fill vertices vectors (trackVertices,showerVertices) taking all the start and end points for tracks and start points for showers. The last two numbers (j,j) indicate only their idx in the xxxVertices vector (which is reduntant at this stage). However, later on, it will be used to identify the vertices uniquely, and the idxs from two vertices will be used to define the parent idx of the mean vertex created between them.

//     // For tracks
//     int j = 0;
//     for(std::vector<int>::size_type i=0; i!=tracks.size(); i++){
//       ScanRecoSelectionParameters::DecayVertex tempV1(tracks[i]->Vertex().X(),tracks[i]->Vertex().Y(),tracks[i]->Vertex().Z(),j,j,"t","t");
//       trackVertices.push_back(tempV1);
//       j++;
//       if (fEndVerticesAlso){
//         ScanRecoSelectionParameters::DecayVertex tempV2(tracks[i]->End().X(),tracks[i]->End().Y(),tracks[i]->End().Z(),j,j,"t","t");
//         trackVertices.push_back(tempV2);
//         j++;
//       }
//     }

//     // And for showers
//     for(std::vector<int>::size_type i=0; i!=showers.size(); i++){
//       ScanRecoSelectionParameters::DecayVertex tempV(showers[i]->ShowerStart().X(),showers[i]->ShowerStart().Y(),showers[i]->ShowerStart().Z(),i,i,"s","s");
//       showerVertices.push_back(tempV);
//     }
//   return;
// } // END function GetOriginVertices

// void ScanRecoSelectionParameters::GetDecayVertices(const std::vector<ScanRecoSelectionParameters::DecayVertex>& trackVertices, const std::vector<ScanRecoSelectionParameters::DecayVertex>& showerVertices, std::vector<ScanRecoSelectionParameters::DecayVertex>& potVertices, std::vector<ScanRecoSelectionParameters::DecayVertex>& cleanVertices)
// {
//   // Given the reco vertices establishing the origin of tracks and showers (trackVertices,showerVertices), take all the possible mean vertices by paring them two by two (potVertices), then impose that in radius (fDistanceCut) only two original vertices are present (cleanVertices). Now the previously used (j,j) index are used to define the parent vertices (which track or shower originating vertex was used to determine the mean vertex between them).
//   // For track-track
//   if (fAnaType == "tt"){
//     for(std::vector<int>::size_type i=0; i!=trackVertices.size(); i++){
//       for(std::vector<int>::size_type j=i+1; j!=trackVertices.size(); j++){
//         ScanRecoSelectionParameters::DecayVertex v1 = trackVertices[i];
//         ScanRecoSelectionParameters::DecayVertex v2 = trackVertices[j];
//         float distance = Distance(v1,v2);
//         pairDistance.push_back(distance);
//         bool isInRadius = (distance<fDistanceCut);
//         if (isInRadius) {
//           ScanRecoSelectionParameters::DecayVertex v3 = MeanVertex(v1, v2);
//           potVertices.push_back(v3);
//           potPairDistance.push_back(distance);
//         }
//       }
//     }
//     // Make sure the potential mean vertices are clean (two decay vertices only in radius)
//     // (start with assumption that isGoodVertex == true, then if you find extra particles in radius, turn it in bad vertex. Only good vertices are saved)
//     for(std::vector<int>::size_type i=0; i!=potVertices.size(); i++){
//       ScanRecoSelectionParameters::DecayVertex mv = potVertices[i];
//       bool isGoodVertex = true;
//       for(std::vector<int>::size_type j=0; j!=trackVertices.size(); j++){
//         ScanRecoSelectionParameters::DecayVertex v1 = trackVertices[j];
//         bool isParent1 = (mv.PIDX1() == v1.PIDX1());
//         bool isParent2 = (mv.PIDX2() == v1.PIDX2());
//         bool notParent = !(isParent1 || isParent2);
//         bool isInRadius = (Distance(mv,v1)<fDistanceCut);
//         if (isInRadius && notParent) isGoodVertex = false;
//       }
//       if (isGoodVertex) cleanVertices.push_back(mv);
//     }
//   }
//   // And for track-shower
//   else if (fAnaType == "ts"){
//     for(std::vector<int>::size_type i=0; i!=trackVertices.size(); i++){
//       for(std::vector<int>::size_type j=0; j!=showerVertices.size(); j++){
//         ScanRecoSelectionParameters::DecayVertex v1 = trackVertices[i];
//         ScanRecoSelectionParameters::DecayVertex v2 = showerVertices[j];
//         float distance = Distance(v1,v2);
//         pairDistance.push_back(distance);
//         bool isInRadius = (distance<fDistanceCut);
//         if (isInRadius) {
//           ScanRecoSelectionParameters::DecayVertex v3 = MeanVertex(v1, v2);
//           potVertices.push_back(v3);
//           potPairDistance.push_back(distance);
//         }
//       }
//     }
//     // Make sure the potential mean vertices are clean (two decay vertices only in radius)
//     // (start with assumption that isGoodVertex == true, then if you find extra particles in radius, turn it in bad vertex. Only good vertices are saved)
//     for(std::vector<int>::size_type i=0; i!=potVertices.size(); i++){
//       ScanRecoSelectionParameters::DecayVertex mv = potVertices[i];
//       bool isGoodVertex = true;
//       for(std::vector<int>::size_type j=0; j!=trackVertices.size(); j++){
//         ScanRecoSelectionParameters::DecayVertex v1 = trackVertices[j];
//         bool notParent = (mv.PIDX1() != v1.PIDX1());
//         bool isInRadius = (Distance(mv,v1)<fDistanceCut);
//         if (isInRadius && notParent) isGoodVertex = false;
//       }
//       for(std::vector<int>::size_type j=0; j!=showerVertices.size(); j++){
//         ScanRecoSelectionParameters::DecayVertex v2 = showerVertices[j];
//         bool notParent = (mv.PIDX2() != v2.PIDX2());
//         bool isInRadius = (Distance(mv,v2)<fDistanceCut);
//         if (isInRadius && notParent) isGoodVertex = false;
//       }
//       if (isGoodVertex) cleanVertices.push_back(mv);
//     }
//   }
//   // Throw error if fAnaType wasn't right.
//   else {
//     throw std::invalid_argument("Invalid fAnaType. Must be 'tt' or 'ts'!");
//   }
//   return;
// } // END function GetDecayVertices

// void ScanRecoSelectionParameters::PerformPandoraAnalysis(art::Event const & evt, const std::vector<recob::PFParticle>& pandora_primaryPFP)
// {
//   // Repeat same analysis as before, this time, loop through each pfparticle which is a primary, find daughters, apply geo cut
//   // Prepare handle labels
//   art::InputTag pfpTag {fPfpLabel};
//   art::InputTag trackTag {fTrackLabel};

//   // Find associations from PFP to tracks and showers
//   const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
//   art::FindMany<recob::Vertex> pfp_to_vertex_assns(pfpHandle,evt,trackTag);
//   art::FindMany<recob::Track> pfp_to_track_assns(pfpHandle,evt,trackTag);
//   art::FindMany<recob::Shower> pfp_to_shower_assns(pfpHandle,evt,trackTag);

//   // Loop through each primary pfp
//   for(std::vector<int>::size_type i=0; i!=pandora_primaryPFP.size(); i++){
//     int temp_nTracks = 0, temp_nShowers = 0, temp_nNearTracks = 0, temp_nNearShowers = 0;
//     // Get primary pfp vertex, (j,j) indices don't matter anymore, so just leave them by default to 0.
//     recob::PFParticle pfp = pandora_primaryPFP[i];
//     pandora_nDaughters.push_back(pfp.NumDaughters());
//     pandora_primaryVertexPDG.push_back(pfp.PdgCode());
//     auto const& self_idx = pfp.Self();
//     std::vector<recob::Vertex const*> tempVertexVector;
//     pfp_to_vertex_assns.get(self_idx,tempVertexVector);
//     auto const& pfpVertex = tempVertexVector.at(0);
//     double tempCoords[3];
//     pfpVertex->XYZ(tempCoords);
//     ScanRecoSelectionParameters::DecayVertex v0(tempCoords[0],tempCoords[1],tempCoords[2],0,0,"n","n");

//     if (pfp.NumDaughters()!=0){
//       for(size_t k=0; k!=(unsigned int)pfp.NumDaughters(); k++){
//         size_t daughter_idx = pfp.Daughters().at(k);
//         std::vector<recob::Track const*> tempTrackVector;
//         std::vector<recob::Shower const*> tempShowerVector;
//         pfp_to_track_assns.get(daughter_idx,tempTrackVector);
//         pfp_to_shower_assns.get(daughter_idx,tempShowerVector);
//         if (tempTrackVector.size()==1){
//           temp_nTracks += 1;
//           auto const& track = tempTrackVector.at(0);
//           ScanRecoSelectionParameters::DecayVertex v1(track->Vertex().X(),track->Vertex().Y(),track->Vertex().Z(),0,0,"t","t");
//           if (Distance(v0,v1)<fDistanceCut) temp_nNearTracks += 1;
//           if (fEndVerticesAlso){
//             ScanRecoSelectionParameters::DecayVertex v2(track->End().X(),track->End().Y(),track->End().Z(),0,0,"t","t");
//             if (Distance(v0,v2)<fDistanceCut) temp_nNearTracks += 1;
//           }
//         }
//         else if (tempShowerVector.size()==1){
//           temp_nShowers += 1;
//           auto const& shower = tempShowerVector.at(0);
//           ScanRecoSelectionParameters::DecayVertex v1(shower->ShowerStart().X(),shower->ShowerStart().Y(),shower->ShowerStart().Z(),0,0,"s","s");
//           if (Distance(v0,v1)<fDistanceCut) temp_nNearShowers += 1;
//         }
//         else{
//           std::cout << "Something's wrong..." << std::endl;
//           std::cout << "A not-primary has PDG: " << (*pfpHandle).at(daughter_idx).PdgCode() << " and nTracks: " << tempTrackVector.size() << "and nShowers: " <<  tempShowerVector.size() << std::endl;
//         }
//       }
//     }
//     if (fAnaType == "tt" && temp_nNearTracks==2) pandora_nCleanVertices += 1;
//     if (fAnaType == "ts" && temp_nNearTracks==1 && temp_nNearShowers==1) pandora_nCleanVertices += 1;
//     pandora_nTracks.push_back(temp_nTracks);
//     pandora_nShowers.push_back(temp_nShowers);
//     pandora_nNearTracks.push_back(temp_nNearTracks);
//     pandora_nNearShowers.push_back(temp_nNearShowers);
//   }// End of primary loop
//   return;
// } // END function PerformPandoraAnalysis

void ScanRecoSelectionParameters::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.

  // Start by clearing all the vectors.
  ClearData();

  // AuxVertex::DecayVertex v1(10.,-52.,1322.,0,0,"t","t");
  // std::cout << AuxVertex::als << std::endl;


  // dv::DecayVertex v2(45.,32.,789.,1,1,"s","s");

  // dv::DecayVertex v0 = AuxVertex::MeanVertex(v1,v2);
  // printf("Distance: %.1f\n", AuxVertex::Distance(v1,v2));
  // printf("Vertex coordinates: [%.1f,%.1f,%.1f]\n", v0.GetX(), v0.GetY(), v0.GetZ());
  // printf("Vertex parent types: [%s,%s]\n",v0.GetParType1().c_str(),v0.GetParType2().c_str());




  // // Determine event ID
  // run = evt.id().run();
  // subrun = evt.id().subRun();
  // event = evt.id().event();

  // // Prepare track, shower and Pandora primary vectors
  // std::vector<recob::Track const*> tracks;
  // std::vector<recob::Shower const*> showers;
  // std::vector<recob::PFParticle> pandora_primaryPFP;
  // GetTrackShowerVectors(evt, pandora_primaryPFP, tracks, showers);

  // // Determine origin vertices for tracks and showers
  // std::vector<ScanRecoSelectionParameters::DecayVertex> trackVertices;
  // std::vector<ScanRecoSelectionParameters::DecayVertex> showerVertices;
  // GetOriginVertices(tracks, showers, trackVertices, showerVertices);

  // // Determine potential and clean decay vertices satisfying selection
  // std::vector<ScanRecoSelectionParameters::DecayVertex> potVertices;
  // std::vector<ScanRecoSelectionParameters::DecayVertex> cleanVertices;
  // GetDecayVertices(trackVertices, showerVertices, potVertices, cleanVertices);

  // // Perform Pandora analysis (similar to all that was before but separate, uses only primary PFP information to find vertices)
  // PerformPandoraAnalysis(evt, pandora_primaryPFP);

  // // Perform calorimetry analysis
  // // Wire and Tick coordinates vector (first dimension = vertex, second dimension = plane)
  // std::vector<std::vector<int>> wire0;
  // std::vector<std::vector<double>> tick0;
  // // Check if each vertex doesn't lie outside the TPC. If not, convert xyz coordinates to wire-tick coordinates.
  // for(auto v: cleanVertices){
  //   // Temporary filling arrays
  //   int temp_channels[3] = {0,0,0};
  //   double temp_ticks[3] = {0.,0.,0.};

  //   double xyz[3] = {v.X(),v.Y(),v.Z()};
  //   bool isVertexInside = IsInsideTPC(xyz);
  //   if (isVertexInside){
  //     CoordinatesToChannels(xyz,temp_channels);
  //     XToTicks(xyz[0],temp_ticks);
  //   }
  //   else nCleanVerticesOutsideTPC += 1;
  //   wire0.push_back({temp_channels[0],temp_channels[1],temp_channels[2]});
  //   tick0.push_back({temp_ticks[0],temp_ticks[1],temp_ticks[2]});
  // }






  // // Calculate final tree quantities
  // nTracks = tracks.size();
  // nShowers = showers.size();
  // nPairs = pairDistance.size();
  // nTrackVertices = trackVertices.size();
  // nShowerVertices = showerVertices.size();
  // nPotVertices = potVertices.size();
  // nCleanVertices = cleanVertices.size();
  // pandora_nPrimaryVertices = pandora_primaryPFP.size();

  // tTree->Fill();
} // END function analyze

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ScanRecoSelectionParameters)

#endif // END def ScanRecoSelectionParameters_module