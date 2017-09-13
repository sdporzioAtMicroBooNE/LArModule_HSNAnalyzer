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
#include "canvas/Persistency/Common/FindOneP.h"


// larsoft object includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackingTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Auxiliary objects includes
#include "DataObjects/DecayVertex.h"


// Analyzer class
class ScanRecoSelectionParameters : public art::EDAnalyzer
{
public:
  explicit ScanRecoSelectionParameters(fhicl::ParameterSet const & pset);
  virtual ~ScanRecoSelectionParameters();
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
  std::string fPfpLabel;
  std::string fHitLabel;
  int fDecayChannel;
  double fSterileMass;
  double fDistanceCut;
  double fCaloCut;
  double fChannelNorm;
  double fTickNorm;
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

  // Geometry functions
  void SetDetectorCoordinates(AuxVertex::DecayVertex& vert);

  void ClearData();
  void GetTrackShowerVectors(art::Event const & evt,
                             std::vector<recob::PFParticle const*>& pandora_primaryPFP,
                             std::vector<recob::PFParticle const*>& tracks,
                             std::vector<recob::PFParticle const*>& showers);
  void GetOriginVertices(art::Event const & evt,
                         const std::vector<recob::PFParticle const*>& tracks,
                         const std::vector<recob::PFParticle const*>& showers,
                         std::vector<AuxVertex::DecayVertex>& trackVertices,
                         std::vector<AuxVertex::DecayVertex>& showerVertices);
  void GetDecayVertices(const std::vector<AuxVertex::DecayVertex>& trackVertices,
                        const std::vector<AuxVertex::DecayVertex>& showerVertices,
                        std::vector<AuxVertex::DecayVertex>& potVertices,
                        std::vector<AuxVertex::DecayVertex>& cleanVertices);
  void GetHitVectors(art::Event const & evt,
                     const std::vector<recob::PFParticle const*>& tracks,
                     const std::vector<recob::PFParticle const*>& showers,
                     std::vector<recob::Hit const*>& totHits,
                     std::vector<std::vector<recob::Hit const*>>& trackHits,
                     std::vector<std::vector<recob::Hit const*>>& showerHits,
                     int& nTrackHits,
                     int& nShowerHits);
}; // End class ScanRecoSelectionParameters

ScanRecoSelectionParameters::ScanRecoSelectionParameters(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fInstanceName(pset.get<std::string>("InstanceName")),
    fIteration(pset.get<int>("Iteration")),
    fMinTpcBound(pset.get<std::vector<double>>("MinTpcBound")),
    fMaxTpcBound(pset.get<std::vector<double>>("MaxTpcBound")),
    fDataType(pset.get<std::string>("DataType")),
    fPfpLabel(pset.get<std::string>("PfpLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fDecayChannel(pset.get<int>("DecayChannel")),
    fSterileMass(pset.get<double>("SterileMass")),
    fDistanceCut(pset.get<double>("DistanceCut")),
    fCaloCut(pset.get<double>("CaloCut")),
    fChannelNorm(pset.get<double>("ChannelNorm")),
    fTickNorm(pset.get<double>("TickNorm")),
    fPrimaryOnly(pset.get<bool>("PrimaryOnly")),
    fEndVerticesAlso(pset.get<bool>("EndVerticesAlso")),
    fAnaType(pset.get<std::string>("AnalysisType"))
{} // END constructor ScanRecoSelectionParameters

ScanRecoSelectionParameters::~ScanRecoSelectionParameters()
{} // END destructor ScanRecoSelectionParameters

void ScanRecoSelectionParameters::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  metaTree = tfs->make<TTree>("Metadata","");
  metaTree->Branch("instanceName",&fInstanceName);
  metaTree->Branch("iteration",&fIteration,"iteration/I");  
  metaTree->Branch("dataType",&fDataType);
  metaTree->Branch("pfpLabel",&fPfpLabel);
  metaTree->Branch("hitLabel",&fHitLabel);
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

  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
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

void ScanRecoSelectionParameters::GetTrackShowerVectors(art::Event const & evt, std::vector<recob::PFParticle const*>& pandora_primaryPFP, std::vector<recob::PFParticle const*>&tracks, std::vector<recob::PFParticle const*>& showers)
{
  // Take an event (evt), then fill the (empty!) track, shower and primaries pfp vectors (tracks,showers,pandora_primaryPFP) with all the PFParticles that are tracks, showers and primaries pfps respectively in the event. These vectors will then be used to determine all the possible vertices.

  // Prepare handle labels
  art::InputTag pfpTag {fPfpLabel};

  // Find associations from PFP to tracks and showers
  const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
  for (auto const& pfp : (*pfpHandle)){
    if (pfp.IsPrimary()) {pandora_primaryPFP.push_back(&pfp);}
    else
    {
      // Find current pfp parent and check if parent is primary (so that pfp are only secondaries, and not tertiaries, etc.), then fill the correspondending vector with the type of particle found
      if (fPrimaryOnly){
        auto const& parent = (*pfpHandle).at(pfp.Parent());
        if (parent.IsPrimary()){
          if (pfp.PdgCode()==13) {tracks.push_back(&pfp);}
          if (pfp.PdgCode()==11) {showers.push_back(&pfp);}
        }
      }
      // Fill the track and shower vectors with the products associated with the pfp (in any case, regardless of their "genealogy")
      else {
        if (pfp.PdgCode()==13) {tracks.push_back(&pfp);}
        if (pfp.PdgCode()==11) {showers.push_back(&pfp);}
      }
    }
  } // End of pfp loop
  return;
} // END function GetTrackShowerVectors

void ScanRecoSelectionParameters::GetOriginVertices(art::Event const & evt, const std::vector<recob::PFParticle const*>& tracks, const std::vector<recob::PFParticle const*>& showers, std::vector<AuxVertex::DecayVertex>& trackVertices, std::vector<AuxVertex::DecayVertex>& showerVertices)
{
  // Given track and shower vectors (tracks,showers), fill vertices vectors (trackVertices,showerVertices) taking all the start and end points for tracks and start points for showers. The last two numbers (j,j) indicate only their idx in the xxxVertices vector (which is reduntant at this stage). However, later on, it will be used to identify the vertices uniquely, and the idxs from two vertices will be used to define the parent idx of the mean vertex created between them.

  // Initialize associations
  art::InputTag pfpTag {fPfpLabel};
  art::FindMany<recob::Track> pta(tracks,evt,pfpTag);
  art::FindMany<recob::Shower> psa(showers,evt,pfpTag);

  // Perform for tracks
  for(std::vector<int>::size_type i=0; i!=tracks.size(); i++)
  {
    // Annoying stupid way of getting associated object
    std::vector<recob::Track const*> vTrack;
    pta.get(i,vTrack);
    if (vTrack.size()==1)
    {
      recob::Track const* track = vTrack[0];
      // Create vertex
      AuxVertex::DecayVertex tempV1(track->Vertex().X(),track->Vertex().Y(),track->Vertex().Z(),i,i,"t","t","front","front");
      SetDetectorCoordinates(tempV1);
      trackVertices.push_back(tempV1);
      if (fEndVerticesAlso)
      {
        AuxVertex::DecayVertex tempV2(track->End().X(),track->End().Y(),track->End().Z(),i,i,"t","t","end","end");
        SetDetectorCoordinates(tempV2);
        trackVertices.push_back(tempV2);
      }
    }
    else {printf("WHAT THE HELL! Why does the association doesn't lead to an object?\nIt looks like one track/object might be missing, don't trust this event.\n");}
  }


  // And for showers
  for(std::vector<int>::size_type i=0; i!=showers.size(); i++)
  {
    // Annoying stupid way of getting associated object
    std::vector<recob::Shower const*> vShower;
    psa.get(i,vShower);
    if (vShower.size()==1)
    {
      recob::Shower const* shower = vShower[0];
      // Create vertex
      AuxVertex::DecayVertex tempV(shower->ShowerStart().X(),shower->ShowerStart().Y(),shower->ShowerStart().Z(),i,i,"s","s","front","front");
      SetDetectorCoordinates(tempV);
      showerVertices.push_back(tempV);
    }
    else {printf("WHAT THE HELL! Why does the association doesn't lead to an object?\nIt looks like one track/object might be missing, don't trust this event.\n");}
  }
  return;
} // END function GetOriginVertices

void ScanRecoSelectionParameters::GetDecayVertices(const std::vector<AuxVertex::DecayVertex>& trackVertices, const std::vector<AuxVertex::DecayVertex>& showerVertices, std::vector<AuxVertex::DecayVertex>& potVertices, std::vector<AuxVertex::DecayVertex>& cleanVertices)
{
  // Given the reco vertices establishing the origin of tracks and showers (trackVertices,showerVertices), take all the possible mean vertices by paring them two by two (potVertices), then impose that in radius (fDistanceCut) only two original vertices are present (cleanVertices). Now the previously used (j,j) index are used to define the parent vertices (which track or shower originating vertex was used to determine the mean vertex between them).
  // For track-track
  if (fAnaType == "tt")
  {
    for(std::vector<int>::size_type i=0; i!=trackVertices.size(); i++)
    {
      for(std::vector<int>::size_type j=i+1; j!=trackVertices.size(); j++)
      {
        AuxVertex::DecayVertex v1 = trackVertices[i];
        AuxVertex::DecayVertex v2 = trackVertices[j];
        float distance = Distance(v1,v2);
        pairDistance.push_back(distance);
        bool isInRadius = (distance<fDistanceCut);
        if (isInRadius)
        {
          AuxVertex::DecayVertex v3 = MeanVertex(v1, v2);
          SetDetectorCoordinates(v3);
          potVertices.push_back(v3);
          potPairDistance.push_back(distance);
        }
      }
    }
    // Make sure the potential mean vertices are clean (two decay vertices only in radius)
    // (start with assumption that isGoodVertex == true, then if you find extra particles in radius, turn it in bad vertex. Only good vertices are saved)
    for(std::vector<int>::size_type i=0; i!=potVertices.size(); i++)
    {
      AuxVertex::DecayVertex mv = potVertices[i];
      bool isGoodVertex = true;
      for(std::vector<int>::size_type j=0; j!=trackVertices.size(); j++)
      {
        AuxVertex::DecayVertex v1 = trackVertices[j];
        bool isParent1 = (mv.GetParIdx1() == v1.GetParIdx1());
        bool isParent2 = (mv.GetParIdx2() == v1.GetParIdx2());
        bool notParent = !(isParent1 || isParent2);
        bool isInRadius = (Distance(mv,v1)<fDistanceCut);
        if (isInRadius && notParent) isGoodVertex = false;
      }
      if (isGoodVertex) cleanVertices.push_back(mv);
    }
  }

  // And for track-shower
  else if (fAnaType == "ts"){
    for(std::vector<int>::size_type i=0; i!=trackVertices.size(); i++)
    {
      for(std::vector<int>::size_type j=0; j!=showerVertices.size(); j++)
      {
        AuxVertex::DecayVertex v1 = trackVertices[i];
        AuxVertex::DecayVertex v2 = showerVertices[j];
        float distance = Distance(v1,v2);
        pairDistance.push_back(distance);
        bool isInRadius = (distance<fDistanceCut);
        if (isInRadius)
        {
          AuxVertex::DecayVertex v3 = MeanVertex(v1, v2);
          SetDetectorCoordinates(v3);
          potVertices.push_back(v3);
          potPairDistance.push_back(distance);
        }
      }
    }
    // Make sure the potential mean vertices are clean (two decay vertices only in radius)
    // (start with assumption that isGoodVertex == true, then if you find extra particles in radius, turn it in bad vertex. Only good vertices are saved)
    for(std::vector<int>::size_type i=0; i!=potVertices.size(); i++)
    {
      AuxVertex::DecayVertex mv = potVertices[i];
      bool isGoodVertex = true;
      for(std::vector<int>::size_type j=0; j!=trackVertices.size(); j++)
      {
        AuxVertex::DecayVertex v1 = trackVertices[j];
        bool notParent = (mv.GetParIdx1() != v1.GetParIdx1());
        bool isInRadius = (Distance(mv,v1)<fDistanceCut);
        if (isInRadius && notParent) isGoodVertex = false;
      }
      for(std::vector<int>::size_type j=0; j!=showerVertices.size(); j++)
      {
        AuxVertex::DecayVertex v2 = showerVertices[j];
        bool notParent = (mv.GetParIdx2() != v2.GetParIdx2());
        bool isInRadius = (Distance(mv,v2)<fDistanceCut);
        if (isInRadius && notParent) isGoodVertex = false;
      }
      if (isGoodVertex) cleanVertices.push_back(mv);
    }
  }
  // Throw error if fAnaType wasn't right.
  else {
    throw std::invalid_argument("Invalid fAnaType. Must be 'tt' or 'ts'!");
  }
  return;
} // END function GetDecayVertices

void ScanRecoSelectionParameters::GetHitVectors(art::Event const & evt, const std::vector<recob::PFParticle const*>& tracks, const std::vector<recob::PFParticle const*>& showers, std::vector<recob::Hit const*>& totHits, std::vector<std::vector<recob::Hit const*>>& trackHits, std::vector<std::vector<recob::Hit const*>>& showerHits, int& nTrackHits, int& nShowerHits)
{
  // Prepare pfpTag
  art::InputTag pfpTag {fPfpLabel};
  art::FindMany<recob::Track> pta(tracks,evt,pfpTag);
  art::FindMany<recob::Shower> psa(showers,evt,pfpTag);
  std::vector<recob::Track const*> realTracks;
  std::vector<recob::Shower const*> realShowers;

  // Perform for tracks
  for(std::vector<int>::size_type i=0; i!=tracks.size(); i++)
  {
    // Annoying stupid way of getting associated object
    std::vector<recob::Track const*> vTrack;
    pta.get(i,vTrack);
    if (vTrack.size()==1)
    {
      recob::Track const* track = vTrack[0];
      realTracks.push_back(track);
    }
  }
  // Perform for showers
  for(std::vector<int>::size_type i=0; i!=showers.size(); i++)
  {
    // Annoying stupid way of getting associated object
    std::vector<recob::Shower const*> vShower;
    psa.get(i,vShower);
    if (vShower.size()==1){
      recob::Shower const* shower = vShower[0];
      realShowers.push_back(shower);
    }
  }

  // Find all hits
  art::Handle<std::vector<recob::Hit>> hitHandle;
  evt.getByLabel(fHitLabel, hitHandle);
  std::vector<recob::Hit> const& tempTotHits(*hitHandle);
  for (auto const& tempTotHit: tempTotHits) {totHits.push_back(&tempTotHit);}
  
  // Find hits associated with tracks
  art::FindMany<recob::Hit> tha(realTracks,evt,pfpTag);
  int tempNTrackHits = 0;
  for (std::vector<int>::size_type i=0; i!=realTracks.size(); i++)
  {
    std::vector<recob::Hit const*> tempHitVector;
    tha.get(i,tempHitVector);
    trackHits.push_back(tempHitVector);
    tempNTrackHits += tempHitVector.size();
  }

  // Find hits associated with showers
  art::FindMany<recob::Hit> sha(realShowers,evt,pfpTag);
  int tempNShowerHits = 0;
  for (std::vector<int>::size_type i=0; i!=realShowers.size(); i++)
  {
    std::vector<recob::Hit const*> tempHitVector;
    sha.get(i,tempHitVector);
    showerHits.push_back(tempHitVector);
    tempNShowerHits += tempHitVector.size();
  }
  nTrackHits = tempNTrackHits;
  nShowerHits = tempNShowerHits;
  return;
} // END function GetHitVectors


void ScanRecoSelectionParameters::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.

  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  printf("------------------------------------------------\n");
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

  // Prepare track, shower and Pandora primary vectors
  std::vector<recob::PFParticle const*> tracks;
  std::vector<recob::PFParticle const*> showers;
  std::vector<recob::PFParticle const*> pandora_primaryPFP;
  GetTrackShowerVectors(evt, pandora_primaryPFP, tracks, showers);
  printf("Loading %lu secondary tracks and %lu secondary showers.\n", tracks.size(), showers.size());

  // Determine origin vertices for tracks and showers
  std::vector<AuxVertex::DecayVertex> trackVertices;
  std::vector<AuxVertex::DecayVertex> showerVertices;
  GetOriginVertices(evt, tracks, showers, trackVertices, showerVertices);
  printf("Processing %lu origin vertices (%lu from tracks, %lu from showers).\n", trackVertices.size()+showerVertices.size(), trackVertices.size(), showerVertices.size());

  // Determine potential and clean decay vertices satisfying selection
  std::vector<AuxVertex::DecayVertex> potVertices;
  std::vector<AuxVertex::DecayVertex> cleanVertices;
  GetDecayVertices(trackVertices, showerVertices, potVertices, cleanVertices);
  printf("Selecting %lu potential vertices, %lu of which are clean.\n", potVertices.size(), cleanVertices.size());

  std::vector<recob::Hit const*> totHits;
  std::vector<std::vector<recob::Hit const*>> trackHits;
  std::vector<std::vector<recob::Hit const*>> showerHits;
  int nTrackHits;
  int nShowerHits;
  GetHitVectors(evt, tracks, showers, totHits, trackHits, showerHits, nTrackHits, nShowerHits);
  printf("Loading %lu hits (%i from %lu secondary tracks, %i from %lu secondary showers).\n", totHits.size(), nTrackHits, tracks.size(), nShowerHits, showers.size());


  for (std::vector<int>::size_type i=0; i!=cleanVertices.size(); i++)
  {
    int startWire[3] = {0,2399,4798};
    int channel0[3] = {cleanVertices[i].GetChannelLoc(0),cleanVertices[i].GetChannelLoc(1),cleanVertices[i].GetChannelLoc(2)};
    double tick0[3] = {cleanVertices[i].GetTickLoc(0),cleanVertices[i].GetTickLoc(1),cleanVertices[i].GetTickLoc(2)};
    int parIdx1 = cleanVertices[i].GetParIdx1();
    int parIdx2 = cleanVertices[i].GetParIdx2();
    cleanVertices[i].PrintInformation();

    printf("|_Clean vertex %lu\n", i);
    printf("  |_Hits from track 1\n");
    for (auto hit : trackHits[parIdx1])
    {
      int hitChannel = hit->Channel();
      double hitTick = (hit->EndTick() + hit->StartTick())/2.;
      int hitPlane = hit->View();
      bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(fCaloCut,2.));
      if (isInsideRadius)
      {
        printf("  | |_P%i Inside: %i -> Center: [%i,%.1f]; Hit:[%i,%.1f]; Distance: %.1f\n", hitPlane, isInsideRadius, channel0[hitPlane]-startWire[hitPlane], tick0[hitPlane], hitChannel-startWire[hitPlane], hitTick, sqrt(pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.)));
      }
    }
    printf("  |_Hits from track 2\n");
    for (auto hit : trackHits[parIdx2])
    {
      int hitChannel = hit->Channel();
      double hitTick = (hit->EndTick() + hit->StartTick())/2.;
      int hitPlane = hit->View();
      bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(fCaloCut,2.));
      if (isInsideRadius)
      {
        printf("    |_P%i Inside: %i -> Center: [%i,%.1f]; Hit:[%i,%.1f]; Distance: %.1f\n", hitPlane, isInsideRadius, channel0[hitPlane]-startWire[hitPlane], tick0[hitPlane], hitChannel-startWire[hitPlane], hitTick, sqrt(pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.)));
      }
    }

  }

  // Calculate final tree quantities
  nTracks = tracks.size();
  nShowers = showers.size();
  nPairs = pairDistance.size();
  nTrackVertices = trackVertices.size();
  nShowerVertices = showerVertices.size();
  nPotVertices = potVertices.size();
  nCleanVertices = cleanVertices.size();
  pandora_nPrimaryVertices = pandora_primaryPFP.size();

  tTree->Fill();
  printf("------------------------------------------------\n");
} // END function analyze


// Geometry functions
void ScanRecoSelectionParameters::SetDetectorCoordinates(AuxVertex::DecayVertex& vert)
{
  // Get spatial coordinates and mark vertex as assigned
  double xyz[3] = {vert.GetX(),vert.GetY(),vert.GetZ()};
  vert.SetIsDetLocAssigned(true);

  // Check whether coordinates are inside TPC 
  bool isInsideX = (xyz[0]>fMinTpcBound[0]+fDistanceCut &&
    xyz[0]<fMaxTpcBound[0]-fDistanceCut);
  bool isInsideY = (xyz[1]>fMinTpcBound[1]+fDistanceCut &&
    xyz[1]<fMaxTpcBound[1]-fDistanceCut);
  bool isInsideZ = (xyz[2]>fMinTpcBound[2]+fDistanceCut &&
    xyz[2]<fMaxTpcBound[2]-fDistanceCut);

  // If vertex is inside TPC, determine channel/tick coordinates and assign them
  if (isInsideX && isInsideY && isInsideZ)
  {
    vert.SetIsInsideTPC(true);
    raw::ChannelID_t channel0 = fGeometry->NearestChannel(xyz,0);
    raw::ChannelID_t channel1 = fGeometry->NearestChannel(xyz,1);
    raw::ChannelID_t channel2 = fGeometry->NearestChannel(xyz,2);
    double tick0 = fDetectorProperties->ConvertXToTicks(xyz[0], 0, 0, 0);
    double tick1 = fDetectorProperties->ConvertXToTicks(xyz[0], 1, 0, 0);
    double tick2 = fDetectorProperties->ConvertXToTicks(xyz[0], 2, 0, 0);
    vert.SetChannelLoc(channel0, channel1, channel2);
    vert.SetTickLoc(tick0, tick1, tick2);
    return;
  }

  // Else flag it as outside the TPC and exit function
  else {vert.SetIsInsideTPC(false); return;}
} // END function SetDetectorCoordinates

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ScanRecoSelectionParameters)

#endif // END def ScanRecoSelectionParameters_module