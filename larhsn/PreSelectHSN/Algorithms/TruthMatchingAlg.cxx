#include "TruthMatchingAlg.h"

namespace TruthMatching
{
  // Constructor/destructor
  TruthMatchingAlg::TruthMatchingAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  TruthMatchingAlg::~TruthMatchingAlg()
  {}

  void TruthMatchingAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fPfpLabel = pset.get<std::string>("PfpLabel");
    fHitLabel = pset.get<std::string>("HitLabel");
    fVerbose = pset.get<bool>("VerboseMode");
  }

  void TruthMatchingAlg::PerformTruthMatching(
          art::Event const & evt,
          AuxEvent::EventDescriptor & evd,
          std::vector<AuxVertex::DecayVertex>& decayVertices)
  {
    // // Clear vectors that will be returned by function
    // evd.match_pdgCode.clear();
    // evd.match_mass.clear();
    // evd.match_energy.clear();
    // evd.match_prong1StartPosition.clear();
    // evd.match_prong2StartPosition.clear();
    // evd.match_prong1Momentum.clear();
    // evd.match_prong2Momentum.clear();

    // auto const& hitHandle = evt.getValidHandle<std::vector<recob::Hit>>(fHitLabel);
    // auto const& trackHandle = evt.getValidHandle<std::vector<recob::Track>>(fPfpLabel);
    // art::FindManyP<recob::Hit> hta(trackHandle, evt, fPfpLabel);

    // // Loop through each decay vertex
    // for (std::vector<int>::size_type i=0; i!=decayVertices.size(); i++)
    // {
    //   // Get current decay vertex and associated hits
    //   AuxVertex::DecayVertex currentVertex = decayVertices[i];
    //   std::vector<recob::Hit const*> trackHits1 = currentVertex.GetProngHits(0);
    //   std::vector<recob::Hit const*> trackHits2 = currentVertex.GetProngHits(1);
    //   // MCTruth information about this decay vertex
    //   std::vector<int> thisNu_pdgCode = {-1,-1};
    //   std::vector<float> thisNu_mass = {-1,-1};
    //   std::vector<float> thisNu_energy = {-1,-1};
    //   std::vector<float> thisNu_prong1StartPosition = {-1,-1,-1};
    //   std::vector<float> thisNu_prong2StartPosition = {-1,-1,-1};
    //   std::vector<float> thisNu_prong1Momentum = {-1,-1};
    //   std::vector<float> thisNu_prong2Momentum = {-1,-1};

    //   // Loop over any recob::track particle in event
    //   for(size_t j = 0; j<hta.size() ; j++)
    //   {
    //     /* The aim of this loop is to create a map for this particular recob::track
    //         We are going to look at every hit from this recob::track, and for each hit,
    //         find the deposited energy that each MCParticle might have contributed.
    //         We create a map [MCParticle track id, energy deposited by this MCParticle towards the track energy]
    //         The element of the map with the highest contributing energy must be the MCParticle associated with the track.
    //     */
    //     // Get vector of recob::hits for this specific recob::track
    //     std::vector<art::Ptr<recob::Hit>> hitsInTrack_v = hta.at(j);

    //     // Total energy deposited for each mcParticle id (id <-> energy contributed)
    //     std::unordered_map<int,double> trackidEnergy_map;
    //     double maximumDepositedEnergy = -1, totalDepositedEnergy = 0;

    //     // Empty pointer to the particle that MC particle that will be matched with the recob::track
    //     simb::MCParticle const* matchedMc = nullptr; 

    //     // Association of each hits to the vectors of mcparticles and matching information
    //     art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> hmba(hitHandle, evt, "gaushitTruthMatch");
    //     // Loop over each recob::hit in the current recob::track
    //     for(size_t i_h=0; i_h<hitsInTrack_v.size(); ++i_h)
    //     {
    //       // Vectors of all the mcparticles and matching information associated with this hit
    //       std::vector<simb::MCParticle const*> mcParticles_v;
    //       std::vector<anab::BackTrackerHitMatchingData const*> matchingStruct_v;

    //       // Fill the mcParticles and matchingStruct vector with the ones associated with this hit
    //       hmba.get(hitsInTrack_v[i_h].key(),mcParticles_v,matchingStruct_v);

    //       // Loop over MCParticle/matchingStruct associated with the current hit
    //       for(size_t i_p=0; i_p<mcParticles_v.size(); ++i_p){
    //         // Get track id of this specific mcparticle
    //         int trackID = mcParticles_v[i_p]->TrackId();
    //         // Add to the (trackID,energy) map the energy we are currently finding
    //         trackidEnergy_map[trackID] += matchingStruct_v[i_p]->energy;
    //         // Sum the contribution to the energy of this recob::track from all the MCParticles
    //         totalDepositedEnergy += matchingStruct_v[i_p]->energy;
    //         // Keep track of maximum value of energy deposited
    //         if( trackidEnergy_map[trackID] > maximumDepositedEnergy )
    //         {
    //           maximumDepositedEnergy = trackidEnergy_map[trackID];
    //           matchedMc = mcParticles_v[i_p];
    //         }
    //       } // End loop over each MCParticle
    //     } // End loop over each reconstructed hits

    //     std::cout << "Final Match (from my loop) is " << matchedMc->TrackId() << " with energy " << maximumDepositedEnergy << " over " << totalDepositedEnergy << " (" << maximumDepositedEnergy/totalDepositedEnergy << ")"
    //                 << " \n\tpdg=" << matchedMc->PdgCode()
    //                 << " trkid=" << matchedMc->TrackId()
    //                 << " ke=" << matchedMc->E()-matchedMc->Mass()
    //                 << "\n\tstart (x,y,z)=(" << matchedMc->Vx()
    //                 << "," << matchedMc->Vy()
    //                 << "," << matchedMc->Vz()
    //                 << ")\tend (x,y,z)=(" << matchedMc->EndX()
    //                 << "," << matchedMc->EndY()
    //                 << "," << matchedMc->EndZ() << ")" << std::endl;

    //     if (hitsInTrack_v.size()==trackHits1.size())
    //     {
    //       thisNu_pdgCode[0] = matchedMc->PdgCode();
    //       thisNu_mass[0] = matchedMc->Mass();
    //       thisNu_energy[0] = matchedMc->E();
    //       thisNu_prong1StartPosition = 
    //       {
    //         (float) matchedMc->Vx(),
    //         (float) matchedMc->Vy(),
    //         (float) matchedMc->Vz()
    //       };
    //       thisNu_prong1Momentum = 
    //       {
    //         (float) matchedMc->Px(),
    //         (float) matchedMc->Py(),
    //         (float) matchedMc->Pz()
    //       };
    //     }

    //     if (hitsInTrack_v.size()==trackHits2.size())
    //     {
    //       thisNu_pdgCode[1] = matchedMc->PdgCode();
    //       thisNu_mass[1] = matchedMc->Mass();
    //       thisNu_energy[1] = matchedMc->E();
    //       thisNu_prong2StartPosition = 
    //       {
    //         (float) matchedMc->Vx(),
    //         (float) matchedMc->Vy(),
    //         (float) matchedMc->Vz()
    //       };
    //       thisNu_prong2Momentum = 
    //       {
    //         (float) matchedMc->Px(),
    //         (float) matchedMc->Py(),
    //         (float) matchedMc->Pz()
    //       };
    //     }
    //   } //end loop over tracks
    //   evd.match_pdgCode.push_back(thisNu_pdgCode);
    //   evd.match_mass.push_back(thisNu_mass);
    //   evd.match_energy.push_back(thisNu_energy);
    //   evd.match_prong1StartPosition.push_back(thisNu_prong1Momentum);
    //   evd.match_prong2StartPosition.push_back(thisNu_prong2StartPosition);
    //   evd.match_prong1Momentum.push_back(thisNu_prong1Momentum);
    //   evd.match_prong2Momentum.push_back(thisNu_prong2Momentum);   
    // } // END loop for each decay vertex
  return;
  } // END function TruthMatching
} // END namespace
