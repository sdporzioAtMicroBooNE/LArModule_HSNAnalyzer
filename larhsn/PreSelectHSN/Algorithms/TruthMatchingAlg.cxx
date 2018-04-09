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

    // Clean vectors that will be returned by function
    evd.prong_matchedPDG.clear();




    //auto const& hit_handle = evt.getValidHandle<std::vector<recob::Hit>>("gaushit");
    auto const& hit_handle = evt.getValidHandle<std::vector<recob::Hit>>( fHitLabel);
    auto const& trk_handle = evt.getValidHandle<std::vector<recob::Track>>("pandoraNu");
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, evt, "pandoraNu"); // Track


    // Loop through each decay vertex
    for (std::vector<int>::size_type i=0; i!=decayVertices.size(); i++)
    {
      AuxVertex::DecayVertex currentVertex = decayVertices[i];
      // Get useful quantities about the decay vertex currently being analyzed, like coordinates and parent indices
      //int channel0[3] = {currentVertex.GetChannelLoc(0),currentVertex.GetChannelLoc(1),currentVertex.GetChannelLoc(2)};
      //float tick0[3] = {currentVertex.GetTickLoc(0),currentVertex.GetTickLoc(1),currentVertex.GetTickLoc(2)};

      // Initialize the vector of hits that will be pushed back to the vector of vector of hits and returned by the function
      //std::vector<recob::Hit const*> totHitsInMaxRadius;
      // Get the vector of hits for the prongs
      std::vector<recob::Hit const*> trackHits1 = currentVertex.GetProngHits(0);
      std::vector<recob::Hit const*> trackHits2 = currentVertex.GetProngHits(1);
      std::vector<int> pdgs (2);

      for(size_t j = 0; j<hits_per_track.size() ; j++){
        std::unordered_map<int,double> trkide;

        double maxe=-1, tote=0;
        simb::MCParticle const* maxp_me = nullptr; //pointer for the particle match we will calculate

        std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(j);

        art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle, evt,"gaushitTruthMatch");
        //art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit1(trackHits1, evt,"gaushitTruthMatch");

        std::vector<simb::MCParticle const*> particle_vec;

        std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
        for(size_t i_h=0; i_h<trk_hits_ptrs.size(); ++i_h){


        particle_vec.clear(); match_vec.clear();

        particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);

        //the .key() gives us the index in the original collection

        //std::cout << "\t\tThere are " << particle_vec.size() << " particles matched to hit " << i_h << std::endl;

        //loop over particles
        for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){

          trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
          tote += match_vec[i_p]->energy; //calculate total energy deposited

          if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
            maxe = trkide[ particle_vec[i_p]->TrackId() ];
            maxp_me = particle_vec[i_p];
          }
        } //end loop over particles per hit

      } // End loop over hits
      std::cout << "Final Match (from my loop) is " << maxp_me->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
                  << " \n\tpdg=" << maxp_me->PdgCode()
                  << " trkid=" << maxp_me->TrackId()
                  << " ke=" << maxp_me->E()-maxp_me->Mass()
                  << "\n\tstart (x,y,z)=(" << maxp_me->Vx()
                  << "," << maxp_me->Vy()
                  << "," << maxp_me->Vz()
                  << ")\tend (x,y,z)=(" << maxp_me->EndX()
                  << "," << maxp_me->EndY()
                  << "," << maxp_me->EndZ() << ")" << std::endl;


      if (trk_hits_ptrs.size()==trackHits1.size()) pdgs.at(0) = maxp_me->PdgCode();
      if (trk_hits_ptrs.size()==trackHits2.size()) pdgs.at(1) = maxp_me->PdgCode();
    } //end loop over tracsk

      evd.prong_matchedPDG.push_back(pdgs);
      pdgs.clear();
    } // END loop for each decay vertex
    return;
  } // END function TruthMatching
} // END namespace
