#include "FindPandoraVertexAlg.h"

namespace FindPandoraVertex
{
  // Constructor/destructor
  FindPandoraVertexAlg::FindPandoraVertexAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  FindPandoraVertexAlg::~FindPandoraVertexAlg()
  {}

  void FindPandoraVertexAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fMinTpcBound = pset.get<std::vector<double>>("MinTpcBound");
    fMaxTpcBound = pset.get<std::vector<double>>("MaxTpcBound");
    fPfpLabel = pset.get<std::string>("PfpLabel");
    fVerbose = pset.get<bool>("VerboseMode");
  }


  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void FindPandoraVertexAlg::GetPotentialNeutrinoVertices(
            art::Event const & evt,
            AuxEvent::EventDescriptor & evd,
            std::vector<AuxVertex::DecayVertex> & ana_decayVertices)
  {
    if (fVerbose) printf("\n--- GetPotentialNeutrinoVertices message ---\n");

    // Clear stuff that will be modified by function
    evd.nNeutrinos = 0;
    evd.nTwoProngedNeutrinos = 0;
    evd.nContainedTwoProngedNeutrinos = 0;
    evd.neutrinoPdgCode.clear();
    evd.neutrinoNumDaughters.clear();
    evd.neutrinoNumTracks.clear();
    evd.neutrinoNumShowers.clear();
    evd.diag_nuWithMissingAssociatedVertex = 0;
    evd.diag_nuWithMissingAssociatedTrack = 0;
    evd.diag_nuProngWithMissingAssociatedHits = 0;


    auto const& trk_handle = evt.getValidHandle<std::vector<recob::Track>>("pandoraNu");
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, evt, "pandoraNu"); // Track
    //Prepare the pfp handle and the vertex/tracks associations
    art::InputTag pfpTag {fPfpLabel};
    const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
    // art::FindMany<recob::Vertex> pva(pfpHandle,evt,pfpTag);
    // art::FindMany<recob::Track> pta(pfpHandle,evt,pfpTag);

    // Loop through each pfp
    for(std::vector<int>::size_type i=0; i!=(*pfpHandle).size(); i++)
    {
      recob::PFParticle pfp = (*pfpHandle)[i];
      // Find out if this pfp is a neutrino
      if (pfp.IsPrimary())
      {
        // Fill useful variables for the tree
        evd.nNeutrinos += 1;
        evd.neutrinoPdgCode.push_back(pfp.PdgCode());
        evd.neutrinoNumDaughters.push_back(pfp.NumDaughters());
        int thisNeutrino_numTracks = 0;
        int thisNeutrino_numShowers = 0;
        // ID of current neutrino and vectors where we will fill pointers to daughters (only tracks for now)
        size_t nuID = pfp.Self();
        std::vector<recob::PFParticle const *> thisNeutrino_pfpTrackPointers;

        // Diagnostic message
        if (fVerbose)
        {
          printf("Neutrino %i (ID: %i, PDG: %i)\n", evd.nNeutrinos, (int) nuID, pfp.PdgCode());
          printf("|_Number of daughters: %i (ID:", pfp.NumDaughters());

          // Prepare vector of ID of neutrino daughters
          auto nuDaughtersID = pfp.Daughters();
          // Loop through each daughter and print their ID
          for (std::vector<int>::size_type j=0; j!=nuDaughtersID.size(); j++)
          {
            printf(" %i", (int) nuDaughtersID[j]);
          }
          printf(" )\n");
        }

        // Loop through each pfp and check if their parent was the neutrino we are currently looping through
        for (auto const& daughter_pfp : (*pfpHandle))
        {
          // Find correct pfp corresponding to daughter we are analyzing
          if (daughter_pfp.Parent()==nuID)
          {
            // Separate in track and shower pfps and save their pointers to corresponding vectors
            if (daughter_pfp.PdgCode()==13)
            {
              if (fVerbose) printf("| |_Found track with ID: %i\n", (int) daughter_pfp.Self());
              thisNeutrino_numTracks += 1;
              thisNeutrino_pfpTrackPointers.push_back(&daughter_pfp);
            }
            if (daughter_pfp.PdgCode()==11)
            {
              if (fVerbose) printf("| |_Found shower with ID: %i\n", (int) daughter_pfp.Self());
              thisNeutrino_numShowers += 1;
            }
          }
        }
        // Save to vector number of tracks and showers found
        evd.neutrinoNumTracks.push_back(thisNeutrino_numTracks);
        evd.neutrinoNumShowers.push_back(thisNeutrino_numShowers);

        // Diagnostic message
        if (fVerbose)
        {
          // Cross check, loop through each pointer, make sure their ID is correct and their parent is as well
          for (auto const& pfpTrack : (thisNeutrino_pfpTrackPointers))
          {
            printf("| |_Checking saved daughter with ID %i and parent ID %i\n", (int) pfpTrack->Self(), (int) pfpTrack->Parent());
          }
        }
        // Diagnostic message
        if (fVerbose)
        {
          printf("|_Summary: %i daughters, %i tracks and %i showers.\n", pfp.NumDaughters(),thisNeutrino_numTracks, thisNeutrino_numShowers);
        }

        // If this neutrino contains two and only two tracks we can create a specific decay vertex for it (to use later for calorimetry), but first we have to make sure we have all the associations we need.
        if (thisNeutrino_numTracks==2)
        {
          evd.nTwoProngedNeutrinos += 1;
          if (fVerbose) printf("|_Neutrino is potential candidate n. %i in event.\n", evd.nTwoProngedNeutrinos);

          // Creating a stupid vector with a pointer to the current neutrino, all other methods won't work so I have to go through this stupid way of retrieving associations
          std::vector<recob::PFParticle const *> thisNeutrino_pfpNeutrinoPointer;
          thisNeutrino_pfpNeutrinoPointer.push_back(&pfp);

          // Association to vertex and tracks objects of our pfpTracks
          art::FindMany<recob::Vertex> nu_pva(thisNeutrino_pfpNeutrinoPointer,evt,pfpTag);
          art::FindMany<recob::Vertex> pva(thisNeutrino_pfpTrackPointers,evt,pfpTag);
          art::FindMany<recob::Track> pta(thisNeutrino_pfpTrackPointers,evt,pfpTag);


          // Make sure we have the necessary vertices associated with the pfp
          std::vector<recob::Vertex const*> nuVertices, t1Vertices, t2Vertices;
          pva.get(0,nuVertices);
          pva.get(0,t1Vertices);
          pva.get(1,t2Vertices);
          bool rightNumVertices = (nuVertices.size()==1 && t1Vertices.size()==1 && t2Vertices.size()==1);

          // Make sure we have the necessary tracks associated with the pfps
          std::vector<recob::Track const*> t1Tracks, t2Tracks;
          pta.get(0,t1Tracks);
          pta.get(1,t2Tracks);
          bool rightNumTracks = (t1Tracks.size()==1 && t2Tracks.size()==1);

          if (!rightNumVertices) evd.diag_nuWithMissingAssociatedVertex += 1;
          if (!rightNumTracks) evd.diag_nuWithMissingAssociatedTrack += 1;
          if (rightNumVertices && rightNumTracks)
          {
            if (fVerbose) printf("| | |_Neutrino has correct number of vertex and tracks associated to PFP.\n");
            // Grab all the pointers to vertices/tracks
            recob::Vertex const * nuVertex = nuVertices[0];
            recob::Vertex const * t1Vertex = t1Vertices[0];
            recob::Vertex const * t2Vertex = t2Vertices[0];
            recob::Track const * t1Track = t1Tracks[0];
            recob::Track const * t2Track = t2Tracks[0];
            std::vector<recob::Track const*> thisNu_tracks = {t1Track, t2Track};

            // Make sure also we have the necessary hits associated to tracks
            art::FindMany<recob::Hit> tha(thisNu_tracks,evt,pfpTag);
            std::vector<recob::Hit const *> t1Hits, t2Hits;
            tha.get(0,t1Hits);
            tha.get(1,t2Hits);
            bool rightNumHits = (t1Hits.size()>1 && t2Hits.size()>1);
            std::cout << "| |_Track 1: There are " << t1Hits.size() << " associated hits." << std::endl;
            std::cout << "| |_Track 2: There are " << t2Hits.size() << " associated hits." << std::endl;

            if (!rightNumHits) evd.diag_nuProngWithMissingAssociatedHits += 1;
            else
            {
              if (fVerbose) printf("| | |_Neutrino has correct number of hits vectors associated to tracks.\n");
              // Time to dump all associations in the neutrino vertex
              AuxVertex::DecayVertex nuV(nuVertex,t1Vertex,t2Vertex,t1Track,t2Track,t1Hits,t2Hits);
              nuV.SetDetectorCoordinates(fMinTpcBound,fMaxTpcBound,fGeometry,fDetectorProperties);
              nuV.PrintInformation();
              if (nuV.IsInsideTPC())
              {
                evd.nContainedTwoProngedNeutrinos += 1;
                ana_decayVertices.push_back(nuV);
              }
            } // END if right number of hits associations
          } // END if right number of vertices / tracks associations
        } // END if neutrino has 2 tracks
      } // END if pfp is a neutrino
    } // END loop for each pfp
  } // END function GetOrderedPFParticles

} // END namespace FindPandoraVertex
