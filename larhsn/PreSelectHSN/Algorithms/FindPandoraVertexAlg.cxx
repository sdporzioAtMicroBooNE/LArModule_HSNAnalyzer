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
            AuxEvent::EventDescriptor & evd,
            art::Event const & evt)
  {
    // Clear vectors that will be returned by function
    ana_pandora_neutrinos.clear();
    ana_pandora_tracks.clear();
    ana_pandora_showers.clear();
    ana_pandora_decayVertices.clear();
    evd.pandora_nNeutrinos = 0;
    evd.pandora_nTwoProngedNeutrinos = 0;
    evd.pandora_nContainedTwoProngedNeutrinos = 0;
    evd.pandora_neutrinoPdgCode.clear();
    evd.pandora_neutrinoNumDaughters.clear();
    evd.pandora_neutrinoNumTracks.clear();
    evd.pandora_neutrinoNumShowers.clear();
    evd.pandora_diag_potentialPfpsWithMissingVertex = 0;

    // Prepare the pfp handle and the vertex associations
    art::InputTag pfpTag {fPfpLabel};
    const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
    art::FindMany<recob::Vertex> pva(pfpHandle,evt,pfpTag);

    // Loop through each pfp
    for(std::vector<int>::size_type i=0; i!=(*pfpHandle).size(); i++)
    {
      recob::PFParticle pfp = (*pfpHandle)[i];
      // Find out if this pfp is a neutrino
      if (pfp.IsPrimary())
      {
        // Save pointer to pfp neutrino
        ana_pandora_neutrinos.push_back(&pfp);
        // Fill useful variables for the tree
        evd.pandora_nNeutrinos += 1;
        evd.pandora_neutrinoPdgCode.push_back(pfp.PdgCode());
        evd.pandora_neutrinoNumDaughters.push_back(pfp.NumDaughters());
        int thisNeutrino_numTracks = 0;
        int thisNeutrino_numShowers = 0;
        // The indices of the tracks belonging to this neutrino in the newly forming tracks vector
        std::vector<int> trackIndices;
        // The indices of the tracks belonging to this neutrino in the pfp original vector
        std::vector<int> trackIndicesInPfpHandle;

        // Prepare vector of ID of neutrino daughters
        auto nuDaughtersID = pfp.Daughters();
        // Loop through each daughter ID
        for (std::vector<int>::size_type j=0; j!=nuDaughtersID.size(); j++)
        {
          size_t id = nuDaughtersID[j];
          // Loop through each pfp
          for (auto const& daughter_pfp : (*pfpHandle))
          {
            // Find correct pfp corresponding to daughter we are analyzing
            if (daughter_pfp.Self()==id)
            {
              // Separate in track and shower pfps and save their pointers to corresponding vectors
              if (daughter_pfp.PdgCode()==13)
              {
                thisNeutrino_numTracks += 1;
                ana_pandora_tracks.push_back(&daughter_pfp);
                trackIndices.push_back(ana_pandora_tracks.size()-1);
                trackIndicesInPfpHandle.push_back(j);
              }
              if (daughter_pfp.PdgCode()==11)
              {
                thisNeutrino_numShowers += 1;
                ana_pandora_showers.push_back(&daughter_pfp);
              }
            }
          }
        }
        evd.pandora_neutrinoNumTracks.push_back(thisNeutrino_numTracks);
        evd.pandora_neutrinoNumShowers.push_back(thisNeutrino_numShowers);

        // Diagnostic message
        if (fVerbose)
        {
          printf("\n--- GetPotentialNeutrinoVertices message ---\n");
          printf("Found neutrino %i in this event with %i daughters: %i tracks and %i showers.\n", evd.pandora_nNeutrinos, pfp.NumDaughters(),thisNeutrino_numTracks, thisNeutrino_numShowers);
        }


        // If this neutrino contains two and only two tracks we can create a specific decay vertex for it (to use later for calorimetry), but first we have to find the vertex data product associated with it.
        if (thisNeutrino_numTracks==2)
        {
          evd.pandora_nTwoProngedNeutrinos += 1;
          if (fVerbose) printf("Neutrino is potential candidate n. %i in event.\n", evd.pandora_nTwoProngedNeutrinos);
          // Get all the necessary associated vertices via the stupid way
          std::vector<recob::Vertex const*> nuVertices, t1Vertices, t2Vertices;
          pva.get(i,nuVertices);
          pva.get(trackIndicesInPfpHandle[0],t1Vertices);
          pva.get(trackIndicesInPfpHandle[1],t2Vertices);
          std::vector<int> thisNu_nVerticesInPfp = {(int) nuVertices.size(),(int) t1Vertices.size(),(int) t2Vertices.size()};
          bool rightNumVertices = (thisNu_nVerticesInPfp[0]==1 &&
                                  thisNu_nVerticesInPfp[1]==1 &&
                                  thisNu_nVerticesInPfp[2]==1);
          if (rightNumVertices)
          {
            if (fVerbose) printf("Neutrino has correct number of vertex associations.\n");
            recob::Vertex const * nuVertex = nuVertices[0];
            recob::Vertex const * t1Vertex = t1Vertices[0];
            recob::Vertex const * t2Vertex = t2Vertices[0];

            // Get vertex and daughters coordinates
            double nuVertexPosition[3], t1VertexPosition[3], t2VertexPosition[3];
            nuVertex->XYZ(nuVertexPosition);
            t1Vertex->XYZ(t1VertexPosition);
            t2Vertex->XYZ(t2VertexPosition);

            AuxVertex::DecayVertex nuV(nuVertexPosition[0],nuVertexPosition[1],nuVertexPosition[2],trackIndices[0],trackIndices[1],"t","t","front","front");
            nuV.SetParXYZ(0, t1VertexPosition[0], t1VertexPosition[1], t1VertexPosition[2]);
            nuV.SetParXYZ(1, t2VertexPosition[0], t2VertexPosition[1], t2VertexPosition[2]) ;
            nuV.SetDetectorCoordinates(fMinTpcBound,fMaxTpcBound,fGeometry,fDetectorProperties);
            if (nuV.IsInsideTPC())
            {
              evd.pandora_nContainedTwoProngedNeutrinos += 1;
              ana_pandora_decayVertices.push_back(nuV);
            }
          }
          else
          {
            if (fVerbose) printf("ERROR! Neutrino has not correct number of vertex associations.\n");
            evd.pandora_diag_potentialPfpsWithMissingAssociatedVertex += 1;
          } // END if all vertices are present
        } // END if neutrino has 2 tracks
      } // END if pfp is a neutrino
    } // END loop for each pfp
  } // END function GetOrderedPFParticles

} // END namespace FindPandoraVertex