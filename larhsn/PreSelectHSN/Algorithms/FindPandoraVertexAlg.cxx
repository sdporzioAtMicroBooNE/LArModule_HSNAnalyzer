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
    fDistanceCut = pset.get<double>("DistanceCut");
    fPrimaryOnly = pset.get<bool>("PrimaryOnly");
    fEndVerticesAlso = pset.get<bool>("EndVerticesAlso");
    fAnaType = pset.get<std::string>("AnalysisType");
    fVerbose = pset.get<bool>("VerboseMode");
  }



  // Find each neutrino, and associated daughter. For each neutrino, fill every vector with the pfp_neutrino and vectors of pfp_tracks and pfp_showers that are its daughters.
  void FindPandoraVertexAlg::GetPotentialNeutrinoVertices(art::Event const & evt)
  {
    // Clear vectors that will be returned by function
    ana_pandora_neutrinos.clear();
    ana_pandora_tracks.clear();
    ana_pandora_showers.clear();
    ana_pandora_decayVertices.clear();
    tree_pandora_nNeutrinos = 0;
    tree_pandora_nTwoProngedNeutrinos = 0;
    tree_pandora_nInsideTwoProngedNeutrinos = 0;
    tree_pandora_neutrinoPdgCode.clear();
    tree_pandora_neutrinoNumDaughters.clear();
    tree_pandora_neutrinoNumTracks.clear();
    tree_pandora_neutrinoNumShowers.clear();
    tree_pandora_neutrinoInTPC.clear();

    // Prepare the pfp handle and the vertex associations
    art::InputTag pfpTag {fPfpLabel};
    const auto& pfpHandle = evt.getValidHandle< std::vector<recob::PFParticle> >(pfpTag);
    art::FindMany<recob::Vertex> pva(pfpHandle,evt,pfpTag);

    // Loop through each pfp

    for(std::vector<int>::size_type i=0; i!=(*pfpHandle).size(); i++)
    {
      recob::PFParticle pfp = (*pfpHandle)[i];

      // Find neutrino
      if (pfp.IsPrimary())
      {
        // Save pointer to pfp neutrino
        ana_pandora_neutrinos.push_back(&pfp);
        // Useful variables for the tree
        tree_pandora_nNeutrinos += 1;
        tree_pandora_neutrinoPdgCode.push_back(pfp.PdgCode());
        tree_pandora_neutrinoNumDaughters.push_back(pfp.NumDaughters());
        int thisNeutrino_numTracks = 0;
        int thisNeutrino_numShowers = 0;
        // The indices of the tracks belonging to this neutrino in the newly forming tracks vector
        std::vector<int> trackIndices;
        // The indices of the tracks belonging to this neutrino in the pfp original vector
        std::vector<int> trackIndicesInPfpHandle;

        // Diagnostic message
        if (fVerbose) printf("\n Found neutrino %i in this event.\n", tree_pandora_nNeutrinos); 

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
        tree_pandora_neutrinoNumTracks.push_back(thisNeutrino_numTracks);
        tree_pandora_neutrinoNumShowers.push_back(thisNeutrino_numShowers);
        // Diagnostic message
        if (fVerbose) printf(" with %i daughters: %i tracks and %i showers.\n", pfp.NumDaughters(),thisNeutrino_numTracks, thisNeutrino_numShowers); 

        // If this neutrino contains two and only two tracks we can create a specific decay vertex for it (to use later for calorimetry), but first we have to find the vertex data product associated with it.
        if (thisNeutrino_numTracks==2)
        {
          tree_pandora_nTwoProngedNeutrinos += 1;
          // Get all the necessary associated vertices via the stupid way
          std::vector<recob::Vertex const*> nuVertices, t1Vertices, t2Vertices;
          pva.get(i,nuVertices);
          pva.get(trackIndicesInPfpHandle[0],t1Vertices);
          pva.get(trackIndicesInPfpHandle[1],t2Vertices);
          printf(" Neutrino has %i vertices.\n", (int) nuVertices.size());
          recob::Vertex const * nuVertex = nuVertices[0];
          recob::Vertex const * t1Vertex = t1Vertices[0];
          recob::Vertex const * t2Vertex = t2Vertices[0];

          // Get vertex and daughters coordinates
          double nuVertexPosition[3], t1VertexPosition[3], t2VertexPosition[3];
          nuVertex->XYZ(nuVertexPosition);
          t1Vertex->XYZ(t1VertexPosition);
          t2Vertex->XYZ(t2VertexPosition);
          // double x1 = ana_pandora_tracks[trackIndices[0]].

          AuxVertex::DecayVertex nuV(nuVertexPosition[0],nuVertexPosition[1],nuVertexPosition[2],trackIndices[0],trackIndices[1],"t","t","front","front");
          nuV.SetParXYZ(0, t1VertexPosition[0], t1VertexPosition[1], t1VertexPosition[2]);
          nuV.SetParXYZ(1, t2VertexPosition[0], t2VertexPosition[1], t2VertexPosition[2]) ;
          nuV.SetDetectorCoordinates(fMinTpcBound,fMaxTpcBound,fDistanceCut,fGeometry,fDetectorProperties);
          if (nuV.IsInsideTPC())
          {
            tree_pandora_neutrinoInTPC.push_back(true);
            tree_pandora_nInsideTwoProngedNeutrinos += 1;
            ana_pandora_decayVertices.push_back(nuV);
          }
          else tree_pandora_neutrinoInTPC.push_back(false);
          if (fVerbose) printf(" |_ Vertex Coordinates: [%.1f,%.1f,%.1f] [%i,%i,%i,%.1f]\n",nuV.GetX(),nuV.GetY(),nuV.GetZ(),nuV.GetChannelLoc(0),nuV.GetChannelLoc(1),nuV.GetChannelLoc(2),nuV.GetTickLoc(0));
        } // END if neutrino has 2 tracks
      } // END if pfp is a neutrino
    } // END loop for each pfp
  } // END function GetOrderedPFParticles

} // END namespace FindPandoraVertex