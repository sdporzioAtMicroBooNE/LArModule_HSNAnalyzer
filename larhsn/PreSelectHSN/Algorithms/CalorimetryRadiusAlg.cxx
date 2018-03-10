#include "CalorimetryRadiusAlg.h"

namespace CalorimetryRadius
{
  // Constructor/destructor
  CalorimetryRadiusAlg::CalorimetryRadiusAlg(fhicl::ParameterSet const & pset)
  {
    reconfigure(pset);
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  }
  CalorimetryRadiusAlg::~CalorimetryRadiusAlg()
  {}

  void CalorimetryRadiusAlg::reconfigure(fhicl::ParameterSet const & pset)
  {
    fPfpLabel = pset.get<std::string>("PfpLabel");
    fHitLabel = pset.get<std::string>("HitLabel");
    fVerbose = pset.get<bool>("VerboseMode");
    fRadiusProfileLimits = pset.get<std::vector<double>>("RadiusProfileLimits");
    fRadiusProfileBins = pset.get<int>("RadiusProfileBins");
    fChannelNorm = pset.get<double>("ChannelNorm");
    fTickNorm = pset.get<double>("TickNorm");

    // Determine profile ticks
    double profileStep = (fRadiusProfileLimits[1] - fRadiusProfileLimits[0]) / float(fRadiusProfileBins);
    double currTick = fRadiusProfileLimits[0];
    for (int i=0; i<fRadiusProfileBins; i++)
    {
      currTick += profileStep;
      profileTicks.push_back(currTick);
    }
  }

  // Fill three vectors, one containing all the hits in the event (totHits), one containing vectors of all the hits pertaining to tracks (trackHits), one for each track, and finally one containing vectors of all the hits from showers (showerHits), one for each shower. These hits will be used to perform calorimetry analysis.
  void CalorimetryRadiusAlg::GetHitVectors(
          AuxEvent::EventDescriptor & evd,
          std::vector<AuxVertex::DecayVertex>& cleanVertices,
          art::Event const & evt,
          const std::vector<recob::PFParticle const*>& tracks,
          const std::vector<recob::PFParticle const*>& showers)
  {

    // Clear vectors that will be returned by function
    ana_calo_totHits.clear();
    ana_calo_trackHits.clear();
    ana_calo_showerHits.clear();
    tree_calo_NumTotHits = 0;

    for (std::vector<int>::size_type i=0; i!=cleanVertices.size(); i++)
    {
      AuxVertex::DecayVertex currentVertex = cleanVertices[i];
      int parIdx1 = cleanVertices[i].GetParIdx1();
      int parIdx2 = cleanVertices[i].GetParIdx2();

    }


    // Prepare pfpTag
    art::InputTag pfpTag {fPfpLabel};
    art::FindMany<recob::Track> pta(tracks,evt,pfpTag);
    art::FindMany<recob::Shower> psa(showers,evt,pfpTag);
    std::vector<recob::Track const*> realTracks;
    std::vector<recob::Shower const*> realShowers;

    // Get the actual recob::Track and recob::Shower objects from the tracks and shower arrays (which instead contain only recob::PFParticles).
    // Hits are associated with tracks/showers objects, not with PFParticles, which is why we must follow this annoying chain of associations (recob::PFParticle->recob::Track/Shower->recob::Hit)

    // Find recob::Tracks
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
      else
      {
        recob::Track const* track;
        realTracks.push_back(track);
      }
      
    }
    // Find recob::Showers
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
    for (auto const& tempTotHit: tempTotHits) {ana_calo_totHits.push_back(&tempTotHit);}
    
    // Find hits associated with recob::Tracks
    art::FindMany<recob::Hit> tha(realTracks,evt,pfpTag);
    for (std::vector<int>::size_type i=0; i!=realTracks.size(); i++)
    {
      std::vector<recob::Hit const*> tempHitVector;
      tha.get(i,tempHitVector);
      ana_calo_trackHits.push_back(tempHitVector);
      tree_calo_NumTotHits += tempHitVector.size();
    }

    // Find hits associated with recob::Showers
    art::FindMany<recob::Hit> sha(realShowers,evt,pfpTag);
    for (std::vector<int>::size_type i=0; i!=realShowers.size(); i++)
    {
      std::vector<recob::Hit const*> tempHitVector;
      sha.get(i,tempHitVector);
      ana_calo_showerHits.push_back(tempHitVector);
      tree_calo_NumTotHits += tempHitVector.size();
    }

    // Diagnostic message
    if (fVerbose) {printf("Loading %i hits.\n",tree_calo_NumTotHits);}
    return;
  } // END function GetHitVectors




  // Perform calorimetry analysis. At this stage we finally calculate all the charge deposited by the hits of track1 and track2 (or shower) within a radius from the assumed HSN decay vertex (cleanVertices[i]), for each decay vertex.
  // The second step looks at all the charge deposited by any hit in radius (which may come from hadronic interaction, in the case of background). And we finally calculate the ratio between the two (caloRatio). We would expect this ratio to be closer to 1 for signal, since HSN decaying in the detector don't interact with any particle, and we'd expect charge deposited by the two decay products to be the only charge within a certain radius from the decay point.
  // Now, we actually repeat this step for different radia in order to build up a profile. The width of the profile is given by fRadiusProfileLimits and the number of bins by fRadiusProfileBin.
  void CalorimetryRadiusAlg::PerformCalorimetry(
          AuxEvent::EventDescriptor & evd,
          std::vector<AuxVertex::DecayVertex>& cleanVertices,
          const std::vector<recob::Hit const*>& totHits,
          const std::vector<std::vector<recob::Hit const*>>& trackHits,
          const std::vector<std::vector<recob::Hit const*>>& showerHits)
  {

    // Clean vectors that will be returned by function
    ana_calo_totHitsInMaxRadius.clear();
    tree_calo_par1ChargeInRadius.clear();
    tree_calo_par2ChargeInRadius.clear();
    tree_calo_totChargeInRadius.clear();
    tree_calo_caloRatio.clear();

    // Loop through each clean vertex
    int vertInd = 0;
    for (std::vector<int>::size_type i=0; i!=cleanVertices.size(); i++)
    {
      // Get useful quantities about the clean vertex currently being analyzed, like coordinates and parent indices
      int channel0[3] = {cleanVertices[i].GetChannelLoc(0),cleanVertices[i].GetChannelLoc(1),cleanVertices[i].GetChannelLoc(2)};
      float tick0[3] = {cleanVertices[i].GetTickLoc(0),cleanVertices[i].GetTickLoc(1),cleanVertices[i].GetTickLoc(2)};
      int parIdx1 = cleanVertices[i].GetParIdx1();
      int parIdx2 = cleanVertices[i].GetParIdx2();
      if (fVerbose) {cleanVertices[i].PrintInformation();}

      // Make sure that you can work with this neutrino and it has all the hits you need, otherwise declare it as pathological
      int nTrackHits1 = (int) trackHits[parIdx1].size();      
      int nTrackHits2 = (int) trackHits[parIdx2].size();
      bool enoughHits = (nTrackHits1>0 && nTrackHits2>0);
      if (!enoughHits)
        {
          printf("ERROR!\n");
          printf("Track indices: [%i, %i]\n", parIdx1, parIdx2);
          printf("Number hits: [%i, %i]\n", nTrackHits1, nTrackHits2);
          cleanVertices[i].SetIsPathological(true,4);
        }
      else
      {
        // Initialize the vector of hits that will be pushed back to the vector of vector of hits and returned by the function
        std::vector<recob::Hit const*> totHitsInMaxRadius_thisDv;

        // Calculate calorimetry for track 1 within radius
        // parCharge1 is a vector, which contains all the charge due to particle1 in a circle of radius r around the vertex.
        // Each element of the vector is that integrated charge in increasing value of r
        // It starts out as a vector of size equal to the number of bins in radius profile, each element is equal to 0.
        // A loop goes then through each hit and for each radius size asks whether the hit is in it. If it is, the charge gets added to the total.
        // Now declare parCharge1 and fill it with zeros
        std::vector<float> parCharge1;
        for (int j=0; j<fRadiusProfileBins; j++) parCharge1.push_back(0.);
        printf("len trackHits1 %i\n", (int) trackHits[parIdx1].size());
        for (auto hit : trackHits[parIdx1])
        {
          int hitChannel = hit->Channel();
          double hitTick = (hit->EndTick() + hit->StartTick())/2.;
          int hitPlane = hit->View();
          for (int j=0; j<fRadiusProfileBins; j++)
          {
            double caloCut = profileTicks[j];
            bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(caloCut,2.));
            if (isInsideRadius)
            {
              double hitCharge = hit->Integral();
              parCharge1[j] += hitCharge;
            }
          }
        }
        // Calculate calorimetry for track 2 within radius
        std::vector<float> parCharge2;
        for (int j=0; j<fRadiusProfileBins; j++) parCharge2.push_back(0.);
        printf("Track indices: [%i, %i]\n", parIdx1, parIdx2);
        printf("Number hits: [%i, %i]\n", nTrackHits1, nTrackHits2);
        printf("len trackHits2 %lu\n", trackHits[parIdx2].size());
        for (auto hit : trackHits[parIdx2])
        {
          int hitChannel = hit->Channel();
          double hitTick = (hit->EndTick() + hit->StartTick())/2.;
          int hitPlane = hit->View();
          for (int j=0; j<fRadiusProfileBins; j++)
          {
            double caloCut = profileTicks[j];
            bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(caloCut,2.));
            if (isInsideRadius)
            {
              double hitCharge = hit->Integral();
              parCharge2[j] += hitCharge;
            }
          }
        }

        // Calculate total calorimetry within radius
        std::vector<float> totCharge;
        for (int j=0; j<fRadiusProfileBins; j++) totCharge.push_back(0.);
        printf("len totHits %i\n", (int) totHits.size());
        for (auto hit : totHits)
        {
          int hitChannel = hit->Channel();
          double hitTick = (hit->EndTick() + hit->StartTick())/2.;
          int hitPlane = hit->View();
          for (int j=0; j<fRadiusProfileBins; j++)
          {
            double caloCut = profileTicks[j];
            bool isInsideRadius = (pow(((hitChannel-channel0[hitPlane])/fChannelNorm),2.) + pow(((hitTick-tick0[hitPlane])/fTickNorm),2.) < pow(caloCut,2.));
            if (isInsideRadius)
            {
              double hitCharge = hit->Integral();
              totCharge[j] += hitCharge;

              // totHitsInMaxRadius are used to draw the evd, you need to do that only for the largest radius
              if (j == fRadiusProfileBins-1) totHitsInMaxRadius_thisDv.push_back(hit);
            }
          } 
        }
        ana_calo_totHitsInMaxRadius.push_back(totHitsInMaxRadius_thisDv);
        totHitsInMaxRadius_thisDv.clear();
        std::vector<float> thisCaloRatio;
        for (int j=0; j<fRadiusProfileBins; j++) thisCaloRatio.push_back((parCharge1[j]+parCharge2[j])/float(totCharge[j]));

        // Calculate tree quantities
        tree_calo_par1ChargeInRadius.push_back(parCharge1);
        tree_calo_par2ChargeInRadius.push_back(parCharge2);
        tree_calo_totChargeInRadius.push_back(totCharge);
        tree_calo_caloRatio.push_back(thisCaloRatio);

        vertInd++;
      } // END if enough hits
    } // END clean vertex loop
    return;
  } // END function PerformCalorimetry



} // END namespace CalorimetryRadius