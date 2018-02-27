#ifndef PreSelectHSN_module
#define PreSelectHSN_module

#include "PreSelectHSN.h"

PreSelectHSN::PreSelectHSN(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fFindDecayVertexAlg(pset),
    fFindPandoraVertexAlg(pset),
    fCalorimetryRadiusAlg(pset),
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
    fRadiusProfileLimits(pset.get<std::vector<double>>("RadiusProfileLimits")),
    fRadiusProfileBins(pset.get<int>("RadiusProfileBins")),
    fChannelNorm(pset.get<double>("ChannelNorm")),
    fTickNorm(pset.get<double>("TickNorm")),
    fPrimaryOnly(pset.get<bool>("PrimaryOnly")),
    fEndVerticesAlso(pset.get<bool>("EndVerticesAlso")),
    fAnaType(pset.get<std::string>("AnalysisType")),
    fManualSearch(pset.get<bool>("ManualSearch")),
    fPandoraSearch(pset.get<bool>("PandoraSearch")),
    fVerbose(pset.get<bool>("VerboseMode")),
    fSaveDrawTree(pset.get<bool>("SaveDrawTree"))
{
  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

    // Determine profile ticks
    double profileStep = (fRadiusProfileLimits[1] - fRadiusProfileLimits[0]) / float(fRadiusProfileBins);
    double currTick = fRadiusProfileLimits[0];
    for (int i=0; i<fRadiusProfileBins; i++)
    {
      currTick += profileStep;
      profileTicks.push_back(currTick);
    }

} // END constructor PreSelectHSN

PreSelectHSN::~PreSelectHSN()
{} // END destructor PreSelectHSN

void PreSelectHSN::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  // Meta tree containing fcl file parameters
  metaTree = tfs->make<TTree>("Metadata","");
  metaTree->Branch("instanceName",&fInstanceName);
  metaTree->Branch("iteration",&fIteration,"iteration/I"); 
  metaTree->Branch("minTpcBound",&fMinTpcBound);
  metaTree->Branch("maxTpcBound",&fMaxTpcBound);  
  metaTree->Branch("dataType",&fDataType);
  metaTree->Branch("pfpLabel",&fPfpLabel);
  metaTree->Branch("hitLabel",&fHitLabel);
  metaTree->Branch("decayChannel",&fDecayChannel,"decayChannel/I");
  metaTree->Branch("distanceCut",&fDistanceCut,"distanceCut/D");
  metaTree->Branch("radiusProfileLimits",&fRadiusProfileLimits);
  metaTree->Branch("radiusProfileBins",&fRadiusProfileBins);
  metaTree->Branch("channelNorm",&fChannelNorm,"channelNorm/D");
  metaTree->Branch("tickNorm",&fTickNorm,"tickNorm/D");
  metaTree->Branch("sterileMass",&fSterileMass,"sterileMass/D");
  metaTree->Branch("primaryOnly",&fPrimaryOnly,"primaryOnly/O");
  metaTree->Branch("endVerticesAlso",&fEndVerticesAlso,"endVerticesAlso/O");
  metaTree->Branch("anaType",&fAnaType);
  metaTree->Branch("saveDrawTree",&fSaveDrawTree,"saveDrawTree/O");
  metaTree->Fill();

  // Tree for manual search
  if (fManualSearch)
  {
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
    tTree->Branch("nTotHits",&nTotHits);
    tTree->Branch("nTrackHits",&nTrackHits);
    tTree->Branch("nShowerHits",&nShowerHits);
    tTree->Branch("profileTicks",&profileTicks);
    tTree->Branch("totChargeInRadius",&totChargeInRadius);
    tTree->Branch("par1ChargeInRadius",&par1ChargeInRadius);
    tTree->Branch("par2ChargeInRadius",&par2ChargeInRadius);
    tTree->Branch("caloRatio",&caloRatio);

    if (fSaveDrawTree)
    {
      drawTree = tfs->make<TTree>("DrawTree","");
      drawTree->Branch("dv_xyzCoordinates",&dv_xyzCoordinates);
      drawTree->Branch("dv_wireCoordinates",&dv_wireCoordinates);
      drawTree->Branch("dv_tickCoordinates",&dv_tickCoordinates);
      drawTree->Branch("par1_xyzCoordinates",&par1_xyzCoordinates);
      drawTree->Branch("par1_wireCoordinates",&par1_wireCoordinates);
      drawTree->Branch("par1_tickCoordinates",&par1_tickCoordinates);
      drawTree->Branch("par2_xyzCoordinates",&par2_xyzCoordinates);
      drawTree->Branch("par2_wireCoordinates",&par2_wireCoordinates);
      drawTree->Branch("par2_tickCoordinates",&par2_tickCoordinates);
      drawTree->Branch("par1_hits_p0_wireCoordinates",&par1_hits_p0_wireCoordinates);
      drawTree->Branch("par1_hits_p1_wireCoordinates",&par1_hits_p1_wireCoordinates);
      drawTree->Branch("par1_hits_p2_wireCoordinates",&par1_hits_p2_wireCoordinates);
      drawTree->Branch("par1_hits_p0_tickCoordinates",&par1_hits_p0_tickCoordinates);
      drawTree->Branch("par1_hits_p1_tickCoordinates",&par1_hits_p1_tickCoordinates);
      drawTree->Branch("par1_hits_p2_tickCoordinates",&par1_hits_p2_tickCoordinates);
      drawTree->Branch("par2_hits_p0_wireCoordinates",&par2_hits_p0_wireCoordinates);
      drawTree->Branch("par2_hits_p1_wireCoordinates",&par2_hits_p1_wireCoordinates);
      drawTree->Branch("par2_hits_p2_wireCoordinates",&par2_hits_p2_wireCoordinates);
      drawTree->Branch("par2_hits_p0_tickCoordinates",&par2_hits_p0_tickCoordinates);
      drawTree->Branch("par2_hits_p1_tickCoordinates",&par2_hits_p1_tickCoordinates);
      drawTree->Branch("par2_hits_p2_tickCoordinates",&par2_hits_p2_tickCoordinates);
      drawTree->Branch("tot_hits_p0_wireCoordinates",&tot_hits_p0_wireCoordinates);
      drawTree->Branch("tot_hits_p1_wireCoordinates",&tot_hits_p1_wireCoordinates);
      drawTree->Branch("tot_hits_p2_wireCoordinates",&tot_hits_p2_wireCoordinates);
      drawTree->Branch("tot_hits_p0_tickCoordinates",&tot_hits_p0_tickCoordinates);
      drawTree->Branch("tot_hits_p1_tickCoordinates",&tot_hits_p1_tickCoordinates);
      drawTree->Branch("tot_hits_p2_tickCoordinates",&tot_hits_p2_tickCoordinates);
    }
  }

  // Tree for Pandora search
  if (fPandoraSearch)
  {
    pandoraTree = tfs->make<TTree>("PandoraData","");
    pandoraTree->Branch("run",&run,"run/I");
    pandoraTree->Branch("subrun",&subrun,"subrun/I");
    pandoraTree->Branch("event",&event,"event/I");
    pandoraTree->Branch("nNeutrinos",&tree_pandora_nNeutrinos);
    pandoraTree->Branch("nTwoProngedNeutrinos",&tree_pandora_nTwoProngedNeutrinos);
    pandoraTree->Branch("nInsideTwoProngedNeutrinos",&tree_pandora_nInsideTwoProngedNeutrinos);
    pandoraTree->Branch("neutrinoPdgCode",&tree_pandora_neutrinoPdgCode);
    pandoraTree->Branch("neutrinoInTpc",&tree_pandora_neutrinoInTPC);
    pandoraTree->Branch("neutrinoNumDaughters",&tree_pandora_neutrinoNumDaughters);
    pandoraTree->Branch("neutrinoNumTracks",&tree_pandora_neutrinoNumTracks);
    pandoraTree->Branch("neutrinoNumShowers",&tree_pandora_neutrinoNumShowers);
    pandoraTree->Branch("nTotHits",&tree_pandoraCalo_nTotHits);
    pandoraTree->Branch("nTrackHits",&tree_pandoraCalo_nTrackHits);
    pandoraTree->Branch("nShowerHits",&tree_pandoraCalo_nShowerHits);
    pandoraTree->Branch("profileTicks",&profileTicks);
    pandoraTree->Branch("totChargeInRadius",&tree_pandoraCalo_totChargeInRadius);
    pandoraTree->Branch("par1ChargeInRadius",&tree_pandoraCalo_par1ChargeInRadius);
    pandoraTree->Branch("par2ChargeInRadius",&tree_pandoraCalo_par2ChargeInRadius);
    pandoraTree->Branch("caloRatio",&tree_pandoraCalo_caloRatio);
    if (fSaveDrawTree)
    {
      pandoraDrawTree = tfs->make<TTree>("PandoraDrawTree","");
      pandoraDrawTree->Branch("dv_xyzCoordinates",&dv_xyzCoordinates);
      pandoraDrawTree->Branch("dv_wireCoordinates",&dv_wireCoordinates);
      pandoraDrawTree->Branch("dv_tickCoordinates",&dv_tickCoordinates);
      pandoraDrawTree->Branch("par1_xyzCoordinates",&par1_xyzCoordinates);
      pandoraDrawTree->Branch("par1_wireCoordinates",&par1_wireCoordinates);
      pandoraDrawTree->Branch("par1_tickCoordinates",&par1_tickCoordinates);
      pandoraDrawTree->Branch("par2_xyzCoordinates",&par2_xyzCoordinates);
      pandoraDrawTree->Branch("par2_wireCoordinates",&par2_wireCoordinates);
      pandoraDrawTree->Branch("par2_tickCoordinates",&par2_tickCoordinates);
      pandoraDrawTree->Branch("par1_hits_p0_wireCoordinates",&par1_hits_p0_wireCoordinates);
      pandoraDrawTree->Branch("par1_hits_p1_wireCoordinates",&par1_hits_p1_wireCoordinates);
      pandoraDrawTree->Branch("par1_hits_p2_wireCoordinates",&par1_hits_p2_wireCoordinates);
      pandoraDrawTree->Branch("par1_hits_p0_tickCoordinates",&par1_hits_p0_tickCoordinates);
      pandoraDrawTree->Branch("par1_hits_p1_tickCoordinates",&par1_hits_p1_tickCoordinates);
      pandoraDrawTree->Branch("par1_hits_p2_tickCoordinates",&par1_hits_p2_tickCoordinates);
      pandoraDrawTree->Branch("par2_hits_p0_wireCoordinates",&par2_hits_p0_wireCoordinates);
      pandoraDrawTree->Branch("par2_hits_p1_wireCoordinates",&par2_hits_p1_wireCoordinates);
      pandoraDrawTree->Branch("par2_hits_p2_wireCoordinates",&par2_hits_p2_wireCoordinates);
      pandoraDrawTree->Branch("par2_hits_p0_tickCoordinates",&par2_hits_p0_tickCoordinates);
      pandoraDrawTree->Branch("par2_hits_p1_tickCoordinates",&par2_hits_p1_tickCoordinates);
      pandoraDrawTree->Branch("par2_hits_p2_tickCoordinates",&par2_hits_p2_tickCoordinates);
      pandoraDrawTree->Branch("tot_hits_p0_wireCoordinates",&tot_hits_p0_wireCoordinates);
      pandoraDrawTree->Branch("tot_hits_p1_wireCoordinates",&tot_hits_p1_wireCoordinates);
      pandoraDrawTree->Branch("tot_hits_p2_wireCoordinates",&tot_hits_p2_wireCoordinates);
      pandoraDrawTree->Branch("tot_hits_p0_tickCoordinates",&tot_hits_p0_tickCoordinates);
      pandoraDrawTree->Branch("tot_hits_p1_tickCoordinates",&tot_hits_p1_tickCoordinates);
      pandoraDrawTree->Branch("tot_hits_p2_tickCoordinates",&tot_hits_p2_tickCoordinates);
    }

  }
} // END function beginJob

void PreSelectHSN::endJob()
{} // END function endJob

void PreSelectHSN::ClearData()
{
  // Manual variables
  nTracks = 0;
  nShowers= 0;
  nPairs = 0;
  nTrackVertices = 0;
  nShowerVertices = 0;
  nPotVertices = 0;
  nCleanVertices = 0;
  nCleanVerticesOutsideTPC = 0;
  pairDistance.clear();
  potPairDistance.clear();
  nTrackHits.clear();
  nShowerHits.clear();
  totChargeInRadius.clear();
  par1ChargeInRadius.clear();
  par2ChargeInRadius.clear();
  caloRatio.clear();

  // Pandora variables
  tree_pandoraCalo_nTotHits = 0;
  tree_pandoraCalo_nTrackHits.clear();
  tree_pandoraCalo_nShowerHits.clear();
  tree_pandoraCalo_par1ChargeInRadius.clear();
  tree_pandoraCalo_par2ChargeInRadius.clear();
  tree_pandoraCalo_totChargeInRadius.clear();
  tree_pandoraCalo_caloRatio.clear();

  // Drawing vectors
  dv_xyzCoordinates.clear();
  dv_wireCoordinates.clear();
  dv_tickCoordinates.clear();
  par1_xyzCoordinates.clear();
  par1_wireCoordinates.clear();
  par1_tickCoordinates.clear();
  par2_xyzCoordinates.clear();
  par2_wireCoordinates.clear();
  par2_tickCoordinates.clear();
  par1_hits_p0_wireCoordinates.clear();
  par1_hits_p1_wireCoordinates.clear();
  par1_hits_p2_wireCoordinates.clear();
  par1_hits_p0_tickCoordinates.clear();
  par1_hits_p1_tickCoordinates.clear();
  par1_hits_p2_tickCoordinates.clear();
  par2_hits_p0_wireCoordinates.clear();
  par2_hits_p1_wireCoordinates.clear();
  par2_hits_p2_wireCoordinates.clear();
  par2_hits_p0_tickCoordinates.clear();
  par2_hits_p1_tickCoordinates.clear();
  par2_hits_p2_tickCoordinates.clear();
  tot_hits_p0_wireCoordinates.clear();
  tot_hits_p1_wireCoordinates.clear();
  tot_hits_p2_wireCoordinates.clear();
  tot_hits_p0_tickCoordinates.clear();
  tot_hits_p1_tickCoordinates.clear();
  tot_hits_p2_tickCoordinates.clear();
} // END function ClearData

void PreSelectHSN::FillDrawTree(const std::vector<AuxVertex::DecayVertex>& cleanVertices, const std::vector<std::vector<recob::Hit const*>>& totHitsInMaxRadius, const std::vector<std::vector<recob::Hit const*>>& trackHits, const std::vector<std::vector<recob::Hit const*>>& showerHits)
{
  for (std::vector<int>::size_type i=0; i!=cleanVertices.size(); i++)
  {
    // Get clean vertex
    auto dv = cleanVertices[i];
    int parIdx1 = dv.GetParIdx1();
    int parIdx2 = dv.GetParIdx2();
    std::vector<recob::Hit const*> par1_hits = trackHits[parIdx1];
    std::vector<recob::Hit const*> par2_hits = trackHits[parIdx2];
    std::vector<recob::Hit const*> thisTot_hits = totHitsInMaxRadius[i];
    std::vector<float> thisPar1_hits_p0_tickCoordinates,
      thisPar1_hits_p1_tickCoordinates,
      thisPar1_hits_p2_tickCoordinates,
      thisPar2_hits_p0_tickCoordinates,
      thisPar2_hits_p1_tickCoordinates,
      thisPar2_hits_p2_tickCoordinates,
      thisTot_hits_p0_tickCoordinates,
      thisTot_hits_p1_tickCoordinates,
      thisTot_hits_p2_tickCoordinates;
    std::vector<int> thisPar1_hits_p0_wireCoordinates,
      thisPar1_hits_p1_wireCoordinates,
      thisPar1_hits_p2_wireCoordinates,
      thisPar2_hits_p0_wireCoordinates,
      thisPar2_hits_p1_wireCoordinates,
      thisPar2_hits_p2_wireCoordinates,
      thisTot_hits_p0_wireCoordinates,
      thisTot_hits_p1_wireCoordinates,
      thisTot_hits_p2_wireCoordinates;

    for (auto hit : par1_hits)
    {
      if (hit->View() == 0) {
        thisPar1_hits_p0_wireCoordinates.push_back(hit->Channel());
        thisPar1_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 1) {
        thisPar1_hits_p1_wireCoordinates.push_back(hit->Channel());
        thisPar1_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 2) {
        thisPar1_hits_p2_wireCoordinates.push_back(hit->Channel());
        thisPar1_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
    }

    for (auto hit : par2_hits)
    {
      if (hit->View() == 0) {
        thisPar2_hits_p0_wireCoordinates.push_back(hit->Channel());
        thisPar2_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 1) {
        thisPar2_hits_p1_wireCoordinates.push_back(hit->Channel());
        thisPar2_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 2) {
        thisPar2_hits_p2_wireCoordinates.push_back(hit->Channel());
        thisPar2_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
    }

    for (auto hit : thisTot_hits)
    {
      if (hit->View() == 0) {
        thisTot_hits_p0_wireCoordinates.push_back(hit->Channel());
        thisTot_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 1) {
        thisTot_hits_p1_wireCoordinates.push_back(hit->Channel());
        thisTot_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 2) {
        thisTot_hits_p2_wireCoordinates.push_back(hit->Channel());
        thisTot_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
    }

    // Get coordinates
    dv_xyzCoordinates.push_back({dv.GetX(),dv.GetY(), (float)dv.GetZ()});
    dv_wireCoordinates.push_back({dv.GetChannelLoc(0),dv.GetChannelLoc(1),dv.GetChannelLoc(2)});
    dv_tickCoordinates.push_back({dv.GetTickLoc(0),dv.GetTickLoc(1),dv.GetTickLoc(2)});
    par1_xyzCoordinates.push_back({dv.GetParX(0),dv.GetParY(0),dv.GetParZ(0)});
    par1_wireCoordinates.push_back({dv.GetParChannelLoc(0,0),dv.GetParChannelLoc(0,1),dv.GetParChannelLoc(0,2)});
    par1_tickCoordinates.push_back({dv.GetParTickLoc(0,0),dv.GetParTickLoc(0,1),dv.GetParTickLoc(0,2)});
    par2_xyzCoordinates.push_back({dv.GetParX(1),dv.GetParY(1),dv.GetParZ(1)});
    par2_wireCoordinates.push_back({dv.GetParChannelLoc(1,0),dv.GetParChannelLoc(1,1),dv.GetParChannelLoc(1,2)});
    par2_tickCoordinates.push_back({dv.GetParTickLoc(1,0),dv.GetParTickLoc(1,1),dv.GetParTickLoc(1,2)});

    par1_hits_p0_wireCoordinates.push_back(thisPar1_hits_p0_wireCoordinates);
    par1_hits_p1_wireCoordinates.push_back(thisPar1_hits_p1_wireCoordinates);
    par1_hits_p2_wireCoordinates.push_back(thisPar1_hits_p2_wireCoordinates);
    par1_hits_p0_tickCoordinates.push_back(thisPar1_hits_p0_tickCoordinates);
    par1_hits_p1_tickCoordinates.push_back(thisPar1_hits_p1_tickCoordinates);
    par1_hits_p2_tickCoordinates.push_back(thisPar1_hits_p2_tickCoordinates);

    par2_hits_p0_wireCoordinates.push_back(thisPar2_hits_p0_wireCoordinates);
    par2_hits_p1_wireCoordinates.push_back(thisPar2_hits_p1_wireCoordinates);
    par2_hits_p2_wireCoordinates.push_back(thisPar2_hits_p2_wireCoordinates);
    par2_hits_p0_tickCoordinates.push_back(thisPar2_hits_p0_tickCoordinates);
    par2_hits_p1_tickCoordinates.push_back(thisPar2_hits_p1_tickCoordinates);
    par2_hits_p2_tickCoordinates.push_back(thisPar2_hits_p2_tickCoordinates);

    tot_hits_p0_wireCoordinates.push_back(thisTot_hits_p0_wireCoordinates);
    tot_hits_p1_wireCoordinates.push_back(thisTot_hits_p1_wireCoordinates);
    tot_hits_p2_wireCoordinates.push_back(thisTot_hits_p2_wireCoordinates);
    tot_hits_p0_tickCoordinates.push_back(thisTot_hits_p0_tickCoordinates);
    tot_hits_p1_tickCoordinates.push_back(thisTot_hits_p1_tickCoordinates);
    tot_hits_p2_tickCoordinates.push_back(thisTot_hits_p2_tickCoordinates);
  }
  return;
} // END function FillDrawTree


// Core analysis. This is where all the functions are executed. Gets repeated event by event.
void PreSelectHSN::analyze(art::Event const & evt)
{
  if (fVerbose) {printf("\n------------------------------------------------\n");}

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  if (fVerbose) {printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);}


  // ----------------- MANUAL SEARCH ----------------- //
  if (fManualSearch)
  {
    // Start by clearing vectors from previous searches
    ClearData();
    // Find all PFParticles in the event. Separate them in track pfps, shower pfps, and primary pfps
    fFindDecayVertexAlg.GetTrackShowerVectors(evt);
    std::vector<recob::PFParticle const*> tracks = fFindDecayVertexAlg.ana_tracks;
    std::vector<recob::PFParticle const*> showers = fFindDecayVertexAlg.ana_showers;
    std::vector<recob::PFParticle const*> pandora_primaryPFP = fFindDecayVertexAlg.ana_pandora_primaryPFP;

    // Determine origin points for all tracks and for all showers
    fFindDecayVertexAlg.GetOriginVertices(evt, tracks, showers);
    std::vector<AuxVertex::DecayVertex> trackVertices = fFindDecayVertexAlg.ana_trackVertices;
    std::vector<AuxVertex::DecayVertex> showerVertices = fFindDecayVertexAlg.ana_showerVertices;

    // Use all vertices to determine potential decay vertices
    fFindDecayVertexAlg.GetDecayVertices(trackVertices, showerVertices);
    std::vector<AuxVertex::DecayVertex> potVertices = fFindDecayVertexAlg.ana_potVertices;
    std::vector<AuxVertex::DecayVertex> cleanVertices = fFindDecayVertexAlg.ana_cleanVertices;
    pairDistance = fFindDecayVertexAlg.ana_pairDistance;
    potPairDistance = fFindDecayVertexAlg.ana_potPairDistance;

    // Determine useful additional quantities that will go in the tree
    nTracks = tracks.size();
    nShowers = showers.size();
    nTrackVertices = trackVertices.size();
    nShowerVertices = showerVertices.size();
    nPairs = pairDistance.size();
    nPotVertices = potVertices.size();
    nCleanVertices = cleanVertices.size();

    // Just being smart, if there's no clean vertices, there's no reason to waste time loading hits that won't be used
    if (nCleanVertices==0) {printf("No clean vertex candidates found. Moving to next event...\n");}
    else
    {
      // Get vectors containing hits for each track/shower object in order to perform calorimetry
      fCalorimetryRadiusAlg.GetHitVectors(evt, tracks, showers);
      std::vector<recob::Hit const*> totHits = fCalorimetryRadiusAlg.ana_totHits;
      std::vector<std::vector<recob::Hit const*>> trackHits = fCalorimetryRadiusAlg.ana_trackHits;
      std::vector<std::vector<recob::Hit const*>> showerHits = fCalorimetryRadiusAlg.ana_showerHits;

      // Perform calorimetry analysis
      fCalorimetryRadiusAlg.PerformCalorimetry(cleanVertices, totHits, trackHits, showerHits);
      std::vector<std::vector<recob::Hit const*>> totHitsInMaxRadius = fCalorimetryRadiusAlg.ana_totHitsInMaxRadius;
      par1ChargeInRadius = fCalorimetryRadiusAlg.ana_par1ChargeInRadius;
      par2ChargeInRadius = fCalorimetryRadiusAlg.ana_par2ChargeInRadius;
      totChargeInRadius = fCalorimetryRadiusAlg.ana_totChargeInRadius;
      caloRatio = fCalorimetryRadiusAlg.ana_caloRatio;

      // Determine useful additional quantities that will go in the tree
      nTotHits = totHits.size();
      for (auto th : trackHits) nTrackHits.push_back(th.size());
      for (auto sh : showerHits) nShowerHits.push_back(sh.size());

      // Fill draw tree (if needed)
      if (fSaveDrawTree) FillDrawTree(cleanVertices, totHitsInMaxRadius, trackHits, showerHits);
    }
    // Fill tree
    tTree->Fill();
    if (fSaveDrawTree) drawTree->Fill();
  }


  // ----------------- PANDORA SEARCH ----------------- //  
  if (fPandoraSearch)
  {
    // Start by clearing vectors from previous searches
    ClearData();
    // Get potential neutrino vectors with two tracks
    fFindPandoraVertexAlg.GetPotentialNeutrinoVertices(evt);
    std::vector<recob::PFParticle const*> ana_pandora_neutrinos = fFindPandoraVertexAlg.ana_pandora_neutrinos;
    std::vector<recob::PFParticle const*> ana_pandora_tracks = fFindPandoraVertexAlg.ana_pandora_tracks;
    std::vector<recob::PFParticle const*> ana_pandora_showers = fFindPandoraVertexAlg.ana_pandora_showers;
    std::vector<AuxVertex::DecayVertex> ana_pandora_decayVertices = fFindPandoraVertexAlg.ana_pandora_decayVertices;
    tree_pandora_nNeutrinos = fFindPandoraVertexAlg.tree_pandora_nNeutrinos;
    tree_pandora_nTwoProngedNeutrinos = fFindPandoraVertexAlg.tree_pandora_nTwoProngedNeutrinos;
    tree_pandora_nInsideTwoProngedNeutrinos = fFindPandoraVertexAlg.tree_pandora_nInsideTwoProngedNeutrinos;
    tree_pandora_neutrinoPdgCode = fFindPandoraVertexAlg.tree_pandora_neutrinoPdgCode;
    tree_pandora_neutrinoNumDaughters = fFindPandoraVertexAlg.tree_pandora_neutrinoNumDaughters;
    tree_pandora_neutrinoNumTracks = fFindPandoraVertexAlg.tree_pandora_neutrinoNumTracks;
    tree_pandora_neutrinoNumShowers = fFindPandoraVertexAlg.tree_pandora_neutrinoNumShowers;
    tree_pandora_neutrinoInTPC = fFindPandoraVertexAlg.tree_pandora_neutrinoInTPC;

    // Just being smart, if there's no clean vertices, there's no reason to waste time loading hits that won't be used
    if (tree_pandora_nInsideTwoProngedNeutrinos==0) {printf("No clean vertex candidates found. Moving to next event...\n");}
    else
    {
      // Get vectors containing hits for each track/shower object in order to perform calorimetry
      fCalorimetryRadiusAlg.GetHitVectors(evt, ana_pandora_tracks, ana_pandora_showers);
      std::vector<recob::Hit const*> ana_pandoraCalo_totHits = fCalorimetryRadiusAlg.ana_totHits;
      std::vector<std::vector<recob::Hit const*>> ana_pandoraCalo_trackHits = fCalorimetryRadiusAlg.ana_trackHits;
      std::vector<std::vector<recob::Hit const*>> ana_pandoraCalo_showerHits = fCalorimetryRadiusAlg.ana_showerHits;
      tree_pandoraCalo_nTotHits = ana_pandoraCalo_totHits.size();
      for (auto th : ana_pandoraCalo_trackHits) tree_pandoraCalo_nTrackHits.push_back(th.size());
      for (auto sh : ana_pandoraCalo_showerHits) tree_pandoraCalo_nShowerHits.push_back(sh.size());

      // Perform calorimetry analysis
      fCalorimetryRadiusAlg.PerformCalorimetry(ana_pandora_decayVertices, ana_pandoraCalo_totHits, ana_pandoraCalo_trackHits, ana_pandoraCalo_showerHits);
      std::vector<std::vector<recob::Hit const*>> ana_pandoraCalo_totHitsInMaxRadius = fCalorimetryRadiusAlg.ana_totHitsInMaxRadius;

      for (size_t i = 0; i<ana_pandoraCalo_totHitsInMaxRadius.size(); i++)
      {
        for (size_t j = 0; j<ana_pandoraCalo_totHitsInMaxRadius.size(); j++)
        {
          printf("%i, %i: %i\n", (int) i, (int) j, (int) ana_pandoraCalo_totHitsInMaxRadius[i][j]->Channel());
        }
      }

      tree_pandoraCalo_par1ChargeInRadius = fCalorimetryRadiusAlg.ana_par1ChargeInRadius;
      tree_pandoraCalo_par2ChargeInRadius = fCalorimetryRadiusAlg.ana_par2ChargeInRadius;
      tree_pandoraCalo_totChargeInRadius = fCalorimetryRadiusAlg.ana_totChargeInRadius;
      tree_pandoraCalo_caloRatio = fCalorimetryRadiusAlg.ana_caloRatio;


      // Fill draw tree (if needed)
      if (fSaveDrawTree) FillDrawTree(ana_pandora_decayVertices, ana_pandoraCalo_totHitsInMaxRadius, ana_pandoraCalo_trackHits, ana_pandoraCalo_showerHits);
    }
    // Fill tree
    pandoraTree->Fill();
    if (fSaveDrawTree) pandoraDrawTree->Fill();
  }

  printf("------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(PreSelectHSN)

#endif // END def PreSelectHSN_module