#ifndef PRESELECTPANDORAHSN_MODULE
#define PRESELECTPANDORAHSN_MODULE

#include "PreSelectPandoraHSN.h"

PreSelectHSN::PreSelectHSN(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fFindPandoraVertexAlg(pset),
    fCalorimetryRadiusAlg(pset),
    fInstanceName(pset.get<std::string>("InstanceName")),
    fIteration(pset.get<int>("Iteration")),
    fMinTpcBound(pset.get<std::vector<double>>("MinTpcBound")),
    fMaxTpcBound(pset.get<std::vector<double>>("MaxTpcBound")),
    fPfpLabel(pset.get<std::string>("PfpLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fDistanceCut(pset.get<double>("DistanceCut")),
    fRadiusProfileLimits(pset.get<std::vector<double>>("RadiusProfileLimits")),
    fRadiusProfileBins(pset.get<int>("RadiusProfileBins")),
    fChannelNorm(pset.get<double>("ChannelNorm")),
    fTickNorm(pset.get<double>("TickNorm")),
    fVerbose(pset.get<bool>("VerboseMode")),
    fSaveDrawTree(pset.get<bool>("SaveDrawTree"))
{
  // Get geometry and detector services
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
  // Declare file service handle
  art::ServiceHandle< art::TFileService > tfs;

  // Meta tree containing fcl file parameters
  metaTree = tfs->make<TTree>("MetaData","");
  metaTree->Branch("instanceName",&fInstanceName);
  metaTree->Branch("iteration",&fIteration,"iteration/I"); 
  metaTree->Branch("minTpcBound",&fMinTpcBound);
  metaTree->Branch("maxTpcBound",&fMaxTpcBound);  
  metaTree->Branch("pfpLabel",&fPfpLabel);
  metaTree->Branch("hitLabel",&fHitLabel);
  metaTree->Branch("distanceCut",&fDistanceCut,"distanceCut/D");
  metaTree->Branch("radiusProfileLimits",&fRadiusProfileLimits);
  metaTree->Branch("radiusProfileBins",&fRadiusProfileBins);
  metaTree->Branch("profileTicks",&profileTicks);
  metaTree->Branch("channelNorm",&fChannelNorm,"channelNorm/D");
  metaTree->Branch("tickNorm",&fTickNorm,"tickNorm/D");
  metaTree->Branch("saveDrawTree",&fSaveDrawTree,"saveDrawTree/O");
  metaTree->Fill();

  pandoraTree = tfs->make<TTree>("PandoraData","");
  pandoraTree->Branch("run",&evd.run);
  pandoraTree->Branch("subrun",&evd.subrun);
  pandoraTree->Branch("event",&evd.event);
  pandoraTree->Branch("nNeutrinos",&evd.pandora_nNeutrinos);
  pandoraTree->Branch("nTwoProngedNeutrinos",&evd.pandora_nTwoProngedNeutrinos);
  pandoraTree->Branch("nContainedTwoProngedNeutrinos",&evd.pandora_nContainedTwoProngedNeutrinos);
  pandoraTree->Branch("neutrinoPdgCode",&evd.pandora_neutrinoPdgCode);
  pandoraTree->Branch("neutrinoNumDaughters",&evd.pandora_neutrinoNumDaughters);
  pandoraTree->Branch("neutrinoNumTracks",&evd.pandora_neutrinoNumTracks);
  pandoraTree->Branch("neutrinoNumShowers",&evd.pandora_neutrinoNumShowers);
  pandoraTree->Branch("nTotHits",&evd.pandora_calo_NumTotHits);
  pandoraTree->Branch("totChargeInRadius",&evd.pandora_calo_totChargeInRadius);
  pandoraTree->Branch("par1ChargeInRadius",&evd.pandora_calo_par1ChargeInRadius);
  pandoraTree->Branch("par2ChargeInRadius",&evd.pandora_calo_par2ChargeInRadius);
  pandoraTree->Branch("caloRatio",&evd.pandora_calo_caloRatio);
  pandoraTree->Branch("diag_potentialPfpsWithMissingAssociatedVertex",&evd.pandora_diag_potentialPfpsWithMissingAssociatedVertex);
  pandoraTree->Branch("diag_pandora_diag_dvWithMissingAssociatedHits",&evd.pandora_diag_dvWithMissingAssociatedHits);
  if (fSaveDrawTree)
  {
    pandoraDrawTree = tfs->make<TTree>("PandoraDrawTree","");
    pandoraDrawTree->Branch("dv_xyzCoordinates",&dtd.dv_xyzCoordinates);
    pandoraDrawTree->Branch("dv_wireCoordinates",&dtd.dv_wireCoordinates);
    pandoraDrawTree->Branch("dv_tickCoordinates",&dtd.dv_tickCoordinates);
    pandoraDrawTree->Branch("par1_xyzCoordinates",&dtd.par1_xyzCoordinates);
    pandoraDrawTree->Branch("par1_wireCoordinates",&dtd.par1_wireCoordinates);
    pandoraDrawTree->Branch("par1_tickCoordinates",&dtd.par1_tickCoordinates);
    pandoraDrawTree->Branch("par2_xyzCoordinates",&dtd.par2_xyzCoordinates);
    pandoraDrawTree->Branch("par2_wireCoordinates",&dtd.par2_wireCoordinates);
    pandoraDrawTree->Branch("par2_tickCoordinates",&dtd.par2_tickCoordinates);
    pandoraDrawTree->Branch("par1_hits_p0_wireCoordinates",&dtd.par1_hits_p0_wireCoordinates);
    pandoraDrawTree->Branch("par1_hits_p1_wireCoordinates",&dtd.par1_hits_p1_wireCoordinates);
    pandoraDrawTree->Branch("par1_hits_p2_wireCoordinates",&dtd.par1_hits_p2_wireCoordinates);
    pandoraDrawTree->Branch("par1_hits_p0_tickCoordinates",&dtd.par1_hits_p0_tickCoordinates);
    pandoraDrawTree->Branch("par1_hits_p1_tickCoordinates",&dtd.par1_hits_p1_tickCoordinates);
    pandoraDrawTree->Branch("par1_hits_p2_tickCoordinates",&dtd.par1_hits_p2_tickCoordinates);
    pandoraDrawTree->Branch("par2_hits_p0_wireCoordinates",&dtd.par2_hits_p0_wireCoordinates);
    pandoraDrawTree->Branch("par2_hits_p1_wireCoordinates",&dtd.par2_hits_p1_wireCoordinates);
    pandoraDrawTree->Branch("par2_hits_p2_wireCoordinates",&dtd.par2_hits_p2_wireCoordinates);
    pandoraDrawTree->Branch("par2_hits_p0_tickCoordinates",&dtd.par2_hits_p0_tickCoordinates);
    pandoraDrawTree->Branch("par2_hits_p1_tickCoordinates",&dtd.par2_hits_p1_tickCoordinates);
    pandoraDrawTree->Branch("par2_hits_p2_tickCoordinates",&dtd.par2_hits_p2_tickCoordinates);
    pandoraDrawTree->Branch("tot_hits_p0_wireCoordinates",&dtd.tot_hits_p0_wireCoordinates);
    pandoraDrawTree->Branch("tot_hits_p1_wireCoordinates",&dtd.tot_hits_p1_wireCoordinates);
    pandoraDrawTree->Branch("tot_hits_p2_wireCoordinates",&dtd.tot_hits_p2_wireCoordinates);
    pandoraDrawTree->Branch("tot_hits_p0_tickCoordinates",&dtd.tot_hits_p0_tickCoordinates);
    pandoraDrawTree->Branch("tot_hits_p1_tickCoordinates",&dtd.tot_hits_p1_tickCoordinates);
    pandoraDrawTree->Branch("tot_hits_p2_tickCoordinates",&dtd.tot_hits_p2_tickCoordinates);
  }

} // END function beginJob

void PreSelectHSN::endJob()
{} // END function endJob

void PreSelectHSN::ClearData()
{} // END function ClearData


// Core analysis. This is where all the functions are executed. Gets repeated event by event.
void PreSelectHSN::analyze(art::Event const & evt)
{
  if (fVerbose) {printf("\n------------------------------------------------\n");}

  // Determine event ID and create an event descriptor
  int run = evt.id().run();
  int subrun = evt.id().subRun();
  int event = evt.id().event();
  evd.Initialize(run,subrun,event);
  dtd.Initialize();
  if (fVerbose) {printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", evd.event, evd.run, evd.subrun);}

  // Search among pfparticles and get vector of potential neutrino pfps with only two tracks. Return vectors of pfps for neutrinos, tracks and showers in event and decay vertices, which contain information about neutrino vertices with exctly two tracks.
  fFindPandoraVertexAlg.GetPotentialNeutrinoVertices(evd, evt);
  std::vector<recob::PFParticle const*> ana_pandora_neutrinos = fFindPandoraVertexAlg.ana_pandora_neutrinos;
  std::vector<recob::PFParticle const*> ana_pandora_tracks = fFindPandoraVertexAlg.ana_pandora_tracks;
  std::vector<recob::PFParticle const*> ana_pandora_showers = fFindPandoraVertexAlg.ana_pandora_showers;
  std::vector<AuxVertex::DecayVertex> ana_pandora_decayVertices = fFindPandoraVertexAlg.ana_pandora_decayVertices;

  // Just being smart, if there's no clean vertices, there's no reason to waste time loading hits that won't be used
  if (evd.pandora_nContainedTwoProngedNeutrinos==0)
  {
    printf("No clean vertex candidates found. Moving to next event...\n");
  }
  else
  {
    // Get vectors containing hits for each track/shower object in order to perform calorimetry
    fCalorimetryRadiusAlg.GetHitVectors(evd, ana_pandora_decayVertices, evt, ana_pandora_tracks, ana_pandora_showers);
    std::vector<recob::Hit const*> ana_pandoraCalo_totHits = fCalorimetryRadiusAlg.ana_calo_totHits;
    std::vector<std::vector<recob::Hit const*>> ana_pandoraCalo_trackHits = fCalorimetryRadiusAlg.ana_calo_trackHits;
    std::vector<std::vector<recob::Hit const*>> ana_pandoraCalo_showerHits = fCalorimetryRadiusAlg.ana_calo_showerHits;
    evd.pandora_calo_NumTotHits = fCalorimetryRadiusAlg.tree_calo_NumTotHits;

    // Perform calorimetry analysis
    fCalorimetryRadiusAlg.PerformCalorimetry(evd, ana_pandora_decayVertices, ana_pandoraCalo_totHits, ana_pandoraCalo_trackHits, ana_pandoraCalo_showerHits);
    std::vector<std::vector<recob::Hit const*>> ana_pandoraCalo_totHitsInMaxRadius = fCalorimetryRadiusAlg.ana_calo_totHitsInMaxRadius;
    evd.pandora_calo_totChargeInRadius = fCalorimetryRadiusAlg.tree_calo_totChargeInRadius;
    evd.pandora_calo_par1ChargeInRadius = fCalorimetryRadiusAlg.tree_calo_par1ChargeInRadius;
    evd.pandora_calo_par2ChargeInRadius = fCalorimetryRadiusAlg.tree_calo_par2ChargeInRadius;
    evd.pandora_calo_caloRatio = fCalorimetryRadiusAlg.tree_calo_caloRatio;

    // Fill draw tree (if needed)
    if (fSaveDrawTree) dtd.FillDrawTreeVariables(ana_pandora_decayVertices, ana_pandoraCalo_totHitsInMaxRadius, ana_pandoraCalo_trackHits, ana_pandoraCalo_showerHits);
  }
  // Fill tree
  pandoraTree->Fill();
  if (fSaveDrawTree) pandoraDrawTree->Fill();

  printf("------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(PreSelectHSN)

#endif // END def PreSelectHSN_module