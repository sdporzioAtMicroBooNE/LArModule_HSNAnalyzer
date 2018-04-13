#ifndef PRESELECTPANDORAHSN_MODULE
#define PRESELECTPANDORAHSN_MODULE

#include "PreSelectPandoraHSN.h"

PreSelectHSN::PreSelectHSN(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fFindPandoraVertexAlg(pset),
    fCalorimetryRadiusAlg(pset),
    fTruthMatchingAlg(pset),
    fInstanceName(pset.get<std::string>("InstanceName")),
    fIteration(pset.get<int>("Iteration")),
    fMinTpcBound(pset.get<std::vector<double>>("MinTpcBound")),
    fMaxTpcBound(pset.get<std::vector<double>>("MaxTpcBound")),
    fPfpLabel(pset.get<std::string>("PfpLabel")),
    fHitLabel(pset.get<std::string>("HitLabel")),
    fRadiusProfileLimits(pset.get<std::vector<double>>("RadiusProfileLimits")),
    fRadiusProfileBins(pset.get<int>("RadiusProfileBins")),
    fChannelNorm(pset.get<double>("ChannelNorm")),
    fTickNorm(pset.get<double>("TickNorm")),
    fVerbose(pset.get<bool>("VerboseMode")),
    fSaveDrawTree(pset.get<bool>("SaveDrawTree")),
    fTruthMatching(pset.get<bool>("TruthMatching"))
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
  pandoraTree->Branch("nNeutrinos",&evd.nNeutrinos);
  pandoraTree->Branch("nTwoProngedNeutrinos",&evd.nTwoProngedNeutrinos);
  pandoraTree->Branch("nContainedTwoProngedNeutrinos",&evd.nContainedTwoProngedNeutrinos);
  pandoraTree->Branch("neutrinoPdgCode",&evd.neutrinoPdgCode);
  pandoraTree->Branch("neutrinoNumDaughters",&evd.neutrinoNumDaughters);
  pandoraTree->Branch("neutrinoNumTracks",&evd.neutrinoNumTracks);
  pandoraTree->Branch("neutrinoNumShowers",&evd.neutrinoNumShowers);
  pandoraTree->Branch("calo_totChargeInRadius",&evd.calo_totChargeInRadius);
  pandoraTree->Branch("calo_prong1ChargeInRadius",&evd.calo_prong1ChargeInRadius);
  pandoraTree->Branch("calo_prong2ChargeInRadius",&evd.calo_prong2ChargeInRadius);
  pandoraTree->Branch("calo_caloRatio",&evd.calo_caloRatio);
  pandoraTree->Branch("phys_prongLength",&evd.phys_prongLength);
  pandoraTree->Branch("phys_prongTheta",&evd.phys_prongTheta);
  pandoraTree->Branch("phys_prongPhi",&evd.phys_prongPhi);
  pandoraTree->Branch("phys_prongStartToNeutrinoDistance",&evd.phys_prongStartToNeutrinoDistance);
  pandoraTree->Branch("phys_prongNumHits",&evd.phys_prongNumHits);
  pandoraTree->Branch("phys_openingAngle",&evd.phys_openingAngle);

  pandoraTree->Branch("diag_nuWithMissingAssociatedVertex",&evd.diag_nuWithMissingAssociatedVertex);
  pandoraTree->Branch("diag_nuWithMissingAssociatedTrack",&evd.diag_nuWithMissingAssociatedTrack);
  pandoraTree->Branch("diag_nuProngWithMissingAssociatedHits",&evd.diag_nuProngWithMissingAssociatedHits);
  pandoraTree->Branch("prong_matchedPDG",&evd.prong_matchedPDG);



  if (fSaveDrawTree)
  {
    pandoraDrawTree = tfs->make<TTree>("PandoraDrawTree","");
    pandoraDrawTree->Branch("dv_xyzCoordinates",&dtd.dv_xyzCoordinates);
    pandoraDrawTree->Branch("dv_wireCoordinates",&dtd.dv_wireCoordinates);
    pandoraDrawTree->Branch("dv_tickCoordinates",&dtd.dv_tickCoordinates);
    pandoraDrawTree->Branch("prong1_xyzCoordinates",&dtd.prong1_xyzCoordinates);
    pandoraDrawTree->Branch("prong1_wireCoordinates",&dtd.prong1_wireCoordinates);
    pandoraDrawTree->Branch("prong1_tickCoordinates",&dtd.prong1_tickCoordinates);
    pandoraDrawTree->Branch("prong2_xyzCoordinates",&dtd.prong2_xyzCoordinates);
    pandoraDrawTree->Branch("prong2_wireCoordinates",&dtd.prong2_wireCoordinates);
    pandoraDrawTree->Branch("prong2_tickCoordinates",&dtd.prong2_tickCoordinates);
    pandoraDrawTree->Branch("prong1_hits_p0_wireCoordinates",&dtd.prong1_hits_p0_wireCoordinates);
    pandoraDrawTree->Branch("prong1_hits_p1_wireCoordinates",&dtd.prong1_hits_p1_wireCoordinates);
    pandoraDrawTree->Branch("prong1_hits_p2_wireCoordinates",&dtd.prong1_hits_p2_wireCoordinates);
    pandoraDrawTree->Branch("prong1_hits_p0_tickCoordinates",&dtd.prong1_hits_p0_tickCoordinates);
    pandoraDrawTree->Branch("prong1_hits_p1_tickCoordinates",&dtd.prong1_hits_p1_tickCoordinates);
    pandoraDrawTree->Branch("prong1_hits_p2_tickCoordinates",&dtd.prong1_hits_p2_tickCoordinates);
    pandoraDrawTree->Branch("prong2_hits_p0_wireCoordinates",&dtd.prong2_hits_p0_wireCoordinates);
    pandoraDrawTree->Branch("prong2_hits_p1_wireCoordinates",&dtd.prong2_hits_p1_wireCoordinates);
    pandoraDrawTree->Branch("prong2_hits_p2_wireCoordinates",&dtd.prong2_hits_p2_wireCoordinates);
    pandoraDrawTree->Branch("prong2_hits_p0_tickCoordinates",&dtd.prong2_hits_p0_tickCoordinates);
    pandoraDrawTree->Branch("prong2_hits_p1_tickCoordinates",&dtd.prong2_hits_p1_tickCoordinates);
    pandoraDrawTree->Branch("prong2_hits_p2_tickCoordinates",&dtd.prong2_hits_p2_tickCoordinates);
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
  // The event descriptor is a special class in which we fill all the information we want to know about the current event. At the end of the event loop the information is taken from the event descriptor and filled into the anatree.
  // The draw tree descriptor is the same class specialized for storing information that will allow to make event displays.
  evd.Initialize(run,subrun,event);
  dtd.Initialize();

  if (fVerbose) {printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", evd.event, evd.run, evd.subrun);}
  // Search among pfparticles and get vector of potential neutrino pfps with only two tracks. Return vectors of pfps for neutrinos, tracks and showers in event and decay vertices, which contain information about neutrino vertices with exctly two tracks.
  std::vector<AuxVertex::DecayVertex> ana_decayVertices;
  fFindPandoraVertexAlg.GetPotentialNeutrinoVertices(evt, evd, ana_decayVertices);

  if (ana_decayVertices.size() == 0) printf("No clean vertex candidates found. Moving to next event...\n");
  else
  {
    // Perform calorimetry analysis
    fCalorimetryRadiusAlg.PerformCalorimetry(evt, evd, ana_decayVertices);
    // Perform truth matching (if requested)
    if (fTruthMatching) fTruthMatchingAlg.PerformTruthMatching(evt, evd ,ana_decayVertices);
  }

  // Fill more physics from the decay vertices (ana_decayVertices) into the event descriptor
  evd.ExtractVertexPhysics(ana_decayVertices);
  // Fill draw tree (if requested)
  if (fSaveDrawTree) dtd.FillDrawTreeVariables(ana_decayVertices);

  // Fill tree
  pandoraTree->Fill();
  if (fSaveDrawTree) pandoraDrawTree->Fill();

  printf("------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(PreSelectHSN)

#endif // END def PreSelectHSN_module
