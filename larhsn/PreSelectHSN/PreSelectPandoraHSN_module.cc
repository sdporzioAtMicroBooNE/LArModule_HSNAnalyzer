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
    fMcsLabel(pset.get<std::string>("McsLabel")),
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
  metaTree->Branch("mcsLabel",&fMcsLabel);
  metaTree->Branch("radiusProfileLimits",&fRadiusProfileLimits);
  metaTree->Branch("radiusProfileBins",&fRadiusProfileBins);
  metaTree->Branch("profileTicks",&profileTicks);
  metaTree->Branch("channelNorm",&fChannelNorm,"channelNorm/D");
  metaTree->Branch("tickNorm",&fTickNorm,"tickNorm/D");
  metaTree->Branch("saveDrawTree",&fSaveDrawTree,"saveDrawTree/O");
  metaTree->Fill();

  pandoraTree = tfs->make<TTree>("PandoraData","");
  // Event
  pandoraTree->Branch("run",&evd.run);
  pandoraTree->Branch("subrun",&evd.subrun);
  pandoraTree->Branch("event",&evd.event);
  // Analysis
  pandoraTree->Branch("nNeutrinos",&evd.nNeutrinos);
  pandoraTree->Branch("nTwoProngedNeutrinos",&evd.nTwoProngedNeutrinos);
  pandoraTree->Branch("nContainedTwoProngedNeutrinos",&evd.nContainedTwoProngedNeutrinos);
  pandoraTree->Branch("neutrinoPdgCode",&evd.neutrinoPdgCode);
  pandoraTree->Branch("neutrinoNumDaughters",&evd.neutrinoNumDaughters);
  pandoraTree->Branch("neutrinoNumTracks",&evd.neutrinoNumTracks);
  pandoraTree->Branch("neutrinoNumShowers",&evd.neutrinoNumShowers);
  // Coordinates
  pandoraTree->Branch("phys_nuPositionX",&evd.phys_nuPosX);
  pandoraTree->Branch("phys_nuPositionY",&evd.phys_nuPosY);
  pandoraTree->Branch("phys_nuPositionZ",&evd.phys_nuPosZ);
  pandoraTree->Branch("phys_prongPositionX",&evd.phys_prongPosX);
  pandoraTree->Branch("phys_prongPositionY",&evd.phys_prongPosY);
  pandoraTree->Branch("phys_prongPositionZ",&evd.phys_prongPosZ);
  pandoraTree->Branch("phys_prongStartPositionX",&evd.phys_prongStartPosX);
  pandoraTree->Branch("phys_prongStartPositionY",&evd.phys_prongStartPosY);
  pandoraTree->Branch("phys_prongStartPositionZ",&evd.phys_prongStartPosZ);
  pandoraTree->Branch("phys_prongEndPositionX",&evd.phys_prongEndPosX);
  pandoraTree->Branch("phys_prongEndPositionY",&evd.phys_prongEndPosY);
  pandoraTree->Branch("phys_prongEndPositionZ",&evd.phys_prongEndPosZ);
  // Direction
  pandoraTree->Branch("phys_prongDirectionX",&evd.phys_prongDirX);
  pandoraTree->Branch("phys_prongDirectionY",&evd.phys_prongDirY);
  pandoraTree->Branch("phys_prongDirectionZ",&evd.phys_prongDirZ);
  pandoraTree->Branch("phys_prongTheta",&evd.phys_prongTheta);
  pandoraTree->Branch("phys_prongPhi",&evd.phys_prongPhi);
  // Hypothesis info
  pandoraTree->Branch("phys_prongPdgCode_h1",&evd.phys_prongPdgCode_h1);
  pandoraTree->Branch("phys_prongPdgCode_h2",&evd.phys_prongPdgCode_h2);
  pandoraTree->Branch("phys_prongMass_h1",&evd.phys_prongMass_h1);
  pandoraTree->Branch("phys_prongMass_h2",&evd.phys_prongMass_h2);
  // Prong momentum (by range, assuming h1)
  pandoraTree->Branch("phys_prongEnergy_ByRange_h1",&evd.phys_prongEnergy_ByRange_h1);
  pandoraTree->Branch("phys_prongMomMag_ByRange_h1",&evd.phys_prongMomMag_ByRange_h1);
  pandoraTree->Branch("phys_prongMom_ByRange_h1_X",&evd.phys_prongMom_ByRange_h1_X);
  pandoraTree->Branch("phys_prongMom_ByRange_h1_Y",&evd.phys_prongMom_ByRange_h1_Y);
  pandoraTree->Branch("phys_prongMom_ByRange_h1_Z",&evd.phys_prongMom_ByRange_h1_Z);
  // Tot momentum (by range, assuming h1)
  pandoraTree->Branch("phys_invariantMass_ByRange_h1",&evd.phys_invariantMass_ByRange_h1);
  pandoraTree->Branch("phys_totEnergy_ByRange_h1",&evd.phys_totEnergy_ByRange_h1);
  pandoraTree->Branch("phys_totMomMag_ByRange_h1",&evd.phys_totMomMag_ByRange_h1);
  pandoraTree->Branch("phys_totMom_ByRange_h1_X",&evd.phys_totMom_ByRange_h1_X);
  pandoraTree->Branch("phys_totMom_ByRange_h1_Y",&evd.phys_totMom_ByRange_h1_Y);
  pandoraTree->Branch("phys_totMom_ByRange_h1_Z",&evd.phys_totMom_ByRange_h1_Z);
  // Tot momentum direction (by range, assuming h1)
  pandoraTree->Branch("phys_totDirection_ByRange_h1_X",&evd.phys_totDir_ByRange_h1_X);
  pandoraTree->Branch("phys_totDirection_ByRange_h1_Y",&evd.phys_totDir_ByRange_h1_Y);
  pandoraTree->Branch("phys_totDirection_ByRange_h1_Z",&evd.phys_totDir_ByRange_h1_Z);
  pandoraTree->Branch("phys_totTheta_ByRange_h1",&evd.phys_totTheta_ByRange_h1);
  pandoraTree->Branch("phys_totPhi_ByRange_h1",&evd.phys_totPhi_ByRange_h1);
  // Prong momentum (by range, assuming h2)
  pandoraTree->Branch("phys_prongEnergy_ByRange_h2",&evd.phys_prongEnergy_ByRange_h2);
  pandoraTree->Branch("phys_prongMomMag_ByRange_h2",&evd.phys_prongMomMag_ByRange_h2);
  pandoraTree->Branch("phys_prongMom_ByRange_h2_X",&evd.phys_prongMom_ByRange_h2_X);
  pandoraTree->Branch("phys_prongMom_ByRange_h2_Y",&evd.phys_prongMom_ByRange_h2_Y);
  pandoraTree->Branch("phys_prongMom_ByRange_h2_Z",&evd.phys_prongMom_ByRange_h2_Z);
  // Tot momentum (by range, assuming h2)
  pandoraTree->Branch("phys_invariantMass_ByRange_h2",&evd.phys_invariantMass_ByRange_h2);
  pandoraTree->Branch("phys_totEnergy_ByRange_h2",&evd.phys_totEnergy_ByRange_h2);
  pandoraTree->Branch("phys_totMomMag_ByRange_h2",&evd.phys_totMomMag_ByRange_h2);
  pandoraTree->Branch("phys_totMom_ByRange_h2_X",&evd.phys_totMom_ByRange_h2_X);
  pandoraTree->Branch("phys_totMom_ByRange_h2_Y",&evd.phys_totMom_ByRange_h2_Y);
  pandoraTree->Branch("phys_totMom_ByRange_h2_Z",&evd.phys_totMom_ByRange_h2_Z);
  // Tot momentum direction (by range, assuming h2)
  pandoraTree->Branch("phys_totDirection_ByRange_h2_X",&evd.phys_totDir_ByRange_h2_X);
  pandoraTree->Branch("phys_totDirection_ByRange_h2_Y",&evd.phys_totDir_ByRange_h2_Y);
  pandoraTree->Branch("phys_totDirection_ByRange_h2_Z",&evd.phys_totDir_ByRange_h2_Z);
  pandoraTree->Branch("phys_totTheta_ByRange_h2",&evd.phys_totTheta_ByRange_h2);
  pandoraTree->Branch("phys_totPhi_ByRange_h2",&evd.phys_totPhi_ByRange_h2);
  // Momentum (By Mcs)
  pandoraTree->Branch("phys_prongPdgCodeHypothesis_ByMcs",&evd.phys_prongPdgCodeHypothesis_ByMcs);
  pandoraTree->Branch("phys_prongIsBestFwd_ByMcs",&evd.phys_prongIsBestFwd_ByMcs);
  // Prong Momentum (By Mcs, best)
  pandoraTree->Branch("phys_prongMomMag_ByMcs_best_h1",&evd.phys_prongMomMag_ByMcs_best_h1);
  pandoraTree->Branch("phys_prongEnergy_ByMcs_best_h1",&evd.phys_prongEnergy_ByMcs_best_h1);
  pandoraTree->Branch("phys_prongMom_ByMcs_best_h1_X",&evd.phys_prongMom_ByMcs_best_h1_X);
  pandoraTree->Branch("phys_prongMom_ByMcs_best_h1_Y",&evd.phys_prongMom_ByMcs_best_h1_Y);
  pandoraTree->Branch("phys_prongMom_ByMcs_best_h1_Z",&evd.phys_prongMom_ByMcs_best_h1_Z);
  pandoraTree->Branch("phys_prongMomMag_ByMcs_best_h2",&evd.phys_prongMomMag_ByMcs_best_h2);
  pandoraTree->Branch("phys_prongEnergy_ByMcs_best_h2",&evd.phys_prongEnergy_ByMcs_best_h2);
  pandoraTree->Branch("phys_prongMom_ByMcs_best_h2_X",&evd.phys_prongMom_ByMcs_best_h2_X);
  pandoraTree->Branch("phys_prongMom_ByMcs_best_h2_Y",&evd.phys_prongMom_ByMcs_best_h2_Y);
  pandoraTree->Branch("phys_prongMom_ByMcs_best_h2_Z",&evd.phys_prongMom_ByMcs_best_h2_Z);
  // Tot momentum (by range, assuming both muons, best)
  pandoraTree->Branch("phys_totMomMag_ByMcs_best_h1",&evd.phys_totMomMag_ByMcs_best_h1);
  pandoraTree->Branch("phys_totEnergy_ByMcs_best_h1",&evd.phys_totEnergy_ByMcs_best_h1);
  pandoraTree->Branch("phys_invariantMass_ByMcs_best_h1",&evd.phys_invariantMass_ByMcs_best_h1);
  pandoraTree->Branch("phys_totMom_ByMcs_best_h1_X",&evd.phys_totMom_ByMcs_best_h1_X);
  pandoraTree->Branch("phys_totMom_ByMcs_best_h1_Y",&evd.phys_totMom_ByMcs_best_h1_Y);
  pandoraTree->Branch("phys_totMom_ByMcs_best_h1_Z",&evd.phys_totMom_ByMcs_best_h1_Z);
  pandoraTree->Branch("phys_totMomMag_ByMcs_best_h2",&evd.phys_totMomMag_ByMcs_best_h2);
  pandoraTree->Branch("phys_totEnergy_ByMcs_best_h2",&evd.phys_totEnergy_ByMcs_best_h2);
  pandoraTree->Branch("phys_invariantMass_ByMcs_best_h2",&evd.phys_invariantMass_ByMcs_best_h2);
  pandoraTree->Branch("phys_totMom_ByMcs_best_h2_X",&evd.phys_totMom_ByMcs_best_h2_X);
  pandoraTree->Branch("phys_totMom_ByMcs_best_h2_Y",&evd.phys_totMom_ByMcs_best_h2_Y);
  pandoraTree->Branch("phys_totMom_ByMcs_best_h2_Z",&evd.phys_totMom_ByMcs_best_h2_Z);
  // Tot momentum direction (by range, assuming both muons, best)
  pandoraTree->Branch("phys_totTheta_ByMcs_best_h1",&evd.phys_totTheta_ByMcs_best_h1);
  pandoraTree->Branch("phys_totPhi_ByMcs_best_h1",&evd.phys_totPhi_ByMcs_best_h1);
  pandoraTree->Branch("phys_totDir_ByMcs_best_h1_X",&evd.phys_totDir_ByMcs_best_h1_X);
  pandoraTree->Branch("phys_totDir_ByMcs_best_h1_Y",&evd.phys_totDir_ByMcs_best_h1_Y);
  pandoraTree->Branch("phys_totDir_ByMcs_best_h1_Z",&evd.phys_totDir_ByMcs_best_h1_Z);
  pandoraTree->Branch("phys_totTheta_ByMcs_best_h2",&evd.phys_totTheta_ByMcs_best_h2);
  pandoraTree->Branch("phys_totPhi_ByMcs_best_h2",&evd.phys_totPhi_ByMcs_best_h2);
  pandoraTree->Branch("phys_totDir_ByMcs_best_h2_X",&evd.phys_totDir_ByMcs_best_h2_X);
  pandoraTree->Branch("phys_totDir_ByMcs_best_h2_Y",&evd.phys_totDir_ByMcs_best_h2_Y);
  pandoraTree->Branch("phys_totDir_ByMcs_best_h2_Z",&evd.phys_totDir_ByMcs_best_h2_Z);
  // Others
  pandoraTree->Branch("phys_prongLength",&evd.phys_prongLength);
  pandoraTree->Branch("phys_prongStartToNeutrinoDistance",&evd.phys_prongStartToNeutrinoDistance);
  pandoraTree->Branch("phys_prongNumHits",&evd.phys_prongNumHits);
  pandoraTree->Branch("phys_openingAngle",&evd.phys_openingAngle);
  // Calorimetry
  pandoraTree->Branch("calo_totChargeInRadius",&evd.calo_totChargeInRadius);
  pandoraTree->Branch("calo_prong1ChargeInRadius",&evd.calo_prong1ChargeInRadius);
  pandoraTree->Branch("calo_prong2ChargeInRadius",&evd.calo_prong2ChargeInRadius);
  pandoraTree->Branch("calo_caloRatio",&evd.calo_caloRatio);
  // Status
  pandoraTree->Branch("status_nuWithMissingAssociatedVertex",&evd.status_nuWithMissingAssociatedVertex);
  pandoraTree->Branch("status_nuWithMissingAssociatedTrack",&evd.status_nuWithMissingAssociatedTrack);
  pandoraTree->Branch("status_nuProngWithMissingAssociatedHits",&evd.status_nuProngWithMissingAssociatedHits);


  if (fSaveDrawTree)
  {
    pandoraDrawTree = tfs->make<TTree>("PandoraDrawTree","");
    pandoraDrawTree->Branch("run",&dtd.run);
    pandoraDrawTree->Branch("subrun",&dtd.subrun);
    pandoraDrawTree->Branch("event",&dtd.event);
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

  // Determine event ID
  int run = evt.id().run();
  int subrun = evt.id().subRun();
  int event = evt.id().event();

  // The event descriptor is a special class in which we fill all the information we want to know about the current event. At the end of the event loop the information is taken from the event descriptor and filled into the anatree.
  // The draw tree descriptor is the same class specialized for storing information that will allow to make event displays.
  evd.Initialize(run,subrun,event);
  dtd.Initialize(run,subrun,event);

  if (fVerbose) {printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", evd.event, evd.run, evd.subrun);}

  // Search among pfparticles and get vector of potential neutrino pfps with only two tracks. Return vectors of pfps for neutrinos, tracks and showers in event and decay vertices, which contain information about neutrino vertices with exctly two tracks.
  std::vector<AuxVertex::DecayVertex> ana_decayVertices;
  fFindPandoraVertexAlg.GetPotentialNeutrinoVertices(evt, evd, ana_decayVertices);

  // Execute functions that should be run only if there's available decay vertices
  if (ana_decayVertices.size() == 0) printf("No clean vertex candidates found. Moving to next event...\n");
  else
  {
    // Perform calorimetry analysis
    fCalorimetryRadiusAlg.PerformCalorimetry(evt, evd, ana_decayVertices);
  }

  // Fill more physics from the decay vertices (ana_decayVertices) into the event descriptor and fill tree
  evd.ExtractVertexPhysics(ana_decayVertices);
  pandoraTree->Fill();
  // If requested, do the same for the draw tree
  if (fSaveDrawTree)
  {
    dtd.FillDrawTreeVariables(ana_decayVertices);
    pandoraDrawTree->Fill();
  }

  printf("------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(PreSelectHSN)

#endif // END def PreSelectHSN_module
