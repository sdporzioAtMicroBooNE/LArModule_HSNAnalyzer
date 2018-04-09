#ifndef SHOWKINEMATICDISTRIBUTIONS_MODULE
#define SHOWKINEMATICDISTRIBUTIONS_MODULE

#include "ShowKinematicDistributions.h"

ShowKinematicDistributions::ShowKinematicDistributions(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMCLabel(pset.get<std::string>("mcLabel"))  
{} // END constructor ShowKinematicDistributions

ShowKinematicDistributions::~ShowKinematicDistributions()
{} // END destructor ShowKinematicDistributions

void ShowKinematicDistributions::beginJob()
{
  art::ServiceHandle< art::TFileService > tfs;
  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run);
  tDataTree->Branch("subrun",&subrun);
  tDataTree->Branch("event",&event);
  tDataTree->Branch("Vx",&Vx);
  tDataTree->Branch("Vy",&Vy);
  tDataTree->Branch("Vz",&Vz);
  tDataTree->Branch("T",&T);
  tDataTree->Branch("EndX",&EndX);
  tDataTree->Branch("EndY",&EndY);
  tDataTree->Branch("EndZ",&EndZ);
  tDataTree->Branch("EndT",&EndT);
  tDataTree->Branch("Px",&Px);
  tDataTree->Branch("Py",&Py);
  tDataTree->Branch("Pz",&Pz);
  tDataTree->Branch("E",&E);
  tDataTree->Branch("P",&P);
  tDataTree->Branch("Pt",&Pt);
  tDataTree->Branch("Theta",&Theta);
  tDataTree->Branch("Phi",&Phi);
  tDataTree->Branch("Nu_Px",&Nu_Px);
  tDataTree->Branch("Nu_Py",&Nu_Py);
  tDataTree->Branch("Nu_Pz",&Nu_Pz);
  tDataTree->Branch("Nu_P",&Nu_P);
  tDataTree->Branch("Nu_E",&Nu_E);
  tDataTree->Branch("Nu_Theta",&Nu_Theta);
  tDataTree->Branch("Nu_Phi",&Nu_Phi);
  tDataTree->Branch("OpeningAngle",&OpeningAngle);
  tDataTree->Branch("InvariantMass",&InvariantMass);

} // END function beginJob

void ShowKinematicDistributions::endJob()
{
} // END function endJob

void ShowKinematicDistributions::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  OpeningAngle = -999;
  InvariantMass = -999;
  pdgCode.clear();
  Vx.clear();
  Vy.clear();
  Vz.clear();
  T.clear();
  EndX.clear();
  EndY.clear();
  EndZ.clear();
  EndT.clear();
  Px.clear();
  Py.clear();
  Pz.clear();
  E.clear();
  P.clear();
  Pt.clear();
  Theta.clear();
  Phi.clear();
} // END function ClearData

void ShowKinematicDistributions::GetTruthParticles(art::Event const & evt, std::vector<simb::MCParticle const*>& mcps)
{
  // Prepare handle labels
  art::InputTag mcTag {fMCLabel};

  // Find mcTruth
  const auto& mctHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTag);
  for (auto const& mct : (*mctHandle)){
    int nParticles = mct.NParticles();
    printf("|_Number of particles: %i\n", nParticles);
    for (int i=0; i<nParticles; i++)
    {
      const simb::MCParticle & mcp = mct.GetParticle(i);
      pdgCode.push_back(mcp.PdgCode());
      Vx.push_back(mcp.Vx());
      Vy.push_back(mcp.Vy());
      Vz.push_back(mcp.Vz());
      T.push_back(mcp.T());
      EndX.push_back(mcp.EndX());
      EndY.push_back(mcp.EndY());
      EndZ.push_back(mcp.EndZ());
      EndT.push_back(mcp.EndT());
      Px.push_back(mcp.Px());
      Py.push_back(mcp.Py());
      Pz.push_back(mcp.Pz());
      E.push_back(mcp.E());
      P.push_back(mcp.P());
      Pt.push_back(mcp.Pt());
      // Calculate other quantities
      Theta.push_back(acos(mcp.Pz()/mcp.P()));
      Phi.push_back(atan2(mcp.Py(),mcp.Px()));
      // Save pointer to mcparticle
      mcps.push_back(&mcp);
    }

    // Calculate angles

    // Calculate opening angle and invariant mass
    if (nParticles==2)
    {
      // Calculate opening angle
      double dotProduct = Px[0]*Px[1] + Py[0]*Py[1] + Pz[0]*Pz[1];
      OpeningAngle = dotProduct / float(P[0]*P[1]);
      // Calculate invariant mass
      double eTerm = pow((E[0] + E[1]),2.);
      double pTerm = pow(P[0],2.) + pow(P[1],2.) + 2.*dotProduct;
      InvariantMass = sqrt(eTerm - pTerm);

      // Calculate neutrino quantities
      Nu_Px = Px[0] + Px[1];
      Nu_Py = Py[0] + Py[1];
      Nu_Pz = Pz[0] + Pz[1];
      Nu_P = sqrt(pow(Nu_Px,2.) + pow(Nu_Py,2.) + pow(Nu_Pz,2.));
      Nu_E = sqrt(pow(InvariantMass,2.) + pow(Nu_P,2.));
      Nu_Theta = acos(Nu_Pz/Nu_P);
      Nu_Phi = atan2(Nu_Py,Nu_Px);
    }
    else
    {
      OpeningAngle = -999;
      InvariantMass = -999;
    }

  } // End of pfp loop

  return;
} // END function GetTruthParticles


void ShowKinematicDistributions::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  printf("\n-------------------------------------------------------\n");
  
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

  // Get vector of primaries and secondaries pfps
  std::vector<simb::MCParticle const*> mcps;
  GetTruthParticles(evt, mcps);

  // Fill tree and finish event loop
  tDataTree->Fill();
  printf("-------------------------------------------------------\n\n");
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ShowKinematicDistributions)

#endif // END def ShowKinematicDistributions_module