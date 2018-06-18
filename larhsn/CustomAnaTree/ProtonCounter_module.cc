/*
ProtonCounter module.
Loop through all MCTruth particle and find all protons, collect various information about them.
*/



#ifndef ProtonCounter_module
#define ProtonCounter_module

#include "AnaHelper.h"

// Analyzer class
class ProtonCounter : public art::EDAnalyzer
{
public:
  explicit ProtonCounter(fhicl::ParameterSet const & pset);
  virtual ~ProtonCounter();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Declare fhiclcpp variables
  std::string fG4Label;

  // Declare trees
  TTree *tDataTree;
  int run, subrun, event;
  std::vector<int> pdgCode;

  int nProtons;
  std::vector<double> protonEnergy, protonMomentum;

  void ClearData();
}; // End class ProtonCounter

ProtonCounter::ProtonCounter(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fG4Label(pset.get<std::string>("g4_label"))
{} // END constructor ProtonCounter

ProtonCounter::~ProtonCounter()
{} // END destructor ProtonCounter

void ProtonCounter::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run,"run/I");
  tDataTree->Branch("subrun",&subrun,"subrun/I");
  tDataTree->Branch("event",&event,"event/I");
  tDataTree->Branch("event",&event,"event/I");
  tDataTree->Branch("pdgCode",&pdgCode);
  tDataTree->Branch("nProtons",&nProtons);
  tDataTree->Branch("protonEnergy",&protonEnergy);
  tDataTree->Branch("protonMomentum",&protonMomentum);

} // END function beginJob

void ProtonCounter::endJob()
{
} // END function endJob

void ProtonCounter::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  pdgCode.clear();
  protonEnergy.clear();
  protonMomentum.clear();
} // END function ClearData

void ProtonCounter::analyze(art::Event const & evt)
{
  
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

  art::InputTag g4Tag {fG4Label};
  const auto& mcpHandle = evt.getValidHandle< std::vector<simb::MCParticle> >(g4Tag);
  for(std::vector<int>::size_type i=0; i!=(*mcpHandle).size(); i++)
  {
    art::Ptr<simb::MCParticle> mcp(mcpHandle,i);
    pdgCode.push_back(mcp->PdgCode());
    // Collect informations about protons in event
    if (mcp->PdgCode()==2212)
    {
      protonEnergy.push_back(mcp->E());
      protonMomentum.push_back(mcp->P());
    }
  }
  nProtons = protonEnergy.size();

  // Fill tree and finish event loop
  tDataTree->Fill();
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(ProtonCounter)

#endif // END def ProtonCounter_module