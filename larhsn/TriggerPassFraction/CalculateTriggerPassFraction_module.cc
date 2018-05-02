#ifndef CalculateTriggerPassFraction_module
#define CalculateTriggerPassFraction_module

#include "AnaHelper.h"

// Analyzer class
class CalculateTriggerPassFraction : public art::EDAnalyzer
{
public:
  explicit CalculateTriggerPassFraction(fhicl::ParameterSet const & pset);
  virtual ~CalculateTriggerPassFraction();
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();
private:
  // Declare fhiclcpp variables
  std::string f_swTriggerLabel;

  // Declare trees
  TTree *tDataTree;

  // Declare analysis variables
  // Declare analysis functions
  void ClearData();
}; // End class CalculateTriggerPassFraction

CalculateTriggerPassFraction::CalculateTriggerPassFraction(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fSwTriggerLabel(pset.get<std::string>("SwTriggerLabel")),
{} // END constructor CalculateTriggerPassFraction

CalculateTriggerPassFraction::~CalculateTriggerPassFraction()
{} // END destructor CalculateTriggerPassFraction

void CalculateTriggerPassFraction::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run);
  tDataTree->Branch("subrun",&subrun);
  tDataTree->Branch("event",&event);
} // END function beginJob

void CalculateTriggerPassFraction::endJob()
{
} // END function endJob

void CalculateTriggerPassFraction::ClearData()
{
} // END function ClearData

void CalculateTriggerPassFraction::analyze(art::Event const & evt)
{
  // Start by clearing all data.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();

   art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
   evt.getByLabel(fSwTriggerLabel, softwareTriggerHandle);
  if (softwareTriggerHandle.isValid()) {
   std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();
   for (int i = 0; i < int(algoNames.size()); i++) {
   if (algoNames[i] == "BNB_FEMBeamTriggerAlgo") {
   EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[0]) ? true : false;
   }
   }
   }

} // END function analyze

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(CalculateTriggerPassFraction)

#endif // END def CalculateTriggerPassFraction_module