#ifndef SHOWKINEMATICDISTRIBUTIONS_MODULE
#define SHOWKINEMATICDISTRIBUTIONS_MODULE
#include "CompareTimeVsTrigger.h"

// Creator, destructor, empty functions
CompareTimeVsTrigger::CompareTimeVsTrigger(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMcLabel(pset.get<std::string>("mcLabel")),
    fSwtrigLabel(pset.get<std::string>("swtrigLabel")),
    fTrigNames(pset.get<std::vector<std::string>>("triggerNames"))
{} 
CompareTimeVsTrigger::~CompareTimeVsTrigger() {}
void CompareTimeVsTrigger::endJob() {}

void CompareTimeVsTrigger::beginJob()
{
  art::ServiceHandle< art::TFileService > tfs;
  tMetaTree = tfs->make<TTree>("MetaData","");
  tMetaTree->Branch("triggerName",&fTrigNames);
  tDataTree = tfs->make<TTree>("Data","");
  tDataTree->Branch("run",&run);
  tDataTree->Branch("subrun",&subrun);
  tDataTree->Branch("event",&event);
  tDataTree->Branch("startTime",&startTime);
  tDataTree->Branch("triggerPass",&triggerPass);

} 

void CompareTimeVsTrigger::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
  startTime = -1;
  triggerPass.clear();
} 


void CompareTimeVsTrigger::analyze(art::Event const & evt)
{
  ClearData();
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);

  // Set producer labels
  art::InputTag mcTag {fMcLabel};
  art::InputTag swtrigTag {fSwtrigLabel};

  const auto& mctHandle = evt.getValidHandle< std::vector<simb::MCTruth> >(mcTag);
  const auto& swtrigHandle = evt.getValidHandle< raw::ubdaqSoftwareTriggerData >(swtrigTag);

  // Get MC Truth information
  for (auto const& mct : (*mctHandle)){
    const simb::MCParticle & mcp = mct.GetParticle(0);
    startTime = mcp.T();
    printf("|_Start time: %.2f\n", startTime);
  }

  // Get trigger information
  for (auto const trigName : fTrigNames)
  {
    bool passedTrig = swtrigHandle->passedAlgo(trigName);
    triggerPass.push_back(passedTrig);
    printf("|_Passed trigger %s: %i", trigName.c_str(), (int) passedTrig);
  }

  tDataTree->Fill();
  printf("-------------------------------------------------------\n\n");
} // END function analyze


DEFINE_ART_MODULE(CompareTimeVsTrigger)
#endif