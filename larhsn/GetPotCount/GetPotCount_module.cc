#ifndef GetPotCount_module
#define GetPotCount_module

#include "AnaHelper.h"

// Analyzer class
class GetPotCount : public art::EDAnalyzer
{
public:
  explicit GetPotCount(fhicl::ParameterSet const & pset);
  virtual ~GetPotCount();
  void analyze(art::Event const & evt) ;
  void endSubRun(art::SubRun const & sr);
  void beginJob();
  void endJob();
private:
  // Declare trees
  TTree *tPotCount;
  TTree *tEventCount;
  bool fVerbose;

  // Declare analysis variables
  int run, subrun, event;
  double pot = 0.;

  // Declare analysis functions
  void ClearData();
}; // End class GetPotCount

GetPotCount::GetPotCount(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fVerbose(pset.get<bool>("verbose"))
{} // END constructor GetPotCount

GetPotCount::~GetPotCount()
{} // END destructor GetPotCount

void GetPotCount::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;
  tPotCount = tfs->make<TTree>("PotCount","");
  tPotCount->Branch("pot",&pot);
  tEventCount = tfs->make<TTree>("EventCount","");
  tEventCount->Branch("run",&run,"run/I");
  tEventCount->Branch("subrun",&subrun,"subrun/I");
  tEventCount->Branch("event",&event,"event/I");
} // END function beginJob

void GetPotCount::endJob()
{
} // END function endJob

void GetPotCount::ClearData()
{
  run = -1;
  subrun = -1;
  event = -1;
} // END function ClearData

void GetPotCount::analyze(art::Event const & evt)
{
  // Core analysis. Use all the previously defined functions to determine success rate. This will be repeated event by event.
  if (fVerbose) printf("\n|-----------------------------------------------------|");
  if (fVerbose) printf("\n|  GETPOTCOUNT MODULE                                 |");
  if (fVerbose) printf("\n|-----------------------------------------------------|\n");  
  // Start by clearing all the vectors.
  ClearData();

  // Determine event ID
  run = evt.id().run();
  subrun = evt.id().subRun();
  event = evt.id().event();
  if (fVerbose) printf("||INFORMATION FOR EVENT %i [RUN %i, SUBRUN %i]||\n", event, run, subrun);
  
  // Start performing analysis

  // Fill tree and finish event loop
  tEventCount->Fill();
  if (fVerbose) printf("-------------------------------------------------------\n\n");
} // END function analyze

void GetPotCount::endSubRun(art::SubRun const & sr) 
{ 
  auto const & POTSummaryHandle = sr.getValidHandle < sumdata::POTSummary >("generator");
  auto const & POTSummary(*POTSummaryHandle);
  const double totalPot = POTSummary.totpot;
  pot = totalPot;
  // pot = 5.33292e+15;
  std::cout << "----------------------------" << std::endl;
  std::cout << "Total POT / subRun: " << pot << std::endl;
  std::cout << "----------------------------" << std::endl;
  tPotCount->Fill();
} // END function endSubRun

// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(GetPotCount)

#endif // END def GetPotCount_module