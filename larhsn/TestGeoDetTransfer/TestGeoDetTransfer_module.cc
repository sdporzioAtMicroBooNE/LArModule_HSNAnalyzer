#ifndef TestGeoDetTransfer_module
#define TestGeoDetTransfer_module

#include "TestGeoDetTransfer.h"

TestGeoDetTransfer::TestGeoDetTransfer(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{} // END constructor TestGeoDetTransfer

TestGeoDetTransfer::~TestGeoDetTransfer()
{} // END destructor TestGeoDetTransfer

void TestGeoDetTransfer::beginJob()
{
  // Declare tree variables
  art::ServiceHandle< art::TFileService > tfs;

  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
} // END function beginJob

void TestGeoDetTransfer::endJob()
{
} // END function endJob


void TestGeoDetTransfer::PrintCoordinates(float* xyz)
{
  raw::ChannelID_t channel0 = fGeometry->NearestChannel(xyz,0);
  raw::ChannelID_t channel1 = fGeometry->NearestChannel(xyz,1);
  raw::ChannelID_t channel2 = fGeometry->NearestChannel(xyz,2);
  double tick0 = fDetectorProperties->ConvertXToTicks(xyz[0], 0, 0, 0);
  double tick1 = fDetectorProperties->ConvertXToTicks(xyz[0], 1, 0, 0);
  double tick2 = fDetectorProperties->ConvertXToTicks(xyz[0], 2, 0, 0);
  std::printf("\n\n----> InsideFunction\n");
  std::printf("----> XYZ: [%.1f,%.1f,%.1f]\n", xyz[0], xyz[1], xyz[2]);
  std::printf("----> Channels: [%i,%i,%i]\n", channel0, channel1, channel2);
  std::printf("----> Ticks: [%.1f,%.1f,%.1f]\n", tick0, tick1, tick2);

  return;
} // END function SetDetectorCoordinates


void TestGeoDetTransfer::analyze(art::Event const & evt)
{
  float xyz[3] = {100,0,400};
  PrintCoordinates(xyz);
  AuxVertex::PrintCoordinates(fGeometry, fDetectorProperties, xyz);

  return;
} // END function analyze


// Name that will be used by the .fcl to invoke the module
DEFINE_ART_MODULE(TestGeoDetTransfer)

#endif // END def TestGeoDetTransfer_module