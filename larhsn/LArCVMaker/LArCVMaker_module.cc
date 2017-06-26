// framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// data product includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Wire.h"

// root includes
#include "TFile.h"
#include "TTree.h"

// c++ includes
#include <vector>
#include <iterator>
#include <typeinfo>
#include <memory>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>

// local includes
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/IOManager.h"

namespace larhsn {

class LArCVMaker : public art::EDAnalyzer {
public:
  explicit LArCVMaker(fhicl::ParameterSet const & pset);
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();

private:
  void ClearData();
  void ResetROI();
  void SetROISize();
  bool IsInsideTPC(double xyz[]);
  void BoxCoordinates(const simb::MCParticle& part, double xyz_origin[], double xyz_end[]);
  int FindBestAPA(std::vector<int> apas);
  int FindROI(int apa, int plane);

  larcv::IOManager fMgr;

  std::string fWireModuleLabel;
  std::string fMcTruthModuleLabel;
  int fMaxTick;
  int fADCCut;
  double fBoxSizeX;
  double fBoxSizeY;
  double fBoxSizeZ;
  int fEventType;

  const int fNumberChannels[3] = { 2400, 2400, 3456 };
  const int fFirstChannel[3] = { 0, 2400, 4800 };
  const int fLastChannel[3] = { 2399, 4799, 8255 };
  const double minTpcBoundd[3] = { 10., -105.53, 10.1 };
  const double maxTpcBound[3] = { 246.35, 107.47, 1026.9 };

  int fEvent;
  int fAPA;
  int fNumberWires;
  int fNumberTicks;

  bool fOriginIsInsideTPC = true;
  bool fEndIsInsideTPC = true;

  std::map<int,std::vector<float> > fWireMap;
  std::vector<std::vector<float> > fImage;

  geo::GeometryCore const*             fGeometry;           ///< pointer to the Geometry service
  detinfo::DetectorProperties const* fDetectorProperties; ///< Pointer to the detector properties
}; // EOClass LArCVMaker

LArCVMaker::LArCVMaker(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMgr(larcv::IOManager::kWRITE),
    fWireModuleLabel(pset.get<std::string>("WireModuleLabel")),
    fMcTruthModuleLabel(pset.get<std::string>("McTruthModuleLabel")),
    fMaxTick(pset.get<int>("MaxTick")),
    fADCCut(pset.get<int>("ADCCut")),
    fBoxSizeX(pset.get<double>("BoxSizeX")),
    fBoxSizeY(pset.get<double>("BoxSizeY")),
    fBoxSizeZ(pset.get<double>("BoxSizeZ")),
    fEventType(pset.get<int>("EventType"))

{} // EOFunction LArCVMaker::LArCVMaker

void LArCVMaker::beginJob() {
  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  std::string filename;
  if (std::getenv("PROCESS") != nullptr) filename = "larcv_" + std::string(std::getenv("PROCESS")) + ".root";
  else filename = "larcv.root";
  fMgr.set_out_file(filename);
  fMgr.initialize();
} // EOFunction LArCVMaker::beginJob

void LArCVMaker::endJob() {
  fMgr.finalize();
} // EOFunction LArCVMaker::endJob

void LArCVMaker::ClearData() {
  fWireMap.clear();
} // EOFunction LArCVMaker::ClearData

bool LArCVMaker::IsInsideTPC(double xyz[]) {
  bool isInsideX = (xyz[0]>minTpcBoundd[0] && xyz[0]<maxTpcBound[0]);
  bool isInsideY = (xyz[1]>-minTpcBoundd[1] && xyz[1]<maxTpcBound[1]);
  bool isInsideZ = (xyz[2]>minTpcBoundd[2] && xyz[2]<maxTpcBound[2]);
  if (isInsideX && isInsideY && isInsideZ) return true;
  else return false;
}

void LArCVMaker::BoxCoordinates(const simb::MCParticle& part, double xyz_origin[], double xyz_end[]) {
  if ((part.Vx()-fBoxSizeX/2.) < minTpcBoundd[0]) xyz_origin[0] = minTpcBoundd[0];
  else xyz_origin[0] = part.Vx()-fBoxSizeX/2.;
  if ((part.Vy()-fBoxSizeY/2.) < minTpcBoundd[1]) xyz_origin[1] = minTpcBoundd[1];
  else xyz_origin[1] = part.Vy()-fBoxSizeY/2.;
  if ((part.Vx()-fBoxSizeX/2.) < minTpcBoundd[2]) xyz_origin[2] = minTpcBoundd[2];
  else xyz_origin[2] = part.Vz()-fBoxSizeZ/2.;

  if ((part.Vx()+fBoxSizeX/2.) > maxTpcBound[0]){
    xyz_end[0] = maxTpcBound[0];
    xyz_origin[0] = xyz_end[0] - fBoxSizeX;
  }
  else xyz_end[0] = xyz_origin[0] + fBoxSizeX;

  if ((part.Vy()+fBoxSizeY/2.) > maxTpcBound[1]){
    xyz_end[1] = maxTpcBound[1];
    xyz_origin[1] = xyz_end[1] - fBoxSizeY;
  }
  else xyz_end[1] = xyz_origin[1] + fBoxSizeY;

  if ((part.Vz()+fBoxSizeZ/2.) > maxTpcBound[2]){
    xyz_end[2] = maxTpcBound[2];
    xyz_origin[2] = xyz_end[2] - fBoxSizeZ;
  }
  else xyz_end[2] = xyz_origin[2] + fBoxSizeZ;

  return;
}

void LArCVMaker::analyze(art::Event const & evt) {
  std::cout << "Changed to Thomas!" << std::endl;
  // Clear data before starting analyzing event
  ClearData();

  // Get event number
  fEvent = evt.event();
  // Set larcv manager
  fMgr.set_id(evt.id().run(),evt.id().subRun(),evt.id().event());

  // Get objects from event
  art::Handle<std::vector<recob::Wire>> wireHandle;
  art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
  evt.getByLabel(fWireModuleLabel,wireHandle);
  evt.getByLabel(fMcTruthModuleLabel,mctruthHandle);

  // Code is built on assumption that there's only 1 mcTruth for event.
  // Return error if that's not the case
  if ((*mctruthHandle).size()!=1){
    std::ostringstream oss;
    oss << "Unexpected number of mcTruth objects (" << (*mctruthHandle).size() << ") in event " << fEvent << std::endl;
    std::string ss = oss.str();
    throw std::invalid_argument(ss);
  }

  // Take first MC particle (make sure it's a neutrino) and take coordinates
  for (auto const& mctruth : (*mctruthHandle)){
    const simb::MCParticle& part = mctruth.GetParticle(0);

    // Report error if not neutrino
    if (part.PdgCode()!=14){
      std::ostringstream oss;
      oss << "First McParticle not a neutrino (PDG: " << part.PdgCode() << ") in event " << fEvent << std::endl;
      std::string ss = oss.str();
      throw std::invalid_argument(ss);
    }

    // Store coordinates in vector from first particle
    double xyz_origin[3], xyz_end[3];
    // geo::WireID wires_origin[3], wires_end[3];
    LArCVMaker::BoxCoordinates(part, xyz_origin, xyz_end);
    std::cout << "Event origin: (" << part.Vx() << ", " << part.Vy() << ", " << part.Vz() << ")" << std::endl;
    std::cout << "Box start: (" << xyz_origin[0] << ", " << xyz_origin[1] << ", " << xyz_origin[2] << ")" << std::endl;
    std::cout << "Box end: (" << xyz_end[0] << ", " << xyz_end[1] << ", " << xyz_end[2] << ")" << std::endl;
    fOriginIsInsideTPC = LArCVMaker::IsInsideTPC(xyz_origin);
    fEndIsInsideTPC = LArCVMaker::IsInsideTPC(xyz_end);
    if (!fOriginIsInsideTPC || !fEndIsInsideTPC){
      std::cout << "Event originates or ends outside TPC!" << std::endl;
      std::cout << "Moving on to next event." << std::endl;
    }
    // The actual code happens here!
    // else{
    //   wires_start[0] = fGeometry->NearestWireID(xyz_start,0);
    //   wires_start[1] = fGeometry->NearestWireID(xyz_start,1);
    //   wires_start[2] = fGeometry->NearestWireID(xyz_start,2);
    //   wires_end[0] = fGeometry->NearestWireID(xyz_end,0);
    //   wires_end[1] = fGeometry->NearestWireID(xyz_end,1);
    //   wires_end[2] = fGeometry->NearestWireID(xyz_end,2);
    //   std::cout << "Nearest origin wires: (" << wires_start[0].Wire << ", " << wires_start[1].Wire << ", " << wires_start[2].Wire << ")" << std::endl;
    //   std::cout << "Nearest end wires: (" << wires_end[0].Wire << ", " << wires_end[1].Wire << ", " << wires_end[2].Wire << ")" << std::endl;

    //   // Loop over each wire and add it to the wire map
    //   for (std::vector<recob::Wire>::const_iterator it = wireHandle->begin(); it != wireHandle->end(); ++it) {
    //     const recob::Wire & wire = *it;
    //     fWireMap.insert(std::pair<int,std::vector<float> >(wire.Channel(),std::vector<float>(wire.Signal())));
    //   }

    //   // Get handle on larcv image
    //   auto images = (larcv::EventImage2D*)(fMgr.get_data(larcv::kProductImage2D, "tpc"));

    //   // Create images from the wire map
    //   for (int it_plane = 0; it_plane < 3; ++it_plane) {
    //     larcv::Image2D image(fNumberChannels[it_plane],fMaxTick);
    //     for (int it_channel = 0; it_channel < fNumberChannels[it_plane]; ++it_channel) {
    //       int channel = it_channel + fFirstChannel[it_plane];
    //       for (int it_tick = 0; it_tick < fMaxTick; ++it_tick) {
    //         if (fWireMap.find(channel) != fWireMap.end()) image.set_pixel(it_channel,it_tick,fWireMap[channel][it_tick]);
    //         else image.set_pixel(it_channel,it_tick,0);
    //       }
    //     }
    //     if (it_plane < 2) {
    //       image.compress(600,640);
    //     }
    //     else {
    //       image.compress(576,640);
    //       image.resize(600,640,0);
    //     }
    //     images->Emplace(std::move(image));
    //   }

    //   auto roi = (larcv::EventROI*)(fMgr.get_data(larcv::kProductROI, "tpc"));
    //   roi->Emplace(larcv::ROI((larcv::ROIType_t)fEventType));

    //   fMgr.save_entry();
    // }//EOExecution over MCParticle
  }//EOLoop over MCTruth
} //EOFunction LArCVMaker::analyze

DEFINE_ART_MODULE(LArCVMaker)

} // namespace larhsn

