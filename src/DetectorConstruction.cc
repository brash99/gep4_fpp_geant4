#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include <stdexcept>

namespace {
void placeBox(const char* name, G4Material* material, G4LogicalVolume* world,
              double halfXcm, double halfYcm, double halfZcm,
              double xCm, double zCm) {
  auto* solid = new G4Box(name, halfXcm * cm, halfYcm * cm, halfZcm * cm);
  auto* logical = new G4LogicalVolume(solid, material, name);
  new G4PVPlacement(nullptr, {xCm * cm, 0., zCm * cm}, logical, name,
                    world, false, 0, true);
}
}  // namespace

G4VPhysicalVolume* DetectorConstruction::Construct() {
  constexpr double analyzer1UpstreamFaceCm = 111.4215;
  constexpr double fpp1UpstreamFaceCm = 172.8;
  if (analyzer1ThicknessCm_ <= 0.0 ||
      analyzer1UpstreamFaceCm + analyzer1ThicknessCm_ >= fpp1UpstreamFaceCm) {
    throw std::runtime_error(
        "Analyzer 1 thickness must be positive and leave clearance before FPP1");
  }
  auto* nist = G4NistManager::Instance();
  auto* air = nist->FindOrBuildMaterial("G4_AIR");
  auto* scintillator = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto* polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

  auto* worldSolid = new G4Box("World", 150. * cm, 130. * cm, 350. * cm);
  auto* worldLogical = new G4LogicalVolume(worldSolid, air, "World");
  auto* world = new G4PVPlacement(nullptr, {0., 0., 0.}, worldLogical,
                                  "World", nullptr, false, 0, true);

  // GEp-III baseline positions from legacy src/ugeom.f. Dimensions are half lengths.
  placeBox("S1", scintillator, worldLogical, 60.96, 40.64, 0.50, 0.0, 63.0995);
  placeBox("S2", scintillator, worldLogical, 60.96, 40.64, 0.50, 0.0, 82.6935);
  placeBox("Analyzer1", polyethylene, worldLogical, 71.12, 53.975,
           analyzer1ThicknessCm_ / 2.0, -5.08,
           analyzer1UpstreamFaceCm + analyzer1ThicknessCm_ / 2.0);
  placeBox("FPP1", air, worldLogical, 83.0, 67.0, 0.05, -5.08, 178.4);
  placeBox("FPP2", air, worldLogical, 83.0, 67.0, 0.05, -5.08, 199.7);
  if (includeAnalyzer2_) {
    placeBox("Analyzer2", polyethylene, worldLogical, 71.12, 53.975, 27.94,
             -5.08, 237.694);
  }
  placeBox("FPP3", air, worldLogical, 83.0, 67.0, 0.05, -5.08, 275.898);
  placeBox("FPP4", air, worldLogical, 83.0, 67.0, 0.05, -5.08, 297.198);

  return world;
}
