#pragma once

#include "FPPGeometry.hh"

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

class DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction(double analyzer1ThicknessCm = 55.88,
                       bool includeAnalyzer2 = true,
                       FPPLayout layout = legacyFPPLayout())
      : analyzer1ThicknessCm_(analyzer1ThicknessCm),
        includeAnalyzer2_(includeAnalyzer2), layout_(layout) {}
  G4VPhysicalVolume* Construct() override;

 private:
  double analyzer1ThicknessCm_;
  bool includeAnalyzer2_;
  FPPLayout layout_;
};
