#pragma once

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;

class DetectorConstruction final : public G4VUserDetectorConstruction {
 public:
  DetectorConstruction(double analyzer1ThicknessCm = 55.88,
                       bool includeAnalyzer2 = true)
      : analyzer1ThicknessCm_(analyzer1ThicknessCm),
        includeAnalyzer2_(includeAnalyzer2) {}
  G4VPhysicalVolume* Construct() override;

 private:
  double analyzer1ThicknessCm_;
  bool includeAnalyzer2_;
};
