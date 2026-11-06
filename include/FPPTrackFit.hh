#pragma once

#include "FPPGeometry.hh"

#include <vector>

struct FPPWireTrack {
  bool valid = false;
  double slopeX = 0.0;       // x(z) = slopeX*z + interceptX
  double interceptX = 0.0;
  double slopeY = 0.0;       // y(z) = slopeY*z + interceptY
  double interceptY = 0.0;
  double reducedChi2 = 0.0;
  int hitCount = 0;
};

struct FPPLeftRightTrack {
  FPPWireTrack track;
  std::vector<double> signedDriftCm;  // aligned with the input hit order
  unsigned long permutationsTried = 0;
};

// Reproduce the Hall C rough-track stage using wire centers only.  The
// historical fitter passed the 2 cm wire spacing directly as sigma^2, so the
// default variance preserves that chi-squared convention.
FPPWireTrack fitFPPWireCenters(const std::vector<FPPPlaneHit>& hits,
                               double varianceCm2 = 2.0);

// Try both signs of every unsigned drift distance and retain the track with
// minimum reduced chi squared.  The default 0.005 cm^2 is the numerical Hall C
// HFPP_resolution value used by the legacy left/right fit.
FPPLeftRightTrack fitFPPLeftRightPermutations(
    const std::vector<FPPPlaneHit>& hits, double varianceCm2 = 0.005);
