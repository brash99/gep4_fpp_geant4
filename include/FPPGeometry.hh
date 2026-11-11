#pragma once

#include <array>
#include <vector>

enum class FPPPlaneOrientation { U, X, V };

struct FPPVector3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct FPPPlane {
  int chamber = 0;                 // one-based FPP chamber number
  FPPPlaneOrientation orientation = FPPPlaneOrientation::X;
  double zCm = 0.0;
  int wireCount = 0;
  double wirePitchCm = 0.0;
  double layerOffsetCm = 0.0;      // wire coordinate = offset + pitch*wire
  double chamberXCm = 0.0;
  double chamberYCm = 0.0;
  double measurementX = 1.0;       // unit vector normal to wires in x-y
  double measurementY = 0.0;
};

struct FPPDigitizedHit {
  bool valid = false;
  int wire = 0;                    // one-based legacy wire number
  double driftDistanceCm = 0.0;    // unsigned 3-D line-to-wire distance
  double measurementCoordinateCm = 0.0;
  double wireCoordinateCm = 0.0;
};

struct FPPPlaneHit {
  int planeIndex = -1;             // zero-based index in hallCFPPPlanes()
  int chamber = 0;                 // one-based chamber number
  FPPPlaneOrientation orientation = FPPPlaneOrientation::X;
  int wire = 0;                    // one-based wire number
  double driftDistanceCm = 0.0;    // unsigned; no truth coordinate retained
};

using FPPPlaneArray = std::array<FPPPlane, 12>;

struct FPPLayout {
  std::array<double, 4> chamberCenterZCm{};
};

const FPPLayout& legacyFPPLayout();
const FPPLayout& compactDownstreamFPPLayout();
FPPPlaneArray makeHallCFPPPlanes(const FPPLayout& layout);
const FPPPlaneArray& hallCFPPPlanes();

FPPDigitizedHit digitizeFPPLine(const FPPPlane& plane,
                                const FPPVector3& pointCm,
                                const FPPVector3& direction);

std::vector<FPPPlaneHit> digitizeFPPTrack(const FPPVector3& pointCm,
                                          const FPPVector3& direction,
                                          int chamberMask = 0xf);
std::vector<FPPPlaneHit> digitizeFPPTrack(const FPPVector3& pointCm,
                                          const FPPVector3& direction,
                                          const FPPPlaneArray& planes,
                                          int chamberMask = 0xf);
