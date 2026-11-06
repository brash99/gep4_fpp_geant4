#include "FPPGeometry.hh"

#include <cmath>
#include <limits>

namespace {
constexpr double kInvSqrt2 = 0.70710678118654752440;

FPPPlane makePlane(int chamber, FPPPlaneOrientation orientation, double zCm) {
  FPPPlane plane;
  plane.chamber = chamber;
  plane.orientation = orientation;
  plane.zCm = zCm;
  plane.wirePitchCm = 2.0;
  plane.chamberXCm = -5.08;
  plane.chamberYCm = 0.0;
  switch (orientation) {
    case FPPPlaneOrientation::U:
      plane.wireCount = 104;
      plane.layerOffsetCm = -104.0;
      plane.measurementX = kInvSqrt2;
      plane.measurementY = kInvSqrt2;
      break;
    case FPPPlaneOrientation::X:
      plane.wireCount = 83;
      plane.layerOffsetCm = -84.0;
      plane.measurementX = 1.0;
      plane.measurementY = 0.0;
      break;
    case FPPPlaneOrientation::V:
      plane.wireCount = 104;
      plane.layerOffsetCm = -104.0;
      plane.measurementX = kInvSqrt2;
      plane.measurementY = -kInvSqrt2;
      break;
  }
  return plane;
}

FPPVector3 subtract(const FPPVector3& left, const FPPVector3& right) {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

FPPVector3 cross(const FPPVector3& left, const FPPVector3& right) {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

double dot(const FPPVector3& left, const FPPVector3& right) {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

double magnitude(const FPPVector3& value) {
  return std::sqrt(dot(value, value));
}
}  // namespace

const std::array<FPPPlane, 12>& hallCFPPPlanes() {
  // Hall C GEp-III values from Analyzer/replay/PARAM/hfpp.param.  The array
  // is ordered by physical chamber, then U/X/V measurement plane.
  static const std::array<FPPPlane, 12> planes = {
      makePlane(1, FPPPlaneOrientation::U, 176.8),
      makePlane(1, FPPPlaneOrientation::X, 178.4),
      makePlane(1, FPPPlaneOrientation::V, 180.0),
      makePlane(2, FPPPlaneOrientation::U, 198.1),
      makePlane(2, FPPPlaneOrientation::X, 199.7),
      makePlane(2, FPPPlaneOrientation::V, 201.3),
      makePlane(3, FPPPlaneOrientation::U, 274.3),
      makePlane(3, FPPPlaneOrientation::X, 275.9),
      makePlane(3, FPPPlaneOrientation::V, 277.5),
      makePlane(4, FPPPlaneOrientation::U, 295.6),
      makePlane(4, FPPPlaneOrientation::X, 297.3),
      makePlane(4, FPPPlaneOrientation::V, 298.8),
  };
  return planes;
}

FPPDigitizedHit digitizeFPPLine(const FPPPlane& plane,
                                const FPPVector3& pointCm,
                                const FPPVector3& direction) {
  FPPDigitizedHit hit;
  if (std::abs(direction.z) < 1.0e-12 || plane.wirePitchCm <= 0.0) {
    return hit;
  }

  const double alongTrack = (plane.zCm - pointCm.z) / direction.z;
  const FPPVector3 crossing = {
      pointCm.x + alongTrack * direction.x,
      pointCm.y + alongTrack * direction.y,
      plane.zCm,
  };
  hit.measurementCoordinateCm =
      plane.measurementX * (crossing.x - plane.chamberXCm) +
      plane.measurementY * (crossing.y - plane.chamberYCm);

  const double floatingWire =
      (hit.measurementCoordinateCm - plane.layerOffsetCm) /
      plane.wirePitchCm;
  hit.wire = static_cast<int>(std::floor(floatingWire + 0.5));
  if (hit.wire < 1 || hit.wire > plane.wireCount) {
    hit.wire = 0;
    return hit;
  }

  hit.wireCoordinateCm =
      plane.layerOffsetCm + plane.wirePitchCm * hit.wire;
  const FPPVector3 wirePoint = {
      plane.chamberXCm + plane.measurementX * hit.wireCoordinateCm,
      plane.chamberYCm + plane.measurementY * hit.wireCoordinateCm,
      plane.zCm,
  };
  const FPPVector3 wireDirection = {
      -plane.measurementY,
      plane.measurementX,
      0.0,
  };
  const auto normal = cross(direction, wireDirection);
  const double normalMagnitude = magnitude(normal);
  if (normalMagnitude < 1.0e-12) {
    hit.wire = 0;
    return hit;
  }
  hit.driftDistanceCm =
      std::abs(dot(subtract(pointCm, wirePoint), normal)) / normalMagnitude;
  hit.valid = true;
  return hit;
}

std::vector<FPPPlaneHit> digitizeFPPTrack(const FPPVector3& pointCm,
                                          const FPPVector3& direction,
                                          int chamberMask) {
  std::vector<FPPPlaneHit> hits;
  const auto& planes = hallCFPPPlanes();
  hits.reserve(planes.size());
  for (std::size_t index = 0; index < planes.size(); ++index) {
    const auto& plane = planes[index];
    const int chamberBit = 1 << (plane.chamber - 1);
    if ((chamberMask & chamberBit) == 0) continue;
    const auto digitized = digitizeFPPLine(plane, pointCm, direction);
    if (!digitized.valid) continue;
    hits.push_back({static_cast<int>(index), plane.chamber, plane.orientation,
                    digitized.wire, digitized.driftDistanceCm});
  }
  return hits;
}
