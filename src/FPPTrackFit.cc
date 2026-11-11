#include "FPPTrackFit.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace {
bool solve4x4(std::array<std::array<double, 5>, 4>& augmented,
              std::array<double, 4>& solution) {
  for (int column = 0; column < 4; ++column) {
    int pivot = column;
    for (int row = column + 1; row < 4; ++row) {
      if (std::abs(augmented[row][column]) >
          std::abs(augmented[pivot][column])) {
        pivot = row;
      }
    }
    if (std::abs(augmented[pivot][column]) < 1.0e-12) return false;
    if (pivot != column) std::swap(augmented[pivot], augmented[column]);

    const double divisor = augmented[column][column];
    for (int entry = column; entry < 5; ++entry) {
      augmented[column][entry] /= divisor;
    }
    for (int row = 0; row < 4; ++row) {
      if (row == column) continue;
      const double factor = augmented[row][column];
      for (int entry = column; entry < 5; ++entry) {
        augmented[row][entry] -= factor * augmented[column][entry];
      }
    }
  }
  for (int index = 0; index < 4; ++index) solution[index] = augmented[index][4];
  return true;
}

FPPWireTrack fitCoordinates(const std::vector<FPPPlaneHit>& hits,
                            const std::vector<double>& coordinateOffsetsCm,
                            double varianceCm2,
                            const FPPPlaneArray& planes) {
  FPPWireTrack result;
  result.hitCount = static_cast<int>(hits.size());
  if (hits.size() < 5 || coordinateOffsetsCm.size() != hits.size() ||
      varianceCm2 <= 0.0) {
    return result;
  }

  double referenceZ = 0.0;
  for (const auto& hit : hits) {
    if (hit.planeIndex < 0 ||
        hit.planeIndex >= static_cast<int>(planes.size())) {
      return result;
    }
    referenceZ += planes[hit.planeIndex].zCm;
  }
  referenceZ /= static_cast<double>(hits.size());

  std::array<std::array<double, 4>, 4> normal{};
  std::array<double, 4> right{};
  const double weight = 1.0 / varianceCm2;
  for (std::size_t index = 0; index < hits.size(); ++index) {
    const auto& hit = hits[index];
    const auto& plane = planes[hit.planeIndex];
    if (hit.wire < 1 || hit.wire > plane.wireCount) return result;
    const double wireCoordinate =
        plane.layerOffsetCm + plane.wirePitchCm * hit.wire;
    const double measured =
        plane.measurementX * plane.chamberXCm +
        plane.measurementY * plane.chamberYCm + wireCoordinate +
        coordinateOffsetsCm[index];
    const double deltaZ = plane.zCm - referenceZ;
    const std::array<double, 4> row = {
        plane.measurementX * deltaZ,
        plane.measurementX,
        plane.measurementY * deltaZ,
        plane.measurementY,
    };
    for (int i = 0; i < 4; ++i) {
      right[i] += weight * row[i] * measured;
      for (int j = 0; j < 4; ++j) {
        normal[i][j] += weight * row[i] * row[j];
      }
    }
  }

  std::array<std::array<double, 5>, 4> augmented{};
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) augmented[i][j] = normal[i][j];
    augmented[i][4] = right[i];
  }
  std::array<double, 4> fitted{};
  if (!solve4x4(augmented, fitted)) return result;

  result.slopeX = fitted[0];
  result.interceptX = fitted[1] - fitted[0] * referenceZ;
  result.slopeY = fitted[2];
  result.interceptY = fitted[3] - fitted[2] * referenceZ;

  double chi2 = 0.0;
  for (std::size_t index = 0; index < hits.size(); ++index) {
    const auto& hit = hits[index];
    const auto& plane = planes[hit.planeIndex];
    const double wireCoordinate =
        plane.layerOffsetCm + plane.wirePitchCm * hit.wire;
    const double measured =
        plane.measurementX * plane.chamberXCm +
        plane.measurementY * plane.chamberYCm + wireCoordinate +
        coordinateOffsetsCm[index];
    const double predictedX = result.slopeX * plane.zCm + result.interceptX;
    const double predictedY = result.slopeY * plane.zCm + result.interceptY;
    const double predicted =
        plane.measurementX * predictedX + plane.measurementY * predictedY;
    const double residual = measured - predicted;
    chi2 += residual * residual / varianceCm2;
  }
  result.reducedChi2 = chi2 / static_cast<double>(hits.size() - 4);
  result.valid = true;
  return result;
}
}  // namespace

FPPWireTrack fitFPPWireCenters(const std::vector<FPPPlaneHit>& hits,
                               double varianceCm2) {
  return fitCoordinates(hits, std::vector<double>(hits.size(), 0.0),
                        varianceCm2, hallCFPPPlanes());
}

FPPWireTrack fitFPPWireCenters(const std::vector<FPPPlaneHit>& hits,
                               const FPPPlaneArray& planes,
                               double varianceCm2) {
  return fitCoordinates(hits, std::vector<double>(hits.size(), 0.0),
                        varianceCm2, planes);
}

FPPLeftRightTrack fitFPPLeftRightPermutations(
    const std::vector<FPPPlaneHit>& hits, double varianceCm2) {
  return fitFPPLeftRightPermutations(hits, hallCFPPPlanes(), varianceCm2);
}

FPPLeftRightTrack fitFPPLeftRightPermutations(
    const std::vector<FPPPlaneHit>& hits, const FPPPlaneArray& planes,
    double varianceCm2) {
  FPPLeftRightTrack best;
  if (hits.size() < 5 || hits.size() > 20 || varianceCm2 <= 0.0) return best;

  const unsigned long permutations = 1UL << hits.size();
  best.permutationsTried = permutations;
  double bestChi2 = std::numeric_limits<double>::infinity();
  std::vector<double> signedDrifts(hits.size(), 0.0);
  for (unsigned long permutation = 0; permutation < permutations;
       ++permutation) {
    for (std::size_t index = 0; index < hits.size(); ++index) {
      const double sign = (permutation & (1UL << index)) ? 1.0 : -1.0;
      signedDrifts[index] = sign * hits[index].driftDistanceCm;
    }
    const auto candidate =
        fitCoordinates(hits, signedDrifts, varianceCm2, planes);
    if (candidate.valid && candidate.reducedChi2 < bestChi2) {
      bestChi2 = candidate.reducedChi2;
      best.track = candidate;
      best.signedDriftCm = signedDrifts;
    }
  }
  return best;
}
