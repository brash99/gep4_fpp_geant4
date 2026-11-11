#include "FPPGeometry.hh"
#include "FPPTrackFit.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {
void require(bool condition, const std::string& message) {
  if (!condition) {
    std::cerr << "FAILED: " << message << '\n';
    std::exit(1);
  }
}

void requireClose(double value, double expected, double tolerance,
                  const std::string& message) {
  require(std::abs(value - expected) <= tolerance, message);
}
}  // namespace

int main() {
  const auto& planes = hallCFPPPlanes();
  require(planes.size() == 12, "Hall C FPP must contain 12 measurement planes");
  require(planes[0].chamber == 1 && planes[0].zCm == 176.8,
          "first U-plane geometry");
  require(planes[11].chamber == 4 && planes[11].zCm == 298.8,
          "last V-plane geometry");

  const auto compactPlanes = makeHallCFPPPlanes(compactDownstreamFPPLayout());
  requireClose(compactPlanes[0].zCm, 231.698, 1.0e-12,
               "compact FPP1 U-plane geometry");
  requireClose(compactPlanes[3].zCm, 252.998, 1.0e-12,
               "compact FPP2 U-plane geometry");
  requireClose(compactPlanes[6].zCm, 274.3, 1.0e-12,
               "compact layout retains FPP3");
  requireClose(compactPlanes[11].zCm, 298.8, 1.0e-12,
               "compact layout retains FPP4");

  const FPPVector3 forward{0.0, 0.0, 1.0};
  const FPPVector3 chamberCenter{-5.08, 0.0, 0.0};
  for (const auto& plane : planes) {
    const auto hit = digitizeFPPLine(plane, chamberCenter, forward);
    require(hit.valid, "central track must digitize in every plane");
    const int expectedWire =
        plane.orientation == FPPPlaneOrientation::X ? 42 : 52;
    require(hit.wire == expectedWire, "central legacy wire number");
    requireClose(hit.driftDistanceCm, 0.0, 1.0e-12,
                 "central track drift distance");
  }

  // Direct checks of the legacy get_wire_numbers formulas at global x=y=0:
  // U,V: round((0 + 3.592 + 104)/2)=54; X: round((0 + 5.08 + 84)/2)=45.
  const FPPVector3 globalOrigin{0.0, 0.0, 0.0};
  require(digitizeFPPLine(planes[0], globalOrigin, forward).wire == 54,
          "legacy U mapping at global origin");
  require(digitizeFPPLine(planes[1], globalOrigin, forward).wire == 45,
          "legacy X mapping at global origin");
  require(digitizeFPPLine(planes[2], globalOrigin, forward).wire == 54,
          "legacy V mapping at global origin");

  // A vertical track displaced by 0.25 cm along the X-plane measurement
  // coordinate has a 0.25 cm unsigned closest approach to central wire 42.
  const FPPVector3 displaced{-4.83, 0.0, 0.0};
  const auto displacedHit = digitizeFPPLine(planes[1], displaced, forward);
  require(displacedHit.valid && displacedHit.wire == 42,
          "sub-cell displacement remains on central X wire");
  requireClose(displacedHit.driftDistanceCm, 0.25, 1.0e-12,
               "unsigned line-to-wire drift distance");

  // Verify the one-based endpoints of the authoritative wire tables.
  const auto& uPlane = planes[0];
  FPPVector3 firstU{uPlane.chamberXCm + uPlane.measurementX * -102.0,
                    uPlane.chamberYCm + uPlane.measurementY * -102.0, 0.0};
  FPPVector3 lastU{uPlane.chamberXCm + uPlane.measurementX * 104.0,
                   uPlane.chamberYCm + uPlane.measurementY * 104.0, 0.0};
  require(digitizeFPPLine(uPlane, firstU, forward).wire == 1,
          "first U wire mapping");
  require(digitizeFPPLine(uPlane, lastU, forward).wire == 104,
          "last U wire mapping");

  const FPPVector3 outside{-200.0, 0.0, 0.0};
  require(!digitizeFPPLine(planes[1], outside, forward).valid,
          "track outside wire coverage must not digitize");

  const auto allHits = digitizeFPPTrack(chamberCenter, forward);
  require(allHits.size() == 12, "central track must yield twelve plane hits");
  require(allHits.front().planeIndex == 0 && allHits.back().planeIndex == 11,
          "digitized hit ordering follows the geometry table");
  const auto firstPairHits = digitizeFPPTrack(chamberCenter, forward, 0x3);
  require(firstPairHits.size() == 6,
          "FPP1+FPP2 mask must yield a six-plane reconstruction input");
  for (const auto& hit : firstPairHits) {
    require(hit.chamber <= 2, "first-pair hit must not retain rear chambers");
  }

  const auto centralTrack = fitFPPWireCenters(firstPairHits);
  require(centralTrack.valid && centralTrack.hitCount == 6,
          "six central wire centers must form a rough track");
  requireClose(centralTrack.slopeX, 0.0, 1.0e-12, "central track x slope");
  requireClose(centralTrack.interceptX, -5.08, 1.0e-10,
               "central track x intercept");
  requireClose(centralTrack.slopeY, 0.0, 1.0e-12, "central track y slope");
  requireClose(centralTrack.interceptY, 0.0, 1.0e-10,
               "central track y intercept");
  requireClose(centralTrack.reducedChi2, 0.0, 1.0e-12,
               "central track reduced chi squared");

  auto fiveHits = firstPairHits;
  fiveHits.pop_back();
  require(fitFPPWireCenters(fiveHits).valid,
          "legacy minimum of five planes must remain fit-capable");
  std::vector<FPPPlaneHit> underconstrained(firstPairHits.begin(),
                                            firstPairHits.begin() + 4);
  require(!fitFPPWireCenters(underconstrained).valid,
          "fewer than five planes must not form an accepted rough track");

  // A small-angle displaced track is genuinely left/right ambiguous: the
  // fitter must find a zero-chi2 solution, but it is not required to recover
  // the truth branch when several sign patterns are exactly degenerate.
  const auto displacedPairHits = digitizeFPPTrack(displaced, forward, 0x3);
  const auto ambiguous = fitFPPLeftRightPermutations(displacedPairHits);
  require(ambiguous.track.valid && ambiguous.permutationsTried == 64,
          "six hits must exercise all 64 left/right permutations");
  requireClose(ambiguous.track.reducedChi2, 0.0, 1.0e-10,
               "small-angle ambiguous solution chi squared");

  // At a larger angle the changing wire pattern breaks the degeneracy.  The
  // line is chosen to cross chamber 1 near its center.
  const FPPVector3 angledPoint{-26.68, -10.8, 0.0};
  const FPPVector3 angledDirection{0.12, 0.06, 1.0};
  const auto angledHits = digitizeFPPTrack(angledPoint, angledDirection, 0x3);
  const auto leftRight = fitFPPLeftRightPermutations(angledHits);
  require(leftRight.track.valid && leftRight.permutationsTried == 64,
          "six hits must exercise all 64 left/right permutations");
  requireClose(leftRight.track.slopeX, 0.12, 2.0e-3,
               "left/right fitted x slope");
  requireClose(leftRight.track.interceptX, -26.68, 0.4,
               "left/right fitted x intercept");
  requireClose(leftRight.track.slopeY, 0.06, 2.0e-3,
               "left/right fitted y slope");
  requireClose(leftRight.track.interceptY, -10.8, 0.4,
               "left/right fitted y intercept");
  require(leftRight.track.reducedChi2 < 1.0e-3,
          "correct large-angle solution chi squared");
  require(leftRight.signedDriftCm.size() == angledHits.size(),
          "one selected signed drift per input hit");

  // The four-chamber lever arm must resolve a small projected angle using all
  // twelve U/X/V measurements.  This is the topology needed below 5 degrees.
  const FPPVector3 smallAnglePoint{-7.08, -1.0, 0.0};
  const FPPVector3 smallAngleDirection{0.012, 0.006, 1.0};
  const auto allFourHits =
      digitizeFPPTrack(smallAnglePoint, smallAngleDirection, 0xf);
  const auto allFourFit = fitFPPLeftRightPermutations(allFourHits);
  require(allFourFit.track.valid && allFourFit.permutationsTried == 4096,
          "twelve hits must exercise all 4096 left/right permutations");
  requireClose(allFourFit.track.slopeX, 0.012, 2.0e-3,
               "four-chamber small-angle x slope");
  requireClose(allFourFit.track.slopeY, 0.006, 2.0e-3,
               "four-chamber small-angle y slope");
  require(allFourFit.track.reducedChi2 < 1.0e-3,
          "four-chamber small-angle solution chi squared");

  const auto compactHits = digitizeFPPTrack(
      smallAnglePoint, smallAngleDirection, compactPlanes, 0xf);
  const auto compactFit =
      fitFPPLeftRightPermutations(compactHits, compactPlanes);
  require(compactHits.size() == 12,
          "compact layout must retain all twelve measurements");
  require(compactFit.track.valid && compactFit.permutationsTried == 4096,
          "compact twelve-hit track must exercise all permutations");
  requireClose(compactFit.track.slopeX, 0.012, 2.0e-3,
               "compact-layout small-angle x slope");
  requireClose(compactFit.track.slopeY, 0.006, 2.0e-3,
               "compact-layout small-angle y slope");
  require(compactFit.track.reducedChi2 < 1.0e-3,
          "compact-layout small-angle solution chi squared");

  std::cout << "FPP geometry and digitizer checks passed\n";
  return 0;
}
