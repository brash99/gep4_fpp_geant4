#include "EventAction.hh"

#include "FPPGeometry.hh"
#include "FPPTrackFit.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

namespace {
constexpr double kChamberCenterXCm = -5.08;
constexpr double kChamberHalfXCm = 83.0;
constexpr double kChamberHalfYCm = 67.0;
constexpr double kFpp2ZCm = 199.7;
constexpr double kFpp4ZCm = 297.198;
}  // namespace

void EventAction::BeginOfEventAction(const G4Event*) {
  // GeneratePrimaries supplies the input event before Geant4 invokes this
  // callback.  Preserve those fields here; resetting them made event_id,
  // weight, and incident_p_GeV identically default-valued in the output.
  analyzer_ = 0;
  processCategory_ = 0;
  interactionZCm_ = 0.0;
  thetaDeg_ = 0.0;
  outgoingMomentumGeV_ = 0.0;
  chargedSecondaries_ = 0;
  primaryMask_ = 0;
  secondaryMask_ = 0;
  tracks_.clear();
}

void EventAction::setInput(long id, double weight, double momentumGeV,
                           const G4ThreeVector& positionCm,
                           const G4ThreeVector& direction) {
  eventId_ = id;
  weight_ = weight;
  momentumGeV_ = momentumGeV;
  incidentPositionCm_ = positionCm;
  incidentDirection_ = direction;
}

void EventAction::recordHadronicInteraction(int analyzer, int processCategory,
                                            double zCm, double thetaDeg,
                                            double outgoingMomentumGeV,
                                            int chargedSecondaries) {
  if (analyzer_ != 0) return;
  analyzer_ = analyzer;
  processCategory_ = processCategory;
  interactionZCm_ = zCm;
  thetaDeg_ = thetaDeg;
  outgoingMomentumGeV_ = outgoingMomentumGeV;
  chargedSecondaries_ = chargedSecondaries;
}

void EventAction::recordChamberCrossing(int chamber, int trackId, int parentId,
                                        bool primary,
                                        const G4ThreeVector& direction,
                                        const G4ThreeVector& positionCm,
                                        double momentumGeV) {
  const int bit = 1 << (chamber - 1);
  auto& observation = tracks_[trackId];
  observation.parentId = parentId;
  observation.mask |= bit;
  observation.direction.at(chamber - 1) = direction;
  observation.positionCm.at(chamber - 1) = positionCm;
  observation.momentumGeV.at(chamber - 1) = momentumGeV;
  if (primary) primaryMask_ |= bit;
  else secondaryMask_ |= bit;
}

bool EventAction::legacyConeTest(const G4ThreeVector& incomingPointCm,
                                 const G4ThreeVector& incomingDirection,
                                 double thetaDeg, double rearPlaneZCm) {
  // Port of the GEp-III cone_test acceptance definition.  At fixed polar
  // scattering angle, require the complete 0--2pi azimuthal cone to remain
  // inside the rectangular active area of the last required rear chamber.
  const auto direction = incomingDirection.unit();
  if (direction.z() <= 0.0) return false;
  const double pathLengthCm =
      (rearPlaneZCm - incomingPointCm.z()) / direction.z();
  if (pathLengthCm <= 0.0) return false;
  const auto projected = incomingPointCm + pathLengthCm * direction;
  const double radiusCm = pathLengthCm * std::tan(thetaDeg * deg);
  return projected.x() - radiusCm >= kChamberCenterXCm - kChamberHalfXCm &&
         projected.x() + radiusCm <= kChamberCenterXCm + kChamberHalfXCm &&
         projected.y() - radiusCm >= -kChamberHalfYCm &&
         projected.y() + radiusCm <= kChamberHalfYCm;
}

EventAction::Candidate EventAction::bestCandidate(int analyzer) const {
  const int firstChamber = analyzer == 2 ? 2 : 0;
  const int requiredMask = analyzer == 1 ? 0x3 : (analyzer == 2 ? 0xc : 0xf);
  Candidate best;

  for (const auto& [trackId, track] : tracks_) {
    (void)trackId;
    if ((track.mask & requiredMask) != requiredMask) continue;

    G4ThreeVector incoming = incidentDirection_;
    G4ThreeVector incomingPoint = incidentPositionCm_;
    if (analyzer == 2) {
      // A particle continuing from the first chamber pair carries its own
      // upstream direction.  A new secondary uses its parent's FPP2 direction.
      if ((track.mask & 0x2) != 0) {
        incoming = track.direction[1];
        incomingPoint = track.positionCm[1];
      } else {
        const auto parent = tracks_.find(track.parentId);
        if (parent != tracks_.end() && (parent->second.mask & 0x2) != 0) {
          incoming = parent->second.direction[1];
          incomingPoint = parent->second.positionCm[1];
        }
      }
    }

    const auto outgoing = track.direction[firstChamber];
    const auto outgoingPoint = track.positionCm[firstChamber];
    const double cosine = std::clamp(incoming.unit().dot(outgoing.unit()), -1.0, 1.0);
    const double theta = std::acos(cosine) / deg;
    // Preserve the GEp-III convention from calc_theta_phi: azimuth is
    // atan2(delta horizontal slope, delta vertical slope), wrapped to 0--2pi.
    const double deltaX = outgoing.x() / outgoing.z() - incoming.x() / incoming.z();
    const double deltaY = outgoing.y() / outgoing.z() - incoming.y() / incoming.z();
    double phi = std::atan2(deltaX, deltaY);
    if (phi < 0.0) phi += 2.0 * CLHEP::pi;
    const double momentum = track.momentumGeV[firstChamber];
    const double rearPlaneZ = analyzer == 1 ? kFpp2ZCm : kFpp4ZCm;
    const bool conePass =
        legacyConeTest(incomingPoint, incoming, theta, rearPlaneZ);

    bool reconstructedValid = false;
    double reconstructedTheta = 0.0;
    double reconstructedPhi = 0.0;
    double reconstructedChi2 = 0.0;
    int reconstructionPermutations = 0;
    if (analyzer == 1 || analyzer == 3) {
      const FPPVector3 digitizerPoint{outgoingPoint.x(), outgoingPoint.y(),
                                      outgoingPoint.z()};
      const FPPVector3 digitizerDirection{outgoing.x(), outgoing.y(),
                                          outgoing.z()};
      const int chamberMask = analyzer == 3 ? 0xf : 0x3;
      const auto hits = digitizeFPPTrack(digitizerPoint, digitizerDirection,
                                         chamberMask);
      const auto reconstructed = fitFPPLeftRightPermutations(hits);
      if (reconstructed.track.valid) {
        const G4ThreeVector reconstructedDirection{
            reconstructed.track.slopeX, reconstructed.track.slopeY, 1.0};
        const double reconstructedCosine =
            std::clamp(incoming.unit().dot(reconstructedDirection.unit()),
                       -1.0, 1.0);
        reconstructedTheta = std::acos(reconstructedCosine) / deg;
        const double reconstructedDeltaX =
            reconstructed.track.slopeX - incoming.x() / incoming.z();
        const double reconstructedDeltaY =
            reconstructed.track.slopeY - incoming.y() / incoming.z();
        reconstructedPhi =
            std::atan2(reconstructedDeltaX, reconstructedDeltaY);
        if (reconstructedPhi < 0.0) reconstructedPhi += 2.0 * CLHEP::pi;
        reconstructedValid = true;
        reconstructedChi2 = reconstructed.track.reducedChi2;
        reconstructionPermutations =
            static_cast<int>(reconstructed.permutationsTried);
      }
    }

    double zClose = 0.0;
    double sClose = 0.0;
    const auto u = incoming.unit();
    const auto v = outgoing.unit();
    const auto w0 = incomingPoint - outgoingPoint;
    const double a = u.dot(u);
    const double b = u.dot(v);
    const double c = v.dot(v);
    const double d = u.dot(w0);
    const double e = v.dot(w0);
    const double denominator = a * c - b * b;
    if (std::abs(denominator) > 1.0e-12) {
      const double alongIncoming = (b * e - c * d) / denominator;
      const double alongOutgoing = (a * e - b * d) / denominator;
      const auto closestIncoming = incomingPoint + alongIncoming * u;
      const auto closestOutgoing = outgoingPoint + alongOutgoing * v;
      zClose = 0.5 * (closestIncoming.z() + closestOutgoing.z());
      sClose = (closestIncoming - closestOutgoing).mag();
    }
    // Use the leading reconstructible charged track when multiple tracks
    // cross the same chamber pair.
    if (!best.valid || momentum > best.momentumGeV) {
      best.valid = true;
      best.conePass = conePass;
      best.thetaDeg = theta;
      best.phiRad = phi;
      best.momentumGeV = momentum;
      best.zCloseCm = zClose;
      best.sCloseCm = sClose;
      best.reconstructedValid = reconstructedValid;
      best.reconstructedThetaDeg = reconstructedTheta;
      best.reconstructedPhiRad = reconstructedPhi;
      best.reconstructedReducedChi2 = reconstructedChi2;
      best.reconstructionPermutations = reconstructionPermutations;
    }
  }
  return best;
}

void EventAction::EndOfEventAction(const G4Event*) {
  auto* out = G4AnalysisManager::Instance();
  out->FillNtupleIColumn(0, static_cast<int>(eventId_));
  out->FillNtupleDColumn(1, weight_);
  out->FillNtupleDColumn(2, momentumGeV_);
  out->FillNtupleIColumn(3, analyzer_);
  out->FillNtupleDColumn(4, interactionZCm_);
  out->FillNtupleIColumn(5, chargedSecondaries_);
  out->FillNtupleIColumn(6, primaryMask_);
  out->FillNtupleIColumn(7, secondaryMask_);
  out->FillNtupleIColumn(8, processCategory_);
  out->FillNtupleDColumn(9, thetaDeg_);
  out->FillNtupleDColumn(10, outgoingMomentumGeV_);

  int bestTrackMask = 0;
  for (const auto& [trackId, track] : tracks_) {
    (void)trackId;
    if (__builtin_popcount(static_cast<unsigned>(track.mask)) >
        __builtin_popcount(static_cast<unsigned>(bestTrackMask))) {
      bestTrackMask = track.mask;
    }
  }
  out->FillNtupleIColumn(11, bestTrackMask);
  out->FillNtupleIColumn(12, (bestTrackMask & 0x3) == 0x3 ? 1 : 0);
  out->FillNtupleIColumn(13, (bestTrackMask & 0xc) == 0xc ? 1 : 0);
  const auto front = bestCandidate(1);
  const auto rear = bestCandidate(2);
  out->FillNtupleDColumn(14, front.thetaDeg);
  out->FillNtupleDColumn(15, front.momentumGeV);
  out->FillNtupleIColumn(16, front.valid ? 1 : 0);
  out->FillNtupleDColumn(17, rear.thetaDeg);
  out->FillNtupleDColumn(18, rear.momentumGeV);
  out->FillNtupleIColumn(19, rear.valid ? 1 : 0);
  out->FillNtupleDColumn(20, front.zCloseCm);
  out->FillNtupleDColumn(21, front.sCloseCm);
  out->FillNtupleDColumn(22, rear.zCloseCm);
  out->FillNtupleDColumn(23, rear.sCloseCm);
  out->FillNtupleDColumn(24, front.phiRad);
  out->FillNtupleDColumn(25, rear.phiRad);
  const auto allFour = bestCandidate(3);
  out->FillNtupleIColumn(26, singleAnalyzer_ && allFour.valid ? 1 : 0);
  out->FillNtupleDColumn(27, allFour.thetaDeg);
  out->FillNtupleDColumn(28, allFour.phiRad);
  out->FillNtupleDColumn(29, allFour.momentumGeV);
  out->FillNtupleDColumn(30, allFour.zCloseCm);
  out->FillNtupleDColumn(31, allFour.sCloseCm);
  out->FillNtupleIColumn(32,
                         singleAnalyzer_ && allFour.valid && allFour.conePass
                             ? 1
                             : 0);
  out->FillNtupleIColumn(33, front.valid && front.conePass ? 1 : 0);
  out->FillNtupleIColumn(34, front.reconstructedValid ? 1 : 0);
  out->FillNtupleDColumn(35, front.reconstructedThetaDeg);
  out->FillNtupleDColumn(36, front.reconstructedPhiRad);
  out->FillNtupleDColumn(37, front.reconstructedReducedChi2);
  out->FillNtupleIColumn(38, front.reconstructionPermutations);
  out->FillNtupleIColumn(39, allFour.reconstructedValid ? 1 : 0);
  out->FillNtupleDColumn(40, allFour.reconstructedThetaDeg);
  out->FillNtupleDColumn(41, allFour.reconstructedPhiRad);
  out->FillNtupleDColumn(42, allFour.reconstructedReducedChi2);
  out->FillNtupleIColumn(43, allFour.reconstructionPermutations);
  out->AddNtupleRow();
}
