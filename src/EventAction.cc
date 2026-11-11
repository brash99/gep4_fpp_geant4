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
  trackMetadata_.clear();
  interactionCount_ = 0;
  primaryInteractionCount_ = 0;
  secondaryInteractionCount_ = 0;
  primaryElasticCount_ = 0;
  primaryInelasticCount_ = 0;
  secondaryElasticCount_ = 0;
  secondaryInelasticCount_ = 0;
  primaryExitState_ = 0;
  primaryExitMomentumGeV_ = 0.0;
  primaryExitThetaDeg_ = 0.0;
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

void EventAction::registerTrack(int trackId, int parentId,
                                int creatorProcessCategory) {
  if (trackMetadata_.find(trackId) != trackMetadata_.end()) return;
  TrackMetadata metadata;
  metadata.parentId = parentId;
  metadata.rootTrackId = trackId;
  if (parentId != 0) {
    const auto parent = trackMetadata_.find(parentId);
    if (parent != trackMetadata_.end()) {
      metadata.rootTrackId = parent->second.rootTrackId;
      metadata.generation = parent->second.generation + 1;
      metadata.lineageProcessCategory =
          parent->second.lineageProcessCategory != 0
              ? parent->second.lineageProcessCategory
              : creatorProcessCategory;
    } else {
      metadata.rootTrackId = parentId;
      metadata.generation = 1;
      metadata.lineageProcessCategory = creatorProcessCategory;
    }
  }
  trackMetadata_.emplace(trackId, metadata);
}

void EventAction::recordHadronicInteraction(
    int analyzer, int trackId, int parentId, int processCategory, double zCm,
    double incomingMomentumGeV, double thetaDeg, double outgoingMomentumGeV,
    int chargedSecondaries, int particlePdg, double particleCharge) {
  registerTrack(trackId, parentId, 0);
  auto& metadata = trackMetadata_.at(trackId);
  ++interactionCount_;
  if (metadata.generation == 0) {
    ++primaryInteractionCount_;
    if (processCategory == 1) ++primaryElasticCount_;
    if (processCategory == 2) ++primaryInelasticCount_;
  } else {
    ++secondaryInteractionCount_;
    if (processCategory == 1) ++secondaryElasticCount_;
    if (processCategory == 2) ++secondaryInelasticCount_;
  }
  if (processCategory == 1) ++metadata.elasticInteractions;
  else if (processCategory == 2) ++metadata.inelasticInteractions;
  else ++metadata.otherHadronicInteractions;

  // Preserve the original first-interaction event summary for compatibility.
  if (analyzer_ == 0) {
    analyzer_ = analyzer;
    processCategory_ = processCategory;
    interactionZCm_ = zCm;
    thetaDeg_ = thetaDeg;
    outgoingMomentumGeV_ = outgoingMomentumGeV;
    chargedSecondaries_ = chargedSecondaries;
  }

  auto* out = G4AnalysisManager::Instance();
  out->FillNtupleIColumn(2, 0, static_cast<int>(eventId_));
  out->FillNtupleIColumn(2, 1, interactionCount_ - 1);
  out->FillNtupleIColumn(2, 2, analyzer);
  out->FillNtupleIColumn(2, 3, trackId);
  out->FillNtupleIColumn(2, 4, parentId);
  out->FillNtupleIColumn(2, 5, metadata.rootTrackId);
  out->FillNtupleIColumn(2, 6, metadata.generation);
  out->FillNtupleIColumn(2, 7, metadata.lineageProcessCategory);
  out->FillNtupleIColumn(2, 8, processCategory);
  out->FillNtupleDColumn(2, 9, zCm);
  out->FillNtupleDColumn(2, 10, incomingMomentumGeV);
  out->FillNtupleDColumn(2, 11, thetaDeg);
  out->FillNtupleDColumn(2, 12, outgoingMomentumGeV);
  out->FillNtupleIColumn(2, 13, chargedSecondaries);
  out->FillNtupleIColumn(2, 14, particlePdg);
  out->FillNtupleDColumn(2, 15, particleCharge);
  out->AddNtupleRow(2);
}

void EventAction::recordPrimaryAnalyzerExit(const G4ThreeVector& positionCm,
                                            const G4ThreeVector& direction,
                                            double momentumGeV) {
  constexpr double upstreamFaceCm = 111.4215;
  constexpr double faceToleranceCm = 1.0e-3;
  const double downstreamFaceCm = upstreamFaceCm + analyzer1ThicknessCm_;
  if (positionCm.z() >= downstreamFaceCm - faceToleranceCm) {
    primaryExitState_ = 1;
  } else if (positionCm.z() <= upstreamFaceCm + faceToleranceCm) {
    primaryExitState_ = 2;
  } else {
    primaryExitState_ = 3;
  }
  primaryExitMomentumGeV_ = momentumGeV;
  primaryExitThetaDeg_ =
      std::acos(std::clamp(incidentDirection_.unit().dot(direction.unit()),
                           -1.0, 1.0)) /
      deg;
}

void EventAction::recordPrimaryKilledInAnalyzer(double momentumGeV) {
  primaryExitState_ = 4;
  primaryExitMomentumGeV_ = momentumGeV;
  primaryExitThetaDeg_ = 0.0;
}

void EventAction::recordChamberCrossing(int chamber, int trackId, int parentId,
                                        bool primary,
                                        const G4ThreeVector& direction,
                                        const G4ThreeVector& positionCm,
                                        double momentumGeV) {
  registerTrack(trackId, parentId, 0);
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

EventAction::Candidate EventAction::candidateForTrack(
    const TrackObservation& track, int analyzer, bool reconstruct) const {
  const int firstChamber = analyzer == 2 ? 2 : 0;
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
  const double rearPlaneZ = analyzer == 1 ? layout_.chamberCenterZCm[1]
                                          : layout_.chamberCenterZCm[3];
  const bool conePass =
      legacyConeTest(incomingPoint, incoming, theta, rearPlaneZ);

  bool reconstructedValid = false;
  double reconstructedTheta = 0.0;
  double reconstructedPhi = 0.0;
  double reconstructedChi2 = 0.0;
  int reconstructionPermutations = 0;
  if (reconstruct && (analyzer == 1 || analyzer == 3)) {
    const FPPVector3 digitizerPoint{outgoingPoint.x(), outgoingPoint.y(),
                                    outgoingPoint.z()};
    const FPPVector3 digitizerDirection{outgoing.x(), outgoing.y(),
                                        outgoing.z()};
    const int chamberMask = analyzer == 3 ? 0xf : 0x3;
    const auto hits = digitizeFPPTrack(digitizerPoint, digitizerDirection,
                                       planes_, chamberMask);
    const auto reconstructed = fitFPPLeftRightPermutations(hits, planes_);
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
  Candidate candidate;
  candidate.valid = true;
  candidate.conePass = conePass;
  candidate.thetaDeg = theta;
  candidate.phiRad = phi;
  candidate.momentumGeV = momentum;
  candidate.zCloseCm = zClose;
  candidate.sCloseCm = sClose;
  candidate.reconstructedValid = reconstructedValid;
  candidate.reconstructedThetaDeg = reconstructedTheta;
  candidate.reconstructedPhiRad = reconstructedPhi;
  candidate.reconstructedReducedChi2 = reconstructedChi2;
  candidate.reconstructionPermutations = reconstructionPermutations;
  return candidate;
}

EventAction::Candidate EventAction::bestCandidate(int analyzer) const {
  const int requiredMask = analyzer == 1 ? 0x3 : (analyzer == 2 ? 0xc : 0xf);
  Candidate best;
  for (const auto& [trackId, track] : tracks_) {
    (void)trackId;
    if ((track.mask & requiredMask) != requiredMask) continue;
    const auto candidate = candidateForTrack(track, analyzer, true);
    // Preserve the original leading-momentum policy for baseline regression.
    if (!best.valid || candidate.momentumGeV > best.momentumGeV) best = candidate;
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
  out->FillNtupleIColumn(44, interactionCount_);
  out->FillNtupleIColumn(45, primaryInteractionCount_);
  out->FillNtupleIColumn(46, secondaryInteractionCount_);
  out->FillNtupleIColumn(47, primaryElasticCount_);
  out->FillNtupleIColumn(48, primaryInelasticCount_);
  out->FillNtupleIColumn(49, secondaryElasticCount_);
  out->FillNtupleIColumn(50, secondaryInelasticCount_);
  out->FillNtupleIColumn(51, primaryExitState_);
  out->FillNtupleDColumn(52, primaryExitMomentumGeV_);
  out->FillNtupleDColumn(53, primaryExitThetaDeg_);
  out->AddNtupleRow();

  // Keep every geometrically reconstructible front-pair track, before cuts.
  // This enables exclusive multiplicity and momentum-blind choices offline.
  for (const auto& [trackId, track] : tracks_) {
    if ((track.mask & 0x3) != 0x3) continue;
    const auto candidate = candidateForTrack(track, 1, false);
    out->FillNtupleIColumn(1, 0, static_cast<int>(eventId_));
    out->FillNtupleIColumn(1, 1, trackId);
    out->FillNtupleIColumn(1, 2, track.parentId);
    out->FillNtupleIColumn(1, 3, track.mask);
    out->FillNtupleDColumn(1, 4, candidate.thetaDeg);
    out->FillNtupleDColumn(1, 5, candidate.phiRad);
    out->FillNtupleDColumn(1, 6, candidate.momentumGeV);
    out->FillNtupleDColumn(1, 7, candidate.zCloseCm);
    out->FillNtupleDColumn(1, 8, candidate.sCloseCm);
    out->FillNtupleIColumn(1, 9, candidate.conePass ? 1 : 0);
    const bool allFour = (track.mask & 0xf) == 0xf;
    out->FillNtupleIColumn(1, 10,
        allFour && legacyConeTest(incidentPositionCm_, incidentDirection_,
                                  candidate.thetaDeg,
                                  layout_.chamberCenterZCm[3]) ? 1 : 0);
    const auto metadata = trackMetadata_.find(trackId);
    const TrackMetadata emptyMetadata;
    const auto& trackInfo = metadata != trackMetadata_.end()
                                ? metadata->second
                                : emptyMetadata;
    out->FillNtupleIColumn(1, 11, trackInfo.rootTrackId);
    out->FillNtupleIColumn(1, 12, trackInfo.generation);
    out->FillNtupleIColumn(1, 13, trackInfo.lineageProcessCategory);
    out->FillNtupleIColumn(1, 14, trackInfo.elasticInteractions);
    out->FillNtupleIColumn(1, 15, trackInfo.inelasticInteractions);
    out->FillNtupleIColumn(1, 16, trackInfo.otherHadronicInteractions);
    out->AddNtupleRow(1);
  }
}
