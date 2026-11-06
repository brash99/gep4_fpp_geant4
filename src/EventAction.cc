#include "EventAction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

void EventAction::BeginOfEventAction(const G4Event*) {
  eventId_ = -1;
  weight_ = 1.0;
  momentumGeV_ = 0.0;
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

EventAction::Candidate EventAction::bestCandidate(int analyzer) const {
  const int firstChamber = analyzer == 1 ? 0 : 2;
  const int requiredMask = analyzer == 1 ? 0x3 : 0xc;
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
      best = {true, theta, phi, momentum, zClose, sClose};
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
  out->AddNtupleRow();
}
