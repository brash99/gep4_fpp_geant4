#include "EventAction.hh"

#include "G4AnalysisManager.hh"

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
  trackMasks_.clear();
}

void EventAction::setInput(long id, double weight, double momentumGeV) {
  eventId_ = id;
  weight_ = weight;
  momentumGeV_ = momentumGeV;
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

void EventAction::recordChamberCrossing(int chamber, int trackId, bool primary) {
  const int bit = 1 << (chamber - 1);
  trackMasks_[trackId] |= bit;
  if (primary) primaryMask_ |= bit;
  else secondaryMask_ |= bit;
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
  for (const auto& [trackId, mask] : trackMasks_) {
    (void)trackId;
    if (__builtin_popcount(static_cast<unsigned>(mask)) >
        __builtin_popcount(static_cast<unsigned>(bestTrackMask))) {
      bestTrackMask = mask;
    }
  }
  out->FillNtupleIColumn(11, bestTrackMask);
  out->FillNtupleIColumn(12, (bestTrackMask & 0x3) == 0x3 ? 1 : 0);
  out->FillNtupleIColumn(13, (bestTrackMask & 0xc) == 0xc ? 1 : 0);
  out->AddNtupleRow();
}
