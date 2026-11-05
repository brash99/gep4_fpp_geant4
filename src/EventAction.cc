#include "EventAction.hh"

#include "G4AnalysisManager.hh"

void EventAction::BeginOfEventAction(const G4Event*) {
  eventId_ = -1;
  weight_ = 1.0;
  momentumGeV_ = 0.0;
  analyzer_ = 0;
  interactionZCm_ = 0.0;
  chargedSecondaries_ = 0;
  primaryMask_ = 0;
  secondaryMask_ = 0;
}

void EventAction::setInput(long id, double weight, double momentumGeV) {
  eventId_ = id;
  weight_ = weight;
  momentumGeV_ = momentumGeV;
}

void EventAction::recordHadronicInteraction(int analyzer, double zCm,
                                            int chargedSecondaries) {
  if (analyzer_ != 0) return;
  analyzer_ = analyzer;
  interactionZCm_ = zCm;
  chargedSecondaries_ = chargedSecondaries;
}

void EventAction::recordChamberCrossing(int chamber, int, bool primary) {
  const int bit = 1 << (chamber - 1);
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
  out->AddNtupleRow();
}
