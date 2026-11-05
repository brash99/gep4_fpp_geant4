#pragma once

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <unordered_map>

class G4Event;

class EventAction final : public G4UserEventAction {
 public:
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event*) override;

  void setInput(long id, double weight, double momentumGeV);
  void recordHadronicInteraction(int analyzer, int processCategory, double zCm,
                                 double thetaDeg, double outgoingMomentumGeV,
                                 int chargedSecondaries);
  void recordChamberCrossing(int chamber, int trackId, bool primary);

 private:
  long eventId_ = -1;
  double weight_ = 1.0;
  double momentumGeV_ = 0.0;
  int analyzer_ = 0;
  int processCategory_ = 0;
  double interactionZCm_ = 0.0;
  double thetaDeg_ = 0.0;
  double outgoingMomentumGeV_ = 0.0;
  int chargedSecondaries_ = 0;
  int primaryMask_ = 0;
  int secondaryMask_ = 0;
  std::unordered_map<int, int> trackMasks_;
};
