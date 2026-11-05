#pragma once

#include "G4UserSteppingAction.hh"

class EventAction;
class G4Step;

class SteppingAction final : public G4UserSteppingAction {
 public:
  explicit SteppingAction(EventAction* eventAction) : eventAction_(eventAction) {}
  void UserSteppingAction(const G4Step*) override;

 private:
  EventAction* eventAction_;
};
