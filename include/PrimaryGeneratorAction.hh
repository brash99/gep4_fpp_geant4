#pragma once

#include "G4VUserPrimaryGeneratorAction.hh"

#include <memory>

class EventAction;
class FPPEventReader;
class G4Event;
class G4ParticleGun;

class PrimaryGeneratorAction final : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGeneratorAction(std::shared_ptr<const FPPEventReader> reader,
                         EventAction* eventAction);
  ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event*) override;

 private:
  std::shared_ptr<const FPPEventReader> reader_;
  EventAction* eventAction_;
  G4ParticleGun* gun_;
};
