#include "PrimaryGeneratorAction.hh"

#include "EventAction.hh"
#include "FPPEventReader.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <stdexcept>

PrimaryGeneratorAction::PrimaryGeneratorAction(
    std::shared_ptr<const FPPEventReader> reader, EventAction* eventAction)
    : reader_(std::move(reader)), eventAction_(eventAction), gun_(new G4ParticleGun(1)) {
  gun_->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("proton"));
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete gun_; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* geantEvent) {
  const auto index = static_cast<std::size_t>(geantEvent->GetEventID());
  if (index >= reader_->size()) throw std::runtime_error("Requested more events than input rows");
  const auto& event = reader_->event(index);
  const double p = std::sqrt(event.pxGeV * event.pxGeV + event.pyGeV * event.pyGeV +
                             event.pzGeV * event.pzGeV);
  const G4ThreeVector positionCm(event.xCm, event.yCm, event.zCm);
  gun_->SetParticlePosition(positionCm * cm);
  const G4ThreeVector direction(event.pxGeV, event.pyGeV, event.pzGeV);
  gun_->SetParticleMomentumDirection(direction);
  const double massGeV = gun_->GetParticleDefinition()->GetPDGMass() / GeV;
  const double kineticEnergyGeV = std::sqrt(p * p + massGeV * massGeV) - massGeV;
  gun_->SetParticleEnergy(kineticEnergyGeV * GeV);
  eventAction_->setInput(event.eventId, event.weight, p, positionCm,
                         direction.unit());
  gun_->GeneratePrimaryVertex(geantEvent);
}
