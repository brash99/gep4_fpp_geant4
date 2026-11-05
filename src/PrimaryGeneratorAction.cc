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
  gun_->SetParticlePosition({event.xCm * cm, event.yCm * cm, event.zCm * cm});
  gun_->SetParticleMomentumDirection({event.pxGeV, event.pyGeV, event.pzGeV});
  gun_->SetParticleMomentum(p * GeV);
  eventAction_->setInput(event.eventId, event.weight, p);
  gun_->GeneratePrimaryVertex(geantEvent);
}
