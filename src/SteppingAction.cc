#include "SteppingAction.hh"

#include "EventAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <cmath>
#include <string>

void SteppingAction::UserSteppingAction(const G4Step* step) {
  const auto* track = step->GetTrack();
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto* preVolume = pre->GetPhysicalVolume();
  const auto* postVolume = post->GetPhysicalVolume();
  if (!preVolume) return;

  const std::string preName = preVolume->GetName();
  const auto* process = post->GetProcessDefinedStep();
  if ((preName == "Analyzer1" || preName == "Analyzer2") && process &&
      process->GetProcessType() == fHadronic) {
    int charged = 0;
    for (const auto* secondary : *step->GetSecondaryInCurrentStep()) {
      if (secondary->GetDefinition()->GetPDGCharge() != 0.) ++charged;
    }
    const auto processName = process->GetProcessName();
    int processCategory = 3;
    if (processName == "hadElastic") processCategory = 1;
    else if (processName.find("Inelastic") != std::string::npos) processCategory = 2;

    const auto incoming = pre->GetMomentumDirection();
    const auto outgoing = post->GetMomentumDirection();
    const double cosine = std::max(-1.0, std::min(1.0, incoming.dot(outgoing)));
    const double thetaDeg = std::acos(cosine) / deg;
    eventAction_->recordHadronicInteraction(
        preName == "Analyzer1" ? 1 : 2, processCategory,
        post->GetPosition().z() / cm, thetaDeg, post->GetMomentum().mag() / GeV,
        charged);
  }

  if (!postVolume || post->GetStepStatus() != fGeomBoundary ||
      track->GetDefinition()->GetPDGCharge() == 0.) return;
  const std::string postName = postVolume->GetName();
  for (int chamber = 1; chamber <= 4; ++chamber) {
    if (postName == "FPP" + std::to_string(chamber)) {
      eventAction_->recordChamberCrossing(
          chamber, track->GetTrackID(), track->GetParentID(),
          track->GetParentID() == 0, post->GetMomentumDirection(),
          post->GetPosition() / cm,
          post->GetMomentum().mag() / GeV);
    }
  }
}
