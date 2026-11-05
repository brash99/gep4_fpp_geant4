#include "SteppingAction.hh"

#include "EventAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

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
    eventAction_->recordHadronicInteraction(preName == "Analyzer1" ? 1 : 2,
                                            post->GetPosition().z() / cm, charged);
  }

  if (!postVolume || post->GetStepStatus() != fGeomBoundary ||
      track->GetDefinition()->GetPDGCharge() == 0.) return;
  const std::string postName = postVolume->GetName();
  for (int chamber = 1; chamber <= 4; ++chamber) {
    if (postName == "FPP" + std::to_string(chamber)) {
      eventAction_->recordChamberCrossing(chamber, track->GetTrackID(),
                                          track->GetParentID() == 0);
    }
  }
}
