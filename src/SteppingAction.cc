#include "SteppingAction.hh"

#include "EventAction.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

#include <algorithm>
#include <cmath>
#include <string>

namespace {
int processCategory(const G4VProcess* process) {
  if (!process) return 0;
  const auto name = process->GetProcessName();
  if (name == "hadElastic") return 1;
  if (name.find("Inelastic") != std::string::npos) return 2;
  return process->GetProcessType() == fHadronic ? 3 : 0;
}
}  // namespace

void SteppingAction::UserSteppingAction(const G4Step* step) {
  const auto* track = step->GetTrack();
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto* preVolume = pre->GetPhysicalVolume();
  const auto* postVolume = post->GetPhysicalVolume();
  if (!preVolume) return;

  eventAction_->registerTrack(
      track->GetTrackID(), track->GetParentID(),
      processCategory(track->GetCreatorProcess()));

  const std::string preName = preVolume->GetName();
  const auto* process = post->GetProcessDefinedStep();
  if ((preName == "Analyzer1" || preName == "Analyzer2") && process &&
      process->GetProcessType() == fHadronic) {
    int charged = 0;
    for (const auto* secondary : *step->GetSecondaryInCurrentStep()) {
      if (secondary->GetDefinition()->GetPDGCharge() != 0.) ++charged;
    }
    const auto incoming = pre->GetMomentumDirection();
    const auto outgoing = post->GetMomentumDirection();
    const double cosine = std::max(-1.0, std::min(1.0, incoming.dot(outgoing)));
    const double thetaDeg = std::acos(cosine) / deg;
    eventAction_->recordHadronicInteraction(
        preName == "Analyzer1" ? 1 : 2, track->GetTrackID(),
        track->GetParentID(), processCategory(process),
        post->GetPosition().z() / cm, pre->GetMomentum().mag() / GeV,
        thetaDeg, post->GetMomentum().mag() / GeV,
        charged, track->GetDefinition()->GetPDGEncoding(),
        track->GetDefinition()->GetPDGCharge());
  }

  if (preName == "Analyzer1" && track->GetTrackID() == 1) {
    if (post->GetStepStatus() == fGeomBoundary &&
        (!postVolume || postVolume->GetName() != "Analyzer1")) {
      eventAction_->recordPrimaryAnalyzerExit(
          post->GetPosition() / cm, post->GetMomentumDirection(),
          post->GetMomentum().mag() / GeV);
    } else if (track->GetTrackStatus() == fStopAndKill ||
               track->GetTrackStatus() == fKillTrackAndSecondaries) {
      eventAction_->recordPrimaryKilledInAnalyzer(
          post->GetMomentum().mag() / GeV);
    }
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
