#include "RunAction.hh"

#include "G4AnalysisManager.hh"

RunAction::RunAction(std::string outputPath) : outputPath_(std::move(outputPath)) {
  auto* out = G4AnalysisManager::Instance();
  out->SetVerboseLevel(1);
  out->CreateNtuple("fpp", "GEp-4 FPP transport events");
  out->CreateNtupleIColumn("event_id");
  out->CreateNtupleDColumn("weight");
  out->CreateNtupleDColumn("incident_p_GeV");
  out->CreateNtupleIColumn("first_analyzer");
  out->CreateNtupleDColumn("interaction_z_cm");
  out->CreateNtupleIColumn("charged_secondaries");
  out->CreateNtupleIColumn("primary_chamber_mask");
  out->CreateNtupleIColumn("secondary_chamber_mask");
  out->FinishNtuple();
}

void RunAction::BeginOfRunAction(const G4Run*) {
  G4AnalysisManager::Instance()->OpenFile(outputPath_);
}

void RunAction::EndOfRunAction(const G4Run*) {
  auto* out = G4AnalysisManager::Instance();
  out->Write();
  out->CloseFile();
}
