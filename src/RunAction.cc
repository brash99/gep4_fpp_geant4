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
  out->CreateNtupleIColumn("process_category");
  out->CreateNtupleDColumn("theta_deg");
  out->CreateNtupleDColumn("outgoing_p_GeV");
  out->CreateNtupleIColumn("best_track_chamber_mask");
  out->CreateNtupleIColumn("same_track_front_pair");
  out->CreateNtupleIColumn("same_track_rear_pair");
  out->CreateNtupleDColumn("front_track_theta_deg");
  out->CreateNtupleDColumn("front_track_p_GeV");
  out->CreateNtupleIColumn("front_track_valid");
  out->CreateNtupleDColumn("rear_track_theta_deg");
  out->CreateNtupleDColumn("rear_track_p_GeV");
  out->CreateNtupleIColumn("rear_track_valid");
  out->CreateNtupleDColumn("front_zclose_cm");
  out->CreateNtupleDColumn("front_sclose_cm");
  out->CreateNtupleDColumn("rear_zclose_cm");
  out->CreateNtupleDColumn("rear_sclose_cm");
  out->CreateNtupleDColumn("front_phi_az_rad");
  out->CreateNtupleDColumn("rear_phi_az_rad");
  out->CreateNtupleIColumn("single_all4_valid");
  out->CreateNtupleDColumn("single_all4_theta_deg");
  out->CreateNtupleDColumn("single_all4_phi_az_rad");
  out->CreateNtupleDColumn("single_all4_p_GeV");
  out->CreateNtupleDColumn("single_all4_zclose_cm");
  out->CreateNtupleDColumn("single_all4_sclose_cm");
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
