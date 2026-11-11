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
  out->CreateNtupleIColumn("single_all4_cone_pass");
  out->CreateNtupleIColumn("front_cone_pass");
  out->CreateNtupleIColumn("front_reco_valid");
  out->CreateNtupleDColumn("front_reco_theta_deg");
  out->CreateNtupleDColumn("front_reco_phi_az_rad");
  out->CreateNtupleDColumn("front_reco_reduced_chi2");
  out->CreateNtupleIColumn("front_reco_permutations");
  out->CreateNtupleIColumn("all4_reco_valid");
  out->CreateNtupleDColumn("all4_reco_theta_deg");
  out->CreateNtupleDColumn("all4_reco_phi_az_rad");
  out->CreateNtupleDColumn("all4_reco_reduced_chi2");
  out->CreateNtupleIColumn("all4_reco_permutations");
  out->CreateNtupleIColumn("hadronic_interaction_count");
  out->CreateNtupleIColumn("primary_hadronic_interaction_count");
  out->CreateNtupleIColumn("secondary_hadronic_interaction_count");
  out->CreateNtupleIColumn("primary_elastic_interaction_count");
  out->CreateNtupleIColumn("primary_inelastic_interaction_count");
  out->CreateNtupleIColumn("secondary_elastic_interaction_count");
  out->CreateNtupleIColumn("secondary_inelastic_interaction_count");
  out->CreateNtupleIColumn("primary_analyzer_exit_state");
  out->CreateNtupleDColumn("primary_analyzer_exit_p_GeV");
  out->CreateNtupleDColumn("primary_analyzer_exit_theta_deg");
  out->FinishNtuple();

  out->CreateNtuple("candidates", "All charged tracks crossing FPP1 and FPP2 before cuts");
  out->CreateNtupleIColumn("event_id");
  out->CreateNtupleIColumn("track_id");
  out->CreateNtupleIColumn("parent_id");
  out->CreateNtupleIColumn("chamber_mask");
  out->CreateNtupleDColumn("theta_deg");
  out->CreateNtupleDColumn("phi_az_rad");
  out->CreateNtupleDColumn("outgoing_p_GeV");
  out->CreateNtupleDColumn("zclose_cm");
  out->CreateNtupleDColumn("sclose_cm");
  out->CreateNtupleIColumn("front_cone_pass");
  out->CreateNtupleIColumn("all4_cone_pass");
  out->CreateNtupleIColumn("root_track_id");
  out->CreateNtupleIColumn("generation");
  out->CreateNtupleIColumn("lineage_process_category");
  out->CreateNtupleIColumn("track_elastic_interaction_count");
  out->CreateNtupleIColumn("track_inelastic_interaction_count");
  out->CreateNtupleIColumn("track_other_hadronic_interaction_count");
  out->FinishNtuple();

  out->CreateNtuple("interactions", "Every hadronic interaction in an analyzer");
  out->CreateNtupleIColumn("event_id");
  out->CreateNtupleIColumn("interaction_index");
  out->CreateNtupleIColumn("analyzer");
  out->CreateNtupleIColumn("track_id");
  out->CreateNtupleIColumn("parent_id");
  out->CreateNtupleIColumn("root_track_id");
  out->CreateNtupleIColumn("generation");
  out->CreateNtupleIColumn("lineage_process_category");
  out->CreateNtupleIColumn("process_category");
  out->CreateNtupleDColumn("z_cm");
  out->CreateNtupleDColumn("incoming_p_GeV");
  out->CreateNtupleDColumn("theta_deg");
  out->CreateNtupleDColumn("outgoing_p_GeV");
  out->CreateNtupleIColumn("charged_secondaries");
  out->CreateNtupleIColumn("particle_pdg");
  out->CreateNtupleDColumn("particle_charge");
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
