#pragma once

#include "G4UserRunAction.hh"

#include <string>

class G4Run;

class RunAction final : public G4UserRunAction {
 public:
  explicit RunAction(std::string outputPath);
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

 private:
  std::string outputPath_;
};
