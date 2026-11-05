#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "FPPEventReader.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "SteppingAction.hh"

#include "G4PhysListFactory.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include <exception>
#include <iostream>
#include <memory>
#include <string>

int main(int argc, char** argv) {
  if (argc < 4 || argc > 5) {
    std::cerr << "Usage: " << argv[0]
              << " run.mac events.csv output.root [physics_list]\n";
    return 2;
  }

  try {
    const std::string physicsName = argc == 5 ? argv[4] : "FTFP_BERT";
    auto reader = std::make_shared<FPPEventReader>(argv[2]);
    // Keep the first validation implementation deterministic.  The shared
    // event reader and output schema can be made thread-local after the
    // single-thread physics path has been validated.
    auto* runManager =
        G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
    runManager->SetUserInitialization(new DetectorConstruction());

    G4PhysListFactory factory;
    auto* physics = factory.GetReferencePhysList(physicsName);
    if (!physics) throw std::runtime_error("Unknown Geant4 physics list: " + physicsName);
    runManager->SetUserInitialization(physics);

    auto* eventAction = new EventAction();
    runManager->SetUserAction(new PrimaryGeneratorAction(reader, eventAction));
    runManager->SetUserAction(new RunAction(argv[3]));
    runManager->SetUserAction(eventAction);
    runManager->SetUserAction(new SteppingAction(eventAction));
    runManager->Initialize();

    const std::string command = "/control/execute " + std::string(argv[1]);
    const int status = G4UImanager::GetUIpointer()->ApplyCommand(command);
    delete runManager;
    return status == 0 ? 0 : 1;
  } catch (const std::exception& error) {
    std::cerr << "Fatal: " << error.what() << '\n';
    return 1;
  }
}
