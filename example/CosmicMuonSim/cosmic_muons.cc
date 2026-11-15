// Usage: cosmic_muons <recipe> <nEvents> [outputName] [seed]
//   recipe: Reyna | Blanpied | Chris | Jokisch | Haruo
#include "Sim.hh"
#include "FTFP_BERT.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <recipe> <nEvents> [outputName] [seed]\n";
    return 1;
  }
  const std::string recipe = argv[1];
  const long nEvents = std::atol(argv[2]);
  const std::string out = argc > 3 ? argv[3] : "muons_" + recipe;
  const unsigned long seed = argc > 4 ? std::strtoul(argv[4], nullptr, 10) : 1759319UL;

  // The generator keeps global RNG state and writes cache files, so run sequentially.
  auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  runManager->SetUserInitialization(new DetectorConstruction);
  auto physics = new FTFP_BERT(0);
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(new ActionInitialization(recipe, out, seed));

  auto ui = G4UImanager::GetUIpointer();
  ui->ApplyCommand("/control/verbose 0");
  ui->ApplyCommand("/run/verbose 0");
  ui->ApplyCommand("/event/verbose 0");
  ui->ApplyCommand("/tracking/verbose 0");
  ui->ApplyCommand("/run/printProgress " + std::to_string(nEvents / 10 > 0 ? nEvents / 10 : 1));

  runManager->Initialize();
  runManager->BeamOn(nEvents);
  delete runManager;
  return 0;
}
