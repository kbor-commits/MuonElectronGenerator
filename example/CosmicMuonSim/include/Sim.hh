// Minimal Geant4 application wrapping MuonElectronGenerator.
// Geometry: a large slab of air with a horizontal scoring plane at z = 0.
// Muons are launched 1 m above the plane, travelling downward.
#ifndef Sim_hh
#define Sim_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4ParticleGun.hh"
#include "MT_Random.hh"
#include <string>

class MuonElectronGenerator;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  G4VPhysicalVolume* Construct() override;
};

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGenerator(const std::string& recipe, unsigned long seed);
  ~PrimaryGenerator() override;
  void GeneratePrimaries(G4Event*) override;
private:
  std::string fRecipe;
  MT_Random fRng;
  MuonElectronGenerator* fGen = nullptr;
  G4ParticleGun* fGun = nullptr;
};

class RunAction : public G4UserRunAction {
public:
  explicit RunAction(const std::string& outName) : fOut(outName) {}
  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;
private:
  std::string fOut;
};

class SteppingAction : public G4UserSteppingAction {
public:
  void UserSteppingAction(const G4Step*) override;
};

class ActionInitialization : public G4VUserActionInitialization {
public:
  ActionInitialization(const std::string& recipe, const std::string& out, unsigned long seed)
    : fRecipe(recipe), fOut(out), fSeed(seed) {}
  void Build() const override;
private:
  std::string fRecipe, fOut;
  unsigned long fSeed;
};

#endif
