#include "Sim.hh"
#include "MuonElectronGenerator.hh"

#include "G4AnalysisManager.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

extern int LIST_OUT;  // defined in MuonElectronGenerator.cc

namespace { const G4double kLaunchHeight = 1.0 * m; }

// ---------------------------------------------------------------- geometry
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  auto air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  // Wide enough that muons up to ~89.97 deg from 1 m above still cross the plane
  const G4double halfXY = 2.5 * km, halfZ = 5 * m, halfPlane = 0.5 * mm;

  auto worldS = new G4Box("World", halfXY, halfXY, halfZ);
  auto worldL = new G4LogicalVolume(worldS, air, "World");
  auto worldP = new G4PVPlacement(nullptr, {}, worldL, "World", nullptr, false, 0);

  auto planeS = new G4Box("Plane", halfXY, halfXY, halfPlane);
  auto planeL = new G4LogicalVolume(planeS, air, "Plane");
  new G4PVPlacement(nullptr, {}, planeL, "Plane", worldL, false, 0);
  return worldP;
}

// ---------------------------------------------------------------- primaries
PrimaryGenerator::PrimaryGenerator(const std::string& recipe, unsigned long seed)
  : fRecipe(recipe)
{
  LIST_OUT = 0;  // per-event text dump not needed; we record via G4AnalysisManager
  fRng.init_random(seed, 0);
  fGen = new MuonElectronGenerator(&fRng, fRecipe);
  fGun = new G4ParticleGun(1);
}

PrimaryGenerator::~PrimaryGenerator() { delete fGen; delete fGun; }

void PrimaryGenerator::GeneratePrimaries(G4Event* event)
{
  fGen->GenerateMuon(fRecipe);  // samples type, zenith angle [deg], momentum [GeV/c]

  G4ParticleDefinition* particle = fGen->GetType();
  const G4double p = fGen->GetEnergy() * GeV;
  const G4double mass = particle->GetPDGMass();
  const G4double ekin = std::sqrt(p * p + mass * mass) - mass;
  const G4double theta = fGen->GetZenithAngle() * deg;
  const G4double phi = twopi * fRng.genrand_real();

  // Downward-going direction with zenith angle theta
  G4ThreeVector dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), -std::cos(theta));
  fGun->SetParticleDefinition(particle);
  fGun->SetParticleMomentumDirection(dir);
  fGun->SetParticleEnergy(ekin);
  fGun->SetParticlePosition({0, 0, kLaunchHeight});
  fGun->GeneratePrimaryVertex(event);

  auto am = G4AnalysisManager::Instance();
  am->FillNtupleIColumn(0, 0, particle->GetPDGEncoding());
  am->FillNtupleDColumn(0, 1, p / GeV);
  am->FillNtupleDColumn(0, 2, ekin / GeV);
  am->FillNtupleDColumn(0, 3, theta / deg);
  am->FillNtupleDColumn(0, 4, phi / deg);
  am->AddNtupleRow(0);
}

// ---------------------------------------------------------------- scoring
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  auto pre = step->GetPreStepPoint();
  if (pre->GetStepStatus() != fGeomBoundary) return;
  if (!pre->GetPhysicalVolume() || pre->GetPhysicalVolume()->GetName() != "Plane") return;
  auto track = step->GetTrack();
  if (track->GetParentID() != 0) return;  // primary muon only

  const G4ThreeVector d = pre->GetMomentumDirection();
  if (d.z() >= 0) return;
  auto am = G4AnalysisManager::Instance();
  am->FillNtupleIColumn(1, 0, track->GetDefinition()->GetPDGEncoding());
  am->FillNtupleDColumn(1, 1, pre->GetMomentum().mag() / GeV);
  am->FillNtupleDColumn(1, 2, pre->GetKineticEnergy() / GeV);
  am->FillNtupleDColumn(1, 3, std::acos(-d.z()) / deg);
  am->FillNtupleDColumn(1, 4, pre->GetPosition().x() / m);
  am->FillNtupleDColumn(1, 5, pre->GetPosition().y() / m);
  am->AddNtupleRow(1);
}

// ---------------------------------------------------------------- run / output
void RunAction::BeginOfRunAction(const G4Run*)
{
  auto am = G4AnalysisManager::Instance();
  am->SetDefaultFileType("csv");
  am->SetVerboseLevel(0);
  am->SetNtupleMerging(false);
  if (am->GetNofNtuples() == 0) {
    am->CreateNtuple("generated", "Primary muons from MuonElectronGenerator");
    am->CreateNtupleIColumn("pdg");
    am->CreateNtupleDColumn("p_GeV");
    am->CreateNtupleDColumn("ekin_GeV");
    am->CreateNtupleDColumn("theta_deg");
    am->CreateNtupleDColumn("phi_deg");
    am->FinishNtuple();
    am->CreateNtuple("plane", "Primary muons crossing the scoring plane");
    am->CreateNtupleIColumn("pdg");
    am->CreateNtupleDColumn("p_GeV");
    am->CreateNtupleDColumn("ekin_GeV");
    am->CreateNtupleDColumn("theta_deg");
    am->CreateNtupleDColumn("x_m");
    am->CreateNtupleDColumn("y_m");
    am->FinishNtuple();
  }
  am->OpenFile(fOut);
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto am = G4AnalysisManager::Instance();
  am->Write();
  am->CloseFile();
}

void ActionInitialization::Build() const
{
  SetUserAction(new PrimaryGenerator(fRecipe, fSeed));
  SetUserAction(new RunAction(fOut));
  SetUserAction(new SteppingAction);
}
