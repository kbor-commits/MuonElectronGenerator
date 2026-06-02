#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

// my headers
#include "MuonElectronGenerator.hh"
#include "DetectorConstruction.hh"
#include "container.h"
// GEANT4 headers
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction{
public:
  PrimaryGeneratorAction(int MyPE,container *);
  ~PrimaryGeneratorAction();

  void GeneratePrimaries(G4Event* anEvent);
  void SetDetectorConstruction(DetectorConstruction* detector){ DetCon = detector; } 
    // Function to select a random point within a specified physical volume
  G4ThreeVector SelectRandomPointInVolume(G4VPhysicalVolume* physVolume);
  
private:
  G4ParticleGun *particleGun;
  MuonElectronGenerator *MEGen;
  // Access the DetectorConstruction object
  DetectorConstruction *DetCon;

  G4ParticleTable *particleTable;
  long PCount;
  container *_dramiel; 
  
  G4double Size_y, Size_z; 
  G4double Size_y1, Size_z1; 
  G4double Size_y2, Size_z2; 
  G4double D1_y;
  G4double D1_z; 
  G4double D1_zup; 
  G4double D2_y;
  G4double D2_z, D2_zup;
  G4double P1_x, P2_x;
  G4double deltax;

  G4double m_par;
  G4float mass; 

};

#endif
