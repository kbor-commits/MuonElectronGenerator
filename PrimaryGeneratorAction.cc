// my headers
#include "PrimaryGeneratorAction.hh"
#include "MuonElectronGenerator.hh"
//#include "DetectorConstruction.hh"
// GEANT4 headers
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
// C++ headers
#include "globals.hh"
#include "MT_Random.hh"
#include "math.h"
#include <assert.h>
// ROOT headers
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TLeaf.h"

#define PI 3.1415926535

MT_Random RNGenerator;

PrimaryGeneratorAction::PrimaryGeneratorAction(int MyPE,container *dramiel){
  _dramiel=dramiel;

  // Initialize random numbers generator
  RNGenerator.init_random(1759319, MyPE);
  
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);
  particleTable = G4ParticleTable::GetParticleTable();

  //G4float mass = (G4float)particle->GetPDGMass();
  MEGen = new MuonElectronGenerator(&RNGenerator);

  // Access the DetectorConstruction object
  //  DetCon = static_cast<DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());

  return;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction(){
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = particleTable->FindParticle("mu-");
  G4float i = (float)RNGenerator.genrand_real();
  if ( i < 0.45 ){
    particle = particleTable->FindParticle("mu-");
  } else {
    particle = particleTable->FindParticle("mu+");
  }
  particleGun->SetParticleDefinition(particle);

  G4double momentum;

  G4double R = 1680.2*mm;
  double theta = MEGen->ChooseTheReyna();
  double phi1 = 2.*PI*RNGenerator.genrand_real();
  double phi2 = 2.*PI*RNGenerator.genrand_real();
  G4double Z1 = 1000.+(0.5-RNGenerator.genrand_real())*1500.*mm;
  
  int f = -1;
  while( f < 0 ){
    if( tan(theta) <= R*sqrt(2.*(1-cos(phi2-phi1)))/(Z1+250.*mm) && tan(theta) >= R*sqrt(2.*(1-cos(phi2-phi1)))/(Z1+1750.*mm) ){
      f = 1;
      G4cout << theta << ", R: " << R*sqrt(2.*(1-cos(phi2-phi1))) << ", Z: " << Z1 << G4endl;
    } else {
      theta = MEGen->ChooseTheReyna();      
    }
  }

  double xp, yp, zp;
  xp = (cos(phi2)-cos(phi1))/sqrt(2.*(1-cos(phi2-phi1)))*sin(theta);
  yp = (sin(phi2)-sin(phi1))/sqrt(2.*(1-cos(phi2-phi1)))*sin(theta);
  zp = -cos(theta);
  G4double gamma = 50/(cos(theta))*mm;
  
  G4ThreeVector v(xp, yp, zp), p(0, 2093.75, 576.25);
  G4cout << "T: " << theta << ", phi1: " << phi1 << ", phi2: "<< phi2 << ", R: " << R*sqrt(2.*(1-cos(phi2-phi1))) << ", Z: " << Z1 << G4endl;
  
  MEGen->GenerateMuon("Reyna",(G4float)theta);
  momentum = (G4double)(MEGen->GetEnergy())*GeV;

  G4VPhysicalVolume *Detector1Volume = DetCon->GetDetector1Volume();
  if (Detector1Volume) {
        // Physical volume found, do something with it
	// DetCon->PrintPhysicalVolume(Detector1Volume);
        p = SelectRandomPointInVolume(Detector1Volume);

    } else {
        // Physical volume not found
	G4cout << "Detector1 Volume not found" << G4endl;
    }

//  particleGun->SetParticlePosition(G4ThreeVector( R*cos(phi1)-gamma*xp, R*sin(phi1)-gamma*yp, Z1-gamma*zp));
  G4cout << "Initial position: (" << p.x() << "," << p.y() << "," << p.z() << ")" << G4endl;
  particleGun->SetParticlePosition(p);
  particleGun->SetParticleMomentumDirection(v);
  particleGun->SetParticleEnergy(momentum);
//  G4cout << "momentum = " << momentum*GeV << " GeV" << G4endl;
  //cout << particleGun->GetParticleEnergy() << endl;

  _dramiel->momentum_1 = momentum/MeV;

  _dramiel->posx_1 = p.x()/cm;
  _dramiel->posy_1 = p.y()/cm;
  _dramiel->posz_1 = p.z()/cm;

  _dramiel->vecx_1 = v.x();
  _dramiel->vecy_1 = v.y();
  _dramiel->vecz_1 = v.z();

  _dramiel->hit1 = true;
 
  particleGun->GeneratePrimaryVertex(anEvent);
  return;
}

G4ThreeVector PrimaryGeneratorAction::SelectRandomPointInVolume(G4VPhysicalVolume* physVolume) {
    if (!physVolume) {
        // Return a zero vector if the input pointer is null
        return G4ThreeVector();
    }

    G4AffineTransform transform = physVolume->GetObjectTranslation();

        // Get the transformation matrix
    G4ThreeVector translation = transform.NetTranslation();

        // Get the coordinates
        // G4double x = translation.x();
        // G4double y = translation.y();
        // G4double z = translation.z()

    // Get the bounding box of the physical volume
    G4LogicalVolume* logicalVolume = physVolume->GetLogicalVolume();
    G4Box* boundingBox = dynamic_cast<G4Box*>(logicalVolume->GetSolid());
    if (!boundingBox) {
        // Return a zero vector if the volume is not box-shaped
        return translation;
    }

    // Get the half-lengths of the bounding box
    G4double halfX = boundingBox->GetXHalfLength();
    G4double halfY = boundingBox->GetYHalfLength();
    G4double halfZ = boundingBox->GetZHalfLength();

    // Generate random coordinates within the bounding box
    G4double x = 2.0 * halfX * RNGenerator.genrand_real() - halfX;
    G4double y = 2.0 * halfY * RNGenerator.genrand_real() - halfY;
    G4double z = 2.0 * halfZ * RNGenerator.genrand_real() - halfZ;

    return translation + G4ThreeVector(x, y, z);
}

