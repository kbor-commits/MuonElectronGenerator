#include "PrimaryGeneratorAction.hh"
#include "MuonElectronGenerator.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "MT_Random.hh"
#include "math.h"
#include <assert.h>

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
      //G4cout << theta << ", R: " << R*sqrt(2.*(1-cos(phi2-phi1))) << ", Z: " << Z1 << G4endl;
    } else {
      theta = MEGen->ChooseTheReyna();      
    }
  }

  double xp, yp, zp;
  xp = (cos(phi2)-cos(phi1))/sqrt(2.*(1-cos(phi2-phi1)))*sin(theta);
  yp = (sin(phi2)-sin(phi1))/sqrt(2.*(1-cos(phi2-phi1)))*sin(theta);
  zp = -cos(theta);
  G4double gamma = 50/(cos(theta))*mm;
  
  G4ThreeVector v(xp, yp, zp);
  //G4cout << "T: " << theta << ", phi1: " << phi1 << ", phi2: "<< phi2 << ", R: " << R*sqrt(2.*(1-cos(phi2-phi1))) << ", Z: " << Z1 << G4endl;
  
  MEGen->GenerateMuon("Reyna",(G4float)theta);
  momentum = (G4double)(MEGen->GetEnergy())*GeV;
  
  particleGun->SetParticlePosition(G4ThreeVector( R*cos(phi1)-gamma*xp, R*sin(phi1)-gamma*yp, Z1-gamma*zp));
  particleGun->SetParticleMomentumDirection(v);
  particleGun->SetParticleEnergy(momentum);
  //cout << particleGun->GetParticleEnergy() << endl;

  /*_dramiel->momentum_1 = momentum/MeV;
  _dramiel->posz_1 = R*cos(phi1)/cm;
  _dramiel->posy_1 = R*sin(phi1)/cm;
  _dramiel->posx_1 = -Z1/cm;
  _dramiel->vecz_1 = -v.x();
  _dramiel->vecy_1 = -v.y();
  _dramiel->vecx_1 =  v.z();*/
 
  particleGun->GeneratePrimaryVertex(anEvent);
  return;
}
