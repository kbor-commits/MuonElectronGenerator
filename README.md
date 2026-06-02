# MuonElectronGenerator
Cosmic ray generator for GEANT4
//
// Id: MuonElectronGenerator.cc,v 0.1 2004/10/18 by KBor
//                              v 1.0 2005/04/19 by KBor
//                              v 1.1 2009/02/04 by KBor
//                              v 1.2 2011/06/06 by KBor
//                              v 1.3 2011/07/20 by KBor
//                              v 1.4 2012/08/29 by KBor
//                              v 1.5 2013/04/26 by KBor
//
// Recent versions include multiple recipes to generate muon spectrum
// instead of single recipe by Blanpied
// The otehr recipes include:
// Reyna - Zenith spectrum from PDG modified by Reyna's method
// Jokisch - 75 deg Jokisch et al. (PR D19, 5, 1368, 1979) spectrum
// Haruo - Jokisch (75 deg) spectrum transformed to LA altitude 2200 m
// Chris - @LA high-altitude spectrum (2200m)
//
// based of FORTRAN file mmuon.f by Gary Blanpied
// with multiple changes in implementation by KBor
// aimed to expand energy range to higher energies,
// smooth generated spectrum at lower energies and
// use less memory (data array of 300 logarithmic energy bins
// instead of 10000 linear energy bins)
//
// constructor of the class initializes variables and then
// function GenerateMuon() generates muon
//
// Here is usage example (in another class in GEANT4):
// avkPrimaryGeneratorAction::avkPrimaryGeneratorAction()
 // {
// ...
// MEGen = new MuonElectronGenerator();
//...
//}
//void avkPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
//{
//  MEGen->GenerateMuon();
//
//  G4ParticleDefinition* particle;
//  particle = MEGen->GetType();
//  particleGun->SetParticleDefinition(particle);
//  G4float momentum = (MEGen->GetEnergy())*GeV;
//  G4float mass = (G4float)particle->GetPDGMass();
//  G4float energy = sqrt (mass*mass + momentum*momentum);
//
//  // set particle kinetic energy to the gun/
//  particleGun->SetParticleEnergy(energy-mass);
//  G4double theta = (G4double)MEGen->GetZenithAngle();
//  ...
// }
// --------------------------------------------------------------
// Other comments:
// - iif input file with list of random numbers can be used for debugging
//   purposes. It may be obsolete, but I leave this feature (definetely,
//   not a bug!) commented for now
// - I consider use of single theta-momentum integral instead of
//   separate theta and momentum integrals. It would save one call
//   of random generator, but make the text more cryptic ... 
//
// Here is usage example (in another class in GEANT4):
// avkPrimaryGeneratorAction::avkPrimaryGeneratorAction()
