//
// Id: MuonElectronGenerator.cc,v 0.1 2004/10/18 by KBor
//                              v 1.0 2005/04/19 by KBor
//                              v 1.1 2009/02/04 by KBor
//                              v 1.2 2011/06/06 by KBor
//                              v 1.3 2011/07/20 by KBor
//                              v 1.4 2012/08/29 by KBor
//				v 1.5 2013/04/26 by KBor
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

#include "MuonElectronGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "assert.h"

int LIST_OUT=1;

MuonElectronGenerator::MuonElectronGenerator()
{
  CalculateCoeff();
  SetVariables();
  if(LIST_OUT) oof.open("muel_gen.list");

  // KBor 10/18/2004::
  // list of random numbers can be read instead
  // of using random generator
  // all "iif" lines should be enabled in this case

  //  iif.open("fort.34");
  RNGenerator.init_random(3759321, 22);
}

MuonElectronGenerator::MuonElectronGenerator(int MyPE)
{
  CalculateCoeff();
  SetVariables();
  if(LIST_OUT) oof.open("muel_gen.list");

  // KBor 10/18/2004::
  // list of random numbers can be read instead
  // of using random generator
  // all "iif" lines should be enabled in this case

  //  iif.open("fort.34");
  RNGenerator.init_random(6789417, MyPE);
}

MuonElectronGenerator::MuonElectronGenerator(MT_Random *RG)
{
  CalculateCoeff();
  SetVariables();
  if(LIST_OUT) oof.open("muel_gen.list");

  // KBor 10/18/2004::
  // list of random numbers can be read instead
  // of using random generator
  // all "iif" lines should be enabled in this case

  //  iif.open("fort.34");
  RNGenerator=*RG;
}

MuonElectronGenerator::MuonElectronGenerator(MT_Random *RG, string recipe_selector)
{
  CalculateCoeff();
  SetRecipe(recipe_selector);
  SetVariables(recipe_selector);
  if(LIST_OUT) oof.open("muel_gen.list");

  // KBor 10/18/2004::
  // list of random numbers can be read instead
  // of using random generator
  // all "iif" lines should be enabled in this case

  //  iif.open("fort.34");
  RNGenerator=*RG;
}

MuonElectronGenerator::~MuonElectronGenerator()
{
  if(LIST_OUT) oof.close();
  //  iif.close();
}

void MuonElectronGenerator::RunMuonElectronGen()
{
  // default particle kinematics
  string selector="";
  SetVariables();
  GenerateCosmic(selector);
}

void MuonElectronGenerator::GenerateMuon()
{
  GenerateMuon(recipe);
}

/* modified by KBor March 6, 2012 */
/* generate muon spectrum according to one of the recipes: */
/* "Blanpied" - based on Blanpied's FORTRAN code */
/* "Jokisch" - based on Jokisch et al. Phys. Rev.D 19, 1368 (1979) - muons at 75 deg */
/* "Reyna" - based on Reyna approximation hep-ph/0604145 */
void MuonElectronGenerator::GenerateMuon(string selector)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4float i = (float)RNGenerator.genrand_real();

  // only +/- muons:
  if( i < (0.45) ) {
    muel = particleTable->FindParticle(particleName="mu-");
  }
  else {
    muel = particleTable->FindParticle(particleName="mu+");
  }

  if( selector == "Blanpied" ) {
    theta_mu = ChooseTheta();
    energy_mu = ChooseEnergy(theta_mu);
  }
  else if( selector == "Jokisch" ) {
    theta_mu = ChooseTheta(Tnorm);
    energy_mu = ChooseJokEnergy((G4int)75);
  }
  else if( selector == "Haruo" ) {
    theta_mu = ChooseTheReyna()/Deg_to_Rad;
    energy_mu = ChooseJokEnergy((const G4float)theta_mu);
  }
  else if( selector == "Reyna" ) {
    theta_mu = ChooseTheReyna()/Deg_to_Rad;
    energy_mu = ChooseEnergy((G4float)0.)/cos(theta_mu*Deg_to_Rad);
  }
  else if( selector == "Chris" ) {
    theta_mu = ChooseTheReyna()/Deg_to_Rad;
    energy_mu = ChooseChrisEnergy((G4float)0.)/cos(theta_mu*Deg_to_Rad);
	}
	else { G4cout << "Unknown recipe, I don't know what to do" << G4endl;
		G4cout << "Will throw vertical monoenergetic muons now" << G4endl;
		theta_mu = 0.; 
		energy_mu = 3.*GeV;
	}

  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta_mu << G4endl;
  }
  return;
}

/* modified by KBor March 19, 2012 */
/* generate muon spectrum according to one of the recipes: */
/* "Blanpied" - based on Blanpied's FORTRAN code */
/* "Jokisch" - based on Jokisch et al. Phys. Rev.D 19, 1368 (1979) - muons at 75 deg */
/* "Reyna" - based on Reyna approximation hep-ph/0604145 */
/* this generates spectrum for a particular angle */
void MuonElectronGenerator::GenerateMuon(string selector, G4float theta)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4float i = (float)RNGenerator.genrand_real();

  // only +/- muons:
  if( i < (0.45) ) {
    muel = particleTable->FindParticle(particleName="mu-");
  }
  else {
    muel = particleTable->FindParticle(particleName="mu+");
  }

  theta_mu = theta;
  if( selector == "Blanpied" ) {
    energy_mu = ChooseEnergy(theta);
  }
  else if( selector == "Jokisch" || selector == "Haruo") {
    energy_mu = ChooseJokEnergy((const G4float)theta);
  }
  else if( selector == "Reyna" ) {
    energy_mu = ChooseEnergy((G4float)0.)/cos(theta*Deg_to_Rad);
  }
  else if( selector == "Chris" ) {
    energy_mu = ChooseChrisEnergy((G4float)0.)/cos(theta*Deg_to_Rad);
        }
        else { G4cout << "Unknown recipe, I don't know what to do" << G4endl;
                G4cout << "Will throw monoenergetic muons now" << G4endl;
                energy_mu = 3.*GeV;
        }

  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta << G4endl;
  }
  return;
}

/* modified by KBor July 19, 2011 */
/* generate muon spectrum for 75 degrees zenith angle */
/* Jokisch et al. 1979 */
void MuonElectronGenerator::GenerateMuon75()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4float i = (float)RNGenerator.genrand_real();

  // only +/- muons:
  if( i < (0.5) ) {
    muel = particleTable->FindParticle(particleName="mu-");
  }
  else {
    muel = particleTable->FindParticle(particleName="mu+");
  }
  theta_mu = ChooseTheta(Tnorm);
  // theta_mu = 75.;
  const G4int th0 = 75;
  energy_mu = ChooseEnergy(th0);

  //    theta_mu = ChooseTheta()*Deg_to_Rad;
  //    energy_mu = 1000.*ChooseEnergy(theta_mu);
  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta_mu/Deg_to_Rad << G4endl;
  }
  //    G4float mass = (G4float)muon->GetPDGMass();
  //    energy_mu = energy_mu - mass;
  return;
}

/* modified by KBor March 6, 2012 */
/* generate muon spectrum according to Reyna recipe */
/* hep-ph/0604145 */
void MuonElectronGenerator::GenerateReynaMuon()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4float i = (float)RNGenerator.genrand_real();

  // only +/- muons:
  if( i < (0.45) ) {
    muel = particleTable->FindParticle(particleName="mu-");
  }
  else {
    muel = particleTable->FindParticle(particleName="mu+");
  }
  theta_mu = ChooseTheReyna();
  energy_mu = ChooseEnergy((G4float)0.)/cos(theta_mu);

  //    theta_mu = ChooseTheta()*Deg_to_Rad;
  //    energy_mu = 1000.*ChooseEnergy(theta_mu);
  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta_mu/Deg_to_Rad << G4endl;
  }
  //    G4float mass = (G4float)muon->GetPDGMass();
  //    energy_mu = energy_mu - mass;
  return;
}

/* modified by KBor Feb 4, 2009 */
/* modified/renamed by KBor Apr, 2012 */
void MuonElectronGenerator::GenerateCosmic()
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  const G4float EMu_Fraction = 1.; // fraction of electrons relative to muons

  G4int i = (int)((1.+EMu_Fraction)*RNGenerator.genrand_real());

  // only +/- muons:
  // G4int i = (int)(RNGenerator.genrand_real());

  //only +/- electrons:
  //G4int i = (int)(EMu_Fraction*RNGenerator.genrand_real()+1.);

  if (i < 1.) {
    if( i < 0.5 ) {
      muel = particleTable->FindParticle(particleName="mu-");
    }
    else {
      muel = particleTable->FindParticle(particleName="mu+");
    }
	  theta_mu = ChooseTheta(Tnorm);
	// KBor - July 11, 2011 - this looks like a hack to generate muon spectrum for zenith angle only  
    // energy_mu = ChooseEnergy(&ELN[0][0]);
	  energy_mu = ChooseEnergy( theta_mu);
  }
  else {
    if( i < 1.+(EMu_Fraction/2.) ) {
      muel = particleTable->FindParticle(particleName="e-");
    }
    else {
      muel = particleTable->FindParticle(particleName="e+");
    }
    theta_mu = ChooseTheta(TElnorm);
    energy_mu = ChooseElEnergy();
  }

  //    theta_mu = ChooseTheta()*Deg_to_Rad;
  //    energy_mu = 1000.*ChooseEnergy(theta_mu);
  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta_mu << G4endl;
  }

  //    G4float mass = (G4float)muon->GetPDGMass();
  //    energy_mu = energy_mu - mass;
}

/* modified by KBor Feb 4, 2009 */
/* modified/renamed by KBor Apr, 2012 */
/* modified by JoPerry/KBor Apr 2013 */
void MuonElectronGenerator::GenerateCosmic(string recipe, G4float theta)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  const G4float EMu_Fraction = 1.; // fraction of electrons relative to muons

  G4float i = (G4float)((1.+EMu_Fraction)*RNGenerator.genrand_real());

  if (i < 1.) {
    GenerateMuon(recipe, theta);
  }
  else {
    if( i < 1.+(EMu_Fraction/2.) ) {
      muel = particleTable->FindParticle(particleName="e-");
    }
    else {
      muel = particleTable->FindParticle(particleName="e+");
    }
    theta_mu = ChooseTheta(TElnorm);
    energy_mu = ChooseElEnergy();
  }

  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta_mu << G4endl;
  }
}



/* modified by KBor Feb 4, 2009 */
/* string selector is added by John Perry */
void MuonElectronGenerator::GenerateCosmic(string selector)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  const G4float EMu_Fraction = 0.5; // fraction of electrons relative to muons

  //G4float i = (float)((1.+EMu_Fraction)*RNGenerator.genrand_real());
  G4float i = (float)RNGenerator.genrand_real();

  if (selector=="m")
    {
      // only +/- muons:
      //G4float i = (float)RNGenerator.genrand_real();
      if( i < (0.5) ) {
	muel = particleTable->FindParticle(particleName="mu-");
      }
      else {
	muel = particleTable->FindParticle(particleName="mu+");
      }
		theta_mu = ChooseTheta(Tnorm);
		// KBor - July 11, 2011 - this looks like a hack to generate muon spectrum for zenith angle only  
		// energy_mu = ChooseEnergy(&ELN[0][0]);
		energy_mu = ChooseEnergy( theta_mu);
    }
  else if (selector=="e")
    {
      //only +/- electrons:
      //G4float i = (float)RNGenerator.genrand_real();
      if( i < (0.5 ) ) {
	muel = particleTable->FindParticle(particleName="e-");
      }
      else {
	muel = particleTable->FindParticle(particleName="e+");
      }
      theta_mu = ChooseTheta(TElnorm);
      energy_mu = ChooseElEnergy();

    }

  else if (selector=="b")
    {
      // both +/- muons and electrons:
      if(i>EMu_Fraction)
	{
	  if( i > (0.5*(1.0+EMu_Fraction)) ) {
	    muel = particleTable->FindParticle(particleName="mu-");
	  }
	  else {
	    muel = particleTable->FindParticle(particleName="mu+");
	  }
		theta_mu = ChooseTheta(Tnorm);
		// KBor - July 11, 2011 - this looks like a hack to generate muon spectrum for zenith angle only  
		// energy_mu = ChooseEnergy(&ELN[0][0]);
		energy_mu = ChooseEnergy( theta_mu);
	}
      else
	{
	  if( i < (0.5*EMu_Fraction ) ) {
	    muel = particleTable->FindParticle(particleName="e-");
	  }
	  else {
	    muel = particleTable->FindParticle(particleName="e+");
	  }
	  theta_mu = ChooseTheta(TElnorm);
	  energy_mu = ChooseElEnergy();
	}
    }
  else
    assert(0);

  //    theta_mu = ChooseTheta()*Deg_to_Rad;
  //    energy_mu = 1000.*ChooseEnergy(theta_mu);
  if( LIST_OUT ) {
    oof.setf(ios::fixed);
    oof << setw(10) << setprecision(4) << energy_mu
	<< setw(10) << theta_mu << G4endl;
  }
  //    G4float mass = (G4float)muon->GetPDGMass();
  //    energy_mu = energy_mu - mass;
}

// Added by KBor April 2012
void MuonElectronGenerator::SetRecipe(string selector)
{
  static char r[10];
  strcpy( r, selector.c_str());
  recipe = (string)r;
  cout << "Spectrum according to " << recipe << " recipe" << endl;
}

// KBor 10/18/2004::
// we prefer here to read data from file,
// but if there is no data to read we calculate the spectrum
/* modified by KBor 08/29/2012 */
void MuonElectronGenerator::SetVariables()
{
  energy_mu = 3.*GeV;
  theta_mu = 0.;
  recipe = "Blanpied";
  Nel_bins=CalcElVariables();
  if(!ReadVariables("tnorm.list", "enorm.list")) {
	CalcVariables(); ReadVariables("tnorm.list", "enorm.list");
	}
  SetReyna();
  ReadJoks();
  return;
}

// KBor 10/18/2004::
// we prefer here to read data from file,
// but if there is no data to read we calculate the spectrum
/* modified by KBor 08/29/2012 */
void MuonElectronGenerator::SetVariables(string recipe_selector)
{
  energy_mu = 3.*GeV;
  theta_mu = 0.;
  recipe = recipe_selector;
  Nel_bins=CalcElVariables();
  if(!ReadVariables("tnorm.list", "enorm.list")) {
        CalcVariables(); ReadVariables("tnorm.list", "enorm.list");
        }
  if(recipe=="Reyna")SetReyna();
  else if(recipe=="Jokisch") ReadJoks();
  else if(recipe=="Haruo") ReadHaruo();	
  else if(recipe=="Chris") ReadChris();
  return;
}


// KBor 3/6/2012::
// we prefer here to read data from file,
// but if there is no data to read we calculate the angle distribution
void MuonElectronGenerator::SetReyna()
{
  if(!ReadReyna("Reyna.list")) { ReynaVariables(); ReadReyna("Reyna.list"); }
  return;
}

G4float MuonElectronGenerator::ThetaReyna(G4int i)
{
  return i*ZANG_MAX*Deg_to_Rad/N_ZANG;
}

// KBor 3/5/2012::
// calculating variables for Reyna scaling of the cosmic-rays
// hep-ph/0604145
void MuonElectronGenerator::ReynaVariables()
{
  ofstream of;
  of.open("Reyna.list");
  of.setf(ios::fixed);

  G4double cosine, total=0., the, sine;
  G4double fReyna[N_ZANG+1];

  for(G4int i = 0; i <= N_ZANG; i++) {
    the = ThetaReyna(i);
    cosine = cos(the);
    sine = sin(the);
    fReyna[i] = cosine*cosine*cosine*sine;
    total += fReyna[i];
    //    of << fReyna[i] << endl;
  }
  total-=(fReyna[0]+fReyna[N_ZANG])/2.;
  cosine = 0.;
  for (G4int i = 0; i < N_ZANG; i++) {
    cosine += (fReyna[i]+fReyna[i+1])/2./total;
    of << setprecision(12) << cosine << " ";
  }

  of << endl;
  of.close();
}

// KBor 10/18/2004::
// if there is no data file to read we calculate
// the input spectrum of cosmic ray muons
// this function is actually written for N_ANGLES=90,
// so N_ANGLES should be changed in header file for this
// function to work properly
void MuonElectronGenerator::CalcVariables()
{
  G4cout << "Calculating variables" << G4endl;
  G4int i, j=0, ii;
  // init to zero
  for(i = 0; i <N_ANGLES; i++)
    {
      theta[i]=1.*(i+1);
      SigmaT[i]=0.;
      Tnorm[i]=0.;
    }
  for(i = 0; i < N_ENERGY; i++)
    {
      SigmaE[i]=0.;
      E[i]=0.01*(i+1);
    }

  ELog[0]=.001;
  for(i=1; i<NLog; i++)
    {
      ELog[i]=ELog[i-1]*1.05;
      //      G4cout << i << " " << j << " " << ELog[j] << " " << E[i] << G4endl;
    }

  for(i=0; i<N_ENERGY; i++)
    {
      while(ELog[j]<E[i] && j<NLog)
	{
	  //	  G4cout << j << " " << i << " " << ELog[j] << " " << E[i] << G4endl;
	  for(ii=N_ENERGY-1; ii>i; ii--) E[ii]=E[ii-1];
	  E[i]=ELog[j];
	  j++;
	}
      if(ELog[j]==E[i]) j++;
    }
  //  G4cout << i << " " << j << " " << ELog[j] << " " << E[i] << G4endl;

  for(ii=j; ii<NLog; ii++) E[N_ENERGY-NLog+ii]=ELog[ii];

  for(i=0; i<N_ENERGY; i++) cout << "En.chan = " << i << " : " << E[i] << endl;

  G4cout.setf(ios::fixed);
  for(i = 0; i < N_ENERGY; i++)
    {
      //      G4cout << setprecision(5) << E[i] << " ";
      //      if((i+1)%8 == 0) G4cout << G4endl;
    }

  //KBor 10/18/2004:: I'm not sure this comment is still valid for this version
  //
  // Calculate with model. Pinss are produced and decay at 15km and follow
  // line to surface of Earth. The fraction 1/cosine(theta) measures
  // the amount of atmosphere.  Energy loss computed in 15 steps which are
  // 1 km thick at zero deg.

  G4double fact;

  ofstream of;
  of.open("sigma.list");
  of.setf(ios::fixed);

  //  of.precision(6);

  for(i = 0; i <N_ANGLES; i++)
    {
      //      thet = theta[i];
      //      c = cos(thet/57.3);
      for (int j=0; j<N_ENERGY; j++)
	{
	  Enorm[i][j]=0.;  //init
	  fact=a_coeff[i]*log10(E[j])*log10(E[j]) + b_coeff[i]*log10(E[j])
	    + c_coeff[i];
	  Sigma[i][j] = pow(10.,fact);
	  if(i==60 && j<2000) of << setw(10) << setprecision(3) << E[j]
				 << setprecision(6) << setw(10) << Sigma[0][j]
				 << setw(10) << Sigma[9][j]
				 << setw(10) << Sigma[19][j]
				 << setw(10) << Sigma[29][j]
				 << setw(10) << Sigma[39][j]
				 << setw(10) << Sigma[49][j]
				 << setw(10) << Sigma[59][j] << G4endl;
	}
    }
  of.unsetf(ios::fixed);
  of.close();

  // Integrate.Variable sumT is over dE and dtheta
  // Variable SigmaT is over dE
  //  SigmaT[0]=0.;
  for(int j=1; j<N_ENERGY; j++)
    {
      SigmaT[0]=SigmaT[0]+(Sigma[0][j-1]+Sigma[0][j])*(E[j]-E[j-1])*
	sin(theta[0]/57.3)*theta[0]*Pi;
    }
  for (i=1; i<N_ANGLES; i++)
    {
      for(int j=1; j<N_ENERGY; j++)
	{
	  SigmaT[i]=SigmaT[i]+(Sigma[i][j-1]+Sigma[i][j])*(E[j]-E[j-1])*
	    sin(theta[i]/57.3)*(theta[i]-theta[i-1])*Pi;
	}
    }

  G4double sumT=0.;
  for (i=0; i<N_ANGLES; i++)
    {
      sumT=sumT+SigmaT[i];
    }

  // Normalized SigmaT values
  // Tnorm is the function to choose angle in random. It goes from 0 to 1.
  of.open("tnorm.list");
  of.setf(ios::fixed);

  G4float tint=0.;
  //  Tnorm[0]=0.;
  for (i=0; i<N_ANGLES; i++)
    {
      SigmaT[i]=SigmaT[i]/sumT;
      tint+=SigmaT[i];
      Tnorm[i]=tint;
      //      cout << "Tnorm: " << i << " -- " << Tnorm[i] << endl;
      of << setw(5) << theta[i]
	 << setprecision(6) << setw(10) << Tnorm[i];
      of << G4endl;
    }

  of.close();

  // We will choose an angle with Tnorm,
  // then will have N_ANGLES energy distributions
  // These are created below: For each of N_ANGLES angles there are
  // N_ENERGY energy weights

  for (i=0; i<N_ANGLES; i++)
    {
      // for each angle we compute sumE and SigmaE(jjj)
      // so these are written over 19 times
      G4double sumE = 0.;
      for(int j=1; j<N_ENERGY; j++)
	{
	  SigmaE[j]=Sigma[i][j]*(E[j]-E[j-1]);
	  sumE += SigmaE[j];
	}
      // Enorm(iii,jjj) are used to select energy after angle is selected
      for(int j=1; j<N_ENERGY; j++)
	{
	  Enorm[i][j]=Enorm[i][j-1]+SigmaE[j]/sumE;
	}
      Enorm[i][0]=Enorm[i][1]/Enorm[i][2]*Enorm[i][1];
    }

  of.open("enorm.list");
  of.setf(ios::scientific);
  //  G4int ebin = 1;
  for(int i=0; i<N_ANGLES; i++)
    {
      ii=0;
      G4cout << i << ": " << theta[i] << G4endl;
      for(int j=0; j<N_ENERGY; j++)
	{
	  if( E[j]==ELog[ii] )
	    {
	      //	      ELN[i][ii]=Enorm[i][j];
	      of << setw(7) << setprecision(2) << theta[i]
		 << setw(15) << setprecision(9) << E[j]
	         << setw(15) << setprecision(9) << Enorm[i][j];
	      of << G4endl;
	      //	      if(ii%10==0) of << G4endl;
	      ii++;
	      //	      ebin=1;
	    }
	  //	  else ebin++;
	}
    }

  of.close();

  //  cout << "Variables ready" << endl;
  return;
}

// KBor 12/08/2004::
// if there is no data file to read we calculate
// the input spectrum of cosmic ray electrons
G4int MuonElectronGenerator::CalcElVariables()
{
  G4cout << "Calculating variables for electron distribution" << G4endl;
  G4int i, j=0, ii;

  // deal with angles first
  G4float diff[N_ANGLES];
  G4float sumd = 0.;

  for(i = 0; i <N_ANGLES; i++)
    {
      theta[i]=1.*(i+1);
      diff[i]=(cos((i+1)/57.3)*cos((i+1)/57.3))*sin((i+1)/57.3);
      sumd+=diff[i];
      TElnorm[i]=sumd;
    }
  for(i = 0; i <N_ANGLES; i++)
    {
      TElnorm[i]/=sumd;
    }

  // now deal with energy
  for(i = 0; i < N_ENERGY; i++)
    {
      E[i]=10.+i;    /* energy in MeV */
    }

  EElLog[0]=10.;           /* log energy bins also in MeV */
  for(i=1; i<NLog; i++)
    {
      EElLog[i]=EElLog[i-1]*1.05;
      //      G4cout << i << " " << j << " " << EElLog[j] << " " << E[i] << G4endl;
    }
  for(i=0; i<N_ENERGY; i++)
    {
      while(EElLog[j]<E[i] && j<NLog)
	{
	  //	  G4cout << j << " " << i << " " << ELog[j] << " " << E[i] << G4endl;
	  for(ii=N_ENERGY-1; ii>i; ii--) E[ii]=E[ii-1];
	  E[i]=EElLog[j];
	  j++;
	}
      if ( EElLog[j] == E[i] ) j++;
    }
  //  G4cout << i << " " << j << " " << EElLog[j] << " " << E[i] << G4endl;

  for(ii=j; ii<NLog; ii++) E[N_ENERGY-NLog+ii]=EElLog[ii];

  G4cout.setf(ios::fixed);
  for(i = 0; i < N_ENERGY; i++)
    {
      //      G4cout << setprecision(5) << E[i] << " ";
      //      if((i+1)%8 == 0) G4cout << G4endl;
    }


  G4String file_sig = "sigma_el.list";
  ifstream ifs;
  ofstream of;
  int exist = 1;

  ifs.open(file_sig);
  if(!ifs) {
    ifs.close();
    exist = 0;
    of.open("sigma_el.list");
    of.setf(ios::fixed);
  }

  //  of.precision(6);
  G4double sume=0.;
  G4double SumE[N_ENERGY];
  /*  for (i=0; i<N_ENERGY; i++)
    {
      if( E[i] <= 10. ) SigmaE[i]=0.;
      else if( E[i] <= 600. ) SigmaE[i]=200.*pow((double)E[i], -1.8);
      else SigmaE[i]=160000.*pow((double)E[i], -2.9);
      sume+=SigmaE[i];
      SumE[i]=sume;
    }
  */
  G4float coef=100.;
  G4float E0=10., al0=-1.6;
  G4float E1=100., al1=al0;
  G4float E2=1400.;

  for (i=0; i<N_ENERGY; i++)
    {
      if( E[i] <= E0 ) SigmaE[i]=0.;
      else {
	if( E[i] <= E1 ) {
	  al1=al0;
	}
	else if( E[i] <= E2 ) {
	  al1-=(E[i]-E[i-1])/1000.;
	  coef*=pow((double)(E[i]), (double)(al0-al1));
	  al0=al1;
	}
	SigmaE[i]=coef*pow((double)(E[i]), (double)(al0))*(E[i]-E[i-1]);
      }
      sume+=SigmaE[i];
      SumE[i]=sume;

      if(!exist) of << SigmaE[i] << " " << sume << G4endl;
    }
  if(exist) {
    ifs.close();
  }
  else {
    of.unsetf(ios::fixed);
    of.close();
  }

  G4String file_eel = "eelnorm.list";
  exist = 1;

  ifs.open(file_eel);
  if(!ifs) {
    ifs.close();
    exist = 0;
    of.open(file_eel);
    of.setf(ios::scientific);
  }

  ii=0;
  for(int j=0; j<N_ENERGY; j++)
    {
      if( E[j]==EElLog[ii] )
	{
	  EELN[ii]=SumE[j]/sume;
	  EElLog[ii]/=1000.;
	  if(!exist) {
	    of << ii << " "
	       << setw(15) << setprecision(9) << EElLog[ii] << " "
	       << setw(15) << setprecision(9) << EELN[ii];
	    of << G4endl;
	  }
	  //	      if(ii%10==0) of << G4endl;
	  ii++;
	}
    }

  if(exist) {
    ifs.close();
  }
  else {
    of.unsetf(ios::fixed);
    of.close();
  }

  //  cout << "Variables ready" << endl;
  return ii;
}

// KBor 3/5/2012::
// file1 should contain N_ZANG theta normalization values
//
G4int MuonElectronGenerator::ReadReyna(G4String file1)
{
  ifstream eth_integral;
  //  G4float thet, energy, eth_norm;

  eth_integral.open(file1);
  if(!eth_integral) return 0;

  for(int i = 0; i < N_ZANG; i++)
    {
      eth_integral >> theReyNorm[i];
    }
  eth_integral.close();

  return 1;
}

// KBor 10/18/2004::
// files "file1" and "file2" are supposed to contain necessary data
// in required format
// these files are created by CalcVariables() function
// as "tnorm.list" and "enorm.list"
G4int MuonElectronGenerator::ReadVariables(G4String file1, G4String file2)
{
  ifstream eth_integral;
  G4float thet, energy, eth_norm;

  eth_integral.open(file1);
  if(!eth_integral) return 0;

  for(int i = 0; i < N_ANGLES; i++)
    {
      eth_integral >> thet >> eth_norm;
      theta[i]=thet;
      Tnorm[i]=eth_norm;
    }
  eth_integral.close();

  eth_integral.open(file2);

  if(!eth_integral) return 0;
  G4cout << "Reading Variables" << G4endl;

  // KBor 18/10/2004::
  // only angles 0, 3, 6, 9, 12 ... 90 are used in this version
  // so the angle array length is 30 (N_ENGLES)
  // alternatively, all 90 (N_ANGLES) angles can be used, it will require
  // change of N_ENGLES value and modification of the loop below
  // KBor 7/12/2011:: introduced a new variable, konv, so that the loop would work for
  // any reasonable values of N_ANGLES and N_ENGLES
  int konv = (int)(N_ANGLES/N_ENGLES);   G4cout << "konv = " << konv << endl;
  int ii;
  for(int i = 0; i < N_ANGLES; i++)
    {
      ii=(int)(i/konv);
      for(int j = 0; j < NLog; j++)
	{
	  eth_integral >> thet >> energy >> eth_norm;
	  //	  G4cout << thet << " " << energy << " " << eth_norm << G4endl;
	  if ((int)(thet/konv)*konv == (int)thet)
	    {
	      // KBor 10/18/2004::
	      // if you are curious, ELog is an array of energy bins,
	      // theta is an array of angle bins
	      // and ELN is integral normalized spectrum
	      // (fraction of the muons below energy ELog,
	      // goes from 0 for Elog=0 to 1 for Elog=infinity)
	      // of course, ELN is just calculated approximately
	      // as shown in CalcVariables() function
	      if (ii < N_ENGLES)
		{
		  ELog[j]=energy;
		  Th[ii]=thet;
		  ELN[ii][j]=eth_norm;
		}
	    }
	}
    }

  eth_integral.close();
  return 1;
}

// KBor 10/18/2004::
// I assume here that muon spectrum is a power-law function
// F=A*E^al or logF = al*logE + b (where b = logA)
// we find al from averaged ELN values:
// al = log(ENprev/ENnext)/log(E1/E2)
// and then calculate (randomized) E as
// E = E1*(x/ENprev)^(1/al), where x is a random number
// it seems to work fine with bins chosen logarithmically
G4float MuonElectronGenerator::ChooseEnergy(G4float thet)
{
  int j, jj;
  G4double l1, l2, al;
  // KBor 10/18/2004::
  // instead of using ELN numbers from table we
  // use their weighted average based on chosen value of theta
  // EN(theta) = ELN1 + (theta-theta1)/(theta2-theta1)*(ELN2-ELN1)
  // ENprev and ENnext is for two different energy values
  G4float ENprev, ENnext;

  for (j=1; j<N_ENGLES; j++)
    {
      if( thet < Th[j] ) break;
    }

  float x = RNGenerator.genrand_real();
  G4float enex;
  G4float coef=(thet-Th[j-1])/(Th[j]-Th[j-1]);

  for (jj=1; jj<NLog; jj++)
    {
      ENnext = ELN[j-1][jj] + coef*(ELN[j][jj]-ELN[j-1][jj]);
      if( x < ENnext )
	{
	  ENprev = ELN[j-1][jj-1] + coef*(ELN[j][jj-1]-ELN[j-1][jj-1]);
	  l1=log(ENnext/ENprev);
	  l2=log(ELog[jj]/ELog[jj-1]);
	  al=l2/l1;
	  enex = pow((double)x/ENprev, al) * ELog[jj-1];
	  if ( enex < 0.01 )
	    {
	      G4cout << G4endl << "mu low energy = "
		     <<  enex*1000. << " MeV";
	    }
	return enex;
	}
    }

  G4cout << G4endl << "mu high energy = "
	 <<  ELog[NLog-1] << " GeV"
	 << j << ": " << thet << " " << x;
  //  G4cout <<  jj << ": " << ELog[jj-1] << G4endl;
  //  G4cout <<  ENprev << " < " << ENnext << G4endl;

  return ELog[NLog-1];
}

// KBor 04/26/2013::
// theta is assumed to be in degrees
G4float MuonElectronGenerator::ChooseChrisEnergy(const G4float theta)
{
// I assume here that muon spectrum is a power-law function
// F=A*E^al or logF = al*logE + b (where b = logA)
// we find al from averaged ELN values:
// al = log(ENprev/ENnext)/log(E1/E2)
// and then calculate (randomized) E as
// E = E1*(x/ENprev)^(1/al), where x is a random number
// it seems to work fine with bins chosen logarithmically
  int jj;
  G4double l1, l2, al;
  G4float ENprev=0., ENnext=1.;

  if( EC[0] != (G4float)0.05 ) { 
	G4cout << "EC[0] = " << EC[0] << G4endl; ReadChris();
	}
  float x = RNGenerator.genrand_real();
  G4float enex;

  for (jj=1; jj<N_EC; jj++)
    {
      ENnext = ELNC[jj];
      if( x < ENnext )
        {
          ENprev = ELNC[jj-1];
          l1=log(ENnext/ENprev);
          l2=log(EC[jj]/EC[jj-1]);
          al=l2/l1;
          enex = pow((double)x/ENprev, al) * EC[jj-1];
          if ( enex < 0.001 )
            {
              G4cout << G4endl << "mu low energy = "
                     <<  enex*1000. << " MeV "
                     <<  x;
            }
        return enex;
        }
    }

  G4cout << G4endl << "mu high energy = "
         <<  EC[N_EC-1]*pow((double)(1-ELNC[N_EC-1])/(1-x), (double)(0.57)) << " GeV"
         << " " << x;
  //  G4cout <<  jj << ": " << ELog[jj-1] << G4endl;
  //  G4cout <<  ENprev << " < " << ENnext << G4endl;

  return EC[N_EC-1]*pow((double)(1-ELNC[N_EC-1])/(1-x), (double)(0.57));
}

// KBor 04/24/2013::
// generalization of the function to all angles, not just 75 degrees
// theta is assumed to be in degrees
G4float MuonElectronGenerator::ChooseJokEnergy(const G4float theta)
{
	return(ChooseJokEnergy((const G4int)75)/cos(theta*Deg_to_Rad)*cos(75.*Deg_to_Rad));
}


// KBor 10/18/2004::
// I assume here that muon spectrum is a power-law function
// F=A*E^al or logF = al*logE + b (where b = logA)
// we find al from averaged ELN values:
// al = log(ENprev/ENnext)/log(E1/E2)
// and then calculate (randomized) E as
// E = E1*(x/ENprev)^(1/al), where x is a random number
// it seems to work fine with bins chosen logarithmically
G4float MuonElectronGenerator::ChooseJokEnergy(const G4int itheta75)
{
  int jj;
  G4double l1, l2, al;
  G4float ENprev=0., ENnext=1.;

  if( itheta75 != 75 ) return 0.;  // only works for angle == 75 degrees

  if( E75[0] != 1. ) ReadJoks();

  float x = RNGenerator.genrand_real();
  G4float enex;

  for (jj=1; jj<N_E75; jj++)
    {
      ENnext = ELN75[jj];
      if( x < ENnext )
	{
	  ENprev = ELN75[jj-1];
	  l1=log(ENnext/ENprev);
	  l2=log(E75[jj]/E75[jj-1]);
	  al=l2/l1;
	  enex = pow((double)x/ENprev, al) * E75[jj-1];
	  if ( enex < 0.001 )
	    {
	      G4cout << G4endl << "mu low energy = "
		     <<  enex*1000. << " MeV "
		     <<  x;
	    }
	return enex;
	}
    }

  G4cout << G4endl << "mu high energy = "
	 <<  E75[N_E75-1]*pow((double)(1-ELN75[N_E75-1])/(1-x), (double)(0.57)) << " GeV"
	 << " " << x;
  //  G4cout <<  jj << ": " << ELog[jj-1] << G4endl;
  //  G4cout <<  ENprev << " < " << ENnext << G4endl;

  return E75[N_E75-1]*pow((double)(1-ELN75[N_E75-1])/(1-x), (double)(0.57));
}

// KBor 12/08/2004::
// I assume here that muon spectrum is a power-law function
// F=A*E^al or logF = al*logE + b (where b = logA)
// we find al from averaged ELN values:
// al = log(ENprev/ENnext)/log(E1/E2)
// and then calculate (randomized) E as
// E = E1*(x/ENprev)^(1/al), where x is a random number
// it seems to work fine with bins chosen logarithmically
G4float MuonElectronGenerator::ChooseElEnergy()
{
  int jj;
  G4double l1, l2, al;
  G4float ENprev=0., ENnext=0.;

  //  for (j=1; j<N_ENGLES; j++)
  //    {
  //      if( theta_mu < theta[j] ) break;
  //    }

  float x = RNGenerator.genrand_real();
  G4float enex;

  for (jj=1; jj<Nel_bins; jj++)
    {
      ENnext = EELN[jj];
      if( x < ENnext )
	{
	  ENprev = EELN[jj-1];
	  if (ENprev == 0.)
	    {
	      al = .68;
	      ENprev=ENnext;
	    }
	  else {
	    l1=log(ENnext/ENprev);
	    l2=log(EElLog[jj]/EElLog[jj-1]);
	    al=l2/l1;
	  }
	  //	  if (jj < 3 ) G4cout << "low-energy electron " << jj << ":"
	  //			  << l1 << "," << l2 << "," << al << G4endl;
	  enex = pow((double)x/ENprev, al) * EElLog[jj-1];

	  if ( ENnext == ENprev && enex < EElLog[jj-1]/10.)
	    {
               G4cout << G4endl << "El low energy = " << enex*1000. << " MeV "
		     << x/ENprev << " , " << EElLog[jj-1];
	      //		  G4cout <<  EElLog[jj-1] << " < " << enex 
	      //	  << " < " << EElLog[jj] << G4endl;
	      //	  G4cout <<  ENprev << " < " << ENnext << G4endl;
	    }
	  return enex;
	}
    }

  if( (1.-x) < 1.e-10) enex = 1.e10;
  else enex = pow((double)(1.-ENnext)/(1.-x), 1./1.9) * EElLog[Nel_bins-1];

  G4cout << "El high energy = " << enex << " GeV ";
  G4cout << x << G4endl;
  //  G4cout <<  jj-2 << ": " << EElLog[jj-2] << " , " <<  EELN[jj-2] << G4endl;
  //  G4cout <<  Nel_bins << ": " << EElLog[Nel_bins-1] << " , " <<  ENnext << G4endl;
  //  G4cout <<  jj << ": " << EElLog[jj-1] << " , " <<  ENnext << G4endl;

  return enex;
}

// KBor 10/18/2004::
// I assume here that muon spectrum is a power-law function
// F=A*E^al or logF = al*logE + b (where b = logA)
// we find al from averaged ELN values:
// al = log(ENprev/ENnext)/log(E1/E2)
// and then calculate (randomized) E as
// E = E1*(x/ENprev)^(1/al), where x is a random number
// it seems to work fine with bins chosen logarithmically
G4float MuonElectronGenerator::ChooseEnergy(G4float* eln)
{
  int j, jj;
  G4double l1, l2, al;
  // KBor 10/18/2004::
  // instead of using ELN numbers from table we
  // use their weighted average based on chosen value of theta
  // EN(theta) = ELN1 + (theta-theta1)/(theta2-theta1)*(ELN2-ELN1)
  // ENprev and ENnext is for two different energy values
  G4float ENprev, ENnext;

  for (j=1; j<N_ENGLES; j++)
    {
      if( theta_mu < Th[j] ) break;
    }

  float x = RNGenerator.genrand_real();
  G4float enex;
  G4float coef=(theta_mu - Th[j-1])/(Th[j]-Th[j-1]);

  for (jj=1; jj<NLog; jj++)
    {
      ENnext = (*(eln+(j-1)*NLog+jj)) + coef*

	( (*(eln+j*NLog+jj)) - (*(eln+(j-1)*NLog+jj)) );
      if( x < ENnext )
	{
	  ENprev = (*(eln+(j-1)*NLog+jj-1)) + coef*
	    ( (*(eln+j*NLog+jj-1)) - (*(eln+(j-1)*NLog+jj-1)) );
	  l1=log(ENnext/ENprev);
	  l2=log(ELog[jj]/ELog[jj-1]);
	  al=l2/l1;
	  enex = pow((double)x/ENprev, al) * ELog[jj-1];
	  if ( enex < 0.01 )
	    {
               G4cout << G4endl << "mu low energy: "
		     <<  enex*1000. << " MeV";
	      //	      G4cout <<  ENprev*1000. << " < " << ENnext*1000. << G4endl;
	    }
	  return enex;
	}
    }

  G4cout << j << ": " << theta_mu << " " << x << G4endl;
  G4cout <<  jj << ": " << *(eln+(j-1)*NLog+jj) << " " << ELog[jj-1] << G4endl;
  G4cout <<  "coef = " << coef << " " << ENnext << G4endl;

  return ELog[NLog-1];
}

G4float MuonElectronGenerator::ChooseTheReyna()
{
  G4float x = RNGenerator.genrand_real();
  // iif >> x;

  for (int j=1; j<N_ZANG; j++)
    {
      if(x < theReyNorm[j]) {
	G4float thet=ThetaReyna(j-1)+
	  (ThetaReyna(j)-ThetaReyna(j-1))*(x-theReyNorm[j-1])/(theReyNorm[j]-theReyNorm[j-1]);
	return thet;
      }
    }

  return ThetaReyna(0);
}

G4float MuonElectronGenerator::ChooseTheta()
{
  G4float x = RNGenerator.genrand_real();
  // iif >> x;

  for (int j=1; j<N_ANGLES; j++)
    {
      if(x < Tnorm[j]) {
	G4float thet=theta[j-1]+
	  (theta[j]-theta[j-1])*(x-Tnorm[j-1])/(Tnorm[j]-Tnorm[j-1]);
	return thet;
      }
    }

  return theta[0];
}

G4float MuonElectronGenerator::ChooseTheta(G4float* tnorm)
{
  G4float x = RNGenerator.genrand_real();
  // iif >> x;

  for (int j=1; j<N_ANGLES; j++)
    {
      if(x < *(tnorm+j)) {
	G4float thet=theta[j-1]+
	  (theta[j]-theta[j-1])*(x - (*(tnorm+j-1)))/
	  ( (*(tnorm+j)) - (*(tnorm+j-1)));
	return thet;
      }
    }

  return theta[0];
}

void MuonElectronGenerator::CalculateCoeff()
{
  G4float pa[5] = {-0.8816e-4, -0.1117e-3,  -0.64, -0.1966e-1, 0.2040e-1};
  G4float pb[5] = { 0.4169e-2, -0.9891e-4,    3.7, -4.3118,   -0.9235e-3};
  G4float pc[5] = {-0.3516e-3,  0.8661e-2, -2.5985,-0.8745e-5,-0.1457};

  for (G4int i = 1; i <= N_ANGLES; i++)
    {
      a_coeff[i-1] = pa[0]/((1./i)+pa[1]*i) + pa[2] + pa[3]*exp(-pa[4]*i);
      b_coeff[i-1] = pb[0]/((1./i)+pb[1]*i) + pb[2] + pb[3]*exp(-pb[4]*i);
      c_coeff[i-1] = pc[0]*i*i   + pc[1]*i  + pc[2] + pc[3]*exp(-pc[4]*i);
      //      G4cout << i << ": " << a_coeff[i-1] << " , " << b_coeff[i-1]
      //	   << " , " << c_coeff[i-1] << G4endl;
    }
}

void MuonElectronGenerator::ReadChris()
{
  // E(GeV)
  // Flux per sr*s*GeV
  // Just a number
  G4float data[] = {
0.05,0.5,0,
0.1,1.171161086,1,
0.2,1.058312605,2,
0.3,0.901913544,3,
0.4,0.777019665,4,
0.5,0.678778965,5,
0.6,0.599923779,6,
0.7,0.535281328,7,
0.8,0.481338101,8,
0.9,0.435659887,9,
1,0.396509759,10,
1.1,0.362612898,11,
1.2,0.333010236,12,
1.3,0.30696477,13,
1.4,0.283899767,14,
1.5,0.263356797,15,
1.6,0.244966486,16,
1.7,0.228427669,17,
1.8,0.213492212,18,
1.9,0.199953762,19,
2,0.187639285,20,
2.1,0.176402586,21,
2.2,0.166119314,22,
2.3,0.156683037,23,
2.4,0.148002141,24,
2.5,0.13999736,25,
2.6,0.132599781,26,
2.7,0.125749219,27,
2.8,0.119392901,28,
2.9,0.113484367,29,
3,0.107982575,30,
3.1,0.102851137,31,
3.2,0.098057697,32,
3.3,0.093573392,33,
3.4,0.089372405,34,
3.5,0.085431582,35,
3.6,0.081730102,36,
3.7,0.078249198,37,
3.8,0.074971918,38,
3.9,0.071882913,39,
4,0.068968259,40,
4.1,0.066215296,41,
4.2,0.063612497,42,
4.3,0.061149343,43,
4.4,0.058816219,44,
4.5,0.056604324,45,
4.6,0.054505588,46,
4.7,0.052512597,47,
5,0.047102572,48,
10,0.011861388,49,
15,0.004599915,50,
20,0.002202706,51,
25,0.001198952,52,
30,0.000711986,53,
35,0.000450526,54,
40,0.00029927,55,
45,0.000206596,56,
50,0.000147159,57,
55,0.000107587,58,
60,8.04074E-05,59,
65,6.12407E-05,60,
70,4.74142E-05,61,
75,3.72414E-05,62};
  G4int i;
  // N_EC = 63; // (muonElectronGenerator.hh)
  // G4float EC[N_EC];  // (muonElectronGenerator.hh)
  // G4float ELNC[N_EC];  // (muonElectronGenerator.hh)
  G4float sint=0.;
  EC[0] = data[0];
  for (i=1; i<N_EC-1; i++) {
    EC[i] = sqrt(data[3*i]*data[3*(i-1)]);
    sint += data[3*i - 2]*(EC[i]-EC[i-1]);
    ELNC[i] = sint;
  } EC[N_EC-1] = data[3*(N_EC-1)];
    sint += data[3*N_EC-2]*(EC[N_EC-1]-EC[N_EC-2]);
    ELNC[N_EC-1] = sint;
  for (i=0; i<N_EC; i++) { ELNC[i] /= 1.001*sint;
    G4cout << EC[i] << " " << ELNC[i] << G4endl;
  }
  return;
}

void MuonElectronGenerator::ReadHaruo()
{
  // E(GeV)  - momentum GeV/c at Jokisch et al. PRD, 19 (1979) 1368
  // Flux per sr*sq.sm*s
  //
  G4float data[] = {
	1.000000,2.50E-05,
	1.221688,1.11E-04,
	1.513586,8.90E-05,
	1.865962,7.58E-05,
	2.399631,6.07E-05,
	2.963546,5.18E-05,
	3.688603,4.33E-05,
	4.615081,3.64E-05,
	5.793331,3.05E-05,
	7.263632,2.50E-05,
	9.106549,2.02E-05,
	11.424666,1.57E-05,
	14.336636,1.19E-05,
	17.999116,8.73E-06,
	22.61224,6.10E-06,
	28.406165,4.18E-06,
	35.701106,2.75E-06,
	44.88733,1.77E-06,
	56.44516,1.09E-06,
	70.995018,6.62E-07,
	89.307424,3.82E-07,
	112.358867,2.24E-07,
	141.358581,1.25E-07,
	177.870811,6.69E-08,
	223.826204,3.76E-08,
	281.675582,1.88E-08,
	354.489972,9.73E-09,
	446.150674,5.46E-09,
	625.86087,1.77E-09,
	990.743091,3.99E-10};

  G4int i;
  // N_E75 = 30; // (muonElectronGenerator.hh)
  // G4float E75[N_E75];  // (muonElectronGenerator.hh)
  // G4float ELN75[N_E75];  // (muonElectronGenerator.hh)
  G4float sint=0.;
  E75[0] = data[0];
  for (i=1; i<N_E75-1; i++) {
    E75[i] = sqrt(data[2*i]*data[2*(i-1)]);
    sint += data[2*i - 1]*(E75[i]-E75[i-1]);
    ELN75[i] = sint;
  } E75[N_E75-1] = data[2*(N_E75-2)];
    sint += data[2*N_E75-1]*(E75[N_E75-1]-E75[N_E75-2]);
    ELN75[N_E75-1] = sint;
  for (i=0; i<N_E75; i++) { ELN75[i] /= 1.001*sint;
    G4cout << E75[i] << " " << ELN75[i] << G4endl;
  }
  return;
}

void MuonElectronGenerator::ReadJoks()
{
  // E(GeV)  - momentum GeV/c at Jokisch et al. PRD, 19 (1979) 1368
  // Flux per sr*sq.sm*s

  G4float data[] = {
    1.,   2.5E-05,
    1.12, 2.97E-05,
    1.41, 3.11E-05,
    1.76, 3.26E-05,
    2.29, 3.17E-05,
    2.85, 3.07E-05,
    3.57, 2.85E-05,
    4.49, 2.61E-05,
    5.66, 2.33E-05,
    7.12, 2.02E-05,
    8.95, 1.71E-05,
    11.26, 1.37E-05,
    14.17, 1.07E-05,
    17.83, 8.03E-06,
    22.44, 5.71E-06,
    28.23, 3.96E-06,
    35.52, 2.64E-06,
    44.7,  1.71E-06,
    56.25, 1.06E-06,
    70.79, 6.48E-07,
    89.09, 3.75E-07,
    112.13, 2.21E-07,
    141.12, 1.23E-07,
    177.62, 6.63E-08,
    223.56, 3.74E-08,
    281.39, 1.87E-08,
    354.18, 9.69E-09,
    445.81, 5.44E-09,
    625.46, 1.77E-09,
    990.22, 3.98E-10};

  G4int i;
  // N_E75 = 30; // (muonElectronGenerator.hh)
  // G4float E75[N_E75];  // (muonElectronGenerator.hh)
  // G4float ELN75[N_E75];  // (muonElectronGenerator.hh)
  G4float sint=0.;
  E75[0] = data[0];
  for (i=1; i<N_E75-1; i++) {
    E75[i] = sqrt(data[2*i]*data[2*(i-1)]);
    sint += data[2*i - 1]*(E75[i]-E75[i-1]);
    ELN75[i] = sint;
  } E75[N_E75-1] = data[2*(N_E75-2)];
    sint += data[2*N_E75-1]*(E75[N_E75-1]-E75[N_E75-2]);
    ELN75[N_E75-1] = sint;
  for (i=0; i<N_E75; i++) { ELN75[i] /= 1.001*sint;
    G4cout << E75[i] << " " << ELN75[i] << G4endl;
  }
  return;
}

/*int main()
{
  LIST_OUT = 1;

  cout << "Testing cosmic-ray generator" << endl;
  MuonElectronGenerator *MEGen = new MuonElectronGenerator();
  MEGen->GenerateMuon();
  G4float momentum = (MEGen->GetEnergy())*GeV;
  G4double theta = (G4double)MEGen->GetZenithAngle();

  cout << momentum << " " << theta << endl;
  delete MEGen;
  }*/
