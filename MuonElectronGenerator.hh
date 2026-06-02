//
//
// $Id: MuonElectronGenerator.hh, ver.1.4 08/29/2012 $
// 				  ver.1.5 04/26/2013
// --------------------------------------------------------------
//

#ifndef MuonElectronGenerator_h
#define MuonElectronGenerator_h 1
#include <string>
#include <fstream>
#include "MT_Random.hh"
using namespace std;

//#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

#ifndef Pi
const double Pi=3.1415927;
#endif

#ifndef Deg_to_Rad
const double Deg_to_Rad=0.0174533;      // Pi/180.
#endif

#ifndef N_E
const int N_E = 42;
#endif

#ifndef N_E75
const int N_E75 = 30;
#endif

#ifndef N_EC
const int N_EC = 63;
#endif

const int N_ZANG = 900;   // number of zenith angle bins
const int ZANG_MAX = 90;  // maximal zenith angle
const int N_ANGLES = 90;
const int N_ENGLES = 30;
const int N_ENERGY = 10000;
//const int N_ENERGY = 300;
const G4int NLog = 300;

#ifndef N_theta
const int N_theta = 20;
#endif

class G4ParticleDefinition;

class MuonElectronGenerator
{
  public:
    MuonElectronGenerator();
    MuonElectronGenerator(int);
    MuonElectronGenerator(MT_Random*);
    MuonElectronGenerator(MT_Random*, string recipe_selector);

    virtual ~MuonElectronGenerator();

  //    virtual void GeneratePrimaries(G4Event*);

  private:

  virtual void SetVariables();
  virtual void SetVariables(string);
  virtual void SetReyna();
  virtual G4int ReadVariables(G4String, G4String);
  virtual G4int ReadReyna(G4String);

  virtual G4float ThetaReyna(const G4int);

  virtual void CalcVariables();
  virtual G4int CalcElVariables();

  virtual G4float ChooseEnergy(G4float);
  virtual G4float ChooseJokEnergy(const G4int);
  virtual G4float ChooseJokEnergy(const G4float theta);
  virtual G4float ChooseChrisEnergy(const G4float theta);
  virtual G4float ChooseEnergy(G4float* eln);
  G4float ChooseElEnergy();

  virtual G4float ChooseTheta();
  virtual G4float ChooseTheta(G4float* tnorm);
public:
  virtual G4float ChooseTheReyna();
private:
  virtual void CalculateCoeff();
  virtual void ReadJoks();
  virtual void ReadHaruo();
  virtual void ReadChris();
  virtual void ReynaVariables();

  G4ParticleDefinition* muel;
  G4float energy_mu;
  G4float theta_mu;
  ofstream oof;
  ifstream iif;
  string recipe;

  MT_Random RNGenerator;
  G4int Nel_bins;

  G4float Sigma[N_ANGLES][N_ENERGY], SigmaE[N_ENERGY], SigmaT[N_ANGLES];
  G4float Enorm[N_ANGLES][N_ENERGY], Tnorm[N_ANGLES];
  G4float TElnorm[N_ANGLES];     /* theta normalization for electrons */
  G4double theReyNorm[N_ZANG];   /* theta normalization for Reyna */

  G4float a_coeff[N_ANGLES],  b_coeff[N_ANGLES],  c_coeff[N_ANGLES];
  G4float theta[N_ANGLES];
  G4float E[N_ENERGY];

  G4float Th[N_ENGLES];        /* angle bins for energy spectrum */
  G4float ELog[NLog];          /* log energy bins */
  G4float E75[N_E75];          /* energy bins for 75 deg */
  G4float EElLog[NLog];        /* log energy bins for electrons */
  G4float ELN[N_ENGLES][NLog]; /* norm.cumulative energy spectrum  for muons */
  G4float ELN75[N_E75];        /* norm.cumulative energy spectrum for 75deg muons */
  G4float EELN[NLog]; /* normalized cumulative energy spectrum for electrons */

  G4float EC[N_EC];  // energy bins for Chris's spectrum
  G4float ELNC[N_EC];  // norm.cumulative energy spectrum for Los Alamos (Chris's version)

  public:
    inline void SetType(G4ParticleDefinition* particle) { muel = particle; }
    inline G4ParticleDefinition* GetType() const { return muel; }
    inline void SetEnergy(G4double val) { energy_mu = val; }
    inline G4double GetEnergy() const { return energy_mu; }
    inline void SetZenithAngle(G4double val) { theta_mu = val; }
    inline G4double GetZenithAngle() const { return theta_mu; }
    inline string GetRecipe() const { return recipe; }
    void SetRecipe(string);
    void RunMuonElectronGen();
    void GenerateMuon75();
    void GenerateReynaMuon();
    void GenerateMuon(string selector);
    void GenerateMuon(string selector, G4float theta);
    void GenerateMuon();
    void GenerateCosmic(string selector);
    void GenerateCosmic(string selector, G4float theta);
    void GenerateCosmic();
};

#endif


