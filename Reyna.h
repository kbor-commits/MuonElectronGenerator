//a ROOT macro that randomly samples from the Reyna muon parametrization (hep-ph/0604145)
//to make realistic cosmic ray muon energy and angle distributions
// Matt's Reyna Implementation

#ifndef Reyna_h
#define Reyna_h 1

#include "TF1.h"
#include "TH2.h"

using namespace std;

class Reyna
{
 public:
  Reyna();
  ~Reyna();
  
 public:
  void Generate(double Ysep, double PosZ, double Zdiff, double* theta, double* energy);
    
  
 private:
  TF1* angle_dist;
  TF1* Reyna_energy_dist;
  double this_angle;
  double this_energy;
};

#endif
