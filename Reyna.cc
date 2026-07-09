//a ROOT macro that randomly samples from the Reyna muon parametrization (hep-ph/0604145)
//to make realistic cosmic ray muon energy and angle distributions
// Matt's Reyna Implementation
#include "math.h"
#include "Reyna.h"

Reyna::Reyna():this_angle(0), this_energy(0){;}

Reyna::~Reyna(){;}

void Reyna::Generate(double Ysep, double PosZ, double Zdiff, double* theta, double* momentum)
{
  
  angle_dist = new TF1("angle_dist","[0]*cos(x)*cos(x)",0,acos(-1.0)/2); angle_dist->SetParameter(0,2);//angular distribution to randomly sample from  
  Reyna_energy_dist = new TF1("Reyna_energy_dist","[1]/x^( [2]+[3]*log10([0]*x)+[4]*log10([0]*x)*log10([0]*x) + [5]*log10([0]*x)*log10([0]*x)*log10([0]*x) )",1,100);//Reyna energy distribution
  Reyna_energy_dist->SetParameter(1,0.00253);
  Reyna_energy_dist->SetParameter(2,0.2455);
  Reyna_energy_dist->SetParameter(3,1.288);
  Reyna_energy_dist->SetParameter(4,-0.2555);
  Reyna_energy_dist->SetParameter(5,0.0209);
  
  
  //TH2D *h_angle_energy = new TH2D("h_angle_energy","h_angle_energy",100,0,100,100,0,100);//histogram to check sampled distributions in degrees, GeV/c
  //h_angle_energy->GetXaxis()->SetTitle("Angle from zenith (degrees)"); h_angle_energy->GetYaxis()->SetTitle("Muon momentum (GeV/c)");
  
  this_angle = angle_dist->GetRandom();
  
  while( this_angle > atan(Ysep/PosZ) || this_angle < atan(Ysep/Zdiff) ){
    this_angle = angle_dist->GetRandom();
  }
  
  Reyna_energy_dist->SetParameter(0,cos(this_angle));
  this_energy = pow(cos(this_angle),3)*Reyna_energy_dist->GetRandom();
  
  *theta = this_angle;
  *momentum = this_energy*1000;

  delete angle_dist;
  delete Reyna_energy_dist;

  /*for(int i=0; i<10000; i++)
    {
    if(i%1000==0) cout<<i<<endl;
    
    this_angle = angle_dist->GetRandom();
    
    if(this_angle*180/acos(-1.0)>30 && this_angle*180/acos(-1.0)<70)
    {
    
    //OK, now I have a random angle from 30 to 70 degrees.  Find an energy from Reyna.
    
    Reyna_energy_dist->SetParameter(0,cos(this_angle));
    
    this_energy = pow(cos(this_angle),3)*Reyna_energy_dist->GetRandom();
    
    h_angle_energy->Fill(this_angle*180/acos(-1.0),this_energy);
    
    }
    }*/

}
