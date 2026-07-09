//container is a class for storing all information that gets passed to a rootfile

#ifndef container_h
#define container_h 1

//includes
#include <string> 
#include <vector>
using namespace std;

class container
{
 public:
  container();    
  ~container();
  
 public:
  void clear(); // clears the container class
  // units in MeV, and meters
  double posx_1, posx_2;
  double posy_1, posy_2;
  double posz_1, posz_2;
  double momentum_1, momentum_2;
  double vecx_1, vecx_2;
  double vecy_1, vecy_2;
  double vecz_1, vecz_2;
  //int inacceptance, outacceptance;
  double rl;
  vector<double> flux_step;
};

#endif

