#include "container.h"

container::container():posx_1(0), posx_2(0),
                       posy_1(0), posy_2(0),
                       posz_1(0), posz_2(0),
                       momentum_1(0), momentum_2(0), 
                       vecx_1(0), vecx_2(0), 
                       vecy_1(0), vecy_2(0), 
                       vecz_1(0), vecz_2(0),
		       //inacceptance(0), outacceptance(0),
		       rl(0)
{ flux_step.clear();}

container::~container(){;}

void container::clear(){
  posx_1=999.99;
  posx_2=999.99;
  posy_1=999.99;
  posy_2=999.99;
  posz_1=999.99;
  posz_2=999.99;
  momentum_1=-1;
  momentum_2=-1;
  vecx_1=999.99;
  vecx_2=999.99;
  vecy_1=999.99;
  vecy_2=999.99;
  vecz_1=999.99;
  vecz_2=999.99;
  rl = -1;
  //inacceptance=0;
  //outacceptance=0; 
  flux_step.clear();
}
