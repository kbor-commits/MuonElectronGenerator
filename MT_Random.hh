/*
   Initializer:
   MT_Random.h
 
      Header file for a random number generator. Function 
      genrand_real returns a random number on [0,1] interval. 
 
      Routine init_random initializes generator on all processors 
      and must be called before drawing random numbers (see 
      MT_Random.cpp for explanation on how's that done).

                      Zarija Lukic, February 2009
                          zarija@lanl.gov
 */

#ifndef MT_Random_Header_Included
#define MT_Random_Header_Included

class MT_Random {
   public: 
//      MT_Random();
//      ~MT_Random();

      void init_random(unsigned long, int);
      double genrand_real();

   private:
      void init_genrand(unsigned long);

};

#endif
