#ifndef VLSN_H
#define VLSN_H
#include "GAP.h"

class VeryLarge
{
   public:
      GeneralizedAssignemnt* GAP;

      VeryLarge(GeneralizedAssignemnt*, int&);
      ~VeryLarge();
      int verylarge(int**,int);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;

};

#endif // VLSN_H
