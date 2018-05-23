#ifndef RINS_H
#define RINS_H
#include "GAP.h"
#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

class Rins
{
   public:
      GeneralizedAssignemnt* GAP;

      Rins(GeneralizedAssignemnt*, int&);
      ~Rins();
      int dive(int**,int, int, int*, bool);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // RINS_H
