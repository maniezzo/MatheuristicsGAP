#ifndef CORRIDOR_H
#define CORRIDOR_H
#include "GAP.h"
#include "cplex.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>    // std::find
#include "MIPCplex.h"

class Corridor
{
   public:
      GeneralizedAssignemnt* GAP;

      Corridor(GeneralizedAssignemnt*, int&);
      ~Corridor();
      int solveByCorridor(int** c, int delta, int maxiter, bool fVerbose);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
      bool  isVerbose;

      void updateModelWithCorridor(MIPCplex* CPX, int* solIter, int delta);
};

#endif // CORRIDOR_H
