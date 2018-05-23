#ifndef VLSN_H
#define VLSN_H
#include "GAP.h"
#include "MIPCplex.h"


class VeryLarge
{
   public:
      GeneralizedAssignemnt* GAP;

      VeryLarge(GeneralizedAssignemnt*, int&);
      ~VeryLarge();
      int verylarge(int** c, int delta, int maxiter, bool fVerbose);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;

      void fixVariables(int** c, MIPCplex* CPX, int* solIter, int k);
};

#endif // VLSN_H
