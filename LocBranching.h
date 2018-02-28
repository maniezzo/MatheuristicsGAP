#ifndef LOCBRANCH_H
#define LOCBRANCH_H
#include "GAP.h"
#include "cplex.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>    // std::find
#include "MIPCplex.h"

class LocBranching
{
   public:
      GeneralizedAssignemnt* GAP;

      LocBranching(GeneralizedAssignemnt*, int&);
      ~LocBranching();
      int localBranching(int** c, int k, int maxiter, bool fVerbose);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
      bool  isVerbose;

      void addLBcut(MIPCplex* CPX, int* solIter, int k); // adds a local branching cut
};

#endif // LOCBRANCH_H
