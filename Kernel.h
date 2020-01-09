#ifndef KERNEL_H
#define KERNEL_H
#include "GAP.h"
#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>    // std::find
#include "MIPCplex.h"

class Kernel
{
   public:
      GeneralizedAssignemnt* GAP;

      Kernel(GeneralizedAssignemnt*, int&);
      ~Kernel();
      int solveByKernel(bool fVerbose, int numBuckets);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
      bool  isVerbose;

      void updateModelWithKernel(MIPCplex* CPX, vector<bool> kernel);
};

#endif // KERNEL_H
