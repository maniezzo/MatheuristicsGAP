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
      int solveByKernel(bool fVerbose);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
      bool  isVerbose;

      int Kernel::genCol(MIPCplex* CPX, vector<double> d);
      void Kernel::addColumn(MIPCplex* CPX, int iserv, int* Ksol, double);
};

#endif // KERNEL_H
