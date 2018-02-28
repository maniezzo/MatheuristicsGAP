#ifndef EJECTION_H
#define EJECTION_H
#include "GAP.h"
#include "Rins.h"
#include "MIPCplex.h"

class Ejection
{
   public:
      GeneralizedAssignemnt* GAP;

      Ejection(GeneralizedAssignemnt*, int&);
      ~Ejection();
      int ejectionChain(int**, int);

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;

      // private functions
      int* generateOneSol(int*);
      int* feasibilityPump(int,double*,double);
      int setInitialBasis(MIPCplex*,double* x, int* cstat, int* rstat, double);
};

#endif // EJECTION_H
