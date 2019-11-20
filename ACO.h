#ifndef ACO_H
#define ACO_H
#include "GAP.h"
#include "MIPCplex.h"
#include <ilcplex/cplex.h>

class ACO
{
   public:
      GeneralizedAssignemnt* GAP;

      ACO(GeneralizedAssignemnt*, int&);
      ~ACO();
      int antColony(int**, int maxiter, int numpop, double alpha);

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;

      vector<vector <double>> tau,eta;

      // private functions
      int montecarlo(vector<double>& v);
};

#endif // ACO_H
