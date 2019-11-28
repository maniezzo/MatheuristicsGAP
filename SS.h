#ifndef SS_H
#define SS_H
#include "GAP.h"
#include "MIPCplex.h"
#include <ilcplex/cplex.h>

class ScatterSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      ScatterSearch(GeneralizedAssignemnt*, int&);
      ~ScatterSearch();
      int go_scatter(int**, int maxiter, int numpop, double alpha);

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

#endif // SS_H
