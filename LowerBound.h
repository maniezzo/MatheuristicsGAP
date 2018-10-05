#ifndef LOWERBOUND_H
#define LOWERBOUND_H
#include "GAP.h"
#include "MIPCplex.h"

class LowerBound
{
   public:
      GeneralizedAssignemnt* GAP;
      MIPCplex* CPX;

      LowerBound(GeneralizedAssignemnt*, int&);
      ~LowerBound();

      double trivialBound();
      double linearBound();
      double linearBound(int** c, int n, int m, int** req, int* cap);
      double lagrangianDecomposition(int** c, double alpha, double alphastep, double minAlpha, int innerIter, int maxiter);
      double benders();

   private:
      int* subproblem_LD(int**, double*, double*, int, vector<vector<double>>&, 
                         vector<vector<double>>&,
                         vector<vector<short>>&,
                         vector<vector<short>>&);

      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub;
};

#endif // LOWERBOUND_H
