#ifndef BENDERS_H
#define BENDERS_H
#include "GAP.h"
#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>    // std::find
#include "MIPCplex.h"

class Benders
{
   public:
      GeneralizedAssignemnt* GAP;

      Benders(GeneralizedAssignemnt*, int&);
      ~Benders();
      int bendersHeu(int** c, int maxiter, bool fVerbose);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
      bool  isVerbose;

      void Benders::AllocateMaster(MIPCplex* CPX, int nServKept, bool isVerbose);
      double Benders::SolveMaster(MIPCplex* CPX, int, int, vector<int> &, vector<double>, bool isVerbose);
      double Benders::SolveDLGAP(vector<int>, vector<double> &,int,bool);
      void AddBendersCut(MIPCplex* CPX, vector<double> wbar, int nServKept);
};

#endif // BENDERS_H
