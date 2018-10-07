#ifndef CPLEX_H
#define CPLEX_H
#include "GAP.h"
#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

struct loginfo 
{
   double lastincumbent;
   int*   lastsol;
   int    lastlog;
   int    numcols;
};
typedef struct loginfo LOGINFO, *LOGINFOptr;

class MIPCplex
{
   public:
      GeneralizedAssignemnt* GAP;

      MIPCplex();
      MIPCplex(int, int, int);
      ~MIPCplex();

      int numRows;
      int numCols;
      int numNZ;

      int      status;
      double*  obj;
      double*  lb;
      double*  ub;
      int*     rmatbeg;
      int*     rmatind;
      double*  rmatval;
      double*  rhs;
      char*    sense;
      char*    ctype;
      char*    lptype;
      char**   colname;
      char**   rowname;
      double*  x;         // cplex solution vector
      double   objval;
      double*  pi;
      double*  slack;
      double*  dj;  

      int allocateMIP(bool isVerbose);
      int allocateMIP(int** c, int n, int m, int** req, int* cap, bool isVerbose);
      int allocateDual(int, int, vector<int> xbar, bool isVerbose);
      int freeMIP();
      int solveMIP(bool fMIPint, bool fVerbose);
      CPXENVptr env;
      CPXLPptr  lp;

   private:
      int    solstat;
      int    *xbest;

      //int populatebycolumn (CPXENVptr env, CPXLPptr lp);
      int populatebyrow (CPXENVptr, CPXLPptr);
      int populatebyrow (CPXENVptr env, CPXLPptr lp, int** c, int n, int m, int** req, int* cap);
      int populateDual(CPXENVptr env, CPXLPptr lp, vector<int>, int, int);
};

#endif // CPLEX_H
