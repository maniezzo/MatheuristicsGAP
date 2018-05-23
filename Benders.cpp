#include "Benders.h"
#include <ilcplex/cplex.h>

Benders::Benders(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m   = GAP->m;
   n   = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

Benders::~Benders()
{
   //dtor
}

/*
BEND HEURISTIC()
1 identify a master MP(z,y) and an “easy” subproblem SP(x), set t = 0
2 repeat
3 solve (heuristically) master problem MPt . Solution (zt , yt )
4 if x are requested to be integer
5 then solve (heuristically) problem SP(x), solution (zH , xt )
6 solve problem DP(wt ), solution (zd , wt ), and add to MP the
corresponding cut
7 if no more cuts can be added
8 then STOP else set t = t + 1
9 until (end_condition)
*/
int Benders::bendersHeu(int** c, int maxiter, bool fVerbose)
{  int* solIter = new int[n];
   double zMaster;       // Bender's bound
   double zSubpr;        // subproblem cost

   if(GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }
   isVerbose = fVerbose;

   int numRows, numCols, numNZrow, nServKept;
   nServKept = m / 3;            // keeping 1/3 of capacity constraints
   numRows  = 1 + nServKept; 
   numCols  = 1 + n*nServKept;   // only vars in the capacities kept
   numNZrow = 1 + n;             // max num of nonzero values in a row
   MIPCplex* CPXmaster = new MIPCplex(numRows, numCols, numNZrow);
   CPXmaster->GAP = GAP;

   vector<int> xbar(numCols-1); // client assignmed fixed in master problem
   vector<double> wbar(nServKept*n + (m - nServKept)); // subproblem duals

   cout << "Benders. Keeping first " << nServKept << " capacity constraints in master" << endl;

   try
   {  
      zSubpr = DBL_MAX;
      zMaster = 0;
      AllocateMaster(CPXmaster,nServKept, isVerbose);

      while ((zSubpr-zMaster) > 0.01)
      {
         AddBendersCut(CPXmaster,wbar,nServKept);
         zMaster = SolveMaster(CPXmaster, numRows, nServKept, xbar, wbar, isVerbose);
         zSubpr  = SolveDLGAP(xbar, wbar, nServKept, isVerbose);
      }
   }
   catch(std::exception const& e)
   {  cout << "Error: " << e.what() << endl;
      goto lend;
   }

   cout << "Bender's, bound = " << zMaster << endl;

lend:
   CPXmaster->freeMIP();
   delete(CPXmaster);
   return 0;
}

// Allocates the CPLEX memory structures for the master core
void Benders::AllocateMaster(MIPCplex* CPX, int nServKept, bool isVerbose)
{
   int i, j, status, numNZrow, idRow;

   try
   {
      status = 0;
      CPX->env = NULL;
      CPX->lp = NULL;
      CPX->pi = NULL;
      CPX->dj = NULL;
      CPX->x = NULL;
      CPX->slack = NULL;

      cout << "Allocating CPLEX for master" << endl;
      CPX->env = CPXopenCPLEX(&status);
      if (CPX->env == NULL)
      {
         char  errmsg[1024];
         cerr << "Could not open CPLEX environment.\n";
         CPXgeterrorstring(CPX->env, status, errmsg);
         cerr << errmsg;
         goto lend;
      }

      // Turn on output to the screen 
      if (isVerbose)
         status = CPXsetintparam(CPX->env, CPX_PARAM_SCRIND, CPX_ON);
      else
         status = CPXsetintparam(CPX->env, CPX_PARAM_SCRIND, CPX_OFF);
      if (status)
      {
         cerr << "Failure to turn on screen indicator, error " << status << endl;
         goto lend;
      }

      // Turn on data checking 
      status = CPXsetintparam(CPX->env, CPX_PARAM_DATACHECK, CPX_ON);
      if (status)
      {
         cerr << "Failure to turn on data checking, error " << status << endl;
         goto lend;
      }

      // Create the problem. 
      CPX->lp = CPXcreateprob(CPX->env, &status, "GAPinst");
      if (CPX->lp == NULL)
      {
         cerr << "Failed to create LP.\n";
         goto lend;
      }

      CPXchgobjsen(CPX->env, CPX->lp, CPX_MIN);  // Problem is minimization 

                                                 // Create the new columns
      int ij = 0;

      // z column      
      CPX->obj[ij] = 1.0;
      CPX->lb[ij] = 0;
      CPX->ub[ij] = DBL_MAX;
      CPX->ctype[ij] = 'C';  // 'B', 'I','C' to indicate binary, general integer, continuous 
      CPX->lptype[ij] = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
      CPX->colname[ij] = (char *)malloc(sizeof(char) * (11));   // why not 11?
      sprintf(CPX->colname[ij], "%s", "z");
      ij++;

      for (i = 0; i<nServKept; i++)
         for (j = 0; j<GAP->n; j++)
         {
            CPX->obj[ij] = 0.0;
            CPX->lb[ij] = 0;
            CPX->ub[ij] = 1.0;
            CPX->ctype[ij] = 'B'; // 'B', 'I','C' to indicate binary, general integer, continuous 
            CPX->lptype[ij] = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
            CPX->colname[ij] = (char *)malloc(sizeof(char) * (11));   // why not 11?
            sprintf(CPX->colname[ij], "%s%d", "x", ij);
            ij++;
         }

      status = CPXnewcols(CPX->env, CPX->lp, CPX->numCols, CPX->obj, CPX->lb, CPX->ub, NULL, CPX->colname);
      if (status)  goto TERMINATE;

      idRow = 0;

      // z contraint
      numNZrow = 0;  // number of nonzero element in the row to add
      CPX->rmatbeg[0] = 0;
      CPX->sense[0] = 'G';
      CPX->rhs[0] = 0;
      CPX->rowname[idRow] = (char *)malloc(sizeof(char) * (11));   // why not 11?
      sprintf(CPX->rowname[0], "%s%d", "z", idRow);
      CPX->rmatind[numNZrow] = 0;
      CPX->rmatval[numNZrow] = 1;
      numNZrow++;
      for (i = 0; i<nServKept; i++)
         for (j = 0; j<GAP->n; j++)
         {
            CPX->rmatind[numNZrow] = 1 + j + GAP->n*i;
            CPX->rmatval[numNZrow] = -GAP->c[i][j];
            numNZrow++;
         }
      CPX->rmatbeg[1] = numNZrow;
      status = CPXaddrows(CPX->env, CPX->lp, 0, 1, numNZrow, CPX->rhs, CPX->sense, CPX->rmatbeg, CPX->rmatind, CPX->rmatval, NULL, CPX->rowname);
      if (status)
         goto TERMINATE;
      idRow++;

      // Capacity constraints
      for (i = 0; i<nServKept; i++)
      {
         numNZrow = 0;  // number of nonzero element in the row to add
         CPX->rmatbeg[0] = 0;
         CPX->sense[0] = 'L';
         CPX->rhs[0] = GAP->cap[i];
         CPX->rowname[idRow] = (char *)malloc(sizeof(char) * (11));   // why not 11?
         sprintf(CPX->rowname[0], "%s%d", "c", idRow);
         for (j = 0; j<GAP->n; j++)
         {
            CPX->rmatind[numNZrow] = 1 + j + GAP->n*i;
            CPX->rmatval[numNZrow] = GAP->req[i][j];
            numNZrow++;
         }
         CPX->rmatbeg[1] = numNZrow;
         status = CPXaddrows(CPX->env, CPX->lp, 0, 1, numNZrow, CPX->rhs, CPX->sense, CPX->rmatbeg, CPX->rmatind, CPX->rmatval, NULL, CPX->rowname);
         if (status)
            goto TERMINATE;
         idRow++;
      }
   }
   catch (std::exception const& e)
   {
      cout << "[SolveMaster] Error: " << e.what() << endl;
      goto lend;
   }

TERMINATE:
   if (status)
      cout << "[AllocateMaster] ---------------  Error in constraint definition" << endl;

lend:
   return;
}

// Add a Benedr's cut
void Benders::AddBendersCut(MIPCplex* CPX, vector<double> wbar, int nServKept)
{
   int i,j,ij, numNZrow, status;
   int numVarDL = nServKept*n + (m - nServKept);

   // z contraint
   numNZrow = 0;  // number of nonzero element in the row to add
   CPX->rmatbeg[0] = 0;
   CPX->sense[0] = 'G';
   CPX->rhs[0] = 0;
   for(i=0;i<GAP->m-nServKept;i++)
      CPX->rhs[0] += wbar[nServKept*n+i]*GAP->cap[nServKept+i];
   CPX->rowname[0] = (char *)malloc(sizeof(char) * (11));   // why not 11?
   sprintf(CPX->rowname[0], "%s%d", "B", CPXgetnumrows(CPX->env, CPX->lp));
   CPX->rmatind[numNZrow] = 0;
   CPX->rmatval[numNZrow] = 1;
   numNZrow++;
   ij=0;
   for (i = 0; i<nServKept; i++)
      for (j = 0; j<GAP->n; j++)
      {
         CPX->rmatind[numNZrow] = 1 + j + GAP->n*i;
         CPX->rmatval[numNZrow] = -wbar[ij];
         CPX->rhs[0] += wbar[ij];
         numNZrow++;
         ij++;
      }
   CPX->rmatbeg[1] = numNZrow;
   status = CPXaddrows(CPX->env, CPX->lp, 0, 1, numNZrow, CPX->rhs, CPX->sense, CPX->rmatbeg, CPX->rmatind, CPX->rmatval, NULL, CPX->rowname);
   if (status)
      cout << "[AddBendersCut] Error.";
}

double Benders::SolveMaster(MIPCplex* CPX, int numRows, int nServKept, vector<int> & x1, vector<double> wbar, bool isVerbose)
{  int i, j, status;
   double zMaster=DBL_MIN;

   try
   {
      status = 0;

      // Write a copy of the problem to a file. 
      if (isVerbose)
      {
         status = CPXwriteprob(CPX->env, CPX->lp, "GAPmaster.lp", NULL);
         if (status)
         {  cerr << "[Error "<<status<<"] Failed to write problem to disk.\n";
            goto lend;
         }
      }

      status = CPX->solveMIP(true, false); // integer solution
      if (!status)
      {
         zMaster = CPX->objval;
         for (i = 0; i < nServKept; ++i)  
            for (j = 0; j < n; ++j)
               x1[i*n + j] = (int) CPX->x[i*n + j +1]; // +1 remoeves the z
      }
   }
   catch (std::exception const& e)
   {
      cout << "[SolveMaster] Error: " << e.what() << endl;
      goto lend;
   }

lend:
   return zMaster;
}

// Duale del ril. lineare
double Benders::SolveDLGAP(vector<int> xbar, vector<double> & wbar, int nServKept, bool isVerbose)
{
   int j, numVarDL, status;
   double zDL = -1;

   int numRows, numCols, numNZrow;
   numRows  = n*(m-nServKept);   // num of constraints
   numCols  = numVarDL = nServKept*n + (m - nServKept);   // num of variables
   numNZrow = 2;   // max num of nonzero values in a row

   try
   {
      MIPCplex* CPX = new MIPCplex(numRows, numCols, numNZrow);
      CPX->GAP = GAP;
      CPX->allocateDual(numRows, numCols, xbar, isVerbose);

      status = CPX->solveMIP(false, false); // linear

      cout << "\nLP Solution status = " << status << endl;
      if(!status)
      {  zDL = CPX->objval;
         cout << "LP Solution value  = " << CPX->objval << endl;
         for (j = 0; j < numVarDL; ++j)
            wbar[j] = CPX->x[j];

         if (n*m < 101)                       // print LP if instance is small
         {  cout << " - Solution: " << endl;
            for (j = 0; j < numVarDL; ++j)
               cout << CPX->x[j] << "\t";
            cout << endl;
         }
      }

      CPX->freeMIP();
      delete(CPX);
   }
   catch (std::exception const& e)
   {
      cout << "[SolveDLGAP]: " << e.what() << endl;
   }

   return(zDL);
}