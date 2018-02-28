#include "Kernel.h"
//#include <vld.h>      /* visual memory leak for VS */

Kernel::Kernel(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

Kernel::~Kernel()
{
   //dtor
}

int Kernel::solveByKernel(bool fVerbose)
{  int i,status,numCol,iter;
   double zlb;
   vector<double> d; // master dual vars

   if(GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }
   isVerbose = fVerbose;

   MIPCplex* CPX = new MIPCplex();
   CPX->GAP = GAP;

   cout << "Dantzig-Wolfe, initial solution" << endl;
   allocateInitialSol(CPX, sol);

   try
   {  iter = 0;
      do
      {
         status = CPX->solveMIP(false,false); // LP solution
         if ( status )
         {  cout << "[solveByKernel] Error" << endl;
            goto lend;
         }
         // reads the solution
         zlb = CPX->objval;
         cout << "iter "<< iter << " zlb " << zlb << endl;
         d.clear();
         for(i=0;i<n+m;i++)
            d.push_back(CPX->pi[i]);

         numCol = genCol(CPX,d);
         iter++;
      } while (numCol > 0 && iter < 100);
   }
   catch(std::exception const& e)
   {  cout << "[solveByKernel] Error: " << e.what() << endl;
      goto lend;
   }

lend:
   CPX->freeMIP();
   delete(CPX);
   return 0;
}

// this frees the  variables in the corridor
void Kernel::allocateInitialSol(MIPCplex* CPX, int* solInit)
{
   int i, j, status;
   double kcost;

   try
   {
      status = 0;
      CPX->env = NULL;
      CPX->lp = NULL;
      CPX->pi = NULL;
      CPX->dj = NULL;
      CPX->x  = NULL;
      CPX->slack = NULL;

      cout << "Allocating CPLEX for master" << endl;
      CPX->env = CPXopenCPLEX(&status);
      if (CPX->env == NULL)
      {  char  errmsg[1024];
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
      {  cerr << "Failure to turn on screen indicator, error " << status << endl;
         goto lend;
      }

      // Turn on data checking 
      status = CPXsetintparam(CPX->env, CPX_PARAM_DATACHECK, CPX_ON);
      if (status)
      {  cerr << "Failure to turn on data checking, error " << status << endl;
         goto lend;
      }

      // Create the problem. 
      CPX->lp = CPXcreateprob(CPX->env, &status, "GAP-SPP");
      if (CPX->lp == NULL)
      {  cerr << "Failed to create LP.\n";
         goto lend;
      }

      CPXchgobjsen(CPX->env, CPX->lp, CPX_MIN);  // Problem is minimization 

      // Create the rows
      int numRows = n+m;
      double*  rhs = (double *)   malloc(numRows * sizeof(double));
      char*    sense = (char *)   malloc(numRows * sizeof(char));
      char**   rowname = (char **)malloc(numRows * sizeof(char*));

      // assignment constraints rows
      for(i=0;i<n;i++)
      {
         rowname[i] = (char *)malloc(sizeof(char) * (7));   // why not 7?
         sprintf(rowname[i], "%s%d", "a", i);
         sense[i] = 'E';
         rhs[i] = 1.0;
      }

      // server constraints rows
      for (i = 0; i<m; i++)
      {
         rowname[n+i] = (char *)malloc(sizeof(char) * (7));   // why not 7?
         sprintf(rowname[n+i], "%s%d", "s", i);
         sense[n + i] = 'E';
         rhs[n + i] = 1.0;
      }
      status = CPXnewrows(CPX->env, CPX->lp, numRows, rhs, sense, NULL, rowname);
      if (status) goto TERMINATE;
      for (i = 0; i<numRows; i++)
         free(rowname[i]);
      free(rowname); rowname = NULL;
      free(sense);   sense = NULL;
      free(rhs);     rhs = NULL;

      // ------------------ Add the columns
      int numCols = m;
      int numNZcol= m+n;   // each client and the iserver knapsacks
      int k = 0;

      double* obj     = (double *)malloc(numCols * sizeof(double));
      int*    cmatbeg = (int    *)malloc(numCols * sizeof(int));
      int*    cmatind = (int    *)malloc(numNZcol * sizeof(int));
      double* cmatval = (double *)malloc(numNZcol * sizeof(double));
      double* lb = (double *)malloc(numCols * sizeof(double));
      double* ub = (double *)malloc(numCols * sizeof(double));
      char** colname = (char **)malloc(numCols * sizeof(char*));

      numCols = 0;   // counting them
      numNZcol = 0;

      for (i = 0; i<GAP->m; i++)
      {
         cmatbeg[i] = k;
         kcost = 0;
         for(j=0;j<n;j++)
            if(sol[j]==i)
            {
               kcost += GAP->c[i][j];
               cmatval[k] = 1;         // assignments
               cmatind[k] = j;
               k++;
            }
         cmatval[k] = 1;               // server
         cmatind[k] = n+i;
         k++;

         obj[numCols] = kcost;
         lb[numCols]  = 0;
         ub[numCols]  = 1.0;
         //ctype[numCols]   = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
         //lptype[numCols]  = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
         colname[numCols] = (char *)malloc(sizeof(char) * (6));   // why not 6?
         sprintf(colname[numCols], "%s%d", "x", i);

         numCols++;
      }
      numNZcol = k;

      status = CPXaddcols(CPX->env, CPX->lp, numCols, numNZcol, obj, cmatbeg, cmatind, cmatval, lb, ub, colname);

      // Free up allocated memory 
      if (obj != NULL)     free(obj);     obj = NULL;
      if (cmatbeg != NULL) free(cmatbeg); cmatbeg = NULL;
      if (cmatind != NULL) free(cmatind); cmatind = NULL;
      if (cmatval != NULL) free(cmatval); cmatval = NULL;
      free(lb); lb = NULL;
      free(ub); ub = NULL;
      for (i = 0; i<numCols; i++)
         free(colname[i]);
      free(colname);       colname = NULL;

      if (status)  goto TERMINATE;
   }
   catch (std::exception const& e)
   {
      cout << "[allocateInitialSol] Error: " << e.what() << endl;
      goto lend;
   }

TERMINATE:
   if (status)
      cout << "[allocateInitialSol] ---------------  Error in constraint definition" << endl;

lend:
   return;
}

// generates negative red. cost columns
int Kernel::genCol(MIPCplex* CPX, vector<double> d)
{  int i,j,numCol;
   double sum,cost;
   int*    q    = new int[n];     // knapsack requests
   double* val  = new double[n];  // knapsack profits
   int*    Ksol = new int[n];     // knapsack solution

   numCol = 0;
   for(i=0;i<m;i++)
   {

      for (j = 0; j < n; j++)
      {  q[j] = GAP->req[i][j];         // requests to the i-th wh by the j-th elem to reassign
         Ksol[j] = 0;
         val[j] = -1*(GAP->c[i][j]-d[j]);
      }

      KDynRecur(n, GAP->cap[i], q, val, Ksol); // solve the knapsack
      sum = cost = 0;
      for (j = 0; j<n; j++)
         if (Ksol[j] > 0)
         {
            sum += -1*val[j];
            cost += GAP->c[i][j];
         }

      if(d[n+i] > sum)
      {  cout << "New column for server "<< i << " red. cost " << sum - d[n + i] << endl;
         addColumn(CPX,i,Ksol,cost);
         numCol++;
      }
   }

lend:
   if (q != NULL)    delete(q);    q = NULL;
   if (val != NULL)  delete(val);  val = NULL;
   if (Ksol != NULL) delete(Ksol); Ksol = NULL;
   return numCol;
}

void Kernel::addColumn(MIPCplex* CPX, int iserv, int* Ksol, double c)
{  int i,j,status=-1, numCols=1, numNZcol;

   double* obj     = (double *) malloc(numCols * sizeof(double));
   int*    cmatbeg = (int    *) malloc(numCols * sizeof(int));
   int*    cmatind = (int    *) malloc((GAP->n+1) * sizeof(int));
   double* cmatval = (double *) malloc((GAP->n+1) * sizeof(double));
   double* lb      = (double *) malloc(numCols * sizeof(double));
   double* ub      = (double *) malloc(numCols * sizeof(double));
   char**  colname = (char  **) malloc(numCols * sizeof(char*));

   colname[0] = (char *)malloc(sizeof(char) * (7));   // why not 7?
   int idCol = CPXgetnumcols(CPX->env,CPX->lp);
   sprintf(colname[0], "%s%d", "x", idCol);

   if (obj == NULL || cmatbeg == NULL || cmatind == NULL || cmatval == NULL) 
   {  status = -2;
      goto TERMINATE;
   }
   int nRows = CPXgetnumrows(CPX->env,CPX->lp);

   numNZcol = 0;
   cmatbeg[0] = 0;
   obj[0] = c;
   lb[0] = 0;
   ub[0] = 1;
   for(j=0;j<n;j++)
      if(Ksol[j]>0)
      {
         cmatind[numNZcol] = j;
         cmatval[numNZcol] = 1.0;
         numNZcol++;
      }
   // server constraint
   cmatind[numNZcol] = n+iserv;
   cmatval[numNZcol] = 1.0;
   numNZcol++;
   //cmatbeg[1]= numNZcol;

   status = CPXaddcols(CPX->env, CPX->lp, numCols, numNZcol, obj, cmatbeg, cmatind, cmatval, lb, ub, colname);
   if (status) 
   {
      cout << "[addColumn] Error, status:" << status << endl;
      goto TERMINATE;
   }

TERMINATE:
   // Free up allocated memory 
   if (obj != NULL) free(obj);     obj = NULL;
   if (cmatbeg != NULL) free(cmatbeg); cmatbeg = NULL;
   if (cmatind != NULL) free(cmatind); cmatind = NULL;
   if (cmatval != NULL) free(cmatval); cmatval = NULL;
   free(lb); lb = NULL;
   free(ub); ub = NULL;
   for(i=0;i<numCols;i++)
      free(colname[i]);
   free(colname); colname = NULL;
   return;
}

