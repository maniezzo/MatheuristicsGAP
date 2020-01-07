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
   cout << "Kernel search, initial solution" << endl;
   isVerbose = fVerbose;

   int numRows,numCols,numNZrow;
   numRows = n+m;   // num of constraints
   numCols = n*m;   // num of variables
   numNZrow= n*m;   // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(isVerbose);

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

         // dual variables, assignments (n) and capacity (m)
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
      {  q[j] = GAP->req[i][j];         // requests to the i-th server by the j-th client to reassign
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

