#include "LocBranching.h"

LocBranching::LocBranching(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

LocBranching::~LocBranching()
{
   //dtor
}

int LocBranching::localBranching(int** c, int k, int maxiter, bool fVerbose)
{  int i,j,iter,status;
   int zubIter,zubOld=0;
   int* solIter = new int[n];

   if(GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }
   isVerbose = fVerbose;

   int numRows,numCols,numNZrow;
   numRows = n+m;   // num of constraints
   numCols = n*m;   // num of variables
   numNZrow= n*m;   // max num of nonzero values in a row

   MIPCplex* CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(isVerbose);
   for(j=0;j<n;j++) 
   {  solIter[j] = sol[j];    // initial solution
      zubOld += c[sol[j]][j]; // initial cost
   }
   if(isVerbose)
   {  cout << "Solution: "; for(j=0;j<n;j++) cout << solIter[j] << " "; cout << endl; }
   cout << "Local branching. Starting from z=" << zubOld << endl; 

   try
   {  
      for (iter=0;iter<maxiter;iter++) // this to ensure a maximum time of completion
      {
         numRows = CPXgetnumrows(CPX->env, CPX->lp);
         addLBcut(CPX, solIter, k);

         status = CPX->solveMIP(true,isVerbose); // integer bounded solution, 
         if ( !status || status == 1001 || status == 1002)
         {
            if(isVerbose && n*m < 201)
            {  cout << " - Solution: " << endl; 
               for (i = 0; i < m; ++i)
               {  cout << "   " << i << ": ";
                  for (j = 0; j < n; ++j)
                     cout << CPX->x[i*n+j] << "\t";
                  cout << endl;
               }
            }
            // reads the solution
            zubIter = 0;
            for(j=0;j<n;j++)
            {  for(i=0;i<m;i++)
                  if(CPX->x[i*n+j] > 0.5)
                  {  solIter[j] = i;
                     break;
                  }
               zubIter += c[solIter[j]][j];
            }

            if(abs(zubIter - GAP->checkSol(solIter)) > GAP->EPS)
               cout << "[localBranching] No feasible solution at this iteration" << endl;
            else
            {  cout << "[localBranching] iter "<< iter <<" zubIter "<< zubIter << endl;
               if(isVerbose)
               {  cout << "Solution: "; for(j=0;j<n;j++) cout << solIter[j] << " "; cout << endl; }

               if(zubIter < zub)
               {  zub = zubIter;
                  for(i=0;i<n;i++) solbest[i]=solIter[i];
                  cout << "[localBranching] ************** new zub. iter "<< iter << " zubIter " << zubIter << endl;
               }

               if(zubIter == zubOld)
               {  cout << "[localBranching] No improvement, exiting " << endl;
                  break;
               }
               zubOld = zubIter;
               // remove the local branching constraint
               numRows = CPXgetnumrows(CPX->env, CPX->lp);
               status = CPXdelrows(CPX->env, CPX->lp, numRows-1, numRows-1);
            }
         }
      }
   }
   catch(std::exception const& e)
   {  cout << "Error: " << e.what() << endl;
      goto lend;
   }

lend:
   CPX->freeMIP();
   delete(CPX);
   delete(solIter);
   return zub;
}

// adds a local branching cut
void LocBranching::addLBcut(MIPCplex* CPX, int* solIter, int k)
{  int i,j,isol;
   int idRow=0,numRows=1,numNZrow=n*m;

   int *    rmatbeg = (int *)    malloc ((numRows+1) * sizeof(int));
   int *    rmatind = (int *)    malloc (numNZrow * sizeof(int));
   double * rmatval = (double *) malloc (numNZrow * sizeof(double));
   double * rhs     = (double *) malloc (numRows * sizeof(double));
   char *   sense   = (char *)   malloc (numRows * sizeof(char));
   char**   rowname = (char **)  malloc (numRows * sizeof(char*));
   rowname[0]       = (char *)   malloc(sizeof(char) * (11));   // there will be only one row

   numNZrow   = 0;  // number of nonzero elements in the row to add
   rmatbeg[0] = 0;     
   sense[0] = 'L';
   rhs[0]   = k;
   sprintf(rowname[0],"%s%d","LB",CPX->numRows);
   for(j=0;j<n;j++)
   {
      isol = solIter[j];
      for(i=0;i<m;i++)
         if(i == isol)
         {  rmatind[numNZrow] = j+GAP->n*i;
            rmatval[numNZrow] = -1; // termine (1-x_j)
            numNZrow++;
            rhs[0] -= 1;
         }
         else
         {  rmatind[numNZrow] = j+GAP->n*i;
            rmatval[numNZrow] = 1;  // termine x_j
            numNZrow++;
         }      
   }

   rmatbeg[1] = numNZrow;
   int status = CPXaddrows (CPX->env, CPX->lp, 0, 1, numNZrow, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
   if ( status )  goto TERMINATE;
   idRow++;

TERMINATE:
   if(rmatbeg != NULL) free(rmatbeg);
   if(rmatind != NULL) free(rmatind);
   if(rmatval != NULL) free(rmatval);
   if(rhs != NULL)     free(rhs);
   if(sense != NULL)   free(sense);
   if (rowname!=NULL)
   {  free(*rowname);
      free(rowname);
   }
   return;
}


