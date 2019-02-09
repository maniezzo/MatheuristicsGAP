#include "Rins.h"
#include "MIPCplex.h"
#include <ilcplex/cplex.h>

Rins::Rins(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

Rins::~Rins()
{
   //dtor
}

int Rins::dive(int** c, int maxNodes, int z, int* sol, bool fVerbose)
{  int i,j,iter,status;
   int zubIter;
   int* solIter = new int[n];
   bool fChanged;

   int numRows,numCols,numNZrow;
   numRows = n+m;   // num of constraints
   numCols = n*m;   // num of variables
   numNZrow= n*m;   // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(true); // verbose output
   CPXsetintparam(CPX->env, CPX_PARAM_RINSHEUR, -1);  // disable rins in cplex (not needed)
   zubIter = z;
   for(j=0;j<n;j++) solIter[j] = sol[j];   // local search around solution sol
   cout << "RINS. Starting from z=" << zubIter << endl; 

   try
   {  
      iter = 0;

      fChanged = false;

      status = CPX->solveMIP(false,false);   // LP of whole instance
      if ( !status )
      {
         if(fVerbose && n*m < 101)                       // print LP if instance is small
         {  cout << " - Solution: " << endl; 
            for (i = 0; i < m; ++i)
            {  cout << "   " << i << ": ";
               for (j = 0; j < n; ++j)
                  cout << CPX->x[i*n+j] << "\t";
               cout << endl;
            }
         }

         // heuristic variable fixing
         zubIter = 0;
         for(j=0;j<n;j++)
         {  for(i=0;i<m;i++)
            {  // linear primal sufficiently 1, not yet fixed, compatible with seed -> fix
               if(CPX->x[i*n+j] > 0.99 && CPX->lb[i*n+j] < 0.99 && solIter[j] == i)
               {  cout << "[RINS] Set variable " << i*n+j << ": client " << j << " to " << i << endl;
                  CPX->lb[i*n+j] = 1;
                  int cnt = 1;
                  int* indices = new int[cnt];     // which var
                  indices[0] = i*n+j;
                  char* lu = new char[cnt];        // lower limit
                  lu[0] = 'L';
                  double * bd = new double[cnt];   // new bound
                  bd[0] = 1.0;
                  status = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
                  free(indices);
                  free(lu);
                  free(bd);
                  fChanged = true;
               }

               // linear primal sufficiently 0, not yet fixed, compatible with seed -> fix
               if(CPX->x[i*n+j] < 0.01 && CPX->ub[i*n+j] > 0.01 && solIter[j] != i)
               {  cout << "[RINS] Set variable " << i*n+j << ": client " << j << " to be different from " << i << endl;
                  CPX->ub[i*n+j] = 0;
                  int cnt = 1;
                  int* indices = new int[cnt];     // which var
                  indices[0] = i*n+j;
                  char* lu = new char[cnt];        // upper limit
                  lu[0] = 'U';
                  double * bd = new double[cnt];   // new bound
                  bd[0] = 0.0;
                  status = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
                  free(indices);
                  free(lu);
                  free(bd);
                  fChanged = true;
               }
               //else if(CPX->lb[i*n+j] < 0.99   && 
               //        CPX->x[i*n+j]  > (1-GAP->EPS)  &&
               //        solIter[j]     != i) // this wouldn't be in rins, it helps with infeasibilities
               //{  cout << "Changed " << j << " to " << i << endl;
               //   solIter[j] = i;           // update of incumbent solution, if LP says so
               //   fChanged = true;
               //   break;
               //}
            }

            if(solIter[j] >= 0 && zubIter != z)
               zubIter += c[solIter[j]][j];    // when completing partial solutions
            else
               zubIter = z;
         }

         // subMIP, go for MIP to try to complete LP solution (node bound maxNode)
         cout << "Solving sub-MIP ..." << endl;
         CPXsetintparam(CPX->env, CPX_PARAM_NODELIM, maxNodes);
         status = CPX->solveMIP(true,false);
         if(status==0)
         {  cout << "SubMIP, z=" << CPX->objval << endl;
            if(CPX->objval < zubIter)
            {  zubIter = (int) floor( CPX->objval + GAP->EPS ); // the integer cost
               for(j=0;j<n;j++)
                  for(i=0;i<m;i++)
                     if(CPX->x[i*n+j] > 0.99)
                        solIter[j] = i;
            }
         }

         if(abs(zubIter - GAP->checkSol(solIter)) > GAP->EPS)
            cout << "[solveGAPbyRins] No feasible solution at this try" << endl;
         else
         {  cout << "[Rins] iter "<< iter <<" zubIter "<< zubIter << endl;
            for(j=0;j<n;j++) sol[j]=solIter[j];

            if(zubIter < zub)
            {  zub = zubIter;
               for(j=0;j<n;j++) solbest[j]=solIter[j];
               cout << "[LB] ************** new zub. iter "<< iter << " zubIter " << zubIter << endl;
            }
         }
      }
      else
         goto lend;
   }
   catch(std::exception const& e)
   {  cout << "Error: " << e.what() << endl;
      goto lend;
   }

lend:
   CPX->freeMIP();
   delete(CPX);
   delete(solIter);
   return zubIter;
}

