#include "Corridor.h"

Corridor::Corridor(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

Corridor::~Corridor()
{
   //dtor
}

int Corridor::solveByCorridor(int** c, int delta, int maxiter, bool fVerbose)
{  int i,j,iter,status;
   int zubIter;
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
   zubIter = zub;
   for(j=0;j<n;j++) solIter[j] = sol[j];
   cout << "Solution: "; for(j=0;j<n;j++) cout << solIter[j] << " "; cout << endl;
   cout << "Corridor method. Starting from z=" << zubIter << endl; 

   try
   {  
      for (iter=0;iter<maxiter;iter++)
      {
         updateModelWithCorridor(CPX,solIter,delta); // defines the corridor of free variables

         status = CPX->solveMIP(true,isVerbose); // integer solution, only corridor variables are free
         if ( !status )
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
               cout << "[solveByCorridor] No feasible solution at this iteration" << endl;
            else
            {  cout << "[Corridor] iter "<< iter <<" zubIter "<< zubIter << endl;
               cout << "Solution: "; for(j=0;j<n;j++) cout << solIter[j] << " "; cout << endl;

               if(zubIter < zub)
               {  zub = zubIter;
                  for(i=0;i<n;i++) solbest[i]=solIter[i];
                  cout << "[Corridor] ************** new zub. iter "<< iter << " zubIter " << zubIter << endl;
               }
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

// this frees the  variables in the corridor
void Corridor::updateModelWithCorridor(MIPCplex* CPX, int* solIter, int delta)
{  int i,j,isol,rand,status;
   bool isInCorridor;
   vector<int> lstClient;
   vector<int> lstCorridor;

   for(i=0;i<n;i++)
      lstClient.push_back(i);

   delta = min(delta,n);
   for(i=0;i<delta;i++)
   {  rand = std::rand() % lstClient.size();     // lists could change size at every iteration
      lstCorridor.push_back(lstClient[rand]);    // corridor contains clients
      lstClient.erase(lstClient.begin()+rand);   // remove client from choosable set
   }

   int cnt = n*m;
   int* indices = new int[cnt];    // indices of the columns corresponding to the variables for which bounds are to be changed
   char*   lu   = new char[cnt];   // whether the corresponding entry in the array bd specifies the lower or upper bound on column indices[j]
   double* bd   = new double[cnt]; // new values of the lower or upper bounds of the variables present in indices

   cout << "Corridor: ";
   for(j=0;j<n;j++)
   {
      // is the client in the corridor?
      isInCorridor = false;
      if(std::find(lstCorridor.begin(), lstCorridor.end(), j)!=lstCorridor.end())
      {  isInCorridor = true;
         cout << j << " ";
      }

      for(i=0;i<m;i++)
      {  if(isInCorridor)              // set it free
         {  
            if(CPX->lb[i*n+j] == 0.0)  // lb OK, only ub could be wrong
            {  CPX->ub[i*n+j] = 1.0; 
               bd[i*n+j] = 1.0;
               lu[i*n+j] = 'U';        // bd[j] is an upper bound
            }
            else                       // lb wrong, ub should be OK
            {  CPX->lb[i*n+j] = 0.0; 
               bd[i*n+j] = 0.0;
               lu[i*n+j] = 'L';        // bd[j] is an upper bound
            }
         }
         else                          // not in corridor
         {  isol = solIter[j];
            if(i==isol)
            {  CPX->lb[i*n+j] = 1.0;   // fixing in the solution
               bd[i*n+j] = 1.0;
            }
            else
            {  CPX->lb[i*n+j] = 0.0;   // forbidding in the solution
               bd[i*n+j] = 0.0;
            }
            lu[i*n+j] = 'B';           // bd[j] is the lower and upper bound
         }
         indices[i*n+j] = i*n+j;
      }
   }
   cout << endl;

   status = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
   free(indices);
   free(lu);
   free(bd);

   return;
}


