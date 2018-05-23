#include "VLSN.h"

/////////
/////////
///////// NOTE: TO BE UPDATED !!!!!!!!!!!!!!!!!!!!!!!
/////////
/////////

VeryLarge::VeryLarge(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

VeryLarge::~VeryLarge()
{
   //dtor
}

int VeryLarge::verylarge(int** c, int k, int maxiter, bool fVerbose)
{
   int i, j, iter, status;
   int zubIter;
   int* solIter = new int[n];

   if (GAP->n == NULL)
   {
      cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }

   int numRows, numCols, numNZrow;
   numRows  = n + m;   // num of constraints
   numCols  = n * m;   // num of variables
   numNZrow = n * m;   // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows, numCols, numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(fVerbose);
   zubIter = zub;
   for (j = 0; j<n; j++) solIter[j] = sol[j];
   cout << "Solution: "; for (j = 0; j<n; j++) cout << solIter[j] << " "; cout << endl;
   cout << "VLSN. Starting from z=" << zubIter << endl;

   try
   {
      for (iter = 0; iter<maxiter; iter++)
      {
         fixVariables(c, CPX, solIter, k); // defines the set of free variables

         status = CPX->solveMIP(true, fVerbose); // integer solution, only ejection set variables are free
         if (!status)
         {
            if (fVerbose && n*m < 201)
            {
               cout << " - Solution: " << endl;
               for (i = 0; i < m; ++i)
               {
                  cout << "   " << i << ":  ";
                  for (j = 0; j < n; ++j)
                     cout << CPX->x[i*n + j] << "   ";
                  cout << endl;
               }
            }
            // reads the solution
            zubIter = 0;
            for (j = 0; j<n; j++)
            {
               for (i = 0; i<m; i++)
                  if (CPX->x[i*n + j] > 0.5)
                  {
                     solIter[j] = i;
                     break;
                  }
               zubIter += c[solIter[j]][j];
            }

            if (abs(zubIter - GAP->checkSol(solIter)) > GAP->EPS)
               cout << "[verylarge] No feasible solution at this iteration" << endl;
            else
            {
               cout << "[VLSN] iter " << iter << " zubIter " << zubIter << endl;
               cout << "Solution: "; for (j = 0; j<n; j++) cout << solIter[j] << " "; cout << endl;

               if (zubIter < zub)
               {
                  zub = zubIter;
                  for (i = 0; i<n; i++) solbest[i] = solIter[i];
                  cout << "[VLSN] ************** new zub. iter " << iter << " zubIter " << zubIter << endl;
               }
            }
         }
      }
   }
   catch (std::exception const& e)
   {
      cout << "Error: " << e.what() << endl;
      goto lend;
   }

lend:
   CPX->freeMIP();
   delete(CPX);
   delete(solIter);
   return zub;
}

// this frees the  variables in the ejection set
void VeryLarge::fixVariables(int** c, MIPCplex* CPX, int* solIter, int k)
{
   int i, j, isol, rand, status,numfix;
   bool isInSet;
   double minr;
   vector<int> lstClients;
   vector<int> lstSet;
   vector<int> indRatios(n);
   vector<double> ratios;  // to be used for client ordering
   auto compRatios = [&ratios](double a, double b) { return ratios[a] < ratios[b]; };    // ASC order

   for (j = 0; j<n; j++)
      lstClients.push_back(j);

   for(j=0;j<n;j++)
   {  minr = DBL_MAX;
      for(i=0;i<m;i++)
         if(c[i][j] < minr)
            minr = c[i][j]/(double)req[i][j];
      ratios.push_back(minr);
   }

   // order indRatios by increasing ratios
   for (j = 0; j<n; j++) indRatios[j] = j;
   std::sort(indRatios.begin(), indRatios.end(), compRatios);

   numfix = n - min(k, n);      // numfix is the number of clients to fix
   for (int ii = 0; ii<numfix; ii++)
   {  i = indRatios[ii];
      lstSet.push_back(lstClients[i]);          // set contains numfix best clients
   }

   int  cnt = n * m;
   int* indices = new int[cnt];  // indices of the columns corresponding to the variables for which bounds are to be changed
   char* lu = new char[cnt];     // whether the corresponding entry in the array bd specifies the lower or upper bound on column indices[j]
   double* bd = new double[cnt]; // new values of the lower or upper bounds of the variables present in indices

   cout << "VLSN: binding set ";
   for (j = 0; j<n; j++)
   {
      // is the client in the binding set?
      isInSet = false;
      if (std::find(lstSet.begin(), lstSet.end(), j) != lstSet.end())
      {
         isInSet = true;
         cout << j << " ";
      }

      for (i = 0; i<m; i++)
      {
         if (!isInSet)              // set it free
         {
            if (CPX->lb[i*n + j] == 0.0)  // lb OK, only ub could be wrong
            {
               CPX->ub[i*n + j] = 1.0;
               bd[i*n + j] = 1.0;
               lu[i*n + j] = 'U';        // bd[j] is an upper bound
            }
            else                       // lb wrong, ub should be OK
            {
               CPX->lb[i*n + j] = 0.0;
               bd[i*n + j] = 0.0;
               lu[i*n + j] = 'L';        // bd[j] is an upper bound
            }
         }
         else                          // not in set
         {
            isol = solIter[j];
            if (i == isol)
            {
               CPX->lb[i*n + j] = 1.0;   // fixing in the solution
               bd[i*n + j] = 1.0;
            }
            else
            {
               CPX->lb[i*n + j] = 0.0;   // forbidding in the solution
               bd[i*n + j] = 0.0;
            }
            lu[i*n + j] = 'B';           // bd[j] is the lower and upper bound
         }
         indices[i*n + j] = i * n + j;
      }
   }
   cout << endl;

   status = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
   free(indices);
   free(lu);
   free(bd);

   return;
}
