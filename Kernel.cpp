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

int Kernel::solveByKernel(bool fVerbose, int numBuckets)
{  int i,j,status,numCol,iter,numKer=0,numNotLB,numVarInBucket,lastVar=0,numAdded;
   double zlb;
   vector<double> x,d; // milp primal and dual vars

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
   vector<bool> kernel(n*m);
   MIPCplex* CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(isVerbose);

   try
   {  iter = 0;
      status = CPX->solveMIP(false,false); // LP solution
      if ( status )
      {  cout << "[solveByKernel] Error" << endl;
         goto lend;
      }
      // reads the solution
      zlb = CPX->objval;
      cout << "iter "<< iter << " zlb " << zlb << endl;

      // primal variables (m*n)
      x.clear();
      for(i=0;i<n*m;i++)
      {  x.push_back(CPX->x[i]);
         if(x[i]>0)
         {  kernel[i] = true; // initial kernel
            numKer++;
         }
      }
      numNotLB = n*m-numKer;
      numVarInBucket = (int) ( 1.0*numNotLB / numBuckets + 0.999);

      // dual variables, assignments (n) and capacity (m). Here unused
      d.clear();
      for(i=0;i<n+m;i++) d.push_back(CPX->pi[i]);

      do
      {
         updateModelWithKernel(CPX,kernel);
         status = CPX->solveMIP(true,false); // LP solution
         if ( status )
         {  cout << "[solveByKernel] No solution" << endl;
         }
         else
         {
            // reads the solution
            zub = CPX->objval;
            cout << "iter "<< iter << " zub " << zub << endl;

            // primal variables (m*n)
            x.clear();
            for(i=0;i<n*m;i++)
               x.push_back(CPX->x[i]);
         }

         numAdded = j = 0;
         while (numAdded < numVarInBucket)
         {
            if(!kernel[j])
               kernel[j] = true;
            j++;
            if(j==kernel.size())
               break;
         }
         iter++;
      } while (iter < 100);
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

// defines the vars that can enter the solution
void Kernel::updateModelWithKernel(MIPCplex* CPX, vector<bool> kernel)
{  int i,j,isol,rand,status;
   vector<int> lstKernel;

   for(i=0;i<n*m;i++)
      if(kernel[i])
         lstKernel.push_back(i);    // kernel contains vars

   int  cnt = n*m;
   int* indices = new int[cnt];    // indices of the columns corresponding to the variables for which bounds are to be changed
   char*   lu   = new char[cnt];   // whether the corresponding entry in the array bd specifies the lower or upper bound on column indices[j]
   double* bd   = new double[cnt]; // new values of the lower or upper bounds of the variables present in indices

   cout << "kernel: ";
   for(i=0;i<n*m;i++)
   {
      if(kernel[i])             // set it free
      {  
         if(CPX->lb[i] == 0.0)  // lb OK, only ub could be wrong
         {  CPX->ub[i] = 1.0; 
            bd[i] = 1.0;
            lu[i] = 'U';        // bd[j] is an upper bound
         }
         else                   // lb wrong, ub should be OK
         {  CPX->lb[i] = 0.0; 
            bd[i] = 0.0;
            lu[i] = 'L';        // bd[j] is an upper bound
         }
      }
      else                      // not in kernel, cannot be chosen
      {  
         CPX->lb[i] = 0.0;   // forbidding in the solution
         bd[i] = 0.0;
         lu[i] = 'B';           // bd[j] is the lower and upper bound
      }
      indices[i] = i;
   }
   cout << endl;

   status = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
   delete(indices);
   free(lu);
   free(bd);

   return;
}
