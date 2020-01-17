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

int Kernel::solveByKernel(bool fVerbose, int numVarInBucket)
{  int i,j,ii,status,numCol,iter,nkinit,numKer=0,numNotLB,numBuckets,lastVar=0,numAdded;
   double zlb,z;
   vector<double> x,d,dj; // milp primal, dual vars and red costs
   auto compRedCost = [&dj](double a, double b){ return dj[a] > dj[b]; }; // DESC order

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
   vector<bool> kernel(n*m),lambda_i(n*m);
   vector<int> indDj(n*m);
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

      // primal variables (m*n) and reduced costs
      x.clear();
      for(i=0;i<n*m;i++)
      {  x.push_back(CPX->x[i]);
         if(x[i]>0)
         {  kernel[i] = true; // initial kernel
            lambda_i[i] = true;
            numKer++;
         }
         dj.push_back(CPX->dj[i]);  // reduced costs (CPLEX naming!)
         indDj[i]=i;
      }
      nkinit   = numKer;
      numNotLB = n*m-numKer;
      numBuckets = (int) (1.0*numNotLB/numVarInBucket + 0.999);
      //numVarInBucket = (int) ( 1.0*numNotLB / numBuckets + 0.999);

      // dual variables, assignments (n) and capacity (m). Here unused
      d.clear();
      for(i=0;i<n+m;i++) d.push_back(CPX->pi[i]);
      std::sort(indDj.begin(), indDj.end(), compRedCost); // sort by increasing red costs

      CPXsetintparam(CPX->env,CPXPARAM_MIP_Display,0);  // cplex output to screen
      do
      {
         updateModelWithKernel(CPX,lambda_i);
         CPXsetdblparam(CPX->env,CPXPARAM_MIP_Tolerances_UpperCutoff,zub); // cost cut
         status = CPX->solveMIP(true,false); // LP solution
         if ( status )
         {  cout << "[solveByKernel] No solution" << endl;
            z = DBL_MAX;
         }
         else
         {
            // reads the solution
            z = CPX->objval;
            if(z<zub) zub = z;

            // primal variables (m*n)
            x.clear();
            for(i=0;i<n*m;i++)
            {  x.push_back(CPX->x[i]);
               if(CPX->x[i] && !kernel[i]) 
               {  kernel[i] = true;
                  numKer++;
                  cout << "variable " << i << " entering kernel" << endl;
               }
            }
         }

         for(int jj=0;jj<n*m;jj++)
         {  j = indDj[jj];
            if(!kernel[j])
               lambda_i[j] = false;
         }

         numAdded = 0;
         j = lastVar;
         while(numAdded < numVarInBucket)
         {  ii = indDj[j]; // considering vars in reduced costs order
            if(kernel[ii])
               lambda_i[ii] = true;
            else
               if (numAdded < numVarInBucket)
               {  lambda_i[ii] = true;
                  numAdded++;
               }
            j++;
            if(j==n*m) break;
            lastVar = j;
         }
         cout << "iter "<< iter << " z= " << z << " zub=" << zub << endl;
         iter++;
      } while (iter < 100 && numAdded > 0);
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
void Kernel::updateModelWithKernel(MIPCplex* CPX, vector<bool> lambda_i)
{  int i,j,isol,rand,status;
   vector<int> lstKernel;

   for(i=0;i<n*m;i++)
      if(lambda_i[i])
         lstKernel.push_back(i);    // kernel contains vars

   int  cnt = n*m;
   int* indices = new int[cnt];    // indices of the columns corresponding to the variables for which bounds are to be changed
   char*   lu   = new char[cnt];   // whether the corresponding entry in the array bd specifies the lower or upper bound on column indices[j]
   double* bd   = new double[cnt]; // new values of the lower or upper bounds of the variables present in indices

   cout << "Lambda_i:";
   for(i=0;i<n*m;i++) if(lambda_i[i]) cout << " " << i;
   for(i=0;i<n*m;i++)
   {
      if(lambda_i[i])             // set it free
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
   delete[](indices);
   delete[](lu);
   delete[](bd);

   return;
}
