#include "ACO.h"

ACO::ACO(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

ACO::~ACO()
{
   //dtor
}

int ACO::antColony(int** c, int maxiter, int numpop, double alpha)
{  int z=0;
   int i,j,k,iter,index;
   double lb, lbMove;
   vector< vector <int> > pop(numpop);
   vector<double> moveProb(m);

   // Lower bound computation, linear bound
   int numRows, numCols, numNZrow, statusMIP;
   numRows = n + m;   // num of constraints
   numCols = n * m;   // num of variables
   numNZrow = n * m;   // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows, numCols, numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(true); // verbose output
   statusMIP = CPX->solveMIP(false, false);   // LP of whole instance
   if (statusMIP == 0)
   {
      lb = CPX->objval;
      cout << "Linear bound: " << lb << endl;
      if (n * m < 101)                       // print LP if instance is small
      {
         cout << " - Solution: " << endl;
         for (i = 0; i < m; ++i)
         {
            cout << "   " << i << ": ";
            for (j = 0; j < n; ++j)
               cout << CPX->x[i * n + j] << "\t";
            cout << endl;
         }
      }
   }
   else
      lb = DBL_MAX;

   // trail initialization to lb
   tau.assign(m, vector < double >(n, 0));
   for (i = 0; i < m; ++i)
      for (j = 0; j < n; ++j)
         tau[i][j] = CPX->x[i * n + j];

   // initialize empty solutions
   pop.assign(numpop, vector <int> (n,-1));

   cout << "Population initialized, starting search" << endl;
   for (iter = 0; iter < maxiter; iter++) // termination condition on num iterations
   {  if(iter%1 == 0)
         cout << "====== iteration " << iter << " zub = " << zub << endl;
      int cnt = 1;
      int* indices = new int[cnt];     // which var
      char*   lu   = new char[cnt];    // lower limit
      double* bd   = new double[cnt];  // new bound

      // for each ant
      for (k = 0; k < numpop; k++)
      {  
         vector<int> sol (n,-1);
         for (j = 0; j < n; j++)    // for each client to assign
         {
            if (j > 0 && sol[j-1] == -1)   // infeasible partial solution
            {  sol[j] = -1;
               continue;
            }

            for (int j1 = 0; j1 < j; j1++)   // partial solution
            {
               indices[j1] = sol[j1] * n + j1;
               lu[j1] = 'L';
               bd[j1] = 1.0;
            }

            for (i = 0; i < m; i++) // for each server it can be assigned to
            {
               indices[j] = i * n + j;
               lu[j] = 'L';
               bd[j] = 1.0;
               statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
               statusMIP = CPX->solveMIP(false, false);   // LP of partial solution
               lbMove = CPX->objval;
               moveProb[i] = lbMove - lb;
               bd[j] = 0.0;                               // change back to free status
               statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
            }
            index = montecarlo(moveProb);
            if (index < m)
               sol[j] = index;
            else
               sol[j] = -1;   // infeasible partial solution

            for (int j1 = 0; j1 < j; j1++)   // reset partial solution
            {
               indices[j1] = sol[j1] * n + j1;
               lu[j1] = 'L';
               bd[j1] = 0.0;
            }
            statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
         }
         int z = GAP->checkSol(sol);
         if(z<INT_MAX)
            for(j=0;j<n;j++)
               pop[k][j] = sol[j];
      }
      delete(indices);
      delete(lu);
      delete(bd);
   }

   double zcheck = 0;
   zcheck = GAP->checkSol(solbest);
   if (abs(zcheck - zub) > GAP->EPS)
   {  cout << "[ANTS] Detected inconsistency" << endl;
      z = INT_MAX;
   }
   cout << "ANTS zub= " << zub << endl;

   if(CPX != NULL) delete CPX;
   return zub;
}

int ACO::montecarlo(vector<double>& v)
{  double sum=0;
   unsigned int i;

   for(i=0;i<v.size();i++)
      sum += v[i];

   double f = sum * ( (double)rand() / RAND_MAX );
   sum = 0;

   for(i=0;i<v.size();i++)
   {  sum += v[i];
      if(sum >= f) break;
   }
   return i;
}

// dynamic fitness scaling of cost value
double ACO::dynFitScal(double lb, double tau0, double zcurr, double zbar)
{  double res =0;
   
   return res;
}


