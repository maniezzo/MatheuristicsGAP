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
   int i,j,k,iter;
   double lb;
   vector< vector <int> > pop(numpop);

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
   CPX->freeMIP();

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

      // for each ant
      for(k=0;k<numpop;k++)
         for(j=0;j<n;j++)   // for each client to assign
         {
            for(i=0;i<m;i++)// for each server it can be assigned to
            {
               CPX->lb[i * n + j] = 1;
               int cnt = 1;
               int* indices = new int[cnt];     // which var
               indices[0] = i * n + j;
               char* lu = new char[cnt];        // lower limit
               lu[0] = 'L';
               double* bd = new double[cnt];    // new bound
               bd[0] = 1.0;
               statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
               free(indices);
               free(lu);
               free(bd);
            }
         }
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


