#include "SS.h"
#include <algorithm>

ScatterSearch::ScatterSearch(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

ScatterSearch::~ScatterSearch()
{
   //dtor
}

int ScatterSearch::go_scatter(int** c, int maxiter, int numpop, double alpha)
{  int z=0,zcheck;
   int i,j,k,iter,index,cnt,numsol=0;
   double lb=0, lbMove,avgtau=0,zavg=0;
   vector< vector <int> > pop(numpop);
   vector<vector <double>> deltaTau;
   vector<double> moveProb(m),zpop(numpop);
   vector<int> indCost;                   // for population initialization

   // ------------------------ initialization
   for(j=0;j<n;j++)
      indCost.push_back(j);

   ConstructHeu* CH  = new ConstructHeu(GAP,GAP->zub); 
   int nind = sizeof(indCost) / sizeof(indCost[0]);
   int zsol, cont=0;
   do { 
      zsol = CH->constructWithInd(indCost, false); 
      if(zsol < INT_MAX)
      {  for(j=0;j<n;j++)
            pop[cont].push_back(GAP->sol[j]);
         // printIntArray(&pop[cont][0],n);
         cout << "cont= " << cont << " zsol = " << zsol << endl;
         cont++;
      }
   } while (next_permutation(&indCost[0], &indCost[0]+nind) && cont < numpop); 

   if(CH != NULL) delete CH;
   CH = NULL;

   // ----------------------------------------- solution improvement
   LocalSearch* LS = new LocalSearch(GAP, GAP->zub);
   for(k=0;k<numpop;k++)
   {
      for(j=0;j<n;j++) GAP->sol[j] = pop[k][j];
      zsol = LS->opt10(GAP->c,true);
      for(j=0;j<n;j++) pop[k][j] = GAP->sol[j];
      //zcheck = GAP->checkSol(pop[k]);
      cout << "cont= " << cont << " zsol = " << zsol << endl;
   }
   if (LS != NULL) delete LS;
   LS = NULL;

   // ----------------------------------------- 

   // Lower bound computation, linear bound
   int numRows, numCols, numNZrow, statusMIP;

   numRows = n + m;   // num of constraints
   numCols = n * m;   // num of variables
   numNZrow = n * m;   // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows, numCols, numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(false); // verbose output
   statusMIP = CPXsetintparam(CPX->env, CPXPARAM_ScreenOutput, CPX_OFF);

   statusMIP = CPX->solveMIP(false, false);   // LP of whole instance
   if (statusMIP == 0)
   {
      lb = CPX->objval;
      cout << "Linear bound: " << lb << endl;
      if (n * m < 101)                       // print LP if instance is small
      {
         cout << " - Solution: " << endl;
         for (i = 0; i < m; ++i)
         {  cout << "   " << i << ": ";
            for (j = 0; j < n; ++j)
               cout << CPX->x[i * n + j] << "\t";
            cout << endl;
         }
      }
   }
   else
      lb = 0;

   zcheck = 0;
   zcheck = GAP->checkSol(solbest);
   if (abs(zcheck - zub) > GAP->EPS)
   {  cout << "[SS] Detected inconsistency" << endl;
      z = INT_MAX;
   }
   else
      cout << "zcheck  " << zcheck << endl;
   cout << "SS zub= " << zub << endl;

   if(CPX != NULL) delete CPX;
   return zub;
}

int ScatterSearch::montecarlo(vector<double>& v)
{  double sum=0;
   unsigned int i;

   for(i=0;i<v.size();i++)
      sum += v[i];
   if(sum<=0) return -1;   // no feasible choice

   double f = sum * ( (double)rand() / RAND_MAX );
   sum = 0;
   for(i=0;i<v.size();i++)
   {  sum += v[i];
      if(sum >= f) break;
   }
   return i;
}


