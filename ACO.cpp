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
   int i,j,k,iter,index,cnt,numsol=0;
   double lb=0, lbMove,avgtau=0,zavg=0;
   vector< vector <int> > pop(numpop);
   vector<vector <double>> deltaTau;
   vector<double> moveProb(m),zpop(numpop);

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

   // trail initialization to lb
   tau.assign(m, vector<double>(n, 0));
   for (i = 0; i < m; ++i)
      for (j = 0; j < n; ++j)
      {  tau[i][j] = CPX->x[i * n + j];
         avgtau += tau[i][j]/(m*n);
      }
   vector<double> zlast(n,lb);

   // initialize empty solutions
   pop.assign(numpop, vector <int> (n,-1));

   cout << "Population initialized, starting search" << endl;
   for (iter = 0; iter < maxiter; iter++) // termination condition on num iterations
   {  if(iter%1 == 0)
         cout << "====== iteration " << iter << " zub = " << zub << endl;

      // for each ant
      for (k = 0; k < numpop; k++)
      {  zpop[k] = -1;
         vector<int> sol (n,-1);
         for (j = 0; j < n; j++)    // for each client to assign
         {
            if (j > 0 && sol[j-1] == -1)   // infeasible partial solution
            {  sol[j] = -1;
               continue;
            }

            cnt = j+1;
            vector<int> indices(cnt);
            vector<char> lu(cnt);
            vector<double> bd(cnt);

            for (int j1 = 0; j1 < j; j1++)   // partial solution
            {  indices[j1] = sol[j1] * n + j1;
               lu[j1] = 'B';
               bd[j1] = 1.0;
            }

            for (i = 0; i < m; i++) // for each server it can be assigned to
            {
               indices[j] = i * n + j;
               lu[j] = 'L';
               bd[j] = 1.0;
               statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, &indices[0], &lu[0], &bd[0]);
               statusMIP = CPX->solveMIP(false, false);   // LP of partial solution
               if(statusMIP)
               {  moveProb[i] = 0;
                  cout << "Infeasible choice" << endl;
               }
               else
               {  lbMove = CPX->objval;
                  moveProb[i] = alpha*tau[i][j] + (1-alpha)*1/(1+lbMove - lb);
                  cout << "lbmove = " << lbMove << " moveprob " << moveProb[i] << endl;
               }
               bd[j] = 0.0;                               // change back to free status
               statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, &indices[0], &lu[0], &bd[0]);
            }

            index = montecarlo(moveProb);
            if (index < m)
            {  cout << "chosen index " << index << endl;
               sol[j] = index;
            }
            else
            {  cout << "infeasible partial solution" << endl;
               sol[j] = -1;   // infeasible partial solution
               zpop[k] = -1;
            }

            // reset partial solution
            cnt = j;
            indices.resize(cnt); // which var
            lu.resize(cnt);      // lower limit
            bd.resize(cnt);      // new bound

            for (int j1 = 0; j1 < j; j1++)   
            {  indices[j1] = sol[j1] * n + j1;
               lu[j1] = 'L';
               bd[j1] = 0.0;
            }
            statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, &indices[0], &lu[0], &bd[0]);
         }
         int z = GAP->checkSol(sol);
         if(z<INT_MAX)
         {  if(z<zub) 
            {  zub = z;
               for(i=0;i<n;i++) solbest[i]=sol[i];
               cout << " ***** New zub = " << zub << endl;
            }
            for(j=0;j<n;j++)
               pop[k][j] = sol[j];
            zpop[k] = z;
            numsol++;   // feasible solution counter
            zlast[numsol%n] = z; // to compute the average of the last n sol
         }
         else
            cout << " -- not good --" << endl;
      }

      // ------------------------------------   trail update
      deltaTau.assign(m, vector<double>(n, 0));
      zavg=0;
      for(j=0;j<n;j++) zavg+=zlast[j];
      zavg = zavg/n;
      for(k=0;k<numpop;k++)
         if(zpop[k]>=0)    // if feasible solution
            for(j=0;j<n;j++)
            {  i = pop[k][j];
               deltaTau[i][j] += avgtau*(1-(zpop[k]-lb)/(zavg-lb));
            }
      for(i=0;i<m;i++)
         for(j=0;j<n;j++)
            tau[i][j] += deltaTau[i][j];
   }

   double zcheck = 0;
   zcheck = GAP->checkSol(solbest);
   if (abs(zcheck - zub) > GAP->EPS)
   {  cout << "[ANTS] Detected inconsistency" << endl;
      z = INT_MAX;
   }
   else
      cout << "zcheck  " << zcheck << endl;
   cout << "ANTS zub= " << zub << endl;

   if(CPX != NULL) delete CPX;
   return zub;
}

int ACO::montecarlo(vector<double>& v)
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

// dynamic fitness scaling of cost value
double ACO::dynFitScal(double lb, double tau0, double zcurr, double zbar)
{  double res =0;
   
   return res;
}


