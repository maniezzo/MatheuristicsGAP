#include "TabuSearch.h"

TabuSearch::TabuSearch(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

TabuSearch::~TabuSearch()
{
   //dtor
}

/*
1.	Generate an initial feasible solution S, 
set S* = S and initialize TL=nil.
2.	Find S' \in N(S), such that 
z(S')=min {z(S^), foreach S^\in N(S), S^\notin TL}.
3.	S=S', TL=TL \in {S}, if (z(S*) > z(S)) set S* = S.
4.	If not(end condition) go to step 2.
*/
// This is based on opt11
int TabuSearch::tabuSearch(int** c, int Ttenure, int maxiter, int fGoMath)
{  int z=0,delta,DeltaMax,cap1,cap2;
   int i,j,i1,i2,j1,j2, i1max, j1max, i2max, j2max,iter;
   vector<int> capleft(m);
   vector<vector<int>> TL(m, vector<int>(n,0));   // m*n 2d vector initialized with 0
   double zcheck;

   // no solution to improve
   if(zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }
   Rins* RINS = new Rins(GAP, GAP->zub);   // in case of matheuristic

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
   {  capleft[sol[j]] -= req[sol[j]][j];
      z += c[sol[j]][j];
   }

   iter = 0;
   Ttenure = Ttenure*m;       // very heuristic guess
   cout << "Starting tabu search" << endl;

start:   
   DeltaMax = i1max = j1max = i2max = j2max =INT_MIN;

   iter++;
   for (j1 = 0; j1<n; j1++)            // 1-1 neighborhood exploration
      for (j2 = j1 + 1; j2<n; j2++)
      {
         i1 = sol[j1];
         i2 = sol[j2];
         delta = (c[i1][j1] + c[i2][j2]) - (c[i1][j2] + c[i2][j1]); // ante minus post
         if (delta > DeltaMax && TL[i1][j1]<iter && TL[i2][j2]<iter)
         {
            cap1 = capleft[i1] + req[i1][j1] - req[i1][j2];
            cap2 = capleft[i2] + req[i2][j2] - req[i2][j1];
            if (cap1 >= 0 && cap2 >= 0)
            {
               DeltaMax = delta;
               j1max = j1;
               j2max = j2;
               i1max = i1;
               i2max = i2;
            }
         }
      }

   if(i1max == INT_MIN || j1max == INT_MIN || i2max == INT_MIN || j2max == INT_MIN)
   {  cout << "No feasible neighbor found, exiting" << endl;
      goto lend;
   }

   // new incumbent solution
   capleft[i1max] += req[i1max][j1max] - req[i1max][j2max];
   capleft[i2max] += req[i2max][j2max] - req[i2max][j1max];
   sol[j1max] = i2max;
   sol[j2max] = i1max;
   z -= DeltaMax;
   TL[i1max][j2max] = iter + Ttenure;
   TL[i2max][j1max] = iter + Ttenure;

   if(z < zub)                            // improved solbest
   {  zcheck = 0;
      zcheck = GAP->checkSol(sol);
      if (abs(zcheck - z) > GAP->EPS)
      {  cout << "[TS] Solution structure inconsistency" << endl;
         z = INT_MAX;
         goto lend;
      }
      for(i=0;i<n;i++) solbest[i]=sol[i];
      zub = z;
      cout << "[TS] new zub " << zub << " iter " << iter << endl;

      if(fGoMath)
      {
         cout << "Trying MIP improvement" << endl;
         dive(RINS,solbest);
      }
   }

   if(iter%1000 == 0)
      cout << iter <<")  z " << z << " iter " << iter << " deltamax " << DeltaMax << endl;

   if(iter<maxiter)          // end condition
      goto start;
   else
      cout << "Tabu search: end" << endl;

lend:
   cout << "TS zub= " << zub << endl;
   printIntArray(sol, n);
   if (RINS != NULL) delete RINS;
   RINS = NULL;
   return z;
}

int TabuSearch::dive(Rins* RINS, int* newsol)
{  int res;

   res = RINS->dive(GAP->c,
      GAP->conf->rinsConf->maxnodes,
      zub,
      newsol, false
   );
   return res;
}