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
int TabuSearch::tabuSearch(int** c, int Ttenure, int maxiter)
{  int z=0,DeltaMax;
   int i,isol,j,imax,jmax,iter;
   vector<int> capleft(m);
   vector<vector<int>> TL(m, vector<int>(n,0));   // m*n 2d vector initialized with 0

   // no solution to improve
   if(zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
   {  capleft[sol[j]] -= req[sol[j]][j];
      z += c[sol[j]][j];
   }

   iter = 0;
   Ttenure = Ttenure*m;       // very heuristic guess
   cout << "Starting tabu search" << endl;

start:   
   DeltaMax=imax=jmax=INT_MIN;

   iter++;
   for (j = 0; j < n; j++)
   {
      isol = sol[j];
      for (i=0; i<m; i++)
      {  if (i == isol) continue;

         if( (c[isol][j] - c[i][j]) > DeltaMax && capleft[i] >= req[i][j]  && TL[i][j]<iter)
         {  imax = i;
            jmax = j;
            DeltaMax = c[isol][j] - c[i][j];
         }
      }
   }

   isol = sol[jmax];
   sol[jmax] = imax;
   capleft[imax] -= req[imax][jmax];
   capleft[isol] += req[isol][jmax];
   z -= DeltaMax;
   if(z < zub)                            // improved solbest
   {  for(i=0;i<n;i++) solbest[i]=sol[i];
      zub = z;
   }
   TL[imax][jmax] = iter+Ttenure;
   if(iter%1000 == 0)
      cout << iter <<")  z " << z << " ite " << iter << " deltamax " << DeltaMax << endl;

   if(iter<maxiter)          // end condition
      goto start;
   else
      cout << "Tabu search: fine" << endl;
   double zcheck = 0;
   zcheck = GAP->checkSol(sol);
   if (abs(zcheck - z) > GAP->EPS)
   {  cout << "[1.0opt] Ahi ahi" << endl;
      z = INT_MAX;
   }
   cout << "Opt01 zub= " << zub << endl;
   return z;
}

