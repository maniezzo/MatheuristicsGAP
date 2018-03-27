#include "SimAnnealing.h"

SimAnnealing::SimAnnealing(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

SimAnnealing::~SimAnnealing()
{
   //dtor
}

/*
1.	Generate an initial feasible solution S, 
initialize S* = S and temperature parameter T.
2.	Generate S’\in N(S).
3.	If z(S') < z(S) then  S=S', if (z(S*) > z(S)) S* = S
else accept to set S=S' with probability p = e-(z(S')-z(S))/kT.
4.	If (annealing condition) decrease T.
5.	If not(end condition) go to step 2.
*/

// This assumes to have in sol a solution to improve
int SimAnnealing::simAnnealing(int** c, double k, double maxT, double alpha, int maxiter, int iterAnneal, double coefAnneal, int fGoMath)
{  int z=0;
   int i,isol,j,i1,j1,iter,rand,lastimpr;
   vector<int> capleft(m);
   double T,p;

   // no solution to improve
   if(zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      zub = INT_MAX;
      goto lend;
   }

   if(fGoMath)
   {  zub = mathSimAnnealing(c, k, maxT, alpha, maxiter, iterAnneal, coefAnneal);
      goto lend;
   }

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
   {  capleft[sol[j]] -= req[sol[j]][j];
      z += c[sol[j]][j];
   }

   z    = zub;
   T    = maxT;
   iter = 0;
   lastimpr   = 0;
   cout << "Starting Simulated Annealing, zub " << zub << endl;

start:
   iter++;   
   cout << setw(15);
   if(iter%1000 == 0)
      cout << "Iter "<<iter<<" T "<< T <<" zub "<< zub <<" z "<<z<<endl;

   // neighborhood 01
   j = std::rand() % n;
   isol = sol[j];
   i = std::rand() % m;

   if(i==isol) i = (i+1) % m;

   if (c[i][j] < c[isol][j] && capleft[i] >= req[i][j])
   {  // remove from isol and assign to i
      sol[j] = i;
      capleft[i]    -= req[i][j];
      capleft[isol] += req[isol][j];
      z -= (c[isol][j] - c[i][j]);
      if(z<zub)
      {  zub = z;
         for(int k=0;k<n;k++) solbest[k] = sol[k];
         cout << "[SA] new zub " << zub << endl;
      }
   }
   else if(capleft[i] >= req[i][j])
   {  //p = e-(z(S')-z(S))/kT 
      p    = exp(-(c[i][j] - c[isol][j])/(k*T)); 
      rand = std::rand() % 100;
      if(rand<p*100)
      {  sol[j] = i;
         capleft[i]    -= req[i][j];
         capleft[isol] += req[isol][j];
         z -= (c[isol][j]-c[i][j]);
         //cout << "[SA] "<< iter <<") Worsening move p= "<< p <<" rand= "<< rand <<" diff. "<<(c[isol][j]-c[i][j]) << endl;
      }
   }        

   if((iter%iterAnneal)==0)                        // annealing condition
      T = coefAnneal*T;

   if(T>0.005 || (iter-lastimpr)<(maxiter/10))     // restart condition
      goto start;
   else if(iter < maxiter)                         // restart
   {  T=maxT;
      for(i=0;i<n;i++)
      {  i1 = std::rand() % n;
         j1 = std::rand() % m;
         if(capleft[i1] >= req[i1][j1])
         {  sol[j1] = i1;
            capleft[i1]   -= req[i1][j1];
            capleft[isol] += req[i1][j1];
            z -= (c[isol][j1]-c[i1][j1]);
         }
      }         
      goto start;
   }         
   else                                            // terminate
   {  //Loc2opt();
      cout << "Simulated annealing: ending." << endl;
   }

   double zcheck = 0;
   zcheck = GAP->checkSol(sol);
   if (abs(zcheck - z) > GAP->EPS)
   {  cout << "[SA] Detected inconsistency" << endl;
      z = INT_MAX;
   }

lend:
   cout << "SA zub= " << zub << endl;
   return zub;
}

// This assumes to have in sol a solution to improve
int SimAnnealing::mathSimAnnealing(int** c, double k, double maxT, double alpha, int maxiter, int iterAnneal, double coefAnneal)
{
   int z = 0;
   int i, isol, j, i1, j1, iter, rand, lastimpr, nInfeas;
   vector<int> capleft(m);
   double T, p, zcheck;

   vector<vector<int>> cPrime(m, vector<int>(n));

   for (i = 0; i<m; i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
   {  capleft[sol[j]] -= req[sol[j]][j];
      for (i = 0; i<m; i++) cPrime[i][j] = c[i][j];
      z += cPrime[sol[j]][j];             // initially cPrime = c
   }

   z = zub;
   T = maxT;
   iter = 0;
   lastimpr = 0;
   cout << "Starting matheuristic simulated annealing, zub " << zub << endl;

start:
   iter++;
   cout << setw(15);
   if (iter % 1000 == 0)
      cout << "Iter " << iter << " T " << T << " zub " << zub << " z " << z << endl;

   // neighborhood 01
   j = std::rand() % n;
   isol = sol[j];
   i = std::rand() % m;
   if (i == isol) i = (i + 1) % m;

   if (cPrime[i][j] < cPrime[isol][j])
   {  // remove from isol and assign to i
      sol[j] = i;
      capleft[i] -= req[i][j];
      capleft[isol] += req[isol][j];
      z -= (cPrime[isol][j] - cPrime[i][j]);
   }
   else
   {  //p = e-(z(S')-z(S))/kT 
      p = exp(-(cPrime[i][j] - cPrime[isol][j]) / (k*T));
      rand = std::rand() % 100;
      if (rand<p * 100)
      {
         sol[j] = i;
         capleft[i] -= req[i][j];
         capleft[isol] += req[isol][j];
         z -= (cPrime[isol][j] - cPrime[i][j]);
         if(maxiter < 0)
            cout << "[SA] "<< iter <<") Worsening move p= "<< p <<" rand= "<< rand <<" diff. "<<(c[isol][j]-c[i][j]) << endl;
      }
   }

   nInfeas = 0;
   for(i=0;i<m;i++)
      if(capleft[i] < 0) 
      {  nInfeas++;
         for(j=0;j<n;j++)
            cPrime[i][j] -= capleft[i];
      }

   if(nInfeas==0) // got a feasible solution
   {
      z = 0;
      for(j=0;j<n;j++)
         z += GAP->c[sol[j]][j];

      zcheck = 0;
      zcheck = GAP->checkSol(sol);
      if (abs(zcheck - z) > GAP->EPS)
      {
         cout << "[MSA] Detected inconsistency" << endl;
         z = INT_MAX;
      }
      else
         cout << "iter " << iter << " feasible! z="<< z << endl;

      if (z<zub)
      {  zub = z;
         for (int k = 0; k<n; k++) solbest[k] = sol[k];
         cout << "[MSA] new zub " << zub << endl;
      }
   }

   if ((iter%iterAnneal) == 0)                        // annealing condition
      T = coefAnneal*T;

   if (T>0.005 || (iter - lastimpr)<(maxiter / 10))     // restart condition
      goto start;
   else if (iter < maxiter)                         // restart
   {
      T = maxT;
      for (i = 0; i<n; i++)
      {
         i1 = std::rand() % n;
         j1 = std::rand() % m;
         if (capleft[i1] >= req[i1][j1])
         {
            sol[j1] = i1;
            capleft[i1] -= req[i1][j1];
            capleft[isol] += req[i1][j1];
            z -= (cPrime[isol][j1] - cPrime[i1][j1]);
         }
      }
      goto start;
   }
   else                                            // terminate
   {  //Loc2opt();
      cout << "Matheuristic simulated annealing: ending." << endl;
   }
   return zub;
}