#include "GAP.h"

Config::Config()
{
   //ctor
   SA = new SimAnn();
   TS = new Tabu();
   GA = new GeneticConf();
   EC = new EjectionConf();
   IterLS  = new IteratedLS();
   lagrAss = new LagrAss();
   lagrCap = new LagrCap();
   rinsConf= new RinsConf();
   beam    = new Beam();
   fbConf  = new FBconf();
   corridorConf = new CorridorConf();
   locBranching = new LocBranching();
   bendersConf  = new BendersConf();
}

GeneralizedAssignemnt::GeneralizedAssignemnt()
{
   //ctor
   conf = new Config();
   EPS = 0.0001;
}

GeneralizedAssignemnt::~GeneralizedAssignemnt()
{  int i;
   //dtor
   if (req!=NULL)
   {  for (i = 0; i<m; i++) free(req[i]);
      free(req);
   }
   if (c!=NULL)
   {  for (i = 0; i<m; i++) free(c[i]);
      free(c);
   }
   if (cap!=NULL)
      free(cap);
}

// controllo ammissibilità soluzione
int GeneralizedAssignemnt::checkSol(int* sol)
{  int cost=0;
   int i,j;
   int* capused = new int[m];
   for(i=0;i<m;i++) capused[i] = 0;

   // controllo assegnamenti
   for(j=0;j<n;j++)
      if(sol[j]<0 || sol[j]>=m)
      {  cost = INT_MAX;
         goto lend;
      }
      else
         cost += c[sol[j]][j];

   // controllo capacità
   for(j=0;j<n;j++)
   {  capused[sol[j]] += req[sol[j]][j];
      if(capused[sol[j]] > cap[sol[j]])
      {  cost = INT_MAX;
         goto lend;
      }
   }
   delete capused;
lend:    
   return cost;
}

// recovers feasibility in case of partial or overassigned solution
int GeneralizedAssignemnt::fixSol(int* infeasSol, int* zsol)
{  int i,j,imin=-1;
   int minreq;
   vector<int> capres, sol;

   for(i=0;i<m;i++) capres.push_back(cap[i]);
   for(i=0;i<n;i++) sol.push_back(infeasSol[i]);

   // ricalcolo capacità residue. Se sovrassegnato, metto a sol a -1
   for(j=0;j<n;j++)
      if(sol[j]>=0 && (capres[sol[j]] >= req[sol[j]][j]))
         capres[sol[j]] -= req[sol[j]][j];
      else
         sol[j] = -1;

   *zsol = 0;
   for(j=0;j<n;j++)
   {  if(sol[j]>=0)              // correct, do nothing
      {  *zsol += c[sol[j]][j];
         continue;
      }

      // reassign i -1
      minreq = INT_MAX;
      imin = -1;
      for(i=0;i<m;i++)
         if(capres[i]>=req[i][j] && req[i][j] < minreq)
         {  minreq = req[i][j];
            imin    = i;
         }

      if(imin<0)
      {  *zsol = INT_MAX;
         goto lend;           // could not recover feasibility
      }
      sol[j]=imin;
      capres[imin] -= req[imin][j];
      *zsol += c[imin][j];
   }

   if(*zsol<zub)
   {  for(i=0;i<n;i++) solbest[i]=sol[i];
      zub = zub = *zsol;
      cout << "[fixSol] -------- zub improved! " << zub << endl;
   }
   for(i=0;i<n;i++) infeasSol[i]=sol[i];

lend:
   return *zsol;
}

// recovers feasibility via knapsacks on residual capacities
int GeneralizedAssignemnt::fixSolViaKnap(int* infeasSol, int* zsol)
{  int i,j,imin=-1,minreq;
   int *sol, *minReqFacility;
   int nelem,nreset;                    // number of clients to reallocate and reallocated
   vector<int> capres,indCap;   
   auto compCap = [&capres](int a, int b){ return capres[a] < capres[b]; };  // ompare, ASC order

   *zsol = INT_MAX;
   for(i=0;i<m;i++) capres.push_back( cap[i] );

   sol = (int*) malloc(n * sizeof(int));
   for(j=0;j<n;j++)     // randomly fill unassigned clients
   {  sol[j]=infeasSol[j];
      if(sol[j]<0 || sol[j] >= m)
         sol[j] = rand()%m;
      // residual capacities
      capres[sol[j]] -= req[sol[j]][j];
   }

   // finds the least requiring assignment for each client
   nelem = 0;
   minReqFacility = (int*) malloc(n * sizeof(int));   // best server for each client
   vector<int> whoIs;
   for(j=0;j<n;j++)
   {  minreq = INT_MAX;
      for(i=0;i<m;i++) 
         if(req[i][j] < minreq)
         {  minReqFacility[j]=i;
            minreq = req[i][j];
         }

      // if j assigned to an overloaded facility, and not the least req one, it must be reassigned
      if(capres[sol[j]] < 0 && sol[j] != minReqFacility[j])
      {  capres[sol[j]] += req[sol[j]][j];
         sol[j] = -1;
         nelem++;
         whoIs.push_back(j);
      }
   }

   int*    q   = new int[nelem];     // knapsack requests
   double* val = new double[nelem];  // knapsack profits
   int*    Ksol= new int[nelem];     // knapsack solution

   if(nelem == 0)    // if we got a feasible solution by pure chance, done
      goto lfeas;

   // order facilities by increasing residual capacity
   for(i=0;i<m;i++) indCap.push_back(i);
   std::sort(indCap.begin(), indCap.end(), compCap);

   for(j=0;j<nelem;j++) val[j] = 1; // any element could be chosen (val > 0)
   nreset = 0;

   for(int ii=0;ii<m;ii++)
   {  i = indCap[ii];               // consider first warehouaes with small residual capacity
      if(capres[i]<0)               // should not happen, to be debugged
         continue;

      for (j = 0; j < nelem; j++)
      {  q[j]    = req[i][whoIs[j]];         // requests to the i-th wh by the j-th elem to reassign
         Ksol[j] = 0;
      }

      KDynRecur(nelem,capres[i],q,val,Ksol); // solve the knapsack
      for(j=0;j<nelem;j++)
         if(Ksol[j] > 0)
         {  sol[whoIs[j]] = i;
            val[j] = -1;                     // won't be chosen again
            nreset++;
         }
      if(nreset == nelem)                    // solution complete
         break;
   }
   if(nreset < nelem) goto lend;             // could not recover fesibility

lfeas:
   for(i=0;i<n;i++) infeasSol[i]=sol[i];
   *zsol = checkSol(sol);
   if(*zsol<zub)
   {  for(i=0;i<n;i++) solbest[i]=sol[i];
      zub = *zsol;
      cout << "[fixSol] -------- zub improved! " << zub;
   }

lend:
   if(q!=NULL) delete(q);
   if(val!=NULL) delete(val);
   if(Ksol!=NULL) delete(Ksol);
   free(sol);
   free(minReqFacility);
   return *zsol;
}

// computes the aversion (opposite preference) of assigning j to i
int GeneralizedAssignemnt::aversion(int i, int j)
{  int res = -1,whichFunc,first;

   whichFunc = conf->aversionf;
   switch (whichFunc)
   {
   case 1:     // assignment cost
      res = c[i][j];
      break;
   case 2:     // required resource amount
      res = req[i][j];
      break;
   case 3:     // regret
         first = INT_MAX;
         for (int ii = 0; ii<m; ii++)
            if (c[ii][j] < first)
               first = c[ii][j];
         res = c[i][j] - first;
      break;
   case 4:     // resource unit cost
      res = (100.0 * c[i][j]) / req[i][j];
      break;
   case 5:     // percentage server occupation
      res = (100.0 * req[i][j])/cap[i];
      break;
   default:
         cout << "Aversion function badly specified" << endl;
         break;
   }

   return res;
}
// **************************************************************************** //
// *************************** Free functions ********************************* //


// computes assignment regrets for each client
void computeRegrets(int** c, int n, int m, vector<int> & regrets)
{  int i,j,first,second;

   for(j=0;j<n;j++)
   {  first = second = INT_MAX;
      for(i=0;i<m;i++)
      {  if(c[i][j] < first)
         {  second = first;
            first  = c[i][j];
         }
         else if (c[i][j] < second)
            second  = c[i][j];
      }
      regrets[j] = second - first;
   }
}

// dynamic programming recursion for the knapsack
double KDynRecur(int n, int Kcap, int* Q, double* val, int* Ksol)
{  int q,i,imax;
   double res=0;

   double** f = (double**) calloc((Kcap+1), sizeof(double*)); // init to 0
   for(i=0;i<Kcap+1;i++)
      f[i] = (double*) calloc(n, sizeof(double));             // init to 0

   for (i = 0; i < n; i++)
   {
      //if (val[i] < 0) continue;
      for (q = 0; q <= Kcap; q++)
         switch (i)
         {
            case 0:
               if (q >= Q[i])
                  f[q][i] = max(0.0, val[i]);
               else
                  f[q][i] = 0;
               break;
            default:
               if (q >= Q[i])
                  f[q][i] = max(f[q][i - 1], f[q - Q[i]][i - 1] + val[i]);
               else
                  f[q][i] = f[q][i - 1];
               break;
         }
   }

   imax=0;
   for(i=0;i<n;i++)
      if(f[Kcap][i] > res)
      {  res = f[Kcap][i];
         imax = i;
      }

   KdecodeSol(imax,Kcap,Q,val,n,f,Ksol);

   // deallocations
   if (f!=NULL)
   {  for (i = 0; i<Kcap+1; i++) 
         free(f[i]);
      free(f);
   }

   return res;
}

// Decodes a knapsack DP recursion, given the recursion matrix f
void KdecodeSol(int i, int Kcap, int* Q, double* val, int n, double** f, int* Ksol)
{  int q=Kcap;
   double eps = 0.0001;

   while(q>0 && i>0)
   {  if(abs(f[q][i-1] - f[q][i]) < eps)
      {  i--;
         continue;
      }

      if(abs(f[q-Q[i]][i-1] - (f[q][i]-val[i])) < eps)
      {  q -= Q[i];
         Ksol[i] = 1;
         i--;
         continue;
      }

      cout << "[KP decodeSol] generic error" << endl;
      goto lend;
   }

   if(i==0 && q>0)
      if(f[q][i] == val[i])
         Ksol[i] = 1;
      else
         Ksol[i] = 0;

lend:
   //checkSol();
   return;
}

// print a 1D array of ints
void printIntArray(int* a, int n)
{
   int i;
   for (i = 0; i<n; i++)
      cout << a[i] << " ";
   cout << endl;
}

// print a 1D array of doubles
void printDblArray(double* a, int n)
{
   int i;
   cout.precision(3);
   for (i = 0; i<n; i++)
      cout << a[i] << " ";
   cout << endl;
}