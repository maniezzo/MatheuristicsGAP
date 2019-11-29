#include "SS.h"
#include <algorithm>    // for permutations
#include <functional>   // for vector hashing

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

/*
CAUTION: algorithm incomplete! Here just the code to generate the example in the text !!
*/
int ScatterSearch::go_scatter(int** c, int maxiter, int numpop, double alpha)
{  int z=0,zcheck,zsol,id;
   int j,k,k1,k2,idsol,numsol=0;
   vector< vector <int> > pop(numpop);
   vector<vector <double>> deltaTau;
   vector<double> moveProb(m),zpop(numpop);
   vector<int> indCost;                   // for population initialization

   cout << "Scatter search algorithm incomplete !!" << endl;

   // ------------------------ initialization
   for(j=0;j<n;j++)
      indCost.push_back(j);

   ConstructHeu* CH  = new ConstructHeu(GAP,GAP->zub); 
   int nind = sizeof(indCost) / sizeof(indCost[0]);
   int cont=0;
   do { 
      zsol = CH->constructWithInd(indCost, false); 
      if(zsol < INT_MAX)
      {  id = cont%numpop;
         pop[id].clear();
         for(j=0;j<n;j++)
            pop[id].push_back(GAP->sol[j]);
         printIntArray(&pop[id][0],n);
         cout << "cont= " << cont << "(" << id << ") zsol = " << zsol << endl;
         cont++;
      }
   } while (next_permutation(&indCost[0], &indCost[0]+nind) && cont < 20*numpop); 

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

   // ----------------------------------------- subset generatino (parent set edfinition)
   // (here, simply randomly select subsets of three different cost solutions)

   tuple <int,int> parent; // id in the population and cost of each parent
   k = rand() % numpop;
   cont = 0;
   vector<tuple<int,int>> parents;
   vector<int> lstdiff;
   int idsol2,zsol2,zson;

   // all pairs of different parents
   for(k1 = 0;k1<numpop-1;k1++)
      for(k2 = k1;k2<numpop;k2++)
      {  zsol = 0;
         idsol = (k+k1)%numpop;
         for(j=0;j<n;j++) zsol += c[pop[idsol][j]][j];

         zsol2 = 0;
         idsol2 = (k+k2)%numpop;
         for(j=0;j<n;j++) zsol2 += c[pop[idsol2][j]][j];

         if(zsol != zsol2)
         {  lstdiff = listDiff(pop[idsol],pop[idsol2]);
            if(lstdiff.size() > 3) continue;
            parent = make_tuple(idsol,zsol);
            parents.push_back(parent);
            parent = make_tuple(idsol2,zsol2);
            parents.push_back(parent);

            //int myints[] = {0,0,1,1,2,2,2,0};
            //pop[idsol].assign (myints,myints+8);   // assigning from array.
            //int myints2[] = {0,0,2,1,1,2,0,2};
            //zsol = GAP->checkSol(pop[idsol]);
            //zsol2 = GAP->checkSol(pop[idsol2]);
            //pop[idsol2].assign (myints2,myints2+8);   // assigning from array.
            zson = pathSon(pop[idsol],pop[idsol2]);
            cout << idsol << "-" << idsol2 << " Numdiff= " << lstdiff.size();
            cout << " z1=" << zsol << " z2=" << zsol2 << " zson=" << zson;
            cout << (zson != zsol && zsol2 != zson ? "*************" : "<>") << endl;
         }
      }

   zcheck = 0;
   zcheck = GAP->checkSol(solbest);
   if (abs(zcheck - zub) > GAP->EPS)
   {  cout << "[SS] Detected inconsistency" << endl;
      z = INT_MAX;
   }
   else
      cout << "zcheck  " << zcheck << endl;
   cout << "SS zub= " << zub << endl;

   cout << "Scatter search algorithm incomplete !!" << endl;
   return zub;
}

// recombines two parents
int ScatterSearch::pathSon(vector<int> s1, vector<int> s2)
{  int res,i,j, cnt=0;
   double lb;

   // Lower bound computation, linear bound
   int numRows, numCols, numNZrow, statusMIP;

   numRows = n + m;   // num of constraints
   numCols = n * m;   // num of variables
   numNZrow = n * m;  // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows, numCols, numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(false); // verbose output
   statusMIP = CPXsetintparam(CPX->env, CPXPARAM_ScreenOutput, CPX_OFF);

   cnt = n - listDiff(s1,s2).size();   // number of equal assignments
   vector<int> indices(cnt);
   vector<char> lu(cnt);
   vector<double> bd(cnt);

   cnt=0;
   for (j = 0; j < n; j++)   // common attributes
      if(s1[j]==s2[j])
      {  indices[cnt] = s1[j] * n + j;
         lu[cnt] = 'B';
         bd[cnt] = 1.0;
         cnt++;
      }

   statusMIP = CPXchgbds(CPX->env, CPX->lp, cnt, &indices[0], &lu[0], &bd[0]);
   statusMIP = CPX->solveMIP(true, true);   // LP of free instance
   if (statusMIP == 0)
   {  lb = CPX->objval;
      cout << "Linear bound: " << lb << endl;
   }
   else
      lb = INT_MAX-1;

   if(CPX != NULL) delete CPX;
   res = (int)lb+0.99;
   return res;
}

// list of different assignments in two solutions
vector<int> ScatterSearch::listDiff(vector<int> s1, vector<int> s2)
{  int j,num;
   vector<int> lstdiff;

   num = 0;
   for(j=0;j<n;j++)
      if(s1[j] != s2[j])
      {  lstdiff.push_back(j);
         num++;
      }

   return lstdiff;
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


