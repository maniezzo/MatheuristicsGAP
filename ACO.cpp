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

int ACO::antColony(int** c, int maxiter)
{  int z=0;
   int i,iter;
   vector<double> fitness(3);
   vector<int*> pop(3),newpop(3);

   Rins* RINS = new Rins(GAP,GAP->zub);   // for crossover
   for(i=0;i<3;i++)
      newpop[i] = new int[n];

   cout << "Population initialized, starting search" << endl;
   for (iter = 0; iter < maxiter; iter++)
   {  if(iter%1 == 0)
         cout << "====== iteration " << iter << " zub = " << zub << endl;
      
   }

   double zcheck = 0;
   zcheck = GAP->checkSol(solbest);
   if (abs(zcheck - zub) > GAP->EPS)
   {  cout << "[GA] Detected inconsistency" << endl;
      z = INT_MAX;
   }
   cout << "GA zub= " << zub << endl;

   if(RINS != NULL) delete RINS;
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


