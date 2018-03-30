#include "MetaLocalSearch.h"

MetaLocalSearch::MetaLocalSearch(GeneralizedAssignemnt* GAPinstance, LocalSearch* LSearch, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   LS  = LSearch;

   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

MetaLocalSearch::~MetaLocalSearch()
{
   //dtor
}

int MetaLocalSearch::iteratedLocSearch(int** c, int maxIter, double alpha)
{  int i,iter;
   int z;

   // allocate perturbed matrix
   int** cPert = (int**) malloc(m * sizeof(int *));
   for(i=0;i<m;i++)
      cPert[i] = (int*) malloc(n * sizeof(int));   

   // the algorithm
   iter = 0;
   while(iter<maxIter)
   {  z = LS->opt10(c);
      dataPerturbation(c,cPert,alpha);
      if(iter%1000 == 0)
         cout << "[ILS] iter "<<iter<< " zub "<<zub<<" z "<<z<<endl;
      iter++;
   }

   // release perturbed matrix
   if (cPert!=NULL)
   {  for (i=0; i<m; i++) 
         free(cPert[i]);
      free(cPert);
   }

   return zub;
}

void MetaLocalSearch::dataPerturbation(int** c,int** cPert, double alpha)
{  int i,j;
   double delta;

   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
      {  delta = alpha * (rand()/RAND_MAX);
         cPert[i][j] = (int) round( (1 + delta - alpha/2) * c[i][j] );
      }

   LS->opt10(cPert);
   if(GAP->checkSol(sol) == INT_MAX)
      cout << "[dataPerturbation] error" << endl;
}
