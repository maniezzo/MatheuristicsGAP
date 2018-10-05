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
   int z,zorg,z2=INT32_MAX;

   VeryLarge* VLSN = new VeryLarge(GAP, GAP->zub);

   // allocate perturbed matrix
   int** cPert = (int**) malloc(m * sizeof(int *));
   for(i=0;i<m;i++)
      cPert[i] = (int*) malloc(n * sizeof(int));   

   // the algorithm
   iter = 0;
   while(iter<maxIter)
   {  zorg = zub;
      z = LS->opt11(c, true);
      // one iteration, fix k servers
      if(z<zorg)
         z2 = VLSN->verylarge(GAP->c, GAP->conf->verylargeConf->k, 1, false);
      dataPerturbation(c,cPert,alpha);
      if(iter%100 == 0)
         cout << "[ILS] iter "<<iter<< " zub "<<zub<<" z "<< z << " z2 " << z2 << endl;
      iter++;
   }

   // release perturbed matrix
   if (cPert!=NULL)
   {  for (i=0; i<m; i++) 
         free(cPert[i]);
      free(cPert);
   }

   // release matheuristic module
   if (VLSN != NULL) delete VLSN;
   VLSN = NULL;

   return zub;
}

void MetaLocalSearch::dataPerturbation(int** c,int** cPert, double alpha)
{  int i,j;
   double delta;

   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
      {  delta = alpha * (rand()/RAND_MAX);
         cPert[i][j] = (int) (round( (1 + delta - alpha/2) * c[i][j] ));
      }

   LS->opt11(cPert,false);
   if(GAP->checkSol(sol) == INT_MAX)
      cout << "[dataPerturbation] error" << endl;
}

double MetaLocalSearch::GRASP(int maxIter, int candNum)
{  double z;
   int iter;

   iter = 0;
   while(iter<maxIter)
   {  z = GRASPcontruct(candNum);
      if(z < zub)
      {  GAP->storeBest(sol,z);
         cout << "[GRASP]: New zub: " << zub;
      }

      if(z < INT_MAX)
      {
         z = LS->opt10(GAP->c, true);
         if(z < zub)
         {  GAP->storeBest(sol,z);
            cout << "[GRASP]: New zub: " << zub;
         }
      }
      //z = LS->opt11(GAP->c, true);
      //if(z < zub)
      //{  GAP->storeBest(sol,z);
      //   cout << "[GRASP]: New zub: " << zub;
      //}
      if(iter%200 == 0)
         cout << "[GRASP] iter "<< iter <<" z " << z << " zub "<< zub << endl;
      iter++;
   }
   return zub;
}

// construction GRASP, with candidate list
double MetaLocalSearch::GRASPcontruct(int candNum)
{  int i,ii,j,jj,m,n,icand;
   double z;

   m = GAP->m;
   n = GAP->n;
   z = 0;
   vector<int> cost(m),capleft(m),indReq(m);
   vector<int> regrets(n),indRegr(n);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };           // ASC order
   auto compRegr = [&regrets](int a, int b){ return regrets[a] > regrets[b]; };  // DESC order

   if(GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   computeRegrets(GAP->c,n,m,regrets);

   for(j=0;j<n;j++) indRegr[j] = j; // sort by decreasing regerets
   std::sort(indRegr.begin(), indRegr.end(), compRegr);

   for(jj=0;jj<n;jj++)
   {  j = indRegr[jj];              // client order by regrets
      for(i=0;i<m;i++)
      {  cost[i]= GAP->aversion(i,j);
         indReq[i] = i;
      }

      std::sort(indReq.begin(), indReq.end(), compCost);  // sort by increasing aversion

      icand=rand() % candNum;
      ii=0;
      while(ii<m)                   // going best to worse
      {  i=indReq[(ii+icand)%m];
         if(capleft[i]>=GAP->req[i][j])
         {  GAP->sol[j]=i;
            GAP->solbest[j]=i;
            capleft[i] -= GAP->req[i][j];
            z += GAP->c[i][j];
            break;
         }
         ii++;
      }
      if(ii==m)
      {  //cout << "[GRASPConstruct] Unable to construct feasible. ii="+ii << endl;
         z = INT_MAX;
         goto end;
      }
   }
   if(abs(GAP->checkSol(GAP->sol)-z) > GAP->EPS)
   {  cout << "[GRASPConstruct]: Error" << endl;
      z = INT_MAX;
   }
   //else
   //{  cout << "GRASP construction terminated. z = " << z << endl;
   //   printIntArray(sol,n);
   //}

   if(z < zub)
   {  GAP->storeBest(sol,z);
      cout << "[GRASPConstruct]: New zub: " << zub;
   }
end: return z;
}

