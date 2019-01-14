#include "MetaLocalSearch.h"
#include "MIPCplex.h"      // per lineare in VNS
#include <ilcplex/cplex.h> // per lineare in VNS

MetaLocalSearch::MetaLocalSearch(GeneralizedAssignemnt* GAPinstance, LocalSearch* LSearch, int & zz) : zub(zz)
{  //ctor
   GAP = GAPinstance;
   LS  = LSearch;
   LB = new LowerBound(GAP,zub);

   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

MetaLocalSearch::~MetaLocalSearch()
{  delete LB;
}

int MetaLocalSearch::iteratedLocSearch(int** c, int maxIter, double alpha)
{  int i,iter;
   double z,zorg,z2=INT32_MAX;

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

   if (GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }

   iter = 0;

   while(iter<maxIter)
   {  z = GRASPcontruct(candNum,true);    // CHOICE WHETHER MATHEURISTC
      if(z < zub)
      {  GAP->storeBest(sol,z);
         cout << "[GRASP]: New zub: " << zub << " iter " << iter << endl;
      }

      if(z < INT_MAX)
      {
         //z = LS->opt10(GAP->c, true);
         //if(z < zub)
         //{  GAP->storeBest(sol,z);
         //   cout << "[GRASP]: opt10 new zub: " << zub << " iter " << iter << endl;
         //}
         z = LS->opt11(GAP->c, true);
         if(z < zub)
         {  GAP->storeBest(sol,z);
            cout << "[GRASP]: opt11 new zub: " << zub << " iter " << iter << endl;
         }
      }
      if(iter%200 == 0)
         cout << "[GRASP] iter "<< iter <<" z " << z << " zub "<< zub << endl;
      iter++;
   }
   return zub;
}

// construction GRASP, with candidate list
double MetaLocalSearch::GRASPcontruct(int candNum, bool isMatheuristic)
{  int i,ii,j,jj,m,n,icand;
   double z,lb;

   m = GAP->m;
   n = GAP->n;
   z = 0;
   vector<int> capleft(m),indReq(m);
   vector<int> regrets(n),indRegr(n);
   vector<double> cost(m);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };           // ASC order
   auto compRegr = [&regrets](int a, int b){ return regrets[a] > regrets[b]; };  // DESC order

   int** iterReq = new int*[m]; // requests for the partial solution
   for(i = 0; i < m; ++i)
      iterReq[i] = new int[n];

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   computeRegrets(GAP->c,n,m,regrets);

   for(j=0;j<n;j++) indRegr[j] = j; // sort clients by decreasing regrets
   std::sort(indRegr.begin(), indRegr.end(), compRegr);

   for (int j = 0; j < n; j++) GAP->sol[j] = -1; // no need for a seed solution

   for(jj=0;jj<n;jj++)
   {  
      for(j=0;j<n;j++)
         if(GAP->sol[j] >=0)  // force sol[j] into the bound
         {  for(i=0;i<m;i++)
               if(i != GAP->sol[j])
                  iterReq[i][j] = INT_MAX;
         }
         else
            for(i=0;i<m;i++)
               iterReq[i][j] = req[i][j];
   
      j = indRegr[jj];                // client order by regrets
      for(i=0;i<m;i++)
      {  indReq[i] = i;
         if(isMatheuristic)
         {  for(ii=0;ii<m;ii++)
               if(ii != i)
                  iterReq[ii][j] = INT_MAX;
            cost[i] = LB->linearBound(GAP->req,n,m,iterReq,GAP->cap);
            for(ii=0;ii<m;ii++)
               iterReq[ii][j] = req[ii][j];
         }
         else
            cost[i] = GAP->aversion(i,j); // aversion of client j for each server
      }

      std::sort(indReq.begin(), indReq.end(), compCost);  // sort servers by increasing aversion

      icand=rand() % candNum;
      ii=0;
      while(ii<m)                   // going best to worse
      {  i=indReq[(ii+icand)%m];
         if(capleft[i] >= GAP->req[i][j])
         {  GAP->sol[j] = i;
            capleft[i] -= GAP->req[i][j];
            z += GAP->c[i][j];
            break;
         }
         ii++;
      }
      if(ii==m)
      {  //cout << "[GRASPConstruct] Unable to construct feasible sol. ii=" << ii << endl;
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

   //if(z < zub)
   //{  GAP->storeBest(sol,z);
   //   cout << "[GRASPConstruct]: New zub: " << zub << endl;
   //}
 
end: ;
   // free
   for(i = 0; i < m; ++i)
      delete [] iterReq[i];
   delete [] iterReq;

   return z;
}

// VNS, vicinanze 1-0 e 1-1
double MetaLocalSearch::VNS(int maxIter, bool isMatheuristic)
{
   double z1, z2;
   int iter;

   if (GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }
   
   // linear bound
   z1 = LB->linearBound(GAP->c, n, m, GAP->req, GAP->cap);

   iter = 0;
   while (iter < maxIter)
   {
loop: z1 = LS->opt10(GAP->c, true);
      if (z1 < zub)
      {  GAP->storeBest(sol, z1);
         cout << "[VNS]: New zub: " << zub << " iter " << iter << endl;
      }
      z2 = LS->opt11(GAP->c, true);
      if (z2 < zub)
      {  GAP->storeBest(sol, z2);
         cout << "[VNS]: New zub: " << zub << " iter " << iter << endl;
      }
      if (z2 < z1)
         goto loop;
      else
      {  LS->neigh21();
         iter++;
      }
      if(iter%100 == 0) cout << "[VNS] iter " << iter << " zub " << zub << endl;
   }

   return zub;
}

