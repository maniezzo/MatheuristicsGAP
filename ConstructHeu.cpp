#include "ConstructHeu.h"

ConstructHeu::ConstructHeu(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m   = GAP->m;
   n   = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

ConstructHeu::~ConstructHeu()
{
   //dtor
}

// constructive: each at its least disliked facility
int ConstructHeu::simpleConstruct()
{  int i,ii,j,jj,m,n;

   m = GAP->m;
   n = GAP->n;
   vector<int> cost(m),capleft(m),indReq(m);
   vector<int> regrets(n),indCost(n);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };           // ASC order
   auto compRegr = [&regrets](int a, int b){ return regrets[a] > regrets[b]; };  // DESC order

   if(GAP->n == NULL)
   {  cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }

   zub = INT_MAX;

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   computeRegrets(GAP->c,n,m,regrets);

   for(j=0;j<n;j++) indCost[j] = j;
   std::sort(indCost.begin(), indCost.end(), compRegr);

   zub = 0;
   for(jj=0;jj<n;jj++)
   {  j = indCost[jj];  // client order by decreasing regrets, higher first
      for(i=0;i<m;i++)
      {  cost[i]= GAP->aversion(i,j);
         indReq[i] = i;
      }

      std::sort(indReq.begin(), indReq.end(), compCost); // prefer low cost

      ii=0;
      while(ii<m)
      {  i=indReq[ii];
         if(capleft[i]>=GAP->req[i][j])
         {  GAP->sol[j]=i;
            GAP->solbest[j]=i;
            capleft[i] -= GAP->req[i][j];
            zub += GAP->c[i][j];
            break;
         }
         ii++;
      }
      if(ii==m)
      {  cout << "[SimpleConstruct] Error. ii="+ii << endl;
         zub = INT_MAX;
         break;
      }
   }
   if(abs(GAP->checkSol(GAP->sol)-zub) > GAP->EPS)
   {  cout << "[simpleContruct]: Error" << endl;
      zub = INT_MAX;
   }
   else
   {  cout << "Construction terminated. Zub = " << zub << endl;
      printIntArray(sol,n);
   }
   for(int k=0;k<n;k++) solbest[k] = sol[k];
   return zub;
}

// constructive: each at its least disliked facili, given the ordering
int ConstructHeu::constructWithInd(vector<int> indCost, bool isVerbose)
{  int i,ii,j,jj,m,n,z;

   m = GAP->m;
   n = GAP->n;
   vector<int> cost(m),capleft(m),indReq(m);
   vector<int> regrets(n);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };           // ASC order
   auto compRegr = [&regrets](int a, int b){ return regrets[a] > regrets[b]; };  // DESC order

   if(GAP->n == NULL)
   {  if(isVerbose) cout << "Instance undefined. Exiting" << endl;
      return INT_MAX;
   }

   z = INT_MAX;
   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];

   z = 0;
   for(jj=0;jj<n;jj++)
   {  j = indCost[jj];  // client order as received
      for(i=0;i<m;i++)
      {  cost[i]= GAP->aversion(i,j);
         indReq[i] = i;
      }

      std::sort(indReq.begin(), indReq.end(), compCost);

      ii=0;
      while(ii<m)
      {  i=indReq[ii];
         if(capleft[i]>=GAP->req[i][j])
         {  GAP->sol[j]=i;
            capleft[i] -= GAP->req[i][j];
            z += GAP->c[i][j];
            break;
         }
         ii++;
      }

      if(ii==m)
      {  if(isVerbose) cout << "[Construct] Error. ii="+ii << endl;
         z = INT_MAX;
         break;
      }
   }

   if(abs(GAP->checkSol(GAP->sol)-z) > GAP->EPS)
   {  if(isVerbose) cout << "[Contruct]: Error" << endl;
      z = INT_MAX;
   }
   else
   {  if(isVerbose) cout << "Construction terminated. z = " << z << endl;
      if(isVerbose) printIntArray(sol,n);

      if(z<zub)
      {  zub = z;
         for(int k=0;k<n;k++) solbest[k] = sol[k];
      }
   }

   return z;
}
