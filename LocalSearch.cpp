#include "LocalSearch.h"

LocalSearch::LocalSearch(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

LocalSearch::~LocalSearch()
{
   //dtor
}

// tries each client with each other facility
int LocalSearch::opt10(int** c)
{  int z=0;
   int i, isol, j;
   vector<int> capleft(m);

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
   {  capleft[sol[j]] -= req[sol[j]][j];
      z += c[sol[j]][j];
   }

l0:
   for (j = 0; j < n; j++)
   {
      isol = sol[j];
      for (i = 0; i < m; i++)
      {
         if (i == isol) continue;
         if (c[i][j] < c[isol][j] && capleft[i] >= req[i][j])
         {  // remove from isol and assign to i
            sol[j] = i;
            capleft[i]    -= req[i][j];
            capleft[isol] += req[isol][j];
            z -= (c[isol][j] - c[i][j]);
            if(z<zub)
            {  zub = z;
               for(int k=0;k<n;k++) solbest[k] = sol[k];
               cout << "[1-0 opt] new zub " << zub << endl;
            }
            goto l0;
         }
      }
   }

   return z;
}

// scambio assegnamento fra due clienti
double LocalSearch::opt11(int** c)
{  int i,j,j1,j2,temp,cap1,cap2;
   int delta, z=0, zcheck;
   vector<int> capleft(m);

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
   {  capleft[sol[j]] -= req[sol[j]][j];
      z += c[sol[j]][j];
   }
   zcheck = GAP->checkSol(sol);

l0:
   for(j1=0;j1<n;j1++)
   {  for(j2=j1+1;j2<n;j2++)
      {  delta = (c[sol[j1]][j1] + c[sol[j2]][j2]) - (c[sol[j1]][j2] + c[sol[j2]][j1]);
         if(delta > 0)
         {  cap1 = capleft[sol[j1]] + req[j1] - req[j2];
            cap2 = capleft[sol[j2]] + req[j2] - req[j1];
            if(cap1>=0 && cap2 >=0)
            {  capleft[sol[j1]] += req[j1] - req[j2];
               capleft[sol[j2]] += req[j2] - req[j1];
               temp    = sol[j1];
               sol[j1] = sol[j2];
               sol[j2] = temp;
               z -= delta;
               zcheck = GAP->checkSol(sol);
               if(abs(z-zcheck) > GAP->EPS)
                  cout << "[1-1] ohi" << endl;
               if(z<zub)
               {  zub = z;
                  cout << "[1-1 opt] new zub " << zub << endl;
               }
               goto l0;
            }
         }
      }
   }

   zcheck = 0;
   for(j=0;j<n;j++)
      zcheck += c[sol[j]][j];
   if(abs(zcheck - z) > GAP->EPS)
      cout << "[1.1opt] Ahi ahi" << endl;
   zcheck = GAP->checkSol(sol);
   return z;
}

// un vicino 21 a caso, scambio due (stesso deposito) vs 1 (altro deposito)
void LocalSearch::neigh21()
{  int i,i1,i2,j,j11,j12,j21,iter;
   vector<int> lst1, lst2;
   vector<int> capleft(m);
   int zcheck;


   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];
   for (j = 0; j < n; j++)
      capleft[sol[j]] -= req[sol[j]][j];

   i1 = rand()%m;
   i2 = rand()%m;
   if(i1==i2)
      i2 = (i2+1) % m;

   for(j=0;j<n;j++)
   {  if(sol[j]==i1)
         lst1.push_back(j);
      if(sol[j]==i2)
         lst2.push_back(j);
   }

   // 2 randomly chosen in i1 and one in i2
   iter = 0;
loop:    
   j11 = rand()%lst1.size();   // first indices then elements, in same variables!!
   j12 = rand()%lst1.size();
   j21 = rand()%lst2.size();
   if(j12==j11)
      j12=(j12+1)%lst1.size();
   j11 = lst1[j11];
   j12 = lst1[j12];
   j21 = lst2[j21];

   if( ((capleft[i1]+req[i1][j11]+req[i1][j12]-req[i1][j21]) >= 0) &&
       ((capleft[i2]-req[i2][j11]-req[i2][j12]+req[i2][j21]) >= 0) )
   {  sol[j11]=i2;
      sol[j12]=i2;
      sol[j21]=i1;
      zcheck = GAP->checkSol(sol);
      if(zcheck == INT_MAX)
         cout << "[2-1] ohi" << endl;
   }
   else  // try another, maybe feasible, neighbor
   {  iter++;
      if(iter < 50)
         goto loop;
   }

   return;
}
