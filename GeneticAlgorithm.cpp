#include "GeneticAlgorithm.h"

GenAlgo::GenAlgo(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

GenAlgo::~GenAlgo()
{
   //dtor
}

int GenAlgo::genAlgo(int** c, int maxiter, int numpop,  double pc)
{  int z=0;
   int i,iter;
   vector<double> fitness(numpop);
   vector<int*> pop(numpop),newpop(numpop);

   if( !initPop(pop, numpop, fitness) )
      return NULL;

   Rins* RINS = new Rins(GAP,GAP->zub);   // for crossover
   for(i=0;i<numpop;i++)
      newpop[i] = new int[n];

   cout << "Population initialized, starting search" << endl;
   for (iter = 0; iter < maxiter; iter++)
   {  if(iter%1 == 0)
         cout << "====== iteration " << iter << " zub = " << zub << endl;
      
      for(i=0;i<numpop;i++)
         selCrossover(pop,numpop,fitness,RINS, newpop[i]);

      // this in case mutation fails
      for(i=0;i<numpop;i++)
         pop[i] = newpop[i];

      if( mutate(newpop,numpop) )
         for(i=0;i<numpop;i++)
            pop[i] = newpop[i];
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

// mutation, false if search stuck
bool GenAlgo::mutate(vector<int*>& pop, int numpop)
{  int i,j,z;
   double r;

   for(i=0;i<numpop;i++)      // for each individual
   {
      for(j=0;j<n;j++)        // for each assignment
      {  r = (double)rand() / RAND_MAX;
         if(r < 1.0 / n)      // mutation probability is 1 / n
            pop[i][j] = std::rand() % m;
      }
      GAP->fixSol(pop[i],&z);
      if(z == INT_MAX)
      {  // cout << "[mutation] cannot proceed, aborting .... " << endl;
         return false;        // could not generate a feasible individual
      }
   }

   return true;
}

// selects two parents and generates one offspring
void GenAlgo::selCrossover(vector<int*>& pop, int numpop, vector<double>& fitness, Rins* RINS, int* newsol)
{  double r;
   int p1,p2;

   p1 = montecarlo(fitness);
   p2 = montecarlo(fitness);

   cout << "Mating " << p1 << " and " << p2 << endl;

   r = (double)rand() / RAND_MAX;
   if(r < GAP->conf->GA->pc)
      crossOver(pop,p1,p2,RINS,newsol);
   else
      newsol = pop[p1];

   return;
}

int GenAlgo::montecarlo(vector<double>& v)
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

// initialization of population and fitness
bool GenAlgo::initPop(vector<int*>& pop, int numpop, vector<double>& fitness)
{
   int i,z;
   for (i=0; i<numpop; i++)
   {  pop[i] = new int[n];
      pop[i] = generateOneSol(&z);
      if(pop[i] == nullptr)
         return false;
      else
         fitness[i] = 1.0/z;

      printIntArray(pop[i],n);
   }
   return true;
}

int* GenAlgo::generateOneSol(int* zval)
{  int i,j,k,c1,c2,temp,numRetry=0;
   int* newsol = new int[n];
   int z = 0;                                // soluton cost
   vector<int> cost(m),capleft(m),indReq(m),ind(n);
   auto compCost = [&cost](int a, int b){ return cost[a] < cost[b]; };  // ASC order, predecate for sort

   for(i=0;i<m;i++) capleft[i] = GAP->cap[i];   // residual capacities
   for(j=0;j<n;j++) ind[j] = j;

l0:for(k=0;k<5*n;k++)
   {  c1 = std::rand() % n;
      c2 = std::rand() % n;
      temp = ind[c1];
      ind[c1]=ind[c2];
      ind[c2]= temp;
   }

   for(int jj=0;jj<n;jj++)
   {
      j = ind[jj];
      for(i=0;i<m;i++)
      {  cost[i]= GAP->req[i][j];   // generalized cost for ordering
         indReq[i] = i;             // server index for requests ordering
      }

      std::sort(indReq.begin(), indReq.end(), compCost); // sort by increasing gen cost
      k = std::rand() % 3;          // random among candate set of best 3
      newsol[j] = indReq[k];
      z += GAP->c[newsol[j]][j];
   }

   if(abs(GAP->checkSol(newsol)-z) > GAP->EPS)
   {  GAP->fixSol(newsol,&z);
      if(abs(GAP->checkSol(newsol)-z) > GAP->EPS)
      {  numRetry++;
         if(numRetry < 5)
            goto l0;
         else if (numRetry < 10)
         {  GAP->fixSolViaKnap(newsol,&z);
            if(abs(GAP->checkSol(newsol)-z) > GAP->EPS)
               goto l0;
         }
         else
         {
            cout << "[generateOneSol]: infeasible" << endl;
            newsol = nullptr;
         }
      }
   }
   cout << "[generateOneSol]: z = " << (*zval = z) << endl;

   return newsol;
}

// fixes common assignments and frees the rest
void GenAlgo::crossOver(vector<int*>& pop, int p1, int p2, Rins* RINS, int* newsol)
{  int j,num=0,res;

   for(j=0;j<n;j++)
      if(pop[p1][j] == pop[p2][j])
      {  newsol[j]   = pop[p1][j];
         num++;
      }
      else
         newsol[j] = -1;

   cout << "[GA] Xover fixes: ";
   printIntArray(newsol,n);

   res = RINS->dive( GAP->c,
                     GAP->conf->rinsConf->maxnodes,
                     INT_MAX,
                     newsol, false
                   );

   cout << "Xover, cost " << res << endl;
}

