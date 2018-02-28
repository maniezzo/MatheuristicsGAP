#include "Ejection.h"

Ejection::Ejection(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

Ejection::~Ejection()
{
   //dtor
}

int Ejection::ejectionChain(int** c, int maxiter)
{  int z=0;
   int iter;
   double* x = new double[n*m];  // an LP solution
   int* sol;

   sol = generateOneSol(&z);

   double zlp = DBL_MAX;
   sol = feasibilityPump(1000,x,zlp);
   cout << "[EC] z = " << z << endl;

   delete x;
   return z;
}

int* Ejection::generateOneSol(int* zval)
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

int* Ejection::feasibilityPump(int maxIter, double* x, double zlp)
{  int i,j,iter,status;
   MIPCplex* CPX;
   vector<double> sigma(n*m);
   vector<int> indSigma(n*m);
   int*    indices  = new int[n*m];       // for cost coef update
   double* ofvalues = new double[n*m];    // for cost coef update
   double* xRound   = new double[n*m];
   double  z = DBL_MAX, delta = DBL_MAX;
   bool isDifferent;
   auto compSigma = [&sigma](int a, int b) { return sigma[a] > sigma[b]; };  // DESC order

   int numRows, numCols, numNZrow;
   numRows  = n+m;   // num of constraints
   numCols  = n*m;   // num of variables
   numNZrow = n*m;   // max num of nonzero values in a row
   try
   {
      CPX = new MIPCplex(numRows, numCols, numNZrow);
      CPX->GAP = GAP;
      CPX->allocateMIP(false);
      int cur_numcols = CPXgetnumcols(CPX->env, CPX->lp);
      int cur_numrows = CPXgetnumrows(CPX->env, CPX->lp);

      // set initial LP solution
      int* cstat = (int *)malloc(cur_numcols * sizeof(int));
      int* rstat = (int *)malloc(cur_numrows * sizeof(int));
      status = setInitialBasis(CPX,x,cstat,rstat,zlp);

      status = CPXcopybase(CPX->env, CPX->lp, cstat, rstat);
      if (status) 
      {  fprintf(stderr, "Failed to copy the basis.\n");
         goto lend;
      }

      // here we get the starting LP and rounded solutions
      status = CPX->solveMIP(false, false);
      if (status == 0)
      {
         double lb = CPX->objval;   // useless, but nice to know
         cout << "Linear bound: " << lb << endl;
         for (i = 0; i < m; ++i)
            for (j = 0; j < n; ++j)
            {  x[i*n + j] = CPX->x[i*n + j];
               xRound[i*n + j] = round(x[i*n + j]);
            }
      }
      cout << "Initial LP solution: "; printDblArray(x,n*m);

      for (iter = 0; (iter<maxIter) && (delta > GAP->EPS); iter++)
      {
         //printDblArray(x, n*m);
         //printDblArray(xRound, n*m);
         for (i = 0; i < m; ++i)          // redefinition of of cost coefficients
            for (j = 0; j < n; ++j)
            {  ofvalues[i*n +j] = (xRound[i*n+j] < GAP->EPS ? x[i*n+j] : (1-x[i*n+j]) );
               indices[i*n + j] = i*n + j;
            }
         status = CPXchgobj(CPX->env, CPX->lp, n*m, indices, ofvalues);
         int statusMIP = CPX->solveMIP(false, false);
         if (statusMIP == 0)
         {
            double cst = CPX->objval;     // useless, but nice to know
            for (i = 0; i < m; ++i)
               for (j = 0; j < n; ++j)
                  x[i*n + j] = CPX->x[i*n + j];
         }

         delta = 0;
         isDifferent = false;
         for (i = 0; i < m; ++i)
            for (j = 0; j < n; ++j)
            {  sigma[i*n + j] = abs(x[i*n + j] - xRound[i*n + j]); // could be also computed as above
               delta += sigma[i*n + j];
               if(round(x[i*n + j]) != xRound[i*n + j])
                  isDifferent = true;
            }
         cout << "FP iter: " << iter << " delta: " << delta << endl;

         if(delta>0)
         {  if(isDifferent && ((iter+1)%n > 0))     // heuristic perturbations added every tot iterations
               for (i = 0; i < m; ++i)
                  for (j = 0; j < n; ++j)
                     xRound[i*n + j] = round(x[i*n + j]);
            else
            {
               for (j = 0; j<n*m; j++) indSigma[j] = j;
               std::sort(indSigma.begin(), indSigma.end(), compSigma);

               int T = rand() % (n/4) + (n/2);   // T in the range n/4 - n/2, no frills heuristic 
               for(i=0;i<T;i++)
               {
                  //j = indSigma[i];  // client order by decreasing sigma
                  j = rand() % n;
                  xRound[j] = 1 - xRound[j];  // restart, very much heuristic
               }
            }
         }
      }

      free(cstat); cstat = NULL;
      free(rstat); rstat = NULL;
   }
   catch (std::exception const& e)
   {
      cout << "Error: " << e.what() << endl;
      goto lend;
   }

   if(delta <= GAP->EPS)
   {  cout << "Delta= " << delta << endl;
      z = 0;
      for (i = 0; i < m; ++i)
         for (j = 0; j < n; ++j)
            if(xRound[i*n + j] == 1)
            {  sol[j] = i;
               z += GAP->c[i][j];
            }
      cout << "Final LP      solution: "; printDblArray(x, n*m);
      cout << "Final rounded solution: "; printDblArray(x, n*m);

      zlp = z;
      if (abs(GAP->checkSol(sol) - z) > GAP->EPS)
      {
         cout << "[FP]: Error" << endl;
         zub = INT_MAX;
      }
      else
         cout << "[FP] Construction terminated. z = " << z << endl;
   }

lend:
   CPX->freeMIP();
   delete(CPX);
   delete(xRound);
   delete(ofvalues);
   delete(indices);
   return sol;
}

// starting LP solution for FP. There must be a better way to impement this but this is how far my cplex mastering goes
int Ejection::setInitialBasis(MIPCplex* CPX, double* x, int* cstat, int* rstat, double z)
{  int status;

   if(z==DBL_MAX)
   {
      status = CPXlpopt(CPX->env, CPX->lp);
      if (status) 
      {  fprintf(stderr, "Failed to optimize LP.\n");
         goto lend;
      }

      int solnstat = CPXgetstat(CPX->env, CPX->lp);

      if (solnstat == CPX_STAT_UNBOUNDED) 
      {  printf("Model is unbounded\n");
         goto lend;
      }
      else if (solnstat == CPX_STAT_INFEASIBLE) 
      {  printf("Model is infeasible\n");
         goto lend;
      }
      else if (solnstat == CPX_STAT_INForUNBD) 
      {  printf("Model is infeasible or unbounded\n");
         goto lend;
      }
      status = CPXgetbase(CPX->env, CPX->lp, cstat, rstat);
   }

lend:
   return status;
}