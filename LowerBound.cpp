#include "LowerBound.h"
#include "LocalSearch.h"
#include "Benders.h"

LowerBound::LowerBound(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

LowerBound::~LowerBound()
{
   //dtor
}

// each at its nearest facility
double LowerBound::trivialBound()
{  int i,j,m,n;
   double lb=0,minc,mini;

   m = GAP->m;
   n = GAP->n;

   for(j=0;j<n;j++)
   {  minc = DBL_MAX;
      for(i=0;i<m;i++)
         if(GAP->c[i][j] < minc)
         {  minc = GAP->c[i][j];
            mini = i;
         }
      lb += minc;
   }

   cout << "Trivial bound: "<< lb << endl;
   return lb;
}

double LowerBound::linearBound()
{  double lb;
   int numRows,numCols,numNZrow;
   int i,j,n,m;
   numRows = GAP->n+GAP->m;   // num of constraints
   numCols = GAP->n*GAP->m;   // num of variables
   numNZrow= GAP->n*GAP->m;   // max num of nonzero values in a row
   n = GAP->n;
   m = GAP->m;
   CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(true);    // verbose output
   int statusMIP = CPX->solveMIP(false,false);
   if(statusMIP == 0)
   {  lb = CPX->objval;
      cout << "Linear bound: "<< lb << endl;
      if(n*m < 101)                       // print LP if instance is small
      {  cout << " - Solution: " << endl; 
         for (i = 0; i < m; ++i)
         {  cout << "   " << i << ": ";
            for (j = 0; j < n; ++j)
               cout << CPX->x[i*n+j] << "\t";
            cout << endl;
         }
      }
   }
   else
      lb = DBL_MAX;
   CPX->freeMIP();

   delete CPX;
   return lb;
}

double LowerBound::linearBound(int** c, int n, int m, int** req, int* cap)
{  double lb;
   int numRows,numCols,numNZrow;
   int i,j;
   numRows = n+m;   // num of constraints
   numCols = n*m;   // num of variables
   numNZrow= n*m;   // max num of nonzero values in a row
   CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(c, n, m, req, cap, true);    // linear bound, verbose output
   int statusMIP = CPX->solveMIP(false,false);
   if(statusMIP == 0)
   {  lb = CPX->objval;
      cout << "Linear bound: "<< lb << endl;
      if(n*m < 101)                       // print LP if instance is small
      {  cout << " - Solution: " << endl; 
         for (i = 0; i < m; ++i)
         {  cout << "   " << i << ": ";
            for (j = 0; j < n; ++j)
               cout << CPX->x[i*n+j] << "\t";
            cout << endl;
         }
      }
   }
   else
      lb = DBL_MAX;
   CPX->freeMIP();

   delete CPX;
   return lb;
}

double LowerBound::lagrangianDecomposition(int** c, double alpha, double alphastep, double minAlpha, int innerIter, int maxiter)
{  int     i,j,iter=0,zcurr;
   double  zlb,step=0,zlbBest=0,fakeZub,sumSubGrad2;
   int*    lbsol   = new int[n];

   vector< vector<short> > x,y;
   vector< vector<double> > lambda,subgrad;
   //from empty 2D-matrix of size (0,0) to m*n, initialized with 0s:
   x.resize( m , vector<short>( n , 0 ) );
   y.resize( m , vector<short>( n , 0 ) );
   lambda.resize( m , vector<double>( n , 0 ) );
   subgrad.resize( m , vector<double>( n , 0 ) );

   zcurr     = INT_MAX;       // initial upper bound
   zlb       = DBL_MIN;

   iter  = 0;
   while(alpha > minAlpha && iter < maxiter)
   {  lbsol = subproblem_LD(c, &zlb, &zlbBest, zub, lambda, subgrad, x, y);
      zcurr = GAP->checkSol(lbsol);

      if(zcurr == zlbBest || (zub-zlbBest) < 1.0)                       // -------------------------- Optimum found 
      {  cout << "[lagrCap] Found the optimum!!! zopt="<< zub << " zlb=" << zlbBest<<endl;
         for(i=0;i<n;i++) solbest[i]=sol[i];
         zub = zcurr;
         goto lend;
      }
      else                                       // -------------------------- Heuristic block
      {  if(zcurr == INT_MAX)
         GAP->fixSolViaKnap(lbsol, &zcurr);   // hope to fix infeasibilities
                                              //GAP->fixSol(lbsol, &zcurr);       
         if(zcurr < 10*zlb)
         {  for(int j=0;j<n;j++) sol[j] = lbsol[j];
            LocalSearch* LS = new LocalSearch(GAP, GAP->zub);
            LS->opt10(c, true);
            delete LS;        
         }
         if(zcurr<zub)
         {  cout << "[lagrCap] -------- zub improved! " << zub;
            for(i=0;i<n;i++) solbest[i]=sol[i];
            zub = zcurr;
         }
      }

      // -------------------------- calcolo passo
      sumSubGrad2 = 0;
      for(i=0;i<m;i++)
         for(j=0;j<n;j++)
            sumSubGrad2 += subgrad[i][j]*subgrad[i][j];
      fakeZub = min((double) zcurr,1.2*zlb);
      fakeZub = max(fakeZub,zlb+1);
      step = alpha*(fakeZub-zlb)/sumSubGrad2;

      for(j=0;j<n;j++)                            // -------------------------- penalty update
         for(i=0;i<m;i++)
            lambda[i][j] += step*subgrad[i][j];
      iter++;
      if(iter % innerIter == 0)
      alpha = alphastep*alpha;

      if(iter%100 == 0)                           // -------------------------- logging
         cout << "[lagrCap] iter="<<iter<<" zub="<<zub<<" zlb="<<zlbBest<<" zcurr="<<zcurr<<endl;
   }

   lend:    
   if (lbsol!=NULL) delete(lbsol);
   return zlbBest;
}

// fills servers up to their capacities
int* LowerBound::subproblem_LD(int** c, double *zlb, double *zlbBest, int zub, 
                               vector<vector<double>> &lambda, 
                               vector<vector<double>> &subgrad, 
                               vector<vector<short>> &x,
                               vector<vector<short>> &y)
{  int i,j,mini;
   double mincost,alpha=0.5,beta=0.5;
   int* Q   = new int[n];
   int* sol = new int[n];
   int* Ksol = new int[n];
   double* val = new double[n];

   *zlb = 0;
   for(j=0;j<n;j++)
   {  sol[j] = -1;
      for(i=0;i<m;i++) 
      {  x[i][j] = y[i][j] = 0;
         subgrad[i][j] = 0; // check
      }
   }

   // solve assignment subproblem
   for(j=0;j<n;j++)
   {  mincost = DBL_MAX;
      mini    = INT_MAX;
      for(i=0;i<m;i++)                    // finds the minimum cost assignment
         if( (alpha*c[i][j]+lambda[i][j]) < mincost)
         {  mincost = alpha*c[i][j]+lambda[i][j];
            mini    = i;
         }
      x[mini][j] = 1;
      *zlb += mincost;
   }

   // solve capacity subproblem
   for(i=0;i<m;i++)
   {  for(j=0;j<n;j++)
      {  Q[j]    = req[i][j];
         val[j]  = -(beta*c[i][j]-lambda[i][j]);  // inverted sign to make it minimzation
         Ksol[j] = 0;
      }
      *zlb -= KDynRecur(n,GAP->cap[i],Q,val,Ksol);  // minus because it is minimization
      for(j=0;j<n;j++)
         y[i][j] = (short)Ksol[j];              
   }

   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
      {  subgrad[i][j] = x[i][j] - y[i][j];
         if(y[i][j]==1 && x[i][j]==1)
            sol[j]=i;   // heuristic
      }

   assert(*zlb<=zub);   // aborts if lb > ub
   if(*zlb>*zlbBest) 
      *zlbBest = *zlb;

   delete(Q);
   delete(val);
   delete(Ksol);
   return sol;
}

double LowerBound::benders()
{
   Benders* BEND = new Benders(GAP, GAP->zub);
   int res = BEND->bendersHeu(GAP->c,
      1,
      (GAP->conf->isVerbose ? true : false));
   if (BEND != NULL) delete BEND;
   BEND = NULL;
   return res;
}
