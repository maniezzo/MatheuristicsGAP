#include "Lagrangian.h"
#include "LocalSearch.h"

Lagrangian::Lagrangian(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

Lagrangian::~Lagrangian()
{
   //dtor
}

// Lagrangian, feasible for the assignments, relaxes capacities
int Lagrangian::lagrAss(int** c, double alpha, double alphastep, double minAlpha, int innerIter, int maxiter)
{  int i,sumSubGrad2,iter=0,zcurr;
   double  zlb,step=0,zlbBest,fakeZub;

   double* lambda  = new double[m];
   int*    subgrad = new int[m];
   int*    lbsol   = new int[n];

   ofstream flog;
   flog.open ("lagr.log");
   flog << fixed << setprecision(3);

   zcurr     = INT_MAX;       // initial upper bound
   zlb       = DBL_MIN;
   zlbBest   = DBL_MIN;
   for(i=0;i<m;i++)
      lambda[i] = 0.0;

   iter  = 0;
   if(maxiter < 0) maxiter = INT_MAX;
   while(alpha > minAlpha && iter < maxiter)
   {  lbsol = subproblem_ass(c, &zlb, &zlbBest, zub, lambda, subgrad);
      zcurr = GAP->checkSol(lbsol);

      if(zcurr == zlbBest  || (zub-zlbBest) < 1.0) // -------------------------- Optimum found 
      {  cout << "[lagrAss] Found the optimum !!! zopt="<<zub<<" iter "<< iter << endl;
         for(i=0;i<n;i++) solbest[i]=sol[i];
         zub = zcurr;
         goto lend;
      }
      else                                       // -------------------------- Heuristic block
      {  if(zcurr == INT_MAX)
            //GAP->fixSolViaKnap(lbsol, &zcurr);   // hope to fix infeasibilities
            GAP->fixSol(lbsol, &zcurr);       
         if(zcurr < 10*zlb)
         {  for(int j=0;j<n;j++) sol[j] = lbsol[j];
            LocalSearch* LS = new LocalSearch(GAP, GAP->zub);
            LS->opt10(c);    
            delete LS;        
         }
         if(zcurr<zub)
         {  for(i=0;i<n;i++) solbest[i]=sol[i];
            zub = zcurr;
            cout << "[lagrAss] -------- zub improved! " << zub << endl;
         }
      }

      sumSubGrad2 = 0;                            // ------------------------ step redefinition
      for(i=0;i<m;i++)
         sumSubGrad2 += subgrad[i]*subgrad[i];
      fakeZub = min((double) zcurr,1.2*zlb);
      fakeZub = max(fakeZub,zlb+1);
      step = alpha*(fakeZub-zlb)/sumSubGrad2;

      for(i=0;i<m;i++)                            // -------------------------- penalty update
         lambda[i] = max(0.0,lambda[i]+step*subgrad[i]);
      iter++;
      if(iter % innerIter == 0)
         alpha = alphastep*alpha;

      if(iter%1000 == 0)                          // -------------------------- logging
      {  cout << "[lagrAss] iter="<<iter<<" zub="<<zub<<" zlb="<<zlbBest<<" zcurr="<<zcurr<<endl;
         writeIterData(flog, iter, zlb, zub, zcurr, alpha, lbsol, subgrad, lambda, step);
      }
   }

lend:    
   if (flog.is_open()) flog.close();
   if (lambda!=NULL)  free(lambda);
   if (subgrad!=NULL) free(subgrad);
   if (lbsol!=NULL)   free(lbsol);
   return zcurr;
}

// assigns each client to a server
int* Lagrangian::subproblem_ass(int** c, double *zlb, double *zlbBest, int zub, double* lambda, int* subgrad)
{  int i,j,mini;
   double mincost;
   int* sol = new int[n];

   for(j=0;j<n;j++) sol[j]     = -1;
   for(i=0;i<m;i++) subgrad[i] =  0;

   *zlb = 0;
   for(i=0;i<m;i++)
   {  subgrad[i] -= GAP->cap[i];
      *zlb -= GAP->cap[i]*lambda[i];       // penalties sum in the lagrangian function
   }

   for(j=0;j<n;j++)
   {  mincost = DBL_MAX;
      mini    = INT_MAX;
      for(i=0;i<m;i++)                    // finds the minimum cost assignment
         if( (c[i][j]+(req[i][j]*lambda[i])) < mincost)
         {  mincost = c[i][j]+(req[i][j]*lambda[i]);
            mini    = i;
         }
      sol[j] = mini;
      subgrad[mini] += req[mini][j];
      *zlb += mincost; 
   }
   assert(*zlb<=zub);  // aborts if lb > ub
   if(*zlb>*zlbBest) 
      *zlbBest = *zlb;
   return sol;
}

// just logging on file flog
void Lagrangian::writeIterData(ofstream& flog, int iter, double zlb, int zub, int zcurr, double alpha,
                               int* lbsol, int* subgrad, double* lambda, double step)
{  int i;

   flog << "iter "<< iter <<" zlb = "<< zlb <<" zub = "<< zub <<" zcurr = "<< zcurr <<" alpha = "<< alpha <<" \nlbsol ";
   for(i=0;i<n;i++)
      flog << " "<<lbsol[i];

   flog << "\nsubgr ";
   for(i=0;i<m;i++)
      flog << " "<<subgrad[i];

   flog << "\nlambda ";
   for(i=0;i<m;i++)
      flog << " "<<lambda[i];

   flog << "\nstep "<<step<< endl;
};

// Lagrangian, feasible for the capacities, relaxes assignments
int Lagrangian::lagrCap(int** c, double alpha, double alphastep, double minAlpha, int innerIter, int maxiter)
{  int i,j,sumSubGrad2,iter=0,zcurr;
   double  zlb,step=0,zlbBest=0,fakeZub;

   double* lambda  = new double[n];
   int*    subgrad = new int[n];
   int*    lbsol   = new int[n];

   ofstream flog;
   flog.open ("lagr.log");
   flog << fixed << setprecision(3);

   zcurr     = INT_MAX;       // initial upper bound
   zlb       = DBL_MIN;
   for(i=0;i<n;i++)
      lambda[i] = 0.0;

   iter  = 0;
   while(alpha > minAlpha && iter < maxiter)
   {  subproblem_cap(c, &zlb, &zlbBest, zub, lambda, subgrad, lbsol);
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
            LS->opt10(c);    
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
      for(j=0;j<n;j++)
         sumSubGrad2 += subgrad[j]*subgrad[j];
      fakeZub = min((double) zcurr,1.2*zlb);
      fakeZub = max(fakeZub,zlb+1);
      step = alpha*(fakeZub-zlb)/sumSubGrad2;

      for(j=0;j<n;j++)                            // -------------------------- penalty update
         lambda[j] += step*subgrad[j];
      iter++;
      if(iter % innerIter == 0)
         alpha = alphastep*alpha;

      if(iter%100 == 0)                           // -------------------------- logging
      {  cout << "[lagrCap] iter="<<iter<<" zub="<<zub<<" zlb="<<zlbBest<<" zcurr="<<zcurr<<endl;
         writeIterData(flog, iter, zlb, zub, zcurr, alpha, lbsol, subgrad, lambda, step);
      }
   }

lend:    
   if (flog.is_open()) flog.close();
   if (lambda!=NULL)  delete(lambda);
   if (subgrad!=NULL) delete(subgrad);
   if (lbsol!=NULL)   delete(lbsol);
   return zcurr;
}

// fills servers up to their capacities
void Lagrangian::subproblem_cap(int** c, double *zlb, double *zlbBest, int zub, double* lambda, int* subgrad, int* lbsol)
{  int i,j;

   int* Q   = new int[n];
   double* val = new double[n];
   int* Ksol = new int[n];

   *zlb = 0;
   for(j=0;j<n;j++)
   {  lbsol[j]   = -1;
      subgrad[j] = 1;
      *zlb += lambda[j];               // penalty sum in the lagrangian function
   }

   for(i=0;i<m;i++)
   {  for(j=0;j<n;j++)
      {  Q[j]    = req[i][j];
         val[j]  = -c[i][j]+lambda[j];  // inverted sign to make it minimzation
         Ksol[j] = 0;
      }
      *zlb -= KDynRecur(n,GAP->cap[i],Q,val,Ksol);  // minus because it is minimization
      for(j=0;j<n;j++)
         if(Ksol[j] > 0)
         {  lbsol[j] = i;              // could be a reassignment, obviously. This is heuristic
            subgrad[j] -= 1;
         }
   }

   assert(*zlb<=zub);  // aborts if lb > ub
   if(*zlb>*zlbBest) 
      *zlbBest = *zlb;

   free(Q);
   free(val);
   free(Ksol);
   return;
}
