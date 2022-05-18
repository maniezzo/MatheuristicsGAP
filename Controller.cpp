#include "Controller.h"

Controller::Controller()
{
   GAP = new GeneralizedAssignemnt();
   P  = new Persistence(); P->GAP = GAP;

   //srand((unsigned)time(NULL));
   srand(550);

   P->loadConfig();
   readJSONdata(GAP->conf->datafile);
}

Controller::~Controller()
{
   //dtor
   delete GAP->sol;
   delete GAP->solbest;
   delete GAP;
   delete P;
   cout << "GAP deleted";
}

// legge i dati da file json
void Controller::readJSONdata(string filename)
{
   int res = P->readJSONdata(filename);
   if(res > 0)
   {  GAP->sol     = new int[GAP->n];
      GAP->solbest = new int[GAP->n];
      GAP->zub     = INT_MAX;
   }
   else
      cout << "something wrong with input files" << endl;
}

// simple constructive heuristic
double Controller::simpleConstruct()
{  int res;

   CH   = new ConstructHeu(GAP,GAP->zub); 
   res = CH->simpleConstruct();
   if(CH != NULL) delete CH;
   CH = NULL;
   return res;
}

// Martello and Toth heuristic
int Controller::run_mthg()
{  int res;
   MT   = new MTHG(GAP,GAP->zub);
   res = MT->run_mthg();
   if(MT != NULL) delete MT;
   return res;
}

// lower bounds
void Controller::computeBounds()
{  double bound;
   LB   = new LowerBound(GAP,GAP->zub);  
   bound = LB->trivialBound();
   cout << "Trivial bound " << bound << endl;
   bound = LB->linearBound();
   cout << "Linear bound " << bound << endl;
   bound = LB->lagrangianDecomposition(GAP->c,
                                       GAP->conf->lagrAss->alpha, 
                                       GAP->conf->lagrAss->alphastep, 
                                       GAP->conf->lagrAss->minalpha,
                                       GAP->conf->lagrAss->innerIter, 
                                       GAP->conf->lagrAss->maxiter);
   cout << "Lagrangian Decomposition bound " << bound << endl;
   cout << "Benders bound to be completed" << endl;
   //bound = LB->benders();
   //cout << "Benders Decomposition bound " << bound << endl;

   if(LB != NULL) delete LB;
   LB = NULL;
}

// exact optimum, MIP
double Controller::exactCplex()
{
   double optMIP;
   int statusMIP;
   int numRows,numCols,numNZrow;
   numRows = GAP->n+GAP->m;   // num of constraints
   numCols = GAP->n*GAP->m;   // num of variables
   numNZrow= GAP->n*GAP->m;   // max num of nonzero values in a row
   MIPCplex* CPX = new MIPCplex(numRows,numCols,numNZrow);
   CPX->GAP = GAP;
   CPX->allocateMIP(true);    // verbose output
   statusMIP = CPX->solveMIP(true,false);
   if(statusMIP == 0)
   {  optMIP = CPX->objval;
      cout << "Optimal cost: "<< optMIP << endl;
   }
   else
      optMIP = DBL_MAX;
   CPX->freeMIP();

   delete CPX;
   return optMIP;
}

// simple local search
int Controller::opt10()
{  int z,zcheck;
   LS = new LocalSearch(GAP, GAP->zub);

   // no solution to improve
   if(GAP->zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }

   z = LS->opt10(GAP->c,true);
   zcheck = GAP->checkSol(GAP->sol);
   if (abs(zcheck - z) > GAP->EPS)
   {  cout << "[1.0opt] Ahi ahi" << endl;
      z = INT_MAX;
   }

   z = (int) LS->opt11(GAP->c,true);
   zcheck = GAP->checkSol(GAP->sol);
   if (abs(zcheck - z) > GAP->EPS)
   {
      cout << "[1.0opt] Ahi ahi" << endl;
      z = INT_MAX;
   }

   if (LS != NULL) delete LS;
   LS = NULL;

   cout << "Opt01 zub= " << z << endl;
   return GAP->zub;
}

// simulated annealing
int Controller::run_simAnn()
{
   SA = new SimAnnealing(GAP, GAP->zub);

   SA->simAnnealing( GAP->c,
      GAP->conf->SA->maxT,
      GAP->conf->SA->k,
      GAP->conf->SA->alpha,
      GAP->conf->SA->maxiter,
      GAP->conf->SA->iterAnneal,
      GAP->conf->SA->coefAnneal,
      GAP->conf->SA->fGoMath
   );

   if (SA != NULL) delete SA;
   SA = NULL;

   return GAP->zub;
}

// tabu search
int Controller::run_tabuSearch()
{
   TS = new TabuSearch(GAP, GAP->zub);

   TS->tabuSearch(GAP->c,
      GAP->conf->TS->Ttenure,
      GAP->conf->TS->maxiter,
      GAP->conf->TS->fGoMath
   );
   if (TS != NULL) delete TS;
   TS = NULL;

   return GAP->zub;
}

// genetic algorithm
int Controller::run_genAlgo()
{
   GA = new GenAlgo(GAP, GAP->zub);
   int res = GA->genAlgo(GAP->c,
      GAP->conf->GA->maxiter,
      GAP->conf->GA->numpop,
      GAP->conf->GA->pc
   );
   if (GA != NULL) delete GA;
   GA = NULL;
   return res;
}

// ant colony optimization (ANTS)
int Controller::run_ACO()
{
   ANTS = new ACO(GAP, GAP->zub);
   int res = ANTS->antColony(GAP->c,
      GAP->conf->ACO->maxiter,
      GAP->conf->ACO->numpop,
      GAP->conf->ACO->alpha
      );
   if (ANTS != NULL) delete ANTS;
   ANTS = NULL;
   return res;
}

// Scatter search
int Controller::run_SS()
{
   SS = new ScatterSearch(GAP, GAP->zub);
   int res = SS->go_scatter(GAP->c,
      GAP->conf->SS->maxiter,
      GAP->conf->SS->numpop,
      GAP->conf->SS->alpha
   );
   if (SS != NULL) delete SS;
   SS = NULL;
   return res;
}

// ejection chain
int Controller::run_ejection()
{
   EC = new Ejection(GAP, GAP->zub);
   int res = EC->ejectionChain(GAP->c,
                         GAP->conf->EC->maxiter
                        );
   if (EC != NULL) delete EC;
   EC = NULL;
   return res;
}

// iterated local search
int Controller::run_iteratedLS()
{  
   if(GAP->zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }
   else
   {
      LS = new LocalSearch(GAP, GAP->zub);
      MLS = new MetaLocalSearch(GAP, LS, GAP->zub);
      MLS->iteratedLocSearch(GAP->c,
                                  GAP->conf->IterLS->maxiter,
                                  GAP->conf->IterLS->alpha);
      if (LS != NULL) delete LS;
      LS = NULL;
      if (MLS != NULL) delete MLS;
      MLS = NULL;

      return GAP->zub;
   }
}

// GRASP
int Controller::run_GRASP()
{  
   if(GAP->zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }
   else
   {
      LS = new LocalSearch(GAP, GAP->zub);
      MLS = new MetaLocalSearch(GAP, LS, GAP->zub);
      MLS->GRASP(GAP->conf->GRASP->maxiter,
                 GAP->conf->GRASP->candNum);
      if (LS != NULL) delete LS;
      LS = NULL;
      if (MLS != NULL) delete MLS;
      MLS = NULL;

      return GAP->zub;
   }
}

// GRASP
int Controller::run_VNS()
{
   if (GAP->zub == INT_MAX)
   {
      cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }
   else
   {
      LS = new LocalSearch(GAP, GAP->zub);
      MLS = new MetaLocalSearch(GAP, LS, GAP->zub);
      MLS->VNS(GAP->conf->VNS->maxiter, GAP->conf->VNS->fGoMath);
      if (LS != NULL) delete LS;
      LS = NULL;
      if (MLS != NULL) delete MLS;
      MLS = NULL;

      return GAP->zub;
   }
}

// lagrangean heuristic, feasible assignments, relax capacities
int Controller::run_lagrAss()
{  
   LAGR = new Lagrangian(GAP,GAP->zub); 
   int res = LAGR->lagrAss(GAP->c,
                           GAP->conf->lagrAss->alpha, 
                           GAP->conf->lagrAss->alphastep, 
                           GAP->conf->lagrAss->minalpha,
                           GAP->conf->lagrAss->innerIter, 
                           GAP->conf->lagrAss->maxiter);
   if(LAGR != NULL) delete LAGR;
   LAGR = NULL;
   return res;
}

// lagrangean heuristic, feasible capacities, relax assignments
int Controller::run_lagrCap()
{  
   LAGR = new Lagrangian(GAP,GAP->zub); 
   int res = LAGR->lagrCap(GAP->c,
                        GAP->conf->lagrCap->alpha, 
                        GAP->conf->lagrCap->alphastep, 
                        GAP->conf->lagrAss->minalpha,
                        GAP->conf->lagrCap->innerIter, 
                        GAP->conf->lagrCap->maxiter);
   if(LAGR != NULL) delete LAGR;
   LAGR = NULL;
   return res;
}

// RINS
int Controller::run_rins()
{  int res; 
   RINS = new Rins(GAP,GAP->zub); 
   if(GAP->zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }
   else
      res = RINS->dive( GAP->c,
                        GAP->conf->rinsConf->maxnodes,
                        GAP->zub,
                        GAP->sol,
                        true
                      );
   if(RINS != NULL) delete RINS;
   RINS = NULL;
   return res;
}

// beam search
int Controller::run_beam()
{  
   BS = new BeamSearch(GAP,GAP->zub); 
   int res = BS->beamSearch(GAP->c,
      GAP->conf->beam->delta,
      GAP->conf->beam->maxnodes,
      (GAP->conf->isVerbose ? true : false) );
   if(BS != NULL) delete BS;
   BS = NULL;
   return res;
}

// forward and backward
int Controller::run_FandB()
{  
   FB = new FandB(GAP,GAP->zub); 
   int res = FB->forwardBackward(GAP->c,
      GAP->conf->fbConf->delta,
      GAP->conf->fbConf->maxnodes,
      (GAP->conf->isVerbose ? true : false) );
   if(FB != NULL) delete FB;
   FB = NULL;
   return res;
}

// corridor
int Controller::run_Corridor()
{  
   CORR = new Corridor(GAP,GAP->zub); 
   int res = CORR->solveByCorridor (GAP->c,
                     GAP->conf->corridorConf->delta,
                     GAP->conf->corridorConf->maxiter,
                     (GAP->conf->isVerbose ? true : false) );
   if(CORR != NULL) delete CORR;
   CORR = NULL;
   return res;
}

// local branching
int Controller::run_localBranching()
{  
   LBR = new LocBranching(GAP,GAP->zub); 
   int res = LBR->localBranching(GAP->c,
                     GAP->conf->locBranching->k,
                     GAP->conf->locBranching->maxiter,
                     (GAP->conf->isVerbose ? true : false) );
   if(LBR != NULL) delete LBR;
   LBR = NULL;
   return res;
}

// benders heuristic
int Controller::run_benders()
{
   BEND = new Benders(GAP, GAP->zub);
   int res = BEND->bendersHeu(GAP->c,
      GAP->conf->bendersConf->maxiter,
      (GAP->conf->isVerbose ? true : false));
   if (BEND != NULL) delete BEND;
   BEND = NULL;
   return res;
}

// kernel search
int Controller::run_kernel()
{
   KER = new Kernel(GAP, GAP->zub);
   int numVarInBucket = 6;  // there should be a way to compute this
   int res = KER->solveByKernel( (GAP->conf->isVerbose ? true : false),numVarInBucket);
   if (KER != NULL) delete KER;
   KER = NULL;
   return res;
}

// veri large scale neighborhood search
int Controller::run_VLSN()
{
   VLSN = new VeryLarge(GAP, GAP->zub);
   int res = VLSN->verylarge(GAP->c, 
      GAP->conf->verylargeConf->k,
      GAP->conf->verylargeConf->maxiter, 
      (GAP->conf->isVerbose ? true : false)
   );
   if (VLSN != NULL) delete VLSN;
   VLSN = NULL;
   return res;
}

