﻿#include "Controller.h"

Controller::Controller()
{
   GAP = new GeneralizedAssignemnt();
   P  = new Persistence(); P->GAP = GAP;

   //srand((unsigned)time(NULL));
   srand(550);

   P->loadConfig();
   readJSONdata(GAP->conf->datafile);

   CH   = new ConstructHeu(GAP,GAP->zub); 
   MT   = new MTHG(GAP,GAP->zub);
   SA   = new SimAnnealing(GAP, GAP->zub);
   TS   = new TabuSearch(GAP,GAP->zub);
   MLS  = new MetaLocalSearch(GAP,LS,GAP->zub); 
   RINS = new Rins(GAP,GAP->zub); 
}

Controller::~Controller()
{
   //dtor
   delete GAP->sol;
   delete GAP->solbest;
   delete GAP;
   delete P;
   if(CH != NULL) delete CH;
   if(SA != NULL) delete SA;
   if(TS != NULL) delete TS;
   if(MLS != NULL) delete MLS;
   if(RINS != NULL) delete RINS;
   cout << "GAP deleted";
}

// legge i dati da file json
void Controller::readJSONdata(string filename)
{
   P->readJSONdata(filename);
   setGAPdata();
   //srand (time(NULL));
   srand (550);
}

// dunno, could do without this one, isn't it?
void Controller::setGAPdata()
{
   GAP->sol     = new int[GAP->n];
   GAP->solbest = new int[GAP->n];
   GAP->zub     = INT_MAX;
}

double Controller::simpleContruct()
{  return CH->simpleContruct();
}

// Martello and Toth heuristic
int Controller::run_mthg()
{  return MT->run_mthg();
}

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
   bound = LB->benders();
   cout << "Benders Decomposition bound " << bound << endl;

   if(LB != NULL) delete LB;
   LB = NULL;
}

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

int Controller::opt10()
{  int z,zcheck;
   LS = new LocalSearch(GAP, GAP->zub);

   // no solution to improve
   if(GAP->zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }

   z = LS->opt10(GAP->c);
   zcheck = GAP->checkSol(GAP->sol);
   if (abs(zcheck - z) > GAP->EPS)
   {  cout << "[1.0opt] Ahi ahi" << endl;
      z = INT_MAX;
   }

   z = LS->opt11(GAP->c);
   zcheck = GAP->checkSol(GAP->sol);
   if (abs(zcheck - z) > GAP->EPS)
   {
      cout << "[1.0opt] Ahi ahi" << endl;
      z = INT_MAX;
   }

   if (LS != NULL) delete LS;

   cout << "Opt01 zub= " << z << endl;
   return GAP->zub;
}

int Controller::simAnn()
{
   return SA->simAnneling( GAP->c, 
      GAP->conf->SA->maxT,
      GAP->conf->SA->k,
      GAP->conf->SA->alpha,
      GAP->conf->SA->maxiter,
      GAP->conf->SA->iterAnneal
   );
}

int Controller::tabuSearch()
{
   return TS->tabuSearch(GAP->c,
      GAP->conf->TS->Ttenure,
      GAP->conf->TS->maxiter
   );
}

int Controller::run_genAlgo()
{
   GA = new GenAlgo(GAP, GAP->zub);
   int res = GA->genAlgo(  GAP->c,
                           GAP->conf->GA->maxiter,
                           GAP->conf->GA->numpop,
                           GAP->conf->GA->pc
                        );
   if(GA != NULL) delete GA;
   GA = NULL;
   return res;
}

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

int Controller::iteratedLS()
{  
   if(GAP->zub==INT_MAX) 
   {  cout << "Uninitialized solution" << endl;
      return INT_MAX;
   }
   else
      return MLS->iteratedLocSearch(GAP->c,
                                    GAP->conf->IterLS->maxiter,
                                    GAP->conf->IterLS->alpha);
}

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

int Controller::run_rins()
{  int res; 
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

int Controller::run_kernel()
{
   KER = new Kernel(GAP, GAP->zub);
   int res = KER->solveByKernel( (GAP->conf->isVerbose ? true : false));
   if (KER != NULL) delete KER;
   KER = NULL;
   return res;
}
