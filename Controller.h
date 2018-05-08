#ifndef CONTROLLER_H
#define CONTROLLER_H
#include "GAP.h"
#include "Persistence.h"
#include "ConstructHeu.h"
#include "LocalSearch.h"
#include "MTHG.h"
#include "LowerBound.h"
#include "SimAnnealing.h"
#include "TabuSearch.h"
#include "GeneticAlgorithm.h"
#include "Ejection.h"
#include "MetaLocalSearch.h"
#include "Lagrangian.h"
#include "Rins.h"
#include "BeamSearch.h"
#include "ForwardBackward.h"
#include "Corridor.h"
#include "LocBranching.h"
#include "Benders.h"
#include "Kernel.h"
#include "VLSN.h"

class Controller
{
   public:
      Persistence*      P;             // read config and instance
      GeneralizedAssignemnt* GAP;      // global data structures
      ConstructHeu*     CH;            // simple constructive
      MTHG*             MT;            // Martello and Toth heuristic
      LowerBound*       LB;            // aka cplex
      LocalSearch*      LS;            // local search (10, 11, 21)
      SimAnnealing*     SA;            // Simulated Annealing (bare)
      TabuSearch*       TS;            // Tabu search
      GenAlgo*          GA;            // Genetic algorithm
      Ejection*         EC;            // Ejection chain
      MetaLocalSearch*  MLS;           // ILS, VND, VNS
      Lagrangian*       LAGR;          // lagrangian, assignment and capacities
      Rins*             RINS;          // Rins, would you guess?
      BeamSearch*       BS;            // Beam search
      FandB*            FB;            // Forward and backward
      Corridor*         CORR;          // Corridor
      LocBranching*     LBR;           // Local Branching
      Benders*          BEND;          // Benders
      Kernel*           KER;           // Kernel
      VeryLarge*        VLSN;          // Very large-scale neighborhood search

      Controller();
      ~Controller();
      void   readJSONdata(string filename);   // reads data from json file
      double simpleContruct();
      void   computeBounds();
      double exactCplex();
      int    run_mthg();
      int    opt10();
      int    simAnn();
      int    tabuSearch();
      int    iteratedLS();
      int    run_genAlgo();
      int    run_ejection();
      int    run_lagrAss();
      int    run_lagrCap();
      int    run_rins();
      int    run_beam();
      int    run_FandB();
      int    run_Corridor();
      int    run_localBranching();
      int    run_benders();
      int    run_kernel();
      int    run_VLSN();

   private:
};

#endif // CONTROLLER_H
