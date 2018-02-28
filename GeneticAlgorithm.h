#ifndef GENALG_H
#define GENALG_H
#include "GAP.h"
#include "Rins.h"

class GenAlgo
{
   public:
      GeneralizedAssignemnt* GAP;

      GenAlgo(GeneralizedAssignemnt*, int&);
      ~GenAlgo();
      int genAlgo(int**, int, int, double);

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;

      // private functions
      bool initPop(vector<int*>&,int,vector<double>&);
      int* generateOneSol(int*);
      void GenAlgo::selCrossover(vector<int*>& pop, int numpop, vector<double>&,Rins* RINS, int*);
      bool mutate(vector<int*>& pop, int numpop);
      int GenAlgo::montecarlo(vector<double>&);
      void GenAlgo::crossOver(vector<int*>& pop, int p1, int p2,Rins* RINS, int*);
};

#endif // GENALG_H
