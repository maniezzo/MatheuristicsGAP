#ifndef METALOCALSEARCH_H
#define METALOCALSEARCH_H
#include "GAP.h"
#include "LocalSearch.h"
#include "VLSN.h"

class MetaLocalSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      MetaLocalSearch(GeneralizedAssignemnt*, LocalSearch*, int&);
      ~MetaLocalSearch();
      int iteratedLocSearch(int**,int,double);
      double MetaLocalSearch::GRASP(int maxIter, int candNum);

   private:
      LocalSearch* LS;
      void dataPerturbation(int**,int**,double);
      double MetaLocalSearch::GRASPcontruct(int candNum);

      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // METALOCALSEARCH_H
