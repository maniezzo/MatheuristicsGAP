#ifndef METALOCALSEARCH_H
#define METALOCALSEARCH_H
#include "GAP.h"
#include "LocalSearch.h"
#include "VLSN.h"
#include "LowerBound.h"

class MetaLocalSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      MetaLocalSearch(GeneralizedAssignemnt*, LocalSearch*, int&);
      ~MetaLocalSearch();
      int iteratedLocSearch(int**,int,double);
      double MetaLocalSearch::GRASP(int maxIter, int candNum);
      double MetaLocalSearch::VNS(int maxIter, bool isMatheuristic);
      double MetaLocalSearch::VNSbasic(int maxIter);

   private:
      LocalSearch* LS;
      LowerBound*  LB;
      void   dataPerturbation(int**,int**,double);
      double MetaLocalSearch::GRASPcontruct(int candNum, bool isMatheuristic);

      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // METALOCALSEARCH_H
