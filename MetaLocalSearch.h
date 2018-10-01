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

   private:
      LocalSearch* LS;
      void dataPerturbation(int**,int**,double);

      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // METALOCALSEARCH_H
