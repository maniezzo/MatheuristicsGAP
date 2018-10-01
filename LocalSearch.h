#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H
#include "GAP.h"

class LocalSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      LocalSearch(GeneralizedAssignemnt*, int&);
      ~LocalSearch();
      int opt10(int**,bool);
      double opt11(int**,bool);
      void neigh21();

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // LOCALSEARCH_H
