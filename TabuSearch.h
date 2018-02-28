#ifndef TABUSEARCH_H
#define TABUSEARCH_H
#include "GAP.h"

class TabuSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      TabuSearch(GeneralizedAssignemnt*, int&);
      ~TabuSearch();
      int tabuSearch(int**,int,int);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // TABUSEARCH_H
