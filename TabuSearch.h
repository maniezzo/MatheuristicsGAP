#ifndef TABUSEARCH_H
#define TABUSEARCH_H
#include "GAP.h"
#include "Rins.h"

class TabuSearch
{
   public:
      GeneralizedAssignemnt* GAP;

      TabuSearch(GeneralizedAssignemnt*, int&);
      ~TabuSearch();
      int tabuSearch(int**,int,int,int);

   private:
      // local mirrors
      int   m,n;
      int   *sol,*solbest;
      int** req;
      int & zub,zlb;

      int dive(Rins* RINS, int* newsol);
};

#endif // TABUSEARCH_H
