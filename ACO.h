#ifndef ACO_H
#define ACO_H
#include "GAP.h"
#include "Rins.h"

class ACO
{
   public:
      GeneralizedAssignemnt* GAP;

      ACO(GeneralizedAssignemnt*, int&);
      ~ACO();
      int antColony(int**, int);

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;

      // private functions
      int ACO::montecarlo(vector<double>&);
};

#endif // ACO_H
