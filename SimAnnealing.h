#ifndef SIMANNEALING_H
#define SIMANNEALING_H
#include "GAP.h"

class SimAnnealing
{
   public:
      GeneralizedAssignemnt* GAP;

      SimAnnealing(GeneralizedAssignemnt*, int&);
      ~SimAnnealing();
      int simAnneling(int**, double, double, double, int, int);

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // SIMANNEALING_H
