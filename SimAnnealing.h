#ifndef SIMANNEALING_H
#define SIMANNEALING_H
#include "GAP.h"

class SimAnnealing
{
   public:
      GeneralizedAssignemnt* GAP;

      SimAnnealing(GeneralizedAssignemnt*, int&);
      ~SimAnnealing();
      int simAnnealing(int**, double, double, double, int, int, double, int);
      int mathSimAnnealing(int**, double, double, double, int, int, double);

   private:
      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // SIMANNEALING_H
