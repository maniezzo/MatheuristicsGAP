#ifndef CONSTRUCTHEU_H
#define CONSTRUCTHEU_H
#include "GAP.h"

class ConstructHeu
{
   public:
      GeneralizedAssignemnt* GAP;

      ConstructHeu(GeneralizedAssignemnt*, int&);
      ~ConstructHeu();
      int simpleConstruct();
      int construct(vector<int> indCost);

   private:
      //void computeRegrets(vector<int> &);

      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // CONSTRUCTHEU_H
