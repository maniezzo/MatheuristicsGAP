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
      int constructWithInd(vector<int> indCost, bool isVerbose);

   private:
      //void computeRegrets(vector<int> &);

      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // CONSTRUCTHEU_H
