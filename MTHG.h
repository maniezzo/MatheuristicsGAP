#ifndef MTHG_H
#define MTHG_H
#include "GAP.h"

class MTHG
{
   public:
      GeneralizedAssignemnt* GAP;

      MTHG(GeneralizedAssignemnt*, int&);
      ~MTHG();
      int run_mthg();

   private:
      int mthg_(int *n, int *m, int *p, int *w, int *c__, int minmax, int *z__, int *xstar, int jck);
      int gha_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
      int feas_(int *, int *, int *, int *, int *, int *, int *, int *, int *);
      int trin_(int *, int *, int *, int *, int *, int *, int *);
      int ghbcd_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, double *, int *);
      int chmthg_(int *, int *, int *, int *, int *, int *, int *, int *);

      // local mirrors
      int m,n;
      int *sol,*solbest;
      int** req;
      int & zub,zlb;
};

#endif // MTHG_H
