#include "VLSN.h"

VeryLarge::VeryLarge(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

VeryLarge::~VeryLarge()
{
   //dtor
}

int VeryLarge::verylarge(int** c, int maxiter)
{  int z=0;

   return z;
}
