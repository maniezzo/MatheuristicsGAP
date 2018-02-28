#include "MTHG.h"

MTHG::MTHG(GeneralizedAssignemnt* GAPinstance, int & zz) : zub(zz)
{
   //ctor
   GAP = GAPinstance;
   m = GAP->m;
   n = GAP->n;
   sol = GAP->sol;
   solbest = GAP->solbest;
   req = GAP->req;
}

MTHG::~MTHG()
{
   //dtor
}

int MTHG::run_mthg()
{  int res,i,j;
   int* xstar = (int *) malloc (n * sizeof(int)); // best solution
   int* p = (int *) malloc (n*m * sizeof(int));   // allocaton costs
   int* w = (int *) malloc (n*m * sizeof(int));   // customers requests

   for(j=0;j<n;j++)
      for(i=0;i<m;i++)
      {  p[j*m+i] = GAP->c[i][j];
         w[j*m+i] = GAP->req[i][j];
      }

   res = mthg_(&n, &m, p, w, GAP->cap, 1, &zub, xstar, 1);

   for(i=0;i<n;i++) GAP->sol[i] = xstar[i]-1;
   res = GAP->checkSol(GAP->sol);

   cout << "MTHG: zub = " << zub << " res= " << res << endl;

   if(xstar != NULL) free(xstar);
   if(p != NULL) free(p);
   if(w != NULL) free(w);
   return res;
}

// What follows is mainly a f2c conversion of MTHG.FOR
int MTHG::mthg_(int *n, int *m, int *p, int *w, int *c__, int minmax, int *z__, int *xstar, int jck)
{
   /* System generated locals */
   int i__1, i__2, res;

   /* Local variables */
   static int *a	/* was [50][500] */, i__, j, zm;
   static int jfi, lam, inf, iub;
   static int* best;

   static int kvst, *dmyc1, *dmyc2, *dmyc3, *dmyc4, 
              *dmyr1, *dmyr2, *dmyr3, *dmyr4, *dmyr5;
   static int jdimc, jdimr, imult, invst;
   static double *dmycr1;

   // allocations
   best = (int *) malloc (*n * sizeof(int)); // int BEST(500)
   dmycr1 = (double *) malloc (*n * sizeof(double)); // double    DMYCR1(500)

   //int DMYR1(50),DMYR2(50),DMYR3(50),DMYR4(50),DMYR5(50)
   dmyr1 = (int *) malloc (*m * sizeof(int)); 
   dmyr2 = (int *) malloc (*m * sizeof(int)); 
   dmyr3 = (int *) malloc (*m * sizeof(int)); 
   dmyr4 = (int *) malloc (*m * sizeof(int)); 
   dmyr5 = (int *) malloc (*m * sizeof(int)); 

   // int DMYC1(500),DMYC2(500),DMYC3(500),DMYC4(500)
   dmyc1 = (int *) malloc (*n * sizeof(int)); 
   dmyc2 = (int *) malloc (*n * sizeof(int)); 
   dmyc3 = (int *) malloc (*n * sizeof(int)); 
   dmyc4 = (int *) malloc (*n * sizeof(int)); 

   // int A(50,500)
   a = (int *) malloc ((*m)*(*n) * sizeof(int));

   /* THIS SUBROUTINE HEURISTICALLY SOLVES THE GENERALIZED ASSIGNMENT */
   /* PROBLEM */

   /* OPT Z = P(1,1)*X(1,1) + ... + P(1,N)*X(1,N) + */
   /*                         ...                 + */
   /*         P(M,1)*X(M,1) + ... + P(M,N)*X(M,N) */

   /*     (WHERE  OPT = MIN  IF  MINMAX = 1 ,  OPT = MAX  IF  MINMAX = 2 ) */

   /* SUBJECT TO: */

   /*       W(I,1)*X(I,1) + ... + W(I,N)*X(I,N) .LE. C(I)  FOR I=1,...,M, */
   /*       X(1,J) + ... + X(M,J) = 1                      FOR J=1,...,N, */
   /*       X(I,J) = 0 OR 1                     FOR I=1,...,M, J=1,...,N. */

   /* THE PROGRAM IS INCLUDED IN THE VOLUME */
   /*   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS */
   /*   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990 */
   /* AND IMPLEMENTS THE POLYNOMIAL-TIME ALGORITHMS DESCRIBED */
   /* IN SECTION  7.4 . */

   /* THE INPUT PROBLEM MUST SATISFY THE CONDITIONS */

   /*   1) 2 .LE. M .LE. JDIMR ; */
   /*   2) 2 .LE. N .LE. JDIMC ; */
   /*      ( JDIMR  AND  JDIMC  ARE DEFINED BY THE FIRST TWO EXECUTABLE */
   /*       STATEMENTS;) */
   /*   3) P(I,J), W(I,J) AND C(I) POSITIVE intS; */
   /*   4) W(I,J) .LE. C(I) FOR AT LEAST ONE I, FOR J=1,...,N; */
   /*   5) C(I) .GE. MIN (W(I,J)) FOR I=1,...,M. */

   /* MTHG CALLS 6 PROCEDURES: CHMTHG, FEAS, GHA, GHBCD, GHX AND TRIN. */

   /* THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS */
   /* ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MTHG. */
   /* NO MACHINE-DEPENDENT CONSTANT IS USED. */
   /* THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN */
   /* AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE */
   /* SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY). */
   /* THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P. */
   /* 9000/840. */

   /* MTHG NEEDS */
   /*   6 ARRAYS ( C ,  DMYR1 ,  DMYR2 ,  DMYR3 ,  DMYR4  AND  DMYR5 ) OF */
   /*              LENGTH AT LEAST  JDIMR ; */
   /*   7 ARRAYS ( XSTAR ,  BEST ,  DMYC1 ,  DMYC2 ,  DMYC3 ,  DMYC4  AND */
   /*              DMYCR1 ) OF LENGTH AT LEAST  JDIMC ; */
   /*   3 ARRAYS ( P ,  W  AND  A ) OF LENGTH AT LEAST  JDIMR X JDIMC . */

   /* THE ARRAYS ARE CURRENTLY DIMENSIONED TO ALLOW PROBLEMS FOR WHICH */
   /*       M .LE. 50 , */
   /*       N .LE. 500 */
   /* (SO, IN THE CALLING PROGRAM, ARRAYS  P  AND  W  MUST BE DIMENSIONED */
   /* AT  (50,500) ). CHANGING SUCH LIMITS NECESSITATES CHANGING THE */
   /* DIMENSION OF ALL THE ARRAYS IN SUBROUTINE MTHG, AS WELL AS THE FIRST */
   /* TWO EXECUTABLE STATEMENTS. */

   /* MEANING OF THE INPUT PARAMETERS: */
   /* N        = NUMBER OF ITEMS; */
   /* M        = NUMBER OF KNAPSACKS; */
   /* P(I,J)   = PROFIT OF ITEM J IF ASSIGNED TO KNAPSACK I */
   /*            (I=1,...,M; J=1,...,N); */
   /* W(I,J)   = WEIGHT OF ITEM J IF ASSIGNED TO KNAPSACK I */
   /*            (I=1,...,M; J=1,...,N); */
   /* C(I)     = CAPACITY OF KNAPSACK I (I=1,...,M); */
   /* MINMAX   = 1 IF THE OBJECTIVE FUNCTION MUST BE MINIMIZED, */
   /*          = 2 IF THE OBJECTIVE FUNCTION MUST BE MAXIMIZED; */
   /* JCK      = 1 IF CHECK ON THE INPUT DATA IS DESIRED, */
   /*          = 0 OTHERWISE. */

   /* MEANING OF THE OUTPUT PARAMETERS: */
   /* Z        = VALUE OF THE SOLUTION FOUND IF Z .GT. 0 , */
   /*          = 0 IF NO FEASIBLE SOLUTION IS FOUND, */
   /*          = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI- */
   /*            TION  - Z  IS VIOLATED; */
   /* XSTAR(J) = KNAPSACK WHERE ITEM J IS INSERTED IN THE SOLUTION FOUND. */

   /* ALL THE PARAMETERS ARE int. ON RETURN OF MTHG ALL THE INPUT */
   /* PARAMETERS ARE UNCHANGED, BUT  P(I,J)  IS SET TO  0  FOR ALL PAIRS */
   /* (I,J)  SUCH THAT  W(I,J) .GT. C(I) . */


   /* DEFINITION OF THE INTERNAL PARAMETERS. */

   /* Parameter adjustments */
   --xstar;
   --c__;
   w -= (*m+1);
   p -= (*m+1);

   /* Function Body */
   jdimr = *m;
   jdimc = *n;
   *z__ = 0;
   if (jck == 1) {
      chmthg_(n, m, &p[(*m+1)], &w[(*m+1)], &c__[1], &jdimr, &jdimc, z__);
   }
   if (*z__ < 0) 
   {  cout << "Error in input data, exiting" << endl;
      goto Lend;
   }

   /* INITIALIZE. */

   invst = 0;
   imult = -1;
   if (minmax == 2)  // maximization
      goto L10;

   /* TRANSFORM THE MINIMIZATION PROBLEM INTO A MAXIMIZATION PROBLEM. */
   trin_(&p[(*m+1)], n, m, &invst, &lam, &jdimr, &jdimc);
   imult = 1;

   /* SOLVE THE MAXIMIZATION PROBLEM. */

   /* CHECK FOR CAPACITY INFEASIBILITY. */
L10:
   res = feas_(n, m, &p[(*m+1)], &w[(*m+1)], &c__[1], &xstar[1], &jfi, &jdimr, &jdimc);
   if (jfi == 1) 
      goto L30;

   /* FIRST HEURISTIC SOLUTION. */
   gha_(&p[(*m+1)], &w[(*m+1)], &c__[1], n, m, z__, &xstar[1], &iub, best, &kvst, &inf, 
        &jdimr, &jdimc, dmyr1, dmyr2, dmyc1, dmyc2, dmyc3, dmyc4);
   if (*z__ == iub) 
      goto L20;

   /* SECOND HEURISTIC SOLUTION. */
   ghbcd_(&p[(*m+1)], &w[(*m+1)], &c__[1], n, m, z__, &xstar[1], &inf, &jdimr, &jdimc, 
          dmyc1, dmyr1, dmyr2, dmyr3, dmyr4, dmyr5, dmyc2, dmyc3, dmyc4, dmycr1, a);

   /* TERMINATE. */

L20:
   zm = *z__;
   *z__ = 0;
   if (zm > kvst) {
      *z__ = invst - zm * imult;
   }
L30:
   if (minmax == 2) 
      goto Lend;

   /* RE-STORE THE ORIGINAL MINIMIZATION PROBLEM. */
   i__1 = *m;
   for (i__ = 1; i__ <= i__1; ++i__) {
      i__2 = *n;
      for (j = 1; j <= i__2; ++j) {
         if (p[i__ + j * *m] > 0) {
            p[i__ + j * *m] = lam - p[i__ + j * *m];
         }
         /* L40: */
      }
      /* L50: */
   }

Lend: 
   // deallocations
   if(best != NULL) free(best);
   if(dmycr1 != NULL) free(dmycr1);
   if(dmyr1 != NULL) free(dmyr1);
   if(dmyr2 != NULL) free(dmyr2);
   if(dmyr3 != NULL) free(dmyr3);
   if(dmyr4 != NULL) free(dmyr4);
   if(dmyr5 != NULL) free(dmyr5);
   if(dmyc1 != NULL) free(dmyc1);
   if(dmyc2 != NULL) free(dmyc2);
   if(dmyc3 != NULL) free(dmyc3);
   if(dmyc4 != NULL) free(dmyc4);

   return *z__;
} /* mthg_ */

int  MTHG::chmthg_(int *n, int *m, int *p, int *w, int *c__, int *jdimr, int *jdimc, int *z__)
{
   /* System generated locals */
   int p_dim1, p_offset, w_dim1, w_offset, i__1, i__2;

   /* Local variables */
   static int i__, j, min__;

   /* CHECK THE INPUT DATA. */

   /* Parameter adjustments */
   --c__;
   w_dim1 = *jdimr;
   w_offset = 1 + w_dim1;
   w -= w_offset;
   p_dim1 = *jdimr;
   p_offset = 1 + p_dim1;
   p -= p_offset;

   /* Function Body */
   if (*m <= 1)     *z__ = -1;
   if (*m > *jdimr) *z__ = -1;
   if (*z__ < 0)    goto Lend;
   if (*n <= 1)     *z__ = -2;
   if (*n > *jdimc) *z__ = -2;
   if (*z__ < 0)    goto Lend;

   i__1 = *m;
   for (i__ = 1; i__ <= i__1; ++i__) 
   {
      if (c__[i__] > 0) 
         goto L10;

      *z__ = -3;
      goto Lend;
   L10:
      min__ = c__[i__] + 1;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j) 
      {
         if (p[i__ + j * p_dim1] > 0 && w[i__ + j * w_dim1] > 0) 
            goto L20;

         *z__ = -3;
         goto Lend;
      L20:
         if (w[i__ + j * w_dim1] < min__) 
            min__ = w[i__ + j * w_dim1];
         /* L30: */
      }
      if (c__[i__] < min__) 
         *z__ = -5;
      /* L40: */
   }
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) 
   {
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) 
         if (w[i__ + j * w_dim1] <= c__[i__]) 
            goto L60;
         /* L50: */

      *z__ = -4;
      goto Lend;
   L60:
      ;
   }
Lend:
   return 0;
} /* chmthg_ */

int  MTHG::feas_(int *n, int *m, int *p, int *w, int *c__, int *xstar, int *jfi, int *jdimr, int * jdimc)
{
   /* System generated locals */
   int p_dim1, p_offset, w_dim1, w_offset, i__1, i__2;

   /* Local variables */
   static int i__, j, kinf;

   /* CHECK FOR INFEASIBILITY. */

   /* Parameter adjustments */
   --c__;
   --xstar;
   w_dim1 = *jdimr;
   w_offset = 1 + w_dim1;
   w -= w_offset;
   p_dim1 = *jdimr;
   p_offset = 1 + p_dim1;
   p -= p_offset;

   /* Function Body */
   *jfi = 0;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) 
   {
      xstar[j] = 0;
      kinf = 0;
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) 
      {
         if (w[i__ + j * w_dim1] <= c__[i__]) 
            goto L10;

         ++kinf;
         p[i__ + j * p_dim1] = 0;
      L10:
         ;
      }
      if (kinf == *m) 
         *jfi = 1;
      /* L20: */
   }
   return kinf;
} /* feas_ */

int  MTHG::gha_(int *p, int *w, int *c__, int *n, int *m, int *z__, int *xstar, int *iub, int *best,
   int *kvst, int *inf, int *jdimr, int *jdimc, int *kw, int *mw, int *pen, int *first, int *second, int *bb)
{
   /* System generated locals */
   int p_dim1, p_offset, w_dim1, w_offset, i__1, i__2;

   /* Local variables */
   static int i__, j, if__, nb, jj, io, jo, is, jjm, fmax, smax, index, ipmin, newsec, maxpen;

   /* APPLY THE APPROXIMATE ALGORITHM GH WITH FUNCTION (A) AND */
   /* DEFINE THE INFINITE VALUE  INF . */

   /* IF IUB = Z THE SOLUTION IS OPTIMAL; */
   /* IF Z = KVST NO FEASIBLE SOLUTION WAS FOUND. */

   /* Parameter adjustments */
   --mw;
   --kw;
   --c__;
   --bb;
   --second;
   --first;
   --pen;
   --best;
   --xstar;
   w_dim1   = *jdimr;
   w_offset = 1 + w_dim1;
   w -= w_offset;
   p_dim1   = *jdimr;
   p_offset = 1 + p_dim1;
   p -= p_offset;

   /* Function Body */
   *inf = 0;
   i__1 = *m;
   for (i__ = 1; i__ <= i__1; ++i__) 
   {
      kw[i__] = c__[i__];
      mw[i__] = 0;
      if (c__[i__] > *inf) 
         *inf = c__[i__];

      /* L10: */
   }
   *iub  = 0;
   *z__  = 0;
   *kvst = 0;
   i__1  = *n;
   for (j = 1; j <= i__1; ++j) 
   {
      ipmin = p[j * p_dim1 + 1];    // i of best choice
      fmax = p[j * p_dim1 + 1];     // best choice profit
      if__ = 1;
      smax = 0;                     // second best
      i__2 = *m;
      for (i__ = 2; i__ <= i__2; ++i__) 
      {
         if (p[i__ + j * p_dim1] < ipmin) 
            ipmin = p[i__ + j * p_dim1];

         if (smax >= p[i__ + j * p_dim1]) 
            goto L30;

         if (fmax >= p[i__ + j * p_dim1]) 
            goto L20;

         smax = fmax;
         is = if__;
         fmax = p[i__ + j * p_dim1];
         if__ = i__;
         goto L30;
L20:
         smax = p[i__ + j * p_dim1];
         is = i__;
L30:
         ;
      }

      *kvst += ipmin;
      first[j]  = if__;             // best choice for j
      best[j]   = if__;
      second[j] = is;               // second best for j
      pen[j]    = fmax - smax;      // regret for j (penalty)
      if (smax == 0) 
         pen[j] = -1;

      bb[j] = j;
      *iub += fmax;
      //cout << "iub: " << iub << " fmax " << fmax << endl;
      if (w[if__ + j * w_dim1] > mw[if__]) 
         mw[if__] = w[if__ + j * w_dim1];

      if (w[is + j * w_dim1] > mw[is]) 
         mw[is] = w[is + j * w_dim1];

      /* L40: */
   }

   if (*kvst > 0) 
      --(*kvst);

   if (*iub > *inf) 
      *inf = *iub;

   i__1 = *n;
   for (j = 1; j <= i__1; ++j) 
      if (pen[j] == -1) 
         pen[j] = *inf;
      /* L50: */

   nb = *n;
L60:
   maxpen = -1;
   i__1 = nb;
   for (jj = 1; jj <= i__1; ++jj) 
   {
      j = bb[jj];
      if (pen[j] <= maxpen) 
         goto L70;

      maxpen = pen[j];
      jjm = jj;
L70:
      ;
   }
   jo = bb[jjm];
   io = first[jo];
   *z__ += p[io + jo * p_dim1];
   xstar[jo] = io;
   bb[jjm] = bb[nb];
   --nb;

   if (nb == 0) 
      goto L130;

   kw[io] -= w[io + jo * w_dim1];
   if (mw[io] <= kw[io]) 
      goto L60;

   i__1 = nb;
   for (jj = 1; jj <= i__1; ++jj) 
   {
      j = bb[jj];
      if (w[io + j * w_dim1] <= kw[io]) 
         goto L120;

      if (first[j] != io) 
         goto L80;

      if (pen[j] == *inf) 
         goto L130;

      first[j] = second[j];
      goto L90;
L80:
      if (second[j] != io) 
         goto L120;

L90:
      index = first[j];
      w[index + j * w_dim1] += *inf;
      newsec = 0;
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) 
      {
         if (w[i__ + j * w_dim1] > kw[i__]) 
            goto L100;

         if (p[i__ + j * p_dim1] <= newsec) 
            goto L100;

         newsec = p[i__ + j * p_dim1];
         is = i__;
L100:
         ;
      }
      w[index + j * w_dim1] -= *inf;
      if (newsec == 0) 
         goto L110;

      second[j] = is;
      pen[j] = p[index + j * p_dim1] - newsec;
      if (w[is + j * w_dim1] > mw[is]) 
         mw[is] = w[is + j * w_dim1];

      goto L120;
L110:
      pen[j] = *inf;
L120:
      ;
   }
   goto L60;
L130:
   *z__ = *kvst;
   return 0;
} /* gha_ */

int  MTHG::ghbcd_(int *p, int *w, int *c__, int *n, 
   int *m, int *z__, int *xstar, int *inf, int *
   jdimr, int *jdimc, int *xsp, int *dmyr1, int *dmyr2, 
   int *dmyr3, int *dmyr4, int *dmyr5, int *dmyc2, 
   int *dmyc3, int *dmyc4, double *dmycr1, int *dmya)
{
   /* System generated locals */
   int p_dim1, p_offset, w_dim1, w_offset, dmya_dim1, dmya_offset, i__1;

   /* Local variables */
   static int j;
   static double a1, a2, a3, a4, a5;
   static int jj;
   extern /* Subroutine */ int ghx_(int *, int *, int *, int 
      *, int *, int *, int *, double *, double *, double *, double *
      , double *, int *, int *, int *, int *, int *, 
      int *, int *, int *, int *, int *, int *, 
      double *, int *);
   static int vsp;


   /* APPLY THE APPROXIMATE ALGORITHM GH WITH FUNCTIONS (B), (C) AND (D). */

   /* Parameter adjustments */
   --dmyr5;
   --dmyr4;
   --dmyr3;
   --dmyr2;
   --dmyr1;
   --c__;
   dmya_dim1 = *jdimr;
   dmya_offset = 1 + dmya_dim1;
   dmya -= dmya_offset;
   --dmycr1;
   --dmyc4;
   --dmyc3;
   --dmyc2;
   --xsp;
   --xstar;
   w_dim1 = *jdimr;
   w_offset = 1 + w_dim1;
   w -= w_offset;
   p_dim1 = *jdimr;
   p_offset = 1 + p_dim1;
   p -= p_offset;

   /* Function Body */
   jj = 2;
   a1 = 1.f;
   a2 = 0.f;
   a3 = 0.f;
   a4 = 0.f;
   a5 = 1.f;
L10:
   ghx_(&p[p_offset], &w[w_offset], &c__[1], n, m, &vsp, &xsp[1], &a1, &a2, &
      a3, &a4, &a5, inf, jdimr, jdimc, &dmyr1[1], &dmyr2[1], &dmyr3[1], 
      &dmyr4[1], &dmyr5[1], &dmyc2[1], &dmyc3[1], &dmyc4[1], &dmycr1[1],
      &dmya[dmya_offset]);
   if (vsp <= *z__) {
      goto L30;
   }
   *z__ = vsp;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      xstar[j] = xsp[j];
      /* L20: */
   }
L30:
   if (jj == 3) {
      goto L40;
   }
   if (jj == 4) {
      goto L50;
   }
   jj = 3;
   a1 = 1.f;
   a2 = 0.f;
   a3 = 1.f;
   a4 = 0.f;
   a5 = 0.f;
   goto L10;
L40:
   jj = 4;
   a1 = 0.f;
   a2 = 1.f;
   a3 = 0.f;
   a4 = 1.f;
   a5 = 0.f;
   goto L10;
L50:
   return 0;
} /* ghbcd_ */

int ghx_(int *p, int *w, int *c__, int *n, 
   int *m, int *z__, int *xstar, double *a1, double *a2, double *
   a3, double *a4, double *a5, int *inf, int *jdimr, int *jdimc, 
   int *kw, int *mw, int *minw, int *kchan, int *kwr,
   int *first, int *second, int *bb, double *pen, int *wl)
{
   /* System generated locals */
   int p_dim1, p_offset, w_dim1, w_offset, wl_dim1, wl_offset, i__1, i__2;

   /* Local variables */
   static int i__, j;
   static double s;
   static int if__, nb, jf, jj, kk, io, jo, ip, is, js;
   static double rp;
   static int jjm;
   static double rwl, rkw, fmax;
   static int maxp;
   static double smax, maxpen;

   /* APPLY THE APPROXIMATE ALGORITHM GH WITH FUNCTION (B) OR (C) OR (D). */

   /* Parameter adjustments */
   --kwr;
   --kchan;
   --minw;
   --mw;
   --kw;
   --c__;
   wl_dim1 = *jdimr;
   wl_offset = 1 + wl_dim1;
   wl -= wl_offset;
   --pen;
   --bb;
   --second;
   --first;
   --xstar;
   w_dim1 = *jdimr;
   w_offset = 1 + w_dim1;
   w -= w_offset;
   p_dim1 = *jdimr;
   p_offset = 1 + p_dim1;
   p -= p_offset;

   /* Function Body */
   i__1 = *m;
   for (i__ = 1; i__ <= i__1; ++i__) {
      kw[i__] = c__[i__];
      mw[i__] = 0;
      minw[i__] = *inf;
      i__2 = *n;
      for (j = 1; j <= i__2; ++j) {
         wl[i__ + j * wl_dim1] = w[i__ + j * w_dim1];
         if (wl[i__ + j * wl_dim1] < minw[i__]) {
            minw[i__] = wl[i__ + j * wl_dim1];
         }
         /* L10: */
      }
      kwr[i__] = kw[i__] - minw[i__];
      /* L20: */
   }
   *z__ = 0;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      fmax = (double) (-(*inf));
      if__ = 0;
      smax = (double) (-(*inf));
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
         if (wl[i__ + j * wl_dim1] > kw[i__]) {
            goto L40;
         }
         if (wl[i__ + j * wl_dim1] > kwr[i__]) {
            wl[i__ + j * wl_dim1] = kw[i__];
         }
         rwl = (double) wl[i__ + j * wl_dim1];
         rp = (double) p[i__ + j * p_dim1];
         rkw = (double) kw[i__];
         s = (-(*a1) * rwl + *a2 * rp) / (*a3 * rkw + *a4 * rwl + *a5);
         if (smax >= s) {
            goto L40;
         }
         if (fmax >= s) {
            goto L30;
         }
         smax = fmax;
         is = if__;
         fmax = s;
         if__ = i__;
         goto L40;
      L30:
         smax = s;
         is = i__;
      L40:
         ;
      }
      first[j] = if__;
      second[j] = is;
      pen[j] = fmax - smax;
      bb[j] = j;
      if (wl[if__ + j * wl_dim1] > mw[if__]) {
         mw[if__] = wl[if__ + j * wl_dim1];
      }
      if (smax > (double) (-(*inf))) {
         goto L50;
      }
      pen[j] = (double) (*inf);
      goto L60;
   L50:
      if (wl[is + j * wl_dim1] > mw[is]) {
         mw[is] = wl[is + j * wl_dim1];
      }
   L60:
      ;
   }
   nb = *n;
L70:
   maxpen = -1.f;
   i__1 = nb;
   for (jj = 1; jj <= i__1; ++jj) {
      j = bb[jj];
      if (pen[j] <= maxpen) {
         goto L80;
      }
      maxpen = pen[j];
      jjm = jj;
   L80:
      ;
   }
   jo = bb[jjm];
   io = first[jo];
   *z__ += p[io + jo * p_dim1];
   xstar[jo] = io;
   bb[jjm] = bb[nb];
   --nb;
   kw[io] -= w[io + jo * w_dim1];
   if (nb == 0) {
      goto L210;
   }
   kk = 0;
   i__1 = *m;
   for (i__ = 1; i__ <= i__1; ++i__) {
      kchan[i__] = 0;
      if (wl[i__ + jo * wl_dim1] > minw[i__]) {
         goto L100;
      }
      minw[i__] = *inf;
      i__2 = nb;
      for (jj = 1; jj <= i__2; ++jj) {
         j = bb[jj];
         if (wl[i__ + j * wl_dim1] < minw[i__]) {
            minw[i__] = wl[i__ + j * wl_dim1];
         }
         /* L90: */
      }
      if (minw[i__] + mw[i__] <= kw[i__]) {
         goto L100;
      }
      kk = 1;
      kchan[i__] = 1;
   L100:
      kwr[i__] = kw[i__] - minw[i__];
      /* L110: */
   }
   if (mw[io] <= kw[io]) {
      goto L120;
   }
   kk = 1;
   kchan[io] = 1;
L120:
   if (kk == 0) {
      goto L70;
   }
   i__1 = nb;
   for (jj = 1; jj <= i__1; ++jj) {
      j = bb[jj];
      jf = first[j];
      if (pen[j] < (double) (*inf)) {
         goto L130;
      }
      if (wl[jf + j * wl_dim1] > kw[jf]) {
         goto L200;
      }
      goto L190;
   L130:
      if (kchan[jf] == 0) {
         goto L140;
      }
      if (wl[jf + j * wl_dim1] > kwr[jf]) {
         goto L150;
      }
   L140:
      js = second[j];
      if (kchan[js] == 0) {
         goto L190;
      }
      if (wl[js + j * wl_dim1] <= kwr[js]) {
         goto L190;
      }
   L150:
      fmax = (double) (-(*inf));
      smax = (double) (-(*inf));
      if__ = 0;
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
         if (wl[i__ + j * wl_dim1] > kw[i__]) {
            goto L170;
         }
         if (wl[i__ + j * wl_dim1] > kwr[i__]) {
            wl[i__ + j * wl_dim1] = kw[i__];
         }
         rwl = (double) wl[i__ + j * wl_dim1];
         rp = (double) p[i__ + j * p_dim1];
         rkw = (double) kw[i__];
         s = (-(*a1) * rwl + *a2 * rp) / (*a3 * rkw + *a4 * rwl + *a5);
         if (smax >= s) {
            goto L170;
         }
         if (fmax >= s) {
            goto L160;
         }
         smax = fmax;
         is = if__;
         fmax = s;
         if__ = i__;
         goto L170;
      L160:
         smax = s;
         is = i__;
      L170:
         ;
      }
      first[j] = if__;
      second[j] = is;
      pen[j] = fmax - smax;
      if (wl[if__ + j * wl_dim1] > mw[if__]) {
         mw[if__] = wl[if__ + j * wl_dim1];
      }
      if (smax > (double) (-(*inf))) {
         goto L180;
      }
      pen[j] = (double) (*inf);
      goto L190;
   L180:
      if (wl[is + j * wl_dim1] > mw[is]) {
         mw[is] = wl[is + j * wl_dim1];
      }
   L190:
      ;
   }
   goto L70;
L200:
   *z__ = 0;
   return 0;
   /* TRY TO IMPROVE ON THE CURRENT SOLUTION Z. */
L210:
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      if__ = xstar[j];
      maxp = p[if__ + j * p_dim1];
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
         if (w[i__ + j * w_dim1] > kw[i__]) {
            goto L220;
         }
         if (p[i__ + j * p_dim1] <= maxp) {
            goto L220;
         }
         maxp = p[i__ + j * p_dim1];
         if__ = i__;
      L220:
         ;
      }
      ip = xstar[j];
      if (if__ == ip) {
         goto L230;
      }
      xstar[j] = if__;
      *z__ = *z__ + p[if__ + j * p_dim1] - p[ip + j * p_dim1];
      kw[ip] += w[ip + j * w_dim1];
      kw[if__] -= w[if__ + j * w_dim1];
   L230:
      ;
   }
   return 0;
} /* ghx_ */

int  MTHG::trin_(int *p, int *n, int *m, int *invst, int *lam, int *jdimr, int *jdimc)
{
   /* System generated locals */
   int p_dim1, p_offset, i__1, i__2;

   /* Local variables */
   static int i__, j, max__;


   /* TRANSFORM AN INSTANCE OF GAP IN MINIMIZATION FORM INTO AN */
   /* EQUIVALENT INSTANCE IN MAXIMIZATION FORM. */

   /* Parameter adjustments */
   p_dim1 = *jdimr;
   p_offset = 1 + p_dim1;
   p -= p_offset;

   /* Function Body */
   *invst = 0;
   max__ = 0;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
         if (p[i__ + j * p_dim1] > max__) {
            max__ = p[i__ + j * p_dim1];
         }
         /* L10: */
      }
      /* L20: */
   }
   *lam = max__ + 1;
   i__1 = *n;
   for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
         p[i__ + j * p_dim1] = *lam - p[i__ + j * p_dim1];
         /* L30: */
      }
      *invst += *lam;
      /* L40: */
   }
   return 0;
} /* trin_ */

/*
// THIS SUBROUTINE HEURISTICALLY SOLVES THE GENERALIZED ASSIGNMENT PROBLEM
double SimpleHeu::MTHG(int n, int m, double **p, int **w, int *c, int minmax, double *zub, int* xStar, int jck)
{
// OPT Z = P(1,1)*X(1,1) + ... + P(1,N)*X(1,N) +
//                         ...                 +
//         P(M,1)*X(M,1) + ... + P(M,N)*X(M,N)
//     (WHERE  OPT = MIN  if  MINMAX = 1 ,  OPT = MAX  if  MINMAX = 2 )
//
// SUBJECT TO:
//       W(I,1)*X(I,1) + ... + W(I,N)*X(I,N) <= C(I)  FOR I=1,...,M,
//       X(1,J) + ... + X(M,J) = 1                      FOR J=1,...,N,
//       X(I,J) = 0 OR 1                     FOR I=1,...,M, J=1,...,N.
//
// THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
//   1) 2 <= M <= JDIMR ;
//   2) 2 <= N <= JDIMC ;
//   3) P(I,J), W(I,J) AND C(I) POSITIVE intS;
//   4) W(I,J) <= C(I) FOR AT LEAST ONE I, FOR J=1,...,N;
//   5) C(I) >= MIN (W(I,J)) FOR I=1,...,M.
//
// MTHG CALLS 6 PROCEDURES: CHMTHG, FEAS, GHA, GHBCD, GHX AND TRIN.
// MTHG NEEDS
//   6 ARRAYS ( C ,  DMYR1 ,  DMYR2 ,  DMYR3 ,  DMYR4  AND  DMYR5 ) OF LENGTH AT LEAST  JDIMR ;
//   7 ARRAYS ( XSTAR ,  BEST ,  DMYC1 ,  DMYC2 ,  DMYC3 ,  DMYC4  AND DMYCR1 ) OF LENGTH AT LEAST  JDIMC ;
//   3 ARRAYS ( P ,  W  AND  A ) OF LENGTH AT LEAST  JDIMR X JDIMC .

// MEANING OF THE INPUT PARAMETERS:
// N        = NUMBER OF ITEMS;
// M        = NUMBER OF KNAPSACKS;
// P(I,J)   = PROFIT OF ITEM J if ASSIGNED TO KNAPSACK I (I=1,...,M; J=1,...,N);
// W(I,J)   = WEIGHT OF ITEM J if ASSIGNED TO KNAPSACK I (I=1,...,M; J=1,...,N);
// C(I)     = CAPACITY OF KNAPSACK I (I=1,...,M);
// MINMAX   = 1 if THE OBJECTIVE FUNCTION MUST BE MINIMIZED,
//          = 2 if THE OBJECTIVE FUNCTION MUST BE MAXIMIZED;
// JCK      = 1 if CHECK ON THE INPUT DATA IS DESIRED,
//          = 0 OTHERWISE.
//
// MEANING OF THE OUTPUT PARAMETERS:
// Z        = VALUE OF THE SOLUTION FOUND if Z > 0 ,
//          = 0 if NO FEASIBLE SOLUTION IS FOUND,
//          = ERROR IN THE INPUT DATA (WHEN JCK=1) if Z < 0 : CONDITION  - Z  IS VIOLATED;
// XSTAR(J) = KNAPSACK WHERE ITEM J IS INSERTED IN THE SOLUTION FOUND.
//
// ALL THE PARAMETERS ARE int. ON return OF MTHG ALL THE INPUT
// PARAMETERS ARE UNCHANGED, BUT  P(I,J)  IS SET TO  0  FOR ALL PAIRS
// (I,J)  SUCH THAT  W(I,J) > C(I) .

   int i,j,jfi;
   double lam,kvst,zm;

   // allocazioni
   
   int* best = (int *) malloc (n * sizeof(int)); // int BEST(500)
   double* dmycr1 = (double *) malloc (n * sizeof(double)); // double    DMYCR1(500)

   //int DMYR1(50),DMYR2(50),DMYR3(50),DMYR4(50),DMYR5(50)
   int* dmyr1 = (int *) malloc (m * sizeof(int)); 
   int* dmyr2 = (int *) malloc (m * sizeof(int)); 
   int* dmyr3 = (int *) malloc (m * sizeof(int)); 
   int* dmyr4 = (int *) malloc (m * sizeof(int)); 
   int* dmyr5 = (int *) malloc (m * sizeof(int)); 

   // int DMYC1(500),DMYC2(500),DMYC3(500),DMYC4(500)
   int* dmyc1 = (int *) malloc (n * sizeof(int)); 
   int* dmyc2 = (int *) malloc (n * sizeof(int)); 
   int* dmyc3 = (int *) malloc (n * sizeof(int)); 
   int* dmyc4 = (int *) malloc (n * sizeof(int)); 

   // int A(50,500)
   int** a = (int **) malloc (m * sizeof(int*));
   for(i = 0 ; i < m; ++i) a[i] = (int*) malloc(sizeof(int) * n);

   // DEFINITION OF THE INTERNAL PARAMETERS.
   int JDIMR = m;
   int JDIMC = n;
   double z = 0;
   if ( jck == 1 ) CHMTHG(n,m,p,w,c,JDIMR,JDIMC,&z);
   if ( z < 0 ) goto Lend;

   // INITIALIZE.
   double sumConst = 0;
   int IMULT = -1;
   if ( minmax < 2 ) 
   {
      // TRANSFORM THE MINIMIZATION PROBLEM INTO A MAXIMIZATION PROBLEM.
      TRIN(p,n,m,&sumConst,&lam); 
      int imult = 1;
   }

   // SOLVE THE MAXIMIZATION PROBLEM.

   // CHECK FOR INFEASIBILITY.
   FEAS(n,m,p,w,c,xStar,&jfi,JDIMR,JDIMC);
   if ( jfi == 1 ) goto L30;

   // FIRST HEURISTIC SOLUTION.
   GHA(p,w,c,n,m,&z,xStar,BEST,&kvst,INF,dmyr1,dmyr2,dmyc1,dmyc2,dmyc3,dmyc4);

   // SECOND HEURISTIC SOLUTION.
   //GHBCD(P,W,C,N,M,Z,XSTAR,INF,JDIMR,JDIMC,DMYC1,DMYR1,DMYR2,DMYR3,DMYR4,DMYR5,DMYC2,DMYC3,DMYC4,DMYCR1,A)

   // TERMINATE.
L20: zm = z;
   z = 0;
   if ( zm > kvst ) 
      z = sumConst - zm*IMULT;

L30: if ( minmax < 2 ) 
      // RE-STORE THE ORIGINAL MINIMIZATION PROBLEM.
      for(i=0; i<m; i++)
         for(j=0;j<n;j++)
            if ( p[i][j] > 0 ) 
               p[i][j] = lam - p[i][j];

   // deallocations
   if(best != NULL) free(best);
   if(dmycr1 != NULL) free(dmycr1);
   if(dmyr1 != NULL) free(dmyr1);
   if(dmyr2 != NULL) free(dmyr2);
   if(dmyr3 != NULL) free(dmyr3);
   if(dmyr4 != NULL) free(dmyr4);
   if(dmyr5 != NULL) free(dmyr5);
   if(dmyc1 != NULL) free(dmyc1);
   if(dmyc2 != NULL) free(dmyc2);
   if(dmyc3 != NULL) free(dmyc3);
   if(dmyc4 != NULL) free(dmyc4);

Lend: *zub = z;
   return z;
}

// CHECK THE INPUT DATA.
void SimpleHeu::CHMTHG(int n, int m, double** p, int **w, int *c, int JDIMR, int JDIMC, double *zub)
{  int i,j,min;
   double z;

   if ( m <= 1 )    z = - 1;
   if ( m > JDIMR ) z = - 1;
   if ( z < 0  ) goto Lend;
   if ( n <= 1 )    z = - 2;
   if ( n > JDIMC ) z = - 2;
   if ( z < 0  ) goto Lend;

   for(i=0;i<m;i++)
   {
      if ( c[i] <= 0 )
      {  z = - 3;
         goto Lend;
      }

      min = c[i] + 1;
      for(j=0;j<n;j++)
      {  if ( p[i][j] <= 0 || w[i][j] <= 0 ) 
         {  z = - 3;
            goto Lend;
         }
         if ( w[i][j] < min ) min = w[i][j];
      }
      if(c[i] < min ) z = - 5;
   }

   for(j=0;j<n;j++)
      for(i=0;i<m;i++)
         if ( w[i][j] > c[i] ) 
         {  z = - 4;
            goto Lend;
         }

Lend: *zub = z;
   return;
}

// CHECK FOR INFEASIBILITY.
void SimpleHeu::FEAS(int n, int m, double** p, int **w, int *c, int* xStar, int* jfi, int JDIMR, int JDIMC)
{  int i,j,numInfeas;

   *jfi = 0;
   for(j=0;j<n;j++)
   {
      xStar[j] = 0;
      numInfeas = 0;
      for(i=0;i<m;i++)
         if ( w[i][j] > c[i] ) 
         {  numInfeas = numInfeas + 1;
            p[i][j] = 0;
         }

      if ( numInfeas == m ) *jfi = 1;
   }
   return;
}

// TRANSFORM AN INSTANCE OF GAP IN MINIMIZATION FORM INTO AN EQUIVALENT INSTANCE IN MAXIMIZATION FORM.
void SimpleHeu::TRIN(double** p, int n, int m, double *sumConst, double *lam)
{  int i,j;

   *sumConst = 0;
   double maxVal = 0;
   for(j=0;j<n;j++)
      for(i=0;i<m;i++)
         if ( p[i][j] > maxVal ) maxVal = p[i][j];

   *lam = maxVal + 1;
   for(j=0;j<n;j++)
   {  for(i=0;i<m;i++)
         p[i][j] = *lam - p[i][j];

      *sumConst += *lam;
   }
}

// APPLY THE APPROXIMATE ALGORITHM GH WITH FUNCTION (A) AND DEFINE THE INFINITE VALUE INF
void SimpleHeu::GHA(double** p, int **w, int *c, int n, int m, double *zub, int* XSTAR, best, double* KVST,INF,KW,MW,PEN,FIRST,SECOND,BB)
{
   // if Z = KVST NO FEASIBLE SOLUTION WAS FOUND.

   int* KW = (int *) malloc (m * sizeof(int)); 
   int* MW = (int *) malloc (m * sizeof(int)); 

   int* PEN = (int *) malloc (n * sizeof(int)); 
   int* FIRST = (int *) malloc (n * sizeof(int)); 
   int* SECOND = (int *) malloc (n * sizeof(int)); 
   int* BB = (int *) malloc (n * sizeof(int)); 

   int i,j,jj,iFirst,iSecond,idFirst,newSecond;
   double Z,pMin,maxRegret;

   int pMax,secondMax,NB;
   double INF = 0;
   for(i=0;i<m;i++)
   {
      KW[i] = c[i];
      MW[i] = 0;
      if ( c[i] > INF ) INF = c[i];
   }
   Z = 0;
   KVST = 0;
   for(j=0;j<n;j++)
   {
      pMin = p[0][j];
      pMax = p[0][j];
      iFirst = 1;
      secondMax = 0;
      for(i=1;i<m;i++)
      {
         if ( p[i][j] < pMin ) pMin = p[i][j];
         if ( secondMax >= p[i][j] ) continue;
         if ( pMax < p[i][j] ) 
         {
            secondMax = pMax;
            iSecond = iFirst;
            pMax = p[i][j];
            iFirst = i;
            continue;
         }
         secondMax = p[i][j];
         iSecond = i;
      }

      KVST = KVST + pMin;
      FIRST[j] = iFirst;
      BEST[j] = iFirst;
      SECOND[j] = iSecond;
      PEN[j] = pMax - secondMax;
      if ( secondMax == 0 ) PEN[j] = - 1;
      BB[j] = j;
      if ( w[iFirst][j] > MW[iFirst] ) MW[iFirst] = w[iFirst][j];
      if ( w[iSecond][j] > MW[iSecond] ) MW[iSecond] = w[iSecond][j];
   }

   if ( KVST > 0 ) KVST = KVST - 1;
   for(j=0;j<n;j++)
      if ( PEN[j] == (- 1) ) 
         PEN[j] = INF;

   NB = n;

L60: ;
   do
   {  maxRegret = - 1;
      for(int jj=0;jj<NB;jj++)
      {
         j = BB[jj];
         if ( PEN[j] > maxRegret ) 
         {
            maxRegret = PEN[j];
            JJM = jj;
         }
      }

      JO = BB[JJM];
      IO = FIRST[JO];
      Z = Z + p[IO][JO];
      XSTAR[JO] = IO;
      BB[JJM] = BB[NB]
      NB = NB - 1;
      if ( NB == 0 ) goto Lend;
      KW[IO] = KW[IO] - w[IO][JO];
   }
   while ( MW[IO] <= KW[IO] );

   for(jj=0;jj<NB;jj++)
   {
      j = BB[jj];
      if ( w[IO][j] <= KW[IO] ) continue;
      if ( FIRST[j] == IO ) 
      {  if ( PEN[j] == INF )
         {  Z = KVST;
            goto Lend;                          // exiting
         }
         FIRST[j] = SECOND[j];
      }
      else
         if ( SECOND[j] != IO ) continue;
      idFirst = FIRST[j];
      w[idFirst][j] = w[idFirst][j] + INF;
      newSecond = 0;
      for(i=0;i<m;i++)
      {
         if ( w[i][j] > KW[i]      ) continue;
         if ( p[i][j] <= newSecond ) continue;
         newSecond = p[i][j];
         iSecond = i;
      }

      w[idFirst][j] = w[idFirst][j] - INF;
      if ( newSecond != 0 ) 
      {
         SECOND[j] = iSecond;
         PEN[j] = p[idFirst][j] - newSecond;
         if ( w[iSecond][j] > MW[iSecond] ) MW[iSecond] = w[iSecond][j];
         continue;
      }
      PEN[j] = INF;
   }

   goto L60;
   Z = KVST;

Lend: *zub = Z;
   if(KW != NULL) free(KW);
   if(MW != NULL) free(MW);
   if(PEN != NULL) free(PEN);
   if(FIRST != NULL) free(FIRST);
   if(SECOND != NULL) free(SECOND);
   if(BB != NULL) free(BB);
   return;
}

// APPLY THE APPROXIMATE ALGORITHM GH WITH FUNCTIONS (B), (C) AND (D).
double SimpleHeu::GHBCD(P,W,C,N,M,Z,XSTAR,INF,JDIMR,JDIMC,XSP, DMYR1,DMYR2,DMYR3,DMYR4,DMYR5,DMYC2,DMYC3,DMYC4,DMYCR1,DMYA)
{
   int P(JDIMR,JDIMC),W(JDIMR,JDIMC),C(JDIMR),XSTAR(JDIMC),Z
   int VSP,XSP(JDIMC)
   int DMYR1(JDIMR),DMYR2(JDIMR),DMYR3(JDIMR),DMYR4(JDIMR),
   1        DMYR5(JDIMR),DMYC2(JDIMC),DMYC3(JDIMC),DMYC4(JDIMC),
   2        DMYA(JDIMR,JDIMC)
   double    DMYCR1(JDIMC)
   JJ = 2
   A1 = 1.
   A2 = 0.
   A3 = 0.
   A4 = 0.
   A5 = 1.
   10 CALL GHX(P,W,C,N,M,VSP,XSP,A1,A2,A3,A4,A5,INF,JDIMR,JDIMC,
      1         DMYR1,DMYR2,DMYR3,DMYR4,DMYR5,DMYC2,DMYC3,DMYC4,
      2         DMYCR1,DMYA)
   if ( VSP <= Z ) goto 30
   Z = VSP
   DO 20 J=1,N
   XSTAR(J) = XSP(J)
   20 CONTINUE
   30 if ( JJ == 3 ) goto 40
   if ( JJ == 4 ) goto 50
   JJ = 3
   A1 = 1.
   A2 = 0.
   A3 = 1.
   A4 = 0.
   A5 = 0.
   goto 10
   40 JJ = 4
   A1 = 0.
   A2 = 1.
   A3 = 0.
   A4 = 1.
   A5 = 0.
   goto 10
   50 return
}

// APPLY THE APPROXIMATE ALGORITHM GH WITH FUNCTION (B) OR (C) OR (D).
double SimpleHeu::GHX(P,W,C,N,M,Z,XSTAR,A1,A2,A3,A4,A5,INF,JDIMR,JDIMC,KW,MW,MINW,KCHAN,KWR,FIRST,SECOND,BB,PEN,WL)
{
   int P(JDIMR,JDIMC),W(JDIMR,JDIMC),XSTAR(JDIMC),C(JDIMR),Z
   int KW(JDIMR),MW(JDIMR),MINW(JDIMR),KCHAN(JDIMR),KWR(JDIMR),
   1        FIRST(JDIMC),SECOND(JDIMC),BB(JDIMC),WL(JDIMR,JDIMC)
   double    PEN(JDIMC),MAXPEN
   DO 20 I=1,M
   KW(I) = C(I)
   MW(I) = 0
   MINW(I) = INF
   DO 10 J=1,N
   WL(I,J) = W(I,J)
   if ( WL(I,J) < MINW(I) ) MINW(I) = WL(I,J)
   10   CONTINUE
   KWR(I) = KW(I) - MINW(I)
   20 CONTINUE
   Z = 0
   DO 60 J=1,N
   FMAX = - INF
   if = 0
   SMAX = - INF
   DO 40 I=1,M
   if ( WL(I,J) > KW(I) ) goto 40
   if ( WL(I,J) > KWR(I) ) WL(I,J) = KW(I)
   RWL = WL(I,J)
   RP = P(I,J)
   RKW = KW(I)
   S = (- A1*RWL + A2*RP)/(A3*RKW + A4*RWL + A5)
   if ( SMAX >= S ) goto 40
   if ( FMAX >= S ) goto 30
   SMAX = FMAX
   IS = if
   FMAX = S
   if = I
   goto 40
   30     SMAX = S
   IS = I
   40   CONTINUE
   FIRST(J) = if
   SECOND(J) = IS
   PEN(J) = FMAX - SMAX
   BB(J) = J
   if ( WL(if,J) > MW(if) ) MW(if) = WL(if,J)
   if ( SMAX > FLOAT(- INF) ) goto 50
   PEN(J) = INF
   goto 60
   50   if ( WL(IS,J) > MW(IS) ) MW(IS) = WL(IS,J)
   60 CONTINUE
   NB = N
   70 MAXPEN = - 1
   DO 80 JJ=1,NB
   J = BB(JJ)
   if ( PEN(J) <= MAXPEN ) goto 80
   MAXPEN = PEN(J)
   JJM = JJ
   80 CONTINUE
   JO = BB(JJM)
   IO = FIRST(JO)
   Z = Z + P(IO,JO)
   XSTAR(JO) = IO
   BB(JJM) = BB(NB)
   NB = NB - 1
   KW(IO) = KW(IO) - W(IO,JO)
   if ( NB == 0 ) goto 210
   KK = 0
   DO 110 I=1,M
   KCHAN(I) = 0
   if ( WL(I,JO) > MINW(I) ) goto 100
   MINW(I) = INF
   DO 90 JJ=1,NB
   J = BB(JJ)
   if ( WL(I,J) < MINW(I) ) MINW(I) = WL(I,J)
   90   CONTINUE
   if ( MINW(I) + MW(I) <= KW(I) ) goto 100
   KK = 1
   KCHAN(I) = 1
   100   KWR(I) = KW(I) - MINW(I)
   110 CONTINUE
   if ( MW(IO) <= KW(IO) ) goto 120
   KK = 1
   KCHAN(IO) = 1
   120 if ( KK == 0 ) goto 70
   DO 190 JJ=1,NB
   J = BB(JJ)
   JF = FIRST(J)
   if ( PEN(J) < FLOAT(INF) ) goto 130
   if ( WL(JF,J) > KW(JF) ) goto 200
   goto 190
   130   if ( KCHAN(JF) == 0 ) goto 140
   if ( WL(JF,J) > KWR(JF) ) goto 150
   140   JS = SECOND(J)
   if ( KCHAN(JS) == 0 ) goto 190
   if ( WL(JS,J) <= KWR(JS) ) goto 190
   150   FMAX = - INF
   SMAX = - INF
   if = 0
   DO 170 I=1,M
   if ( WL(I,J) > KW(I) ) goto 170
   if ( WL(I,J) > KWR(I) ) WL(I,J) = KW(I)
   RWL = WL(I,J)
   RP = P(I,J)
   RKW = KW(I)
   S = (- A1*RWL + A2*RP)/(A3*RKW + A4*RWL + A5)
   if ( SMAX >= S ) goto 170
   if ( FMAX >= S ) goto 160
   SMAX = FMAX
   IS = if
   FMAX = S
   if = I
   goto 170
   160     SMAX = S
   IS = I
   170   CONTINUE
   FIRST(J) = if
   SECOND(J) = IS
   PEN(J) = FMAX - SMAX
   if ( WL(if,J) > MW(if) ) MW(if) = WL(if,J)
   if ( SMAX > FLOAT(- INF) ) goto 180
   PEN(J) = INF
   goto 190
   180   if ( WL(IS,J) > MW(IS) ) MW(IS) = WL(IS,J)
   190 CONTINUE
   goto 70
   200 Z = 0
   return
   C TRY TO IMPROVE ON THE CURRENT SOLUTION Z.
   210 DO 230 J=1,N
   if = XSTAR(J)
   MAXP = P(if,J)
   DO 220 I=1,M
   if ( W(I,J) > KW(I) ) goto 220
   if ( P(I,J) <= MAXP ) goto 220
   MAXP = P(I,J)
   if = I
   220   CONTINUE
   IP = XSTAR(J)
   if ( if == IP ) goto 230
   XSTAR(J) = if
   Z = Z + P(if,J) - P(IP,J)
   KW(IP) = KW(IP) + W(IP,J)
   KW(if) = KW(if) - W(if,J)
   230 CONTINUE
}
*/