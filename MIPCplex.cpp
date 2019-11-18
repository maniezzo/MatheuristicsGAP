#include "MIPCplex.h"
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <ilcplex/cplex.h>

// This simple routine frees up the pointer *ptr, and sets *ptr to NULL 
static void free_and_null (char **ptr)
{  if ( *ptr != NULL ) 
   {  try
      { free (*ptr); }
      catch (const std::exception& e) 
      {  cout << e.what(); }
      *ptr = NULL;
   }
} // END free_and_null  

// callback, records the best solution so far
static int CPXPUBLIC logcallback (CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle)
{  int j,status = 0;

   LOGINFOptr info = (LOGINFOptr) cbhandle;
   int        hasincumbent = 0, newincumbent = 0, nodecnt, nodesleft;
   double     objval=INT_MAX;
   double     bound;
   double     *x = NULL;

   status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_MIP_FEAS, &hasincumbent);
   if ( status )  goto TERMINATE;

   if ( hasincumbent ) // a feasible solution exists
   {  // cost so far
      status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &objval);
      if ( status )  goto TERMINATE;
      if ( fabs(info->lastincumbent - objval) > 0.1)  
      {  newincumbent = 1;
         info->lastincumbent = objval;
      }
   }

   status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &nodecnt);
   if ( status )  goto TERMINATE;

   // output to console every 1000 new nodes or when a new incumbent is met
   if ( nodecnt > 0 && (nodecnt >= info->lastlog + 10000  ||  newincumbent) ) 
   {
      if ( !newincumbent )  info->lastlog = nodecnt;

      // these are just for the console output
      status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &bound);
      if ( status )  goto TERMINATE;

      status = CPXgetcallbackinfo (env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODES_LEFT, &nodesleft);
      if ( status )  goto TERMINATE;

      cout << nodecnt <<"("<< nodesleft << ") bound = " << bound << " -> " << (newincumbent ? "NEW" : "" ) << " incumbent objective = " << objval << endl;

      int numcols = info->numcols;
      x = (double *) malloc (numcols*sizeof(double));
      if ( x == NULL ) {
         status = CPXERR_NO_MEMORY;
         goto TERMINATE;
      }
      status = CPXgetcallbackincumbent (env, cbdata, wherefrom, x, 0, numcols-1);
      if ( status )  goto TERMINATE;

      for (j = 0; j < numcols; j++) 
         if ( fabs(x[j]) > 0.5 ) 
            info->lastsol[j] = 1;
         else
            info->lastsol[j] = 0;
   }

TERMINATE:
   free_and_null ((char **) &x);
   return status;
} /* END callback */

MIPCplex::MIPCplex()
{
   //ctor
}

MIPCplex::MIPCplex(int numr, int numc, int numnonz)
{
   //ctor
   status    = 0;
   numCols   = numc;
   numRows   = numr;
   numNZ     = numnonz;
   obj = (double *) malloc (numCols * sizeof(double));
   lb  = (double *) malloc (numCols * sizeof(double));
   ub  = (double *) malloc (numCols * sizeof(double));
   rmatbeg = (int *) malloc (numRows * sizeof(int));
   rmatind = (int *) malloc (numNZ *   sizeof(int));
   rmatval = (double *) malloc (numNZ * sizeof(double));
   rhs     = (double *) malloc (numRows * sizeof(double));
   sense   = (char *)   malloc (numRows * sizeof(char));
   ctype   = (char *)   malloc(numCols * sizeof(char));
   lptype  = (char *)   malloc(numCols * sizeof(char));
   rowname = (char **)  malloc (numRows * sizeof(char*));
   //for(int i = 0 ; i < numRows; ++i) rowname[i] = (char*) malloc(sizeof(char) * sizeOfString);
   colname = (char **)  malloc (numCols * sizeof(char*));
   dj = NULL;
   pi = NULL;
   slack = NULL;
}

MIPCplex::~MIPCplex()
{
   //dtor
   if(obj != NULL)     free(obj);
   if(lb != NULL)      free(lb);
   if(ub != NULL)      free(ub);
   if(rmatbeg != NULL) free(rmatbeg);
   if(rmatind != NULL) free(rmatind);
   if(rmatval != NULL) free(rmatval);
   if(rhs != NULL)     free(rhs);
   if(sense != NULL)   free(sense);
   if(ctype != NULL)   free(ctype);
   if(lptype != NULL)  free(lptype);
   if(colname != NULL) free(colname);
   if(rowname != NULL) free(rowname);
   if(x != NULL)       free(x);
   //if(slack != NULL)   free(slack);
}

// To populate by row, we first create the columns, and then add the rows.  
int MIPCplex::populatebyrow (CPXENVptr env, CPXLPptr lp)
{  int i,j,ij,idRow,numNZrow;

   CPXchgobjsen (env, lp, CPX_MIN);  // Problem is minimization 

   // Create the new columns.
   ij = 0;
   for(i=0;i<GAP->m;i++)
      for(j=0;j<GAP->n;j++)
      {  obj[ij] = GAP->c[i][j];
         lb[ij]  = 0;
         ub[ij]  = 1.0;
         ctype[ij]   = 'B'; // 'B', 'I','C' to indicate binary, general integer, continuous 
         lptype[ij]  = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
         colname[ij] = (char *) malloc(sizeof(char) * (11));   // why not 11?
         sprintf(colname[ij],"%s%d","x",ij);
         ij++;
      }

   status = CPXnewcols (env, lp, numCols, obj, lb, ub, NULL, colname);
   if ( status )  goto TERMINATE;

   // Assignment constraints.  
   idRow = 0;
   for(j=0;j<GAP->n;j++)
   {
      numNZrow = 0;  // number of nonzero element in the row to add
      rmatbeg[0] = 0;     
      sense[0] = 'E';
      rhs[0]   = 1.0;
      rowname[idRow] = (char *) malloc(sizeof(char) * (11));   // why not 11?
      sprintf(rowname[0],"%s%d","c",idRow);
      for(i=0;i<GAP->m;i++)
      {
         rmatind[i] = j+ GAP->n*i;
         rmatval[i] = 1.0;
         numNZrow++;
      }
      rmatbeg[1] = numNZrow;
      status = CPXaddrows (env, lp, 0, 1, numNZrow, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
      if ( status )  goto TERMINATE;
      idRow++;
   }

   // Capacity constraints.  
   for(i=0;i<GAP->m;i++)
   {
      numNZrow = 0;  // number of nonzero element in the row to add
      rmatbeg[0] = 0;     
      sense[0] = 'L';
      rhs[0]   = GAP->cap[i];
      rowname[idRow] = (char *) malloc(sizeof(char) * (11));   // why not 11?
      sprintf(rowname[0],"%s%d","c",idRow);
      for(j=0;j<GAP->n;j++)
      {
         rmatind[numNZrow] = j+ GAP->n*i;
         rmatval[numNZrow] = GAP->req[i][j];
         numNZrow++;
      }
      rmatbeg[1] = numNZrow;
      status = CPXaddrows (env, lp, 0, 1, numNZrow, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
      if ( status )  goto TERMINATE;
      idRow++;
   }

TERMINATE:
   if( status )
      cout << " >>>>>>>>>>>> Error in constraint definition <<<<<<<<<<<<<<<" << endl;
   return (status);

}  // END populatebyrow 

// Populate by row a subproblem
int MIPCplex::populatebyrow (CPXENVptr env, CPXLPptr lp, int** c, int n, int m, int** req, int*  cap)
{  int i,j,ij,idRow,numNZrow;

   CPXchgobjsen (env, lp, CPX_MIN);  // Problem is minimization 

   // Create the new columns.
   ij = 0;
   for(i=0;i<m;i++)
      for(j=0;j<n;j++)
      {  obj[ij] = c[i][j];
         lb[ij]  = 0;
         ub[ij]  = 1.0;
         ctype[ij]   = 'B'; // 'B', 'I','C' to indicate binary, general integer, continuous 
         lptype[ij]  = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
         colname[ij] = (char *) malloc(sizeof(char) * (11));   // why not 11?
         sprintf(colname[ij],"%s%d","x",ij);
         ij++;
      }

   status = CPXnewcols (env, lp, numCols, obj, lb, ub, NULL, colname);
   if ( status )  goto TERMINATE;

   // Assignment constraints.  
   idRow = 0;
   for(j=0;j<n;j++)
   {
      numNZrow = 0;  // number of nonzero element in the row to add
      rmatbeg[0] = 0;     
      sense[0] = 'E';
      rhs[0]   = 1.0;
      rowname[idRow] = (char *) malloc(sizeof(char) * (11));   // why not 11?
      sprintf(rowname[0],"%s%d","c",idRow);
      for(i=0;i<m;i++)
      {
         rmatind[i] = j+ n*i;
         rmatval[i] = 1.0;
         numNZrow++;
      }
      rmatbeg[1] = numNZrow;
      status = CPXaddrows (env, lp, 0, 1, numNZrow, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
      if ( status )  goto TERMINATE;
      idRow++;
   }

   // Capacity constraints.  
   for(i=0;i<m;i++)
   {
      numNZrow = 0;  // number of nonzero element in the row to add
      rmatbeg[0] = 0;     
      sense[0] = 'L';
      rhs[0]   = cap[i];
      rowname[idRow] = (char *) malloc(sizeof(char) * (11));   // why not 11?
      sprintf(rowname[0],"%s%d","c",idRow);
      for(j=0;j<n;j++)
      {
         rmatind[numNZrow] = j+ n*i;
         rmatval[numNZrow] = req[i][j];
         numNZrow++;
      }
      rmatbeg[1] = numNZrow;
      status = CPXaddrows (env, lp, 0, 1, numNZrow, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
      if ( status )  goto TERMINATE;
      idRow++;
   }

TERMINATE:
   if( status )
      cout << " >>>>>>>>>>>> Error in constraint definition <<<<<<<<<<<<<<<" << endl;
   return (status);

}  // END populatebyrow 

// To populate by row, we first create the columns, and then add the rows.  
int MIPCplex::populateDual(CPXENVptr env, CPXLPptr lp, vector<int> xbar, int nRows, int nCols)
{
   int i, j, ij, idRow, numNZrow;

   CPXchgobjsen(env, lp, CPX_MAX);  // Problem is maximization 

   // Create the new columns.
   ij = 0;

   for (j = 0; j<GAP->n; j++)    // columns for assignment duals, w'
   {
      obj[ij] = 1 -xbar[j];
      lb[ij] = DBL_MIN;
      ub[ij] = DBL_MAX;
      ctype[ij]  = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
      lptype[ij] = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
      colname[ij] = (char *)malloc(sizeof(char) * (7));   // why not 7?
      sprintf(colname[ij], "%s%d", "w'", ij);
      ij++;
   }

   int istart = GAP->m - (nCols - GAP->n);

   for (i = istart; i<GAP->m; i++)    // columns for capcity duals, w"
   {
      obj[ij] = GAP->cap[i];
      lb[ij] = DBL_MIN;
      ub[ij] = 0.0;
      ctype[ij]  = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
      lptype[ij] = 'C'; // 'B', 'I','C' to indicate binary, general integer, continuous 
      colname[ij] = (char *)malloc(sizeof(char) * (7));   // why not?
      sprintf(colname[ij], "%s%d", "w''", ij);
      ij++;
   }

   status = CPXnewcols(env, lp, nCols, obj, lb, ub, NULL, colname);
   if (status)  goto TERMINATE;

   // Dual constraints.  
   idRow = ij = 0;
   for (i = istart; i<GAP->m; i++)
   {  
      for (j = 0; j<GAP->n; j++)
      {
         numNZrow = 0;  // number of nonzero element in the row to add
         rmatbeg[0] = 0;
         sense[0] = 'L';
         rhs[0] = GAP->c[i][j];
         rowname[idRow] = (char *)malloc(sizeof(char) * (7));   // why not?
         sprintf(rowname[0], "%s%d", "c", idRow);

         // var w'
         rmatind[numNZrow] = j;
         rmatval[numNZrow] = 1.0;
         numNZrow++;

         // var w''
         rmatind[numNZrow] = (i-istart) + GAP->n;
         rmatval[numNZrow] = GAP->req[i][j];
         numNZrow++;

         rmatbeg[1] = numNZrow;
         status = CPXaddrows(env, lp, 0, 1, numNZrow, rhs, sense, rmatbeg, rmatind, rmatval, NULL, rowname);
         if (status)  goto TERMINATE;
         idRow++;
      }
   }

TERMINATE:
   if (status)
      cout << " >>>>>>>>>>>> Error in constraint definition <<<<<<<<<<<<<<<" << endl;
   return (status);

}  // END populatedual

// Just as AllocateMIP but specialized to instantiate dual (could have done with an if inside)
int MIPCplex::allocateDual(int numRows, int numCols, vector<int> xbar, bool isVerbose)
{
   env = NULL;
   lp  = NULL;
   status = 0;

   pi    = NULL;
   slack = NULL;
   dj = NULL;
   x  = NULL;

   cout << "Allocating CPLEX" << endl;

   // Initialize the CPLEX environment 
   env = CPXopenCPLEX(&status);
   /* A call to CPXgeterrorstring will produce the text of
   the error message.  For other CPLEX routines, the errors will
   be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */
   if (env == NULL)
   {
      char  errmsg[1024];
      cerr << "Could not open CPLEX environment.\n";
      CPXgeterrorstring(env, status, errmsg);
      cerr << errmsg;
      goto lend;
   }

   // Turn on output to the screen 
   if (isVerbose)
      status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
   else
      status = CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
   if (status)
   {
      cerr << "Failure to turn on screen indicator, error " << status << endl;
      goto lend;
   }

   // Turn on data checking 
   status = CPXsetintparam(env, CPX_PARAM_DATACHECK, CPX_ON);
   if (status)
   {
      cerr << "Failure to turn on data checking, error " << status << endl;
      goto lend;
   }

   // Create the problem. 
   lp = CPXcreateprob(env, &status, "GAPinst");
   if (lp == NULL)
   {
      cerr << "Failed to create LP.\n";
      goto lend;
   }

   // Populate the problem with the GAP data.
   status = populateDual(env, lp, xbar, numRows, numCols);
   if (status)
   {
      cerr << "Failed to populate problem.\n";
      goto lend;
   }

   // Write a copy of the problem to a file. 
   if (isVerbose)
   {
      status = CPXwriteprob(env, lp, "GAPdual.lp", NULL);
      if (status)
      {
         cerr << "Failed to write problem to disk.\n";
         goto lend;
      }
   }

   return status;
lend:
   freeMIP();     // releases menory structures
   return status;
}

// allocates and initializes the cplex environment
int MIPCplex::allocateMIP(bool isVerbose)
{  
   env    = NULL;
   lp     = NULL;
   status = 0;

   pi    = NULL;
   slack = NULL;
   dj    = NULL;  
   x     = NULL;

   cout << "Allocating CPLEX" << endl;

   // Initialize the CPLEX environment 
   env = CPXopenCPLEX (&status);
   /* A call to CPXgeterrorstring will produce the text of
   the error message.  For other CPLEX routines, the errors will
   be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */
   if ( env == NULL ) 
   {  char  errmsg[1024];
      cerr << "Could not open CPLEX environment.\n";
      CPXgeterrorstring (env, status, errmsg);
      cerr << errmsg;
      goto lend;
   }

   // Turn on output to the screen 
   if(isVerbose)
      status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
   else
      status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) 
   {  cerr << "Failure to turn on screen indicator, error " << status << endl;
      goto lend;
   }

   // Turn on data checking 
   status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
   if ( status ) 
   {  cerr << "Failure to turn on data checking, error " << status << endl;
      goto lend;
   }

   // Create the problem. 
   lp = CPXcreateprob (env, &status, "GAPinst");
   if ( lp == NULL ) 
   {  cerr << "Failed to create LP.\n";
      goto lend;
   }

   // Populate the problem with the GAP data.
   status = populatebyrow (env, lp);
   if ( status ) 
   {  cerr << "Failed to populate problem.\n";
      goto lend;
   }

   return status;
lend:
   freeMIP();     // releases menory structures
   return status;
}

// allocates and initializes the cplex environment for a subproblem
int MIPCplex::allocateMIP(int** c, int n, int m, int** req, int* cap, bool isVerbose)
{  
   env    = NULL;
   lp     = NULL;
   status = 0;

   pi    = NULL;
   slack = NULL;
   dj    = NULL;  
   x     = NULL;

   if(isVerbose)
      cout << "Allocating CPLEX" << endl;

   // Initialize the CPLEX environment 
   env = CPXopenCPLEX (&status);
   /* A call to CPXgeterrorstring will produce the text of
   the error message.  For other CPLEX routines, the errors will
   be seen if the CPX_PARAM_SCRIND indicator is set to CPX_ON. */
   if ( env == NULL ) 
   {  char  errmsg[1024];
      cerr << "Could not open CPLEX environment.\n";
      CPXgeterrorstring (env, status, errmsg);
      cerr << errmsg;
      goto lend;
   }

   // Turn on output to the screen 
   if(isVerbose)
      status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);
   else
      status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
   if ( status ) 
   {  cerr << "Failure to turn on screen indicator, error " << status << endl;
      goto lend;
   }

   // Turn on data checking 
   status = CPXsetintparam (env, CPX_PARAM_DATACHECK, CPX_ON);
   if ( status ) 
   {  cerr << "Failure to turn on data checking, error " << status << endl;
      goto lend;
   }

   // Create the problem. 
   lp = CPXcreateprob (env, &status, "GAPsubinst");
   if ( lp == NULL ) 
   {  cerr << "Failed to create LP.\n";
      goto lend;
   }

   // Populate the problem with the GAP data.
   status = populatebyrow (env, lp, c, n, m, req, cap);
   if ( status ) 
   {  cerr << "Failed to populate problem.\n";
      goto lend;
   }

   return status;
lend:
   freeMIP();     // releases menory structures
   return status;
}

// Optimize the problem and obtain solution. 
int MIPCplex::solveMIP(bool fMIPint, bool isVerbose)
{  
   int i,j,cur_numrows,cur_numcols;

   // The actual size of the problem is in cur_numrows and cur_numcols
   cur_numrows = CPXgetnumrows(env, lp);
   cur_numcols = CPXgetnumcols(env, lp);

   if (xbest == NULL) xbest = (int *)   malloc(cur_numcols * sizeof(int));
   if (x == NULL)     x     = (double *)malloc(cur_numcols * sizeof(double));
   if (slack == NULL) slack = (double *)malloc(cur_numrows * sizeof(double));
   if (xbest == NULL || x == NULL || slack == NULL)
   {  status = CPXERR_NO_MEMORY;
      cerr << "Could not allocate memory for solution structures.\n";
      goto lend;
   }

   // Write a copy of the problem to a file. 
   if(isVerbose)
   {  status = CPXwriteprob (env, lp, "GAP.lp", NULL);
      if ( status ) 
      {  cerr << "Failed to write problem to disk.\n";
         goto lend;
      }
   }

   // linear case
   status = CPXchgprobtype (env, lp, CPXPROB_LP);
   if ( status ) 
   {  cerr << "Failed to set problem to LP.\n";
      goto lend;
   }

   // solve linear
   status = CPXlpopt (env, lp);
   if ( status ) 
   {  cerr << "Failed to optimize LP.\n";
      goto lend;
   }
   else
   {
      status = CPXsolution(env, lp, &solstat, &objval, x, pi, slack, dj);
      if (status || solstat>1)
      {  cerr << "Failed to obtain LP solution.\n";
         status = max(status,solstat);
         goto lreturn;
      }
      if(isVerbose)
      {  cout << "\nLP Solution status = " << status << endl;
         cout << "LP Solution value  = " << objval << endl;
      }

      if (dj == NULL) dj = (double *)malloc(cur_numcols * sizeof(double));
      if (pi == NULL) pi = (double *)malloc(cur_numrows * sizeof(double));
      if (dj == NULL || pi == NULL)
      {  status = CPXERR_NO_MEMORY;
         cerr << "Could not allocate memory for solution.\n";
         goto lend;
      }
      status = CPXgetpi(env, lp, pi, 0, CPXgetnumrows(env, lp) - 1);

      if(isVerbose)
      {  status = CPXgetlb(env, lp, lb, 0, cur_numcols - 1);
         cout << "LB:  "; for (j = 0; j < cur_numcols; j++) cout << lb[j] << ", "; cout << endl;
         status = CPXgetub(env, lp, ub, 0, cur_numcols - 1);
         cout << "UB:  "; for (j = 0; j < cur_numcols; j++) cout << ub[j] << ", "; cout << endl;
         cout << "sol: "; for (j = 0; j < cur_numcols; j++) cout << x[j] << ", "; cout << endl;
      }
   }

   // go for intergality
   if(fMIPint)
   {  // logging best solutions
      LOGINFO lastSol;
      lastSol.lastincumbent = objval = CPXgetobjsen (env, lp) * 1e+35;
      lastSol.lastsol = (int *) malloc (cur_numcols * sizeof(int));
      lastSol.lastlog = -100000;
      lastSol.numcols = cur_numcols;

      status = CPXsetinfocallbackfunc (env, logcallback, &lastSol); // mind the logcallback
      if ( status ) 
      {  cerr << "Failed to set logging callback function.\n";
         goto lend;
      }

      status = CPXcopyctype(env, lp, ctype);
      if ( status ) 
      {  cerr << "Failed to set integrality on vars.\n";
         goto lend;
      }

      status = CPXmipopt(env, lp);
      if ( status ) 
      {  cerr << "Failed to optimize MIP. zbest = " << lastSol.lastincumbent << endl;
         for(j=0;j<cur_numcols;j++) x[j] = lastSol.lastsol[j];
         goto lreturn;
      }

      status = CPXsolution (env, lp, &solstat, &objval, x, NULL, slack, NULL);
      if ( status ) 
      {  cerr << "Failed to obtain MIP solution.\n";
         goto lreturn;
      }
      cout << "\nMIP Solution status (101 optimal ok) = " << solstat << endl;
      cout << "MIP Solution value  = "   << objval << endl;
   }

   if(isVerbose)
   {
      int solnmethod, solntype, pfeasind, dfeasind;
      status = CPXsolninfo(env, lp, &solnmethod, &solntype, &pfeasind, &dfeasind);

      for (i = 0; i < cur_numrows; i++)
         cout << "Row: "<< i << " Slack = " << slack[i] << endl;
      for (j = 0; j < cur_numcols; j++) 
         if(x[j] > GAP->EPS)
            cout <<"Column " << j << " - " << colname[j] << "  Value = " << x[j] << endl;

      status = CPXgetlb(env, lp, lb, 0, cur_numcols - 1);
      status = CPXgetub(env, lp, ub, 0, cur_numcols - 1);

      cout << "LB:  "; for (j = 0; j < cur_numcols; j++) cout << lb[j] << ", "; cout << endl;
      cout << "UB:  "; for (j = 0; j < cur_numcols; j++) cout << ub[j] << ", "; cout << endl;
      cout << "sol: "; for (j = 0; j < cur_numcols; j++) cout << x[j]  << ", "; cout << endl;
      for(j=0;j<GAP->n;j++)
         for(i=0;i<GAP->m;i++)
            if(x[i*GAP->n + j] > 0.001)
               cout << i << " ";
      cout << endl;
   }

lreturn:
   return status;
lend:
   freeMIP();        // releases memory structures
   return status;
}

int MIPCplex::freeMIP()
{  
   // Free up the solution 
   // free_and_null ((char **) &x);    // in the destructor
   // if(slack != NULL) free_and_null ((char **) &slack);
   if(dj != NULL) free_and_null ((char **) &dj);
   if(pi != NULL) free_and_null ((char **) &pi);

   // Free up the problem as allocated by CPXcreateprob, if necessary 
   if ( lp != NULL ) 
   {  status = CPXfreeprob (env, &lp);
      if ( status ) 
         cerr << "CPXfreeprob failed, error code " << status << endl;
   }

   // Free up the CPLEX environment, if necessary 
   if ( env != NULL ) 
   {  status = CPXcloseCPLEX (&env);
      if ( status ) 
      {  char  errmsg[1024];
         cerr << "Could not close CPLEX environment.\n";
         CPXgeterrorstring (env, status, errmsg);
         cerr << errmsg << endl;
      }
   }

   if(lp != NULL )
      cout << "CPLEX released" << endl;
   return (status);
}

// this frees and fixes variables 
void MIPCplex::fixVariables(MIPCplex* CPX, vector<int> fixVal)
{
   int k, status, numfix;

   numfix = 0;      // numfix is the number of clients to fix
   for(k=0;k<fixVal.size();k++)
      if(fixVal[k] != INT_MAX) // if NAN, variable is free
         numfix++;

   int  cnt = GAP->n * GAP->m;
   int* indices = new int[cnt];  // indices of the columns corresponding to the variables for which bounds are to be changed
   char* lu = new char[cnt];     // whether the corresponding entry in the array bd specifies the lower or upper bound on column indices[j]
   double* bd = new double[cnt]; // new values of the lower or upper bounds of the variables present in indices

   for (k = 0; k<cnt; k++)
   {
      if (fixVal[k]==INT_MAX)           // set it free
      {
         if (CPX->lb[k] == 0.0)     // lb OK, only ub could be wrong
         {  CPX->ub[k] = 1.0;
         bd[k] = 1.0;
         lu[k] = 'U';            // bd[j] is an upper bound
         }
         else                       // lb wrong, ub should be OK
         {  CPX->lb[k] = 0.0;
         bd[k] = 0.0;
         lu[k] = 'L';            // bd[j] is an upper bound
         }
      }
      else                          // fix the var
      {
         CPX->lb[k] = fixVal[k];    // fixing in the solution
         bd[k] = fixVal[k];
         lu[k] = 'B';               // bd[k] is the lower and upper bound
      }
      indices[k] = k;
   }

   // change bounds in the model
   status = CPXchgbds(CPX->env, CPX->lp, cnt, indices, lu, bd);
   delete(indices);
   free(lu);
   free(bd);

   return;
}