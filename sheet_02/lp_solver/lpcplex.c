/****************************************************************************
    This file is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    It is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this file.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

/***  Interface for ILOG CPLEX 12.1  ***/

#include <stdlib.h>
#include <stdio.h>
/* #include <strings.h> */
#include <math.h>
#include <cplex.h>

#include "lp.h"

struct LINGOlp {
    CPXENVptr cplex_env;
    CPXLPptr  cplex_lp;
    int       noptcalls;
    int       ncols;
    int       nintegers;
     double   dbl_cutoff;
};

const double int_tolerance = 0.00001;

int LINGOlp_init (LINGOlp **p, const char *name)
{
   int rval = 0;
   int cpx_dbglvl = CPX_ON;

    (*p) = LINGO_SAFE_MALLOC (1, LINGOlp);
    if ((*p) == (LINGOlp *) NULL) {
        fprintf (stderr, "Out of memory in LINGOlp_init\n");
        rval = 1; goto CLEANUP;
    }

    (*p)->noptcalls = 0;
    (*p)->ncols     = 0;
    (*p)->nintegers = 0;
    (*p)->dbl_cutoff = 0.0;
    (*p)->cplex_env = (CPXENVptr) NULL;
    (*p)->cplex_lp = (CPXLPptr) NULL;

    (*p)->cplex_env = CPXopenCPLEX (&rval);
    if (rval) {
        fprintf (stderr, "CPXopenCPLEX failed, return code %d\n", rval);
        goto CLEANUP;
    }


    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_SCRIND,
                           cpx_dbglvl);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_SCRIND failed");


    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_ADVIND, 1);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_ADVIND failed");

    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_DPRIIND,
                           CPX_DPRIIND_STEEP);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_DPRIIND failed");

    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_PPRIIND,
                           CPX_PPRIIND_STEEP);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_PPRIIND failed");

    /* The following three parameters were set by Bix in TSP code */

    rval = CPXsetdblparam ((*p)->cplex_env, CPX_PARAM_EPPER, 1.0E-6);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_EPPER failed");

    rval = CPXsetdblparam ((*p)->cplex_env, CPX_PARAM_EPOPT, 1.0E-9);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_EPOPT failed");

    rval = CPXsetdblparam ((*p)->cplex_env, CPX_PARAM_EPRHS, 1.0E-9);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_EPRHS failed");

    rval = CPXsetintparam ((*p)->cplex_env, CPX_PARAM_THREADS,1);
    LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_THREADS failed");



    (*p)->cplex_lp = CPXcreateprob ((*p)->cplex_env, &rval, name);
    if (!(*p)->cplex_lp || rval) {
       fprintf (stderr, "CPXcreateprob failed, return code %d\n", rval);
       goto CLEANUP;
    }

CLEANUP:
    return rval;
}

void LINGOlp_free (LINGOlp **p)
{
     if (*p) {
        if ((*p)->cplex_env) {
            if ((*p)->cplex_lp) {
                CPXfreeprob ((*p)->cplex_env, &((*p)->cplex_lp));
            }
            CPXcloseCPLEX (&((*p)->cplex_env));
        }
        LINGO_FREE (*p, LINGOlp);
    }
}

static
int LINGOmip_optimize (LINGOlp *p)
{
   int rval   = 0;
   int status = 0;
   int solstat;

   status = CPXmipopt (p->cplex_env, p->cplex_lp);
   LINGOcheck_rval (status, "CPXmipopt failed");

   solstat = CPXgetstat (p->cplex_env, p->cplex_lp);
   if (solstat == CPX_STAT_INFEASIBLE) {
      printf ("Infeasible MIP\n");
      LINGOlp_write (p, "joethelion.lp");
      rval = 2;  goto CLEANUP;
   } else if (solstat != CPXMIP_OPTIMAL       &&
              solstat != CPXMIP_OPTIMAL_INFEAS  )
   {
      fprintf (stderr, "Cplex optimization status %d\n", solstat);
      goto CLEANUP;
   }

CLEANUP:
   return rval;
}
int LINGOlp_optimize (LINGOlp *p)
{
    int rval = 0;
    int solstat;

    if (p->nintegers) {
       rval = LINGOmip_optimize (p);
       LINGOcheck_rval (rval, "LINGOmip_optimize failed");
    } else {

       /** When restarting coloring with a set of stable sets (-r <filename>).
           The Simplex Method has a very long running time for solving the first LP.
           On C4000.5 the dual simplex would take several hours.
           In such situations we presolve the first LP with the Barrier Method.

           The steepest edge norms are not initialized properly but only
           to 1. Therefore, Barrier should not be used when there are
           only a few rows.
       */
       if (p->noptcalls == 0) {
          int ncols = CPXgetnumcols (p->cplex_env, p->cplex_lp);
          int nrows = CPXgetnumrows (p->cplex_env, p->cplex_lp);
          if (ncols > 5 * nrows) {
             rval = CPXbaropt (p->cplex_env, p->cplex_lp);
             LINGOcheck_rval (rval, "CPXbaropt failed");
          }
       }
       rval = CPXdualopt (p->cplex_env, p->cplex_lp);
       LINGOcheck_rval (rval, "CPXdualopt failed");


       solstat = CPXgetstat (p->cplex_env, p->cplex_lp);
       if (solstat == CPX_STAT_INFEASIBLE) {
          printf ("Infeasible LP\n");
          LINGOlp_write (p, "joethelion.lp");
          rval = 2;  goto CLEANUP;
       } else if (solstat != CPX_STAT_OPTIMAL       &&
                  solstat != CPX_STAT_OPTIMAL_INFEAS  ) {
          fprintf (stderr, "Cplex optimization status %d\n", solstat);
          if (solstat == CPX_STAT_ABORT_IT_LIM) {
             int itlim;
             rval = CPXgetintparam (p->cplex_env, CPX_PARAM_ITLIM, &itlim);
             if (!rval) {
                printf ("cplex iteration limit: %d\n", itlim);
                fflush (stdout);
             }
          }
          rval  = 1;  goto CLEANUP;
       }
    }
    (p->noptcalls)++;

    CLEANUP:
       return rval;
    }

int LINGOlp_objval (LINGOlp *p, double *obj)
{
    int rval = 0;

    rval = CPXgetobjval (p->cplex_env, p->cplex_lp, obj);
    LINGOcheck_rval (rval, "CPXgetobjval failed");

CLEANUP:
    return rval;
}

int LINGOlp_change_objective(LINGOlp *p, int start, int len, double* values)
{
    int i, rval = 0;
    int *indices = (int *) NULL;

    indices = LINGO_SAFE_MALLOC (len, int);
    LINGOcheck_NULL (indices, "out of memory for indices");
    for (i = 0; i < len; i++) {
        indices[i] = start+i;
    }

    rval = CPXchgobj (p->cplex_env, p->cplex_lp, len, indices, values);
    LINGOcheck_rval (rval, "CPXchgobj failed");

CLEANUP:
   LINGO_IFFREE (indices, int);
   return rval;
}

int LINGOlp_addrow (LINGOlp *p, int nzcount, int *cind, double *cval,
       char sense, double rhs, char *name)
{
    int rval = 0;
    char isense[1];
    char *iname[1];
    double irhs[1];
    int matbeg[1];


    switch (sense) {
    case LINGOlp_EQUAL:
        isense[0] = 'E'; break;
    case LINGOlp_LESS_EQUAL:
        isense[0] = 'L'; break;
    case LINGOlp_GREATER_EQUAL:
        isense[0] = 'G'; break;
    default:
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    irhs[0] = rhs;
    iname[0] = name;
    matbeg[0] = 0;

    if (nzcount == 0) {
        rval = CPXnewrows (p->cplex_env, p->cplex_lp, 1, irhs,
                       isense, (double *) NULL, (char **) NULL);
        LINGOcheck_rval (rval, "CPXnewrows failed");
    } else {
        rval = CPXaddrows (p->cplex_env, p->cplex_lp, 0, 1, nzcount, irhs,
                  isense, matbeg, cind, cval, (char **) NULL, (char **) NULL);
        LINGOcheck_rval (rval, "CPXaddrows failed");
    }

CLEANUP:
    return rval;
}

int LINGOlp_deleterows (LINGOlp *p, int first_cind, int last_cind)
{
    int rval = 0;

    rval = CPXdelrows (p->cplex_env, p->cplex_lp, first_cind, last_cind);
    LINGOcheck_rval (rval, "CPXdelrows failed");

CLEANUP:
   return rval;
}


int LINGOlp_addcol (LINGOlp *p, int nzcount, int *cind, double *cval,
       double obj, double lb, double ub, char sense, char *name)
{
    int rval = 0;
    int matbeg[1];
    double iobj[1], ilb[1], iub[1];
    char *iname[1];
    int ncolind = p->ncols;

    if (sense < 0) {
        printf ("bad sense, but cplex does not use this anyway\n");
        rval = 1; goto CLEANUP;
    }

    iobj[0] = obj;
    ilb[0] = lb;
    iub[0] = ub;
    iname[0] = name;
    matbeg[0] = 0;

    rval = CPXaddcols (p->cplex_env, p->cplex_lp, 1, nzcount, iobj, matbeg,
                       cind, cval, ilb, iub, (char **) NULL);
    p->ncols++;
    if (sense == LINGOlp_BINARY) {
       char cpx_binary = CPX_BINARY;
       rval = CPXchgctype (p->cplex_env, p->cplex_lp, 1, &ncolind, &cpx_binary);
       LINGOcheck_rval (rval, "CPXchgctype failed");
       p->nintegers++;
    }
    if (sense == LINGOlp_INTEGER) {
       char cpx_integer = CPX_INTEGER;
       rval = CPXchgctype (p->cplex_env, p->cplex_lp, 1, &ncolind, &cpx_integer);
       LINGOcheck_rval (rval, "CPXchgctype failed");
       p->nintegers++;
    }

    LINGOcheck_rval (rval, "CPXaddcols failed");

CLEANUP:
    return rval;
}

int LINGOlp_deletecols (LINGOlp *p, int first_cind, int last_cind)
{
    int rval = 0;

    rval = CPXdelcols (p->cplex_env, p->cplex_lp, first_cind, last_cind);
    LINGOcheck_rval (rval, "CPXdelcols failed");

CLEANUP:
   return rval;
}


int LINGOlp_pi (LINGOlp *p, double *pi)
{
    int rval = 0;
    int nrows;

    nrows = CPXgetnumrows (p->cplex_env, p->cplex_lp);
    rval = CPXgetpi (p->cplex_env, p->cplex_lp, pi, 0, nrows - 1);
    LINGOcheck_rval (rval, "CPXgetpi failed");

CLEANUP:
    return rval;
}

int LINGOlp_x (LINGOlp *p, double *x)
{
    int rval = 0;
    int ncols;

    ncols = CPXgetnumcols (p->cplex_env, p->cplex_lp);
    rval = CPXgetx (p->cplex_env, p->cplex_lp, x, 0, ncols - 1);
    LINGOcheck_rval (rval, "CPXgetx failed");

CLEANUP:
    return rval;
}

int LINGOlp_basis_cols (LINGOlp *p, int *cstat)
{
   int rval = 0;
   int* rstat = (int*) NULL;

   rval =  CPXgetbase(p->cplex_env, p->cplex_lp, cstat, rstat);
   LINGOcheck_rval (rval, "CPXgetbase failed");

 CLEANUP:
   return rval;
}

int LINGOlp_set_all_coltypes (LINGOlp *p, char sense)
{
   int rval = 0;

   if (!p) {
       printf ("LINGOlp_set_all_coltypes called without an LP\n");
       rval = 1;  goto CLEANUP;
   }

   if (sense != LINGOlp_CONTINUOUS) {
       printf ("Not set up to parse integer variables\n");
       rval = 1;  goto CLEANUP;
   }

CLEANUP:
   return rval;
}

int LINGOlp_objective_sense (LINGOlp *p, int sense)
{
    int rval = 0;
    char isense;

    if (sense == LINGOlp_MIN) isense = CPX_MIN;
    else                      isense = CPX_MAX;

    CPXchgobjsen (p->cplex_env, p->cplex_lp, isense);

    return rval;
}

int LINGOlp_setbound (LINGOlp *p, int col, char lower_or_upper, double bnd)
{
    int rval = 0;
    int cindex[1];
    double bd[1];
    char lu[1];

    cindex[0] = col;
    lu[0] = lower_or_upper;
    bd[0] = bnd;

    rval = CPXchgbds (p->cplex_env, p->cplex_lp, 1, cindex, lu, bd);
    LINGOcheck_rval (rval, "CPXchgbds failed");

CLEANUP:
    return rval;
}

int LINGOlp_setnodelimit (LINGOlp *p, int mip_node_limit)
{
    int rval = 0;

    if (!p || mip_node_limit < 1) {
        printf ("called with empty limit data, not set up in cplex\n");
        rval = 1;
    }
    return rval;
}



static int intercept_cplex_incumbent_cb(CPXCENVptr cpx_env, void *cbdata,
					int where, void *cbhandle)
{
   int rval = 0;

   /* Avoid warning on unused parameter usrdata:*/
   double dbl_cutoff = ((LINGOlp*)cbhandle)->dbl_cutoff;

   if (where ==CPX_CALLBACK_MIP) {
      double objective, objbound;
      int    feas_exists;
      rval = CPXgetcallbackinfo(cpx_env, cbdata, where,CPX_CALLBACK_INFO_MIP_FEAS,
                                (void*) &feas_exists);
      LINGOcheck_rval (rval, "CPXgetcallbackinfo CPX_CALLBACK_INFO_MIP_FEAS failed");

      if (feas_exists) {

         rval = CPXgetcallbackinfo(cpx_env, cbdata, where,CPX_CALLBACK_INFO_BEST_INTEGER,
                                   (void*) &objective);
         LINGOcheck_rval (rval, "CPXgetcallbackinfo CPX_CALLBACK_INFO_BEST_INTEGER failed");

         rval = CPXgetcallbackinfo(cpx_env, cbdata, where,CPX_CALLBACK_INFO_BEST_REMAINING,
                                   (void*) &objbound);
         LINGOcheck_rval (rval, "CPXgetcallbackinfo CPX_CALLBACK_INFO_BEST_REMAINING failed");

         if (objective < objbound && objective > dbl_cutoff + LINGOlp_int_tolerance()) {
	   printf("Terminating gurobicplex based on current objective value %f\n.",
		  objective);
	   
	   rval = 99;
	   goto CLEANUP;
         }
      }
   }
 CLEANUP:
      return rval;
}


int LINGOlp_set_cutoff (LINGOlp *p, double cutoff)
{
   int rval = 0;
   int objsens = -1;

   LINGOcheck_NULL(p, "LINGOlp_set_cutoff called with NULL lp.");


   rval = LINGOlp_objective_sense (p, objsens);
   LINGOcheck_rval (rval, "LINGOlp_objective_sense CPX_PARAM_EPPER failed");

   if (objsens == LINGOlp_MAX) {
      rval = CPXsetdblparam (p->cplex_env, CPX_PARAM_CUTLO, cutoff);
      LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_CUTLO failed");
   }

   if (objsens == LINGOlp_MIN) {
      rval = CPXsetdblparam (p->cplex_env, CPX_PARAM_CUTUP, cutoff);
      LINGOcheck_rval (rval, "CPXsetintparam CPX_PARAM_CUTLO failed");
   }

   if (cutoff > 0) {
      p->dbl_cutoff = cutoff;
      rval = CPXsetmipcallbackfunc(p->cplex_env,intercept_cplex_incumbent_cb,
                                   (void*) p);
      LINGOcheck_rval (rval, "CPXsetmipcallbackfunc  failed");
   }

CLEANUP:
   return rval;
}


int LINGOlp_write (LINGOlp *p, const char *fname)
{
    int rval = 0;

    rval = CPXwriteprob (p->cplex_env, p->cplex_lp, fname, "RLP");
    LINGOcheck_rval (rval, "CPXwriteprob failed");

CLEANUP:

    return rval;
}

void LINGOlp_printerrorcode (int c)
{
    printf ("CPLEX error code: %d\n", c);
    fflush (stdout);
}

double LINGOlp_int_tolerance ()
{
   return int_tolerance;
}


void *LINGOutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void LINGOutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }
    free (p);
}
