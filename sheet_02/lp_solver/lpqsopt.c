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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <qsopt.h>

#include "lp.h"

struct LINGOlp {
    QSprob p;
};

const double int_tolerance = 0.00001;

int LINGOlp_init (LINGOlp **p, const char *name)
{
    int rval = 0;

    (*p) = LINGO_SAFE_MALLOC (1, LINGOlp);
    if ((*p) == (LINGOlp *) NULL) {
        fprintf (stderr, "Out of memory in LINGOlp_init\n");
        rval = 1; goto CLEANUP;
    }

    (*p)->p = QScreate_prob (name, QS_MIN);
    if ((*p)->p == (QSprob) NULL) {
        fprintf (stderr, "QScreate_prob failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param ((*p)->p, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

void LINGOlp_free (LINGOlp **p)
{
    if (*p) {
        QSfree_prob ((*p)->p);
        LINGO_FREE (*p, LINGOlp);
    }
}

int LINGOlp_optimize (LINGOlp *p)
{
    int rval = 0;
    int status;

    rval = QSopt_dual (p->p, &status);
    if (rval) {
        fprintf (stderr, "QSopt_dual failed\n"); goto CLEANUP;
    }

    if (status == QS_LP_ITER_LIMIT) {
        printf ("Dual LP Solver reached iteration limit\n"); fflush (stdout);
        rval = QSwrite_prob (p->p, "iter.lp", "LP");
        if (rval) {
            fprintf (stderr, "QSwrite_prob failed\n");
            goto CLEANUP;
        }
        printf ("Saved LP as iter.lp\n"); fflush (stdout);
        rval = QSwrite_basis (p->p, (QSbas) NULL, "iter.bas");
        if (rval) {
            fprintf (stderr, "QSwrite_basis failed\n"); goto CLEANUP;
        }
        printf ("Saved LP as iter.bas\n"); fflush (stdout);
    } else if (status == QS_LP_TIME_LIMIT) {
        printf ("Dual LP Solver reached time limit\n"); fflush (stdout);
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = 2; goto CLEANUP;
    } else if (status != QS_LP_OPTIMAL) {
        fprintf (stderr, "no optimal LP-solution exists: %d\n", status);
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int LINGOlp_objval (LINGOlp *p, double *obj)
{
    int rval = 0;

    rval = QSget_objval (p->p, obj);
    if (rval) {
        fprintf (stderr, "QSget_objval failed");
        goto CLEANUP;
    }


CLEANUP:
    return rval;

}

int LINGOlp_change_objective(LINGOlp *p, int start, int len, double* values)
{
    int i, rval = 0;

    for (i = 0; i < len; i++) {
        rval = QSchange_objcoef (p->p, start+i, values[i]);
        LINGOcheck_rval (rval, "Failed in QSchange_objcoef");
    }

CLEANUP:
   return rval;
}

int LINGOlp_addrow (LINGOlp *p, int nzcount, int *cind, double *cval,
       char sense, double rhs, char *name)
{
    int rval = 0;
    char isense;

    switch (sense) {
    case LINGOlp_EQUAL:
        isense = 'E'; break;
    case LINGOlp_LESS_EQUAL:
        isense = 'L'; break;
    case LINGOlp_GREATER_EQUAL:
        isense = 'G'; break;
    default:
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    rval = QSadd_row (p->p, nzcount, cind, cval, rhs, isense, name);
    LINGOcheck_rval (rval, "QSadd_row failed");

CLEANUP:
    return rval;
}

int LINGOlp_deleterows (LINGOlp *p, int first_cind, int last_cind)
{
   int rval = 0;
   int* dellist = (int*) NULL;
   int numdel  = last_cind - first_cind + 1;
   int i;

   dellist = LINGO_SAFE_MALLOC(numdel,int);
   LINGOcheck_NULL(dellist, "Failed to allocate dellist");

   for (i = 0; i < numdel; ++i) {
      dellist[i] = first_cind + i;
   }

   rval = QSdelete_rows (p->p, numdel,dellist);
   LINGOcheck_rval (rval, "QSdelete_rows failed");

CLEANUP:
   if (dellist) free(dellist);
   return rval;

}
int LINGOlp_addcol (LINGOlp *p, int nzcount, int *cind, double *cval,
       double obj, double lb, double ub, char sense, char *name)
{
    int rval = 0;

    if (sense < 0) {
        printf ("bad sense, but qsopt does not use this anyway\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSadd_col (p->p, nzcount, cind, cval, obj, lb, ub, name);
    LINGOcheck_rval (rval, "QSadd_col failed");

CLEANUP:
    return rval;
}

int LINGOlp_deletecols (LINGOlp *p, int first_cind, int last_cind)
{
   int rval = 0;
   int* dellist = (int*) NULL;
   int numdel  = last_cind - first_cind + 1;
   int i;

   dellist = LINGO_SAFE_MALLOC(numdel,int);
   LINGOcheck_NULL(dellist, "Failed to allocate dellist");

   for (i = 0; i < numdel; ++i) {
      dellist[i] = first_cind + i;
   }

   rval = QSdelete_cols (p->p, numdel,dellist);
   LINGOcheck_rval (rval, "QSdelete_col failed");

CLEANUP:
   if (dellist) free(dellist);
   return rval;
}


int LINGOlp_pi (LINGOlp *p, double *pi)
{
    int rval = 0;
    int status;

    rval = QSget_status (p->p, &status);
    if (rval) {
        fprintf (stderr, "QSget_status failed"); goto CLEANUP;
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = QSget_infeas_array (p->p, pi);
        LINGOcheck_rval (rval, "QSget_infeas_array failed");
    } else {
        rval = QSget_pi_array (p->p, pi);
        LINGOcheck_rval (rval, "QSget_pi_array failed");
    }

CLEANUP:
    return rval;
}

int LINGOlp_x (LINGOlp *p, double *x)
{
    int rval = 0;

    rval = QSget_x_array (p->p, x);
    LINGOcheck_rval (rval, "QSget_x_array failed");

 CLEANUP:
    return rval;
}

int LINGOlp_basis_cols (LINGOlp *p, int *int_cstat)
{
   int rval = 0;
   char* rstat = (char*) NULL;
   char* cstat = (char*) NULL;
   int   ncols = QSget_colcount (p->p);
   int   nrows = QSget_rowcount (p->p);
   int   i;
   cstat = LINGO_SAFE_MALLOC(ncols, char);
   LINGOcheck_NULL(cstat,"Failed to allocate cstat");

   rstat = LINGO_SAFE_MALLOC(nrows, char);
   LINGOcheck_NULL(rstat,"Failed to allocate rstat");


   rval = QSget_basis_array (p->p, cstat,rstat);
   LINGOcheck_rval (rval, "QSget_basis_array failed");

   for (i = 0; i < ncols; ++i) {
      switch (cstat[i]) {
      case QS_COL_BSTAT_LOWER:
         int_cstat[i] = LINGOlp_LOWER;
         break;
      case QS_COL_BSTAT_UPPER:
         int_cstat[i] = LINGOlp_UPPER;
         break;
      case QS_COL_BSTAT_FREE:
         int_cstat[i] = LINGOlp_FREE;
         break;
      case QS_COL_BSTAT_BASIC:
         int_cstat[i] = LINGOlp_BASIC;
         break;
      default:
         rval = 1;
         LINGOcheck_rval(rval,"ERROR: Received unknown cstat");
      }
   }
 CLEANUP:
   if(cstat) free(cstat);
   if(rstat) free(rstat);

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
       printf ("QSopt does not handled integer variables\n");
       rval = 1;  goto CLEANUP;
   }

CLEANUP:
   return rval;
}

int LINGOlp_objective_sense (LINGOlp *p, int sense)
{
    int rval = 0;
     char isense;

    if (sense == LINGOlp_MIN) isense = QS_MIN;
    else                      isense = QS_MAX;

    rval = QSchange_objsense (p->p, isense);
    LINGOcheck_rval (rval, "QSchange_objsense failed");

CLEANUP:
    return rval;
}

int LINGOlp_setbound (LINGOlp *p, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    rval = QSchange_bound (p->p, col, lower_or_upper, bnd);
    LINGOcheck_rval (rval, "QSchange_bounds failed");

CLEANUP:
    return rval;
}

int LINGOlp_setnodelimit (LINGOlp *p, int mip_node_limit)
{
    int rval = 0;

    if (!p || mip_node_limit < 1) {
        printf ("called with empty limit data, no meaning in QSopt\n");
        rval = 1;
    }
    return rval;
}

int LINGOlp_write (LINGOlp *p, const char *fname)
{
    int rval = 0;

    rval = QSwrite_prob (p->p, fname, "LP");
    LINGOcheck_rval (rval, "QSwrite_prob failed");

CLEANUP:

    return rval;
}

void LINGOlp_printerrorcode (int c)
{
    printf ("QSopt error code: %d\n", c);
    fflush (stdout);
}

double LINGOlp_int_tolerance ()
{
   return int_tolerance;
}

int LINGOlp_set_cutoff (LINGOlp *p, double cutoff)
{
   int rval = 1;

   (void) p;
   (void) cutoff;

   LINGOcheck_rval(rval,"LINGOlp_set_cutoff not yet implemented.");

CLEANUP:
   return rval;
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
