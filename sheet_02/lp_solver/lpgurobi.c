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
#include "lp.h"
#include <gurobi_c.h>

struct LINGOlp {
    GRBenv *env;
    GRBmodel *model;

    double dbl_cutoff;
};

const double int_tolerance = 0.00001;

#define LINGOcheck_rval_grb(rval,msg,env) {                            \
      if ((rval)) {                                                    \
         fprintf (stderr, "%s at %s, line %d: %s\n",                   \
                  (msg), __FILE__, __LINE__,GRBgeterrormsg (env));     \
         goto CLEANUP;                                                 \
    }                                                                  \
}


int LINGOlp_init (LINGOlp **p, const char *name)
{
    int rval = 0;

    (*p) = (LINGOlp *) LINGO_SAFE_MALLOC (1,LINGOlp);
    LINGOcheck_NULL (*p, "out of memory for lp");

    (*p)->env = (GRBenv *) NULL;
    (*p)->model = (GRBmodel *) NULL;

    rval = GRBloadenv (&((*p)->env), NULL);
    LINGOcheck_rval (rval, "GRBloadenv failed");

    /* Set to 1 to turn on Gurobi output, 0 to turn off output */
    rval = GRBsetintparam ((*p)->env, GRB_INT_PAR_OUTPUTFLAG, 1);
    LINGOcheck_rval_grb (rval, "GRBsetintparam OUTPUTFLAG failed",(*p)->env);

    /* Use primal simplex. */
    rval = GRBsetintparam ((*p)->env, GRB_INT_PAR_LPMETHOD, GRB_METHOD_DUAL);
    LINGOcheck_rval_grb (rval, "GRBsetintparam LPMETHOD failed",(*p)->env);


    rval = GRBsetintparam ((*p)->env, GRB_INT_PAR_THREADS , 1);
    LINGOcheck_rval_grb (rval, "GRBsetintparam THREADS failed",(*p)->env);

    rval = GRBnewmodel ((*p)->env, &((*p)->model), name, 0, (double *) NULL,
                    (double *) NULL, (double *) NULL, (char *) NULL, NULL);
    LINGOcheck_rval_grb (rval, "GRBnewmodel failed",(*p)->env);

CLEANUP:
    return rval;
}

void LINGOlp_free (LINGOlp **p)
{
    if (*p) {
        if ((*p)->model) GRBfreemodel ((*p)->model);
        if ((*p)->env) GRBfreeenv ((*p)->env);
        free (*p);
        *p = (LINGOlp *) NULL;
    }
}

int LINGOlp_optimize (LINGOlp *p)
{
    int rval = 0;
    int status;

    rval = GRBoptimize (p->model);
    LINGOcheck_rval_grb (rval, "GRBoptimize failed",p->env);

    rval = GRBgetintattr(p->model,GRB_INT_ATTR_STATUS,&status);
    LINGOcheck_rval_grb (rval, "GRBgetintattr failed",p->env);

    if (status != GRB_OPTIMAL) {
       printf("Failed to solve model to optimality. status = ");
       switch (status) {
       case GRB_LOADED:
          printf("GRB_LOADED ");
          rval = 1;break;
       case GRB_INFEASIBLE:
          printf("GRB_INFEASIBLE ");
          rval = GRBcomputeIIS(p->model);
          LINGOcheck_rval_grb (rval, "GRBcomputeIIS failed",p->env);
          rval = GRBwrite(p->model,"grbinfeas_debug.lp");
          LINGOcheck_rval_grb (rval, "GRBwrite lp failed",p->env);

          rval = 1;break;
       case GRB_INF_OR_UNBD:
          printf("GRB_INF_OR_UNBD ");
          rval = 1;break;
       case GRB_UNBOUNDED:
          printf("GRB_UNBOUNDED ");
          rval = 1;break;
       default:
          printf("%d",status);
       }
       printf("\n");

       if(rval) { goto CLEANUP; }
    }

CLEANUP:
    return rval;
}

int LINGOlp_objval (LINGOlp *p, double *obj)
{
    int rval = 0;

    rval = GRBgetdblattr (p->model, GRB_DBL_ATTR_OBJVAL, obj);
    LINGOcheck_rval_grb (rval, "GRBgetdblattr OBJVAL failed",p->env);

CLEANUP:
    return rval;

}

int LINGOlp_change_objective(LINGOlp *p, int start, int len, double* values)
{
   int rval = 0;
   rval = GRBsetdblattrarray(p->model,GRB_DBL_ATTR_OBJ,
                             start,len,values);
   LINGOcheck_rval_grb(rval,"Failed in GRBsetdblattrarray",p->env);

   rval = GRBupdatemodel (p->model);
   LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

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
        isense = GRB_EQUAL; break;
    case LINGOlp_LESS_EQUAL:
        isense = GRB_LESS_EQUAL; break;
    case LINGOlp_GREATER_EQUAL:
        isense = GRB_GREATER_EQUAL; break;
    default:
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    rval = GRBaddconstr (p->model, nzcount, cind, cval, isense, rhs, name);
    LINGOcheck_rval_grb (rval, "GRBaddconstr failed", p->env);
    rval = GRBupdatemodel (p->model);
    LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

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

   rval = GRBdelconstrs(p->model,numdel,dellist);
   LINGOcheck_rval_grb (rval, "GRBdelconstrs failed", p->env);
   rval = GRBupdatemodel (p->model);
   LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

 CLEANUP:
   if(dellist) {free (dellist);}

   return rval;
}



int LINGOlp_addcol (LINGOlp *p, int nzcount, int *cind, double *cval,
       double obj, double lb, double ub, char sense, char *name)
{
    int rval = 0;
    char isense;

    switch (sense) {
    case LINGOlp_CONTINUOUS:
        isense = GRB_CONTINUOUS; break;
    case LINGOlp_BINARY:
        isense = GRB_BINARY; break;
    case LINGOlp_INTEGER:
        isense = GRB_INTEGER; break;
    default:
        fprintf (stderr, "unknown variable sense: %c\n", sense);
        rval = 1;  goto CLEANUP;
    }

    rval = GRBaddvar (p->model, nzcount, cind, cval, obj, lb, ub, isense,
                      name);
    LINGOcheck_rval_grb (rval, "GRBaddvar failed", p->env);
    rval = GRBupdatemodel (p->model);
    LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

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

   rval = GRBdelvars(p->model,numdel,dellist);
   LINGOcheck_rval_grb (rval, "GRBdelvars failed", p->env);
   rval = GRBupdatemodel (p->model);
   LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

 CLEANUP:
   if(dellist) {free (dellist);}

   return rval;
}


int LINGOlp_pi (LINGOlp *p, double *pi)
{
    int rval = 0;
    int nrows;
    int solstat;

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_STATUS, &solstat);
    LINGOcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_STATUS failed",
                         p->env);
    if (solstat == GRB_INFEASIBLE) {
        fprintf (stderr, "Problem is infeasible\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_NUMCONSTRS, &nrows);
    LINGOcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_NUMCONSTRS failed",
                         p->env);

    if (nrows == 0) {
        fprintf (stderr, "No rows in LP\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetdblattrarray(p->model, GRB_DBL_ATTR_PI, 0, nrows, pi);
    LINGOcheck_rval_grb (rval, "GRBgetdblattrarray GRB_DBL_ATTR_PI failed",
                         p->env);

CLEANUP:
    return rval;
}

int LINGOlp_x (LINGOlp *p, double *x)
{
    int rval = 0;
    int ncols;
    int solstat;

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_STATUS, &solstat);
    LINGOcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_STATUS failed",
                         p->env);
    if (solstat == GRB_INFEASIBLE) {
        fprintf (stderr, "Problem is infeasible\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetintattr(p->model, GRB_INT_ATTR_NUMVARS, &ncols);
    LINGOcheck_rval(rval,"Failed in GRBgetintattr");

    if (ncols == 0) {
        fprintf (stderr, "No columns in LP\n");
        rval = 1;  goto CLEANUP;
    }

    rval = GRBgetdblattrarray(p->model, GRB_DBL_ATTR_X, 0, ncols, x);
    LINGOcheck_rval_grb (rval, "GRBgetdblattrarray GRB_DBL_ATTR_X failed",
                         p->env);

CLEANUP:
    return rval;
}

int LINGOlp_basis_cols (LINGOlp *p, int *cstat)
{
   int rval = 0;
   int ncols,i;

   rval = GRBgetintattr(p->model, GRB_INT_ATTR_NUMVARS, &ncols);
   LINGOcheck_rval(rval,"Failed in GRBgetintattr");

   rval = GRBgetintattrarray(p->model, GRB_INT_ATTR_VBASIS,0,ncols,cstat);
   LINGOcheck_rval(rval,"Failed in GRBgetintattrarray");

   for (i = 0; i < ncols; ++i) {
      switch (cstat[i]) {
      case GRB_BASIC:
         cstat[i] = LINGOlp_BASIC;
         break;
      case GRB_NONBASIC_LOWER:
         cstat[i] = LINGOlp_LOWER;
         break;
      case GRB_NONBASIC_UPPER:
         cstat[i] = LINGOlp_UPPER;
         break;
      case GRB_SUPERBASIC:
         cstat[i] = LINGOlp_FREE;
         break;
      default:
         rval = 1;
         LINGOcheck_rval(rval,"ERROR: Received unknown cstat");
      }
   }
 CLEANUP:
   return rval;
}

int LINGOlp_set_all_coltypes (LINGOlp *p, char sense)
{
   int nvars,i, rval = 0;
   char isense;

   switch (sense) {
   case LINGOlp_CONTINUOUS:
       isense = GRB_CONTINUOUS; break;
   case LINGOlp_BINARY:
       isense = GRB_BINARY; break;
   case LINGOlp_INTEGER:
       isense = GRB_INTEGER; break;
   default:
       fprintf (stderr, "unknown variable sense: %c\n", sense);
       rval = 1;  goto CLEANUP;
   }

   rval= GRBgetintattr(p->model,GRB_INT_ATTR_NUMVARS,&nvars);
   LINGOcheck_rval_grb (rval, "GRBgetintattr GRB_INT_ATTR_NUMVARS failed",
                        p->env);

   for (i = 0; i < nvars; i++) {
      rval = GRBsetcharattrelement(p->model,GRB_CHAR_ATTR_VTYPE,i,isense);
      LINGOcheck_rval_grb (rval, "GRBsetintattrelement GRB_CHAR_ATTR_VTYPE failed",
                           p->env);
   }

   rval = GRBupdatemodel (p->model);
   LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

 CLEANUP:
   return rval;
}

int LINGOlp_objective_sense (LINGOlp *p, int sense)
{
    int rval = 0;

    /* Min = 1   Max = -1 */
    rval = GRBsetintattr (p->model, GRB_INT_ATTR_MODELSENSE, sense);
    LINGOcheck_rval_grb (rval, "GRBsetintattr failed", p->env);

CLEANUP:

    return rval;
}

int LINGOlp_setbound (LINGOlp *p, int col, char lower_or_upper, double bnd)
{
    int rval = 0;

    if (lower_or_upper == 'L') {
        rval = GRBsetdblattrelement (p->model, GRB_DBL_ATTR_LB, col, bnd);
    } else {
        rval = GRBsetdblattrelement (p->model, GRB_DBL_ATTR_UB, col, bnd);
    }
    LINGOcheck_rval (rval, "GRBsetdblattr LB or UB failed");

    rval = GRBupdatemodel (p->model);
    LINGOcheck_rval_grb (rval, "GRBupdatemodel failed", p->env);

CLEANUP:

    return rval;
}

int LINGOlp_setnodelimit (LINGOlp *p, int mip_node_limit)
{
   int rval = GRBsetdblparam (GRBgetenv(p->model), GRB_DBL_PAR_NODELIMIT, mip_node_limit);
   LINGOcheck_rval_grb (rval, "GRBsetdblparam NODELIMIT failed",p->env);
 CLEANUP:
   return rval;
}



static int intercept_grb_cb(GRBmodel *grb_model, void *cbdata, int where, void *usrdata)
{
   int rval = 0;

   /* Avoid warning on unused parameter usrdata:*/

   double dbl_cutoff = ((LINGOlp*)usrdata)->dbl_cutoff;

   if (where ==GRB_CB_MIPSOL) {
      double objective, objbound;

      rval = GRBcbget(cbdata,where,GRB_CB_MIPSOL_OBJBST,(void*) &objective);
      LINGOcheck_rval (rval, "GRBcbget OBJBST failed");

      rval = GRBcbget(cbdata,where,GRB_CB_MIPSOL_OBJBND,(void*) &objbound);
      LINGOcheck_rval (rval, "GRBcbget OBJBND failed");


      if (objective < objbound && objective > dbl_cutoff + LINGOlp_int_tolerance()) {
	 printf("Terminating gurobi based on current objective value %f\n.",
                   objective);
	 GRBterminate(grb_model);
      }
   }

 CLEANUP:
   return rval;
}



int LINGOlp_set_cutoff (LINGOlp *p, double cutoff)
{
   int rval = 0;

   rval = GRBsetdblparam (GRBgetenv(p->model),GRB_DBL_PAR_CUTOFF, cutoff);
   LINGOcheck_rval(rval,"Failed in GRBsetdblparam GRB_DBL_PAR_CUTOFF");

   if (cutoff > 0) {
      p->dbl_cutoff = cutoff;

      rval  = GRBsetcallbackfunc(p->model, intercept_grb_cb, (void*) p);
      LINGOcheck_rval (rval, "GRBsetcallbackfunc failed");
   }

CLEANUP:
   return rval;
}

int LINGOlp_write (LINGOlp *p, const char *fname)
{
    int rval = 0;

    rval = GRBwrite (p->model, fname);
    LINGOcheck_rval_grb (rval, "GRBwrite failed", p->env);

CLEANUP:

    return rval;
}

void LINGOlp_printerrorcode (int c)
{
    switch (c) {
    case GRB_ERROR_OUT_OF_MEMORY:
        printf ("Available memory was exhausted\n");
        break;
    case GRB_ERROR_NULL_ARGUMENT:
        printf ("NULL input value provided for a required argument\n");
        break;
    case GRB_ERROR_INVALID_ARGUMENT:
        printf ("An invalid value was provided for a routine argument\n");
        break;
    case GRB_ERROR_UNKNOWN_ATTRIBUTE:
        printf ("Tried to query or set an unknown attribute\n");
        break;
    case GRB_ERROR_DATA_NOT_AVAILABLE:
        printf ("Attempted to query or set an attribute that could\n");
        printf ("not be accessed at that time\n");
        break;
    case GRB_ERROR_INDEX_OUT_OF_RANGE:
        printf ("Tried to query or set an attribute, but one or more\n");
        printf ("of the provided indices (e.g., constraint index, variable \n");
        printf ("index) was outside the range of valid values\n");
        break;
    case GRB_ERROR_UNKNOWN_PARAMETER:
        printf ("Tried to query or set an unknown parameter\n");
        break;
    case GRB_ERROR_VALUE_OUT_OF_RANGE:
        printf ("Tried to set a parameter to a value that is outside\n");
        printf ("the parameter's valid range\n");
        break;
    case GRB_ERROR_NO_LICENSE:
        printf ("Failed to obtain a valid license\n");
        break;
    case GRB_ERROR_SIZE_LIMIT_EXCEEDED:
        printf ("Attempted to solve a model that is larger than the\n");
        printf ("limit for a demo license\n");
        break;
    case GRB_ERROR_CALLBACK:
        printf ("Problem in callback\n");
        break;
    case GRB_ERROR_FILE_READ:
        printf ("Failed to read the requested file\n");
        break;
    case GRB_ERROR_FILE_WRITE:
        printf ("Failed to write the requested file\n");
        break;
    case GRB_ERROR_NUMERIC:
        printf ("Numerical error during requested operation\n");
        break;
    case GRB_ERROR_IIS_NOT_INFEASIBLE:
        printf ("Attempted to perform infeasibility analysis on a\n");
        printf ("feasible model\n");
        break;
    case GRB_ERROR_NOT_FOR_MIP:
        printf ("Requested operation not valid for a MIP model\n");
        break;
    case GRB_ERROR_OPTIMIZATION_IN_PROGRESS:
        printf ("Tried to query or modify a model while optimization\n");
        printf ("was in progress\n");
        break;
     default:
        printf ("Unknown error code: %d\n", c);
    }
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

