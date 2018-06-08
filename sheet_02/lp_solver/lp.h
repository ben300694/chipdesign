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

#ifndef __LP_H
#define __LP_H


/** lp.h is an interface to various MIP/LP-solvers. There are implementations
    for CPLEX, Gurobi, and QSopt.
*/

typedef struct LINGOlp LINGOlp;

/** Create a MIP/LP-problem with the given name. */
int  LINGOlp_init (LINGOlp **p, const char *name);

/** Free a MIP/LP-problem with the given name. */
void LINGOlp_free (LINGOlp **p);

/** Variables/columns can be continuous (LPs), binary or integer (MIPs).*/
#define LINGOlp_CONTINUOUS 0
#define LINGOlp_BINARY     1
#define LINGOlp_INTEGER    2

/** Constraints/Rows refer to one of the following inequality-types*/
#define LINGOlp_EQUAL         'E'
#define LINGOlp_LESS_EQUAL    'L'
#define LINGOlp_GREATER_EQUAL 'G'

/** A column/variable in the current solution can be 
    tight at the lower/upper bounds, basic, or free (non-basic).
 */
#define LINGOlp_LOWER      0
#define LINGOlp_BASIC      1
#define LINGOlp_UPPER      2
#define LINGOlp_FREE       3

/** The problem can either be a maximization or minimization problem.*/
#define LINGOlp_MIN  1
#define LINGOlp_MAX -1

/** Optimize the given problem. */
int LINGOlp_optimize (LINGOlp *p);
/** Query the optimum value of problem (requires LINGOlp_optimize to be run beforehand).*/
int LINGOlp_objval (LINGOlp *p, double *obj);
/** Query the dual variables, where pi must be an array of the length of the number of rows.*/
int LINGOlp_pi (LINGOlp *p, double *pi);
/** Query the primal variables, where x must be an array of the length of the number of columns.*/
int LINGOlp_x (LINGOlp *p, double *x);

int LINGOlp_basis_cols (LINGOlp *p, int *cstat);

int LINGOlp_change_objective(LINGOlp *p, int start, int len, double* values);

/** Add a row/constraint to the LP (sparsely). Parameters:
    p       : The LP
    nzcount : number of (non-zero) column-entries in the new row.
    cind    : array of length nzcount, specifying the (non-zero) column/variable
              indices of the new row.
    cval    : array of length nzcount, specifying the coefficients of the
              columns/variables specified in cind.
    sense   : type of inequalitym either LINGOlp_EQUAL, LINGOlp_LESS_EQUAL, or 
              LINGOlp_GREATER_EQUAL            
    rhs     : the right hand side (b-component) of the constraint.
    name    : name of the constraint, which can be set to (char*) NULL.   
 */
int LINGOlp_addrow (LINGOlp *p, int nzcount, int *cind, double *cval,
                    char sense, double rhs, char *name);

/** Delete all rows/constraints with index between first_cind and last_cind.*/
int LINGOlp_deleterows (LINGOlp *p, int first_cind, int last_cind);

/**Add a column/variable to the LP.
   p       : The LP
   nzcount : number of non-zero row-entries of the new column.
   rind    : array of length nzcount, specifying the (non-zero) row/constraint
             indices of the new column.
   rvals   : array of length nzcount, specifying the coefficients of the
              rows/contraints specified in rind.
   obj     : coefficient in the objective (c-vector) of the new variable.
   lb      : lower bound for the new variable.
   ub      : upper bound for the new variable.
   vartype : Variable type: one of LINGOlp_CONTINUOUS, LINGOlp_INTEGER and LINGOlp_BINARY.
             The latter two are not valid in this context.
   name    : name of the constraint, which can be set to (char*) NULL.   
*/
int LINGOlp_addcol (LINGOlp *p, int nzcount, int *rind, double *rvals,
                    double obj, double lb, double ub, char vartype, char *name);

/** Delete all columns/variables with index between first_cind and last_cind.*/
int LINGOlp_deletecols (LINGOlp *p, int first_cind, int last_cind);

int LINGOlp_set_all_coltypes (LINGOlp *p, char sense);

int LINGOlp_objective_sense (LINGOlp *p, int sense);
int LINGOlp_setbound (LINGOlp *p, int col, char lower_or_upper, double bnd);
int LINGOlp_setnodelimit (LINGOlp *p, int mip_node_limit);
int LINGOlp_set_cutoff (LINGOlp *p, double cutoff);


int LINGOlp_write (LINGOlp *p, const char *fname);
void LINGOlp_printerrorcode (int c);

double LINGOlp_int_tolerance (void);


void *LINGOutil_allocrus (size_t size);
void  LINGOutil_freerus (void *p);


/** Check whether rval is 0. If not, jump to CLEANUP.*/
#define LINGOcheck_rval(rval,msg) {                                        \
    if ((rval)) {                                                          \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       goto CLEANUP;                                                       \
    }                                                                      \
}

/** Check whether item is NULL. If so, jump to CLEANUP.*/
#define LINGOcheck_NULL(item,msg) {			                   \
    if ((!item)) {                                                         \
       fflush(stdout);                                                     \
       fprintf (stderr, "%s at %s, line %d\n", (msg),__FILE__,__LINE__);   \
       rval = 1;                                                           \
       goto CLEANUP;                                                       \
    }                                                                      \
}

/** Allocate an array of  nnum objects of type type.*/
#define LINGO_SAFE_MALLOC(nnum,type)                                       \
    (type *) LINGOutil_allocrus (((size_t) (nnum)) * sizeof (type))

/** Free the pointer object of type type.*/
#define LINGO_FREE(object,type) {                                          \
    LINGOutil_freerus ((void *) (object));                                 \
    object = (type *) NULL;                                                \
}

/** Free the pointer object of type type if it is not NULL.*/
#define LINGO_IFFREE(object,type) {                                        \
    if ((object)) LINGO_FREE ((object),type);                              \
}


#endif  /* __LP_H */
