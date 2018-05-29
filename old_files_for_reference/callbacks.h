#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>

#include "cplex.h"

extern double PSA_x_val;
extern double PSA_y_val;
extern int changed;
               
int CPXPUBLIC
userincumbent (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p);
               
void chg_coefs(CPXCENVptr env, CPXLPptr prob, int *indices, double slope);

void PSA_full_right(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig);
void PSA_full_left(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig);
void PSA_full(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig);

