#ifndef BB_BICRITERIA_H

/* Input problems should be available at testP1.lp and testP2.lp files */
/* Only difference between two testP.lp files are the objective function coefficients */
/* Note that these problems should be maximization type and equality constraints */
/* 
   /* Max c1x
   /* subject to Ax = b */
/*
  /* All feasible points found are written to output.txt */
/* Before running this code, FindNadir2.c should be run. Nadir points are then written to init_nadir.txt */ 

/* This version does not update Pareto set or nadir set. Initial points are used */


/*The diffence of this version from version ComputeFeasSol4.c is that nadir points are computed by another code
  This program only reads them */

/* The difference of this version from BB-bicriteriaMILP.c  is that 
   y_1^SE = max{c1x: x \in \tilde X}
   y_2^SE = min{c2x: x \in \tilde X}
   y_1^NW = min{c1x: x \in \tilde X}
   y_2^NW = max{c2x: x \in \tilde X}    */ 

/* Preprocessing is not avaialble */
/* Fathoming Rule 2b is closed. It does not work. */
/* The difference of this code from BB-bicriteriaMILP3.c is that I am leaving the CPXptr dlp free now. */



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>
#include "cplex.h"

#ifdef SOL_tree
#include "max_tree.h"
#else
#ifdef SOL_list
#include "max_list.h"
#else
#error "Please define either SOL_tree or SOL_list"
#endif
#endif

#define REALLOC_STEP 1024

extern node *tree;
extern node *leftmost_node;
extern node *rightmost_node;

extern int cnt;
extern double x_ideal;
extern double y_ideal;

extern int obj1_index;
extern int obj2_index;
extern int total_num_integer;
extern int *indices;
extern int *integer_indices;
extern double *weighted_coefs;
extern double *obj_coef1;
extern double *obj_coef2;

extern int insert_counter;
extern int insert_counter2;
extern int rebuild_count;
extern int another_counter;
extern int last_insert_num;

extern double NW_extreme_x;
extern double NW_extreme_y;
extern double SE_extreme_x;
extern double SE_extreme_y;
extern double leftmost_val;
extern double leftmost_val_y;
extern double rightmost_val;
extern double rightmost_val_y;

extern int recursion_count;

extern int insert_level;

extern CPXLPptr global_mip;
extern int *global_beg;
extern int *global_varindices; 
extern double *global_values; 
extern int *global_effortlevel;
extern int global_num_starts;
extern int global_startspace;

extern CPXLPptr   lp1;
extern CPXLPptr   lp2;

extern struct nadir *theta;

extern int insert_at_beginning;


static int CPXPUBLIC 
nodeoperations (CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int *useraction_p);

static int CPXPUBLIC
usersetbranch  (CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int brtype, int sos, int nodes,
		int bdcnt, const double *nodeest, const int *nodebeg,
		const int *indices, const char *lu, const int *bd,
		int *useraction_p);

static int CPXPUBLIC 
userselectnode (CPXCENVptr env, void *cbdata, int wherefrom,
		void *cbhandle, int *nodeid_p,
		int *useraction_p);

//static int CPXPUBLIC 
//userincumbent (CPXCENVptr env, void *cbdata, int wherefrom,
//	       void *cbhandle, double objval,
//	       double *x, int *isfeas_p,
//	       int *useraction_p);


int computeextremes (CPXENVptr  env,
		     CPXLPptr   lp1,
		     CPXLPptr   lp2,
		     double *obj_coef1,
		     double *obj_coef2);

int computefeassol (CPXENVptr  env,
		    CPXLPptr   lp1,
		    CPXLPptr   lp2,
		    double *obj_coef1,
		    double *obj_coef2,
		    int nPoints,
		    int *nNadir);

int parametricsimplex(CPXLPptr lp);
double mini(double, double);

void free_and_null (char **ptr);


struct pointSeg
{double end1_z1; /* end1_z1 < end2_z1   */
  double end1_z2;
  double end2_z1;
  double end2_z2;
  enum {ISOLATED = 1, SEGMENT} type;       /* 1-isolated, 2-segment  */
};

int lp2stdform (CPXENVptr env, CPXLPptr lp);

#endif
