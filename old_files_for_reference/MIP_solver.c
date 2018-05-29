/* File created by Nathan Adelgren, Graduate Assisistant at Clemson University.
Started: 9/1/2014 
Finished: N/A
This work serves as part of my doctoral research under the advisement of Dr. Akshay Gupte. 
The goal of this program is to solve Biobjective Mixed-Integer Linear Programs using an MIP
based technique and employing the Parametric Simplex Algorithm (PSA) to find Pareto optimal
integer feasible line segments. (Small pieces of this code may be copied from a BB based 
code authored by Drs. Banu Soylu and Pietro Belotti. I collaborated on this work.)*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>

#include "cplex.h"
#include "max_tree.h"
#include "bb-bicriteria.h"
#include "MIP_solver.h"
#include "callbacks.h"

/*******************************************************************/
/*              CPLEX VARIABLES                                    */
/*******************************************************************/

CPXENVptr  env=NULL;
CPXLPptr   lp1=NULL;
CPXLPptr   lp2=NULL;

double *obj_coef1=NULL;
double *obj_coef2=NULL;
double *obj_coef3=NULL;
int cur_numcols ;
int cur_numrows ;
int status = 0;
int status2 = 0;
double objval = 0.;
double *x = NULL;
int numsols = 0;
int real_numsols = 0;
int num_mips_solved = 0;

double PSA_x_val = 0.;
double PSA_y_val = 0.;
int changed = 1, sol_added = 0;

clock_t start_time, finish_time;

/*double x_ps[10000];*/
/*double x_ps_old[10000];*/
/*char *ctype;*/
/*double objval_ps;*/
/*int *indices_new = NULL;*/
/*int *indices_still_to_check = NULL;*/
/*int *indices_f_n_l = NULL;*/
/*int a1,b1,c1,d1,p;*/
/*double x_proj, y_proj;*/
/*double y_loc_ideal_p2s, x_loc_ideal_p2s;*/
/*double *infeasout=NULL;*/
/*double infeas_amt;*/
/******************************************************************/

/*******************************************************************/
/*    Other Global Variables 			                   */
/*******************************************************************/

double x_ideal = 0.;
double y_ideal = 0.;
split_pt *first_split_pt = NULL;
split_pt *last_split_pt = NULL;

int obj1_index;
int obj2_index;
int total_num_integer;
int *integer_indices;
double *weighted_coefs;
int *indices;

double NW_extreme_x;
double NW_extreme_y;
double SE_extreme_x;
double SE_extreme_y;

double prev_cheby_sol = 0.;

/*******************************************************************/

void free_and_null (char **ptr);

void provide_xctype(CPXENVptr env, char *types, int num_cols)
{
	integer_indices = (int *) malloc ((num_cols+3)*sizeof(int));

	total_num_integer = 0;
	int i;
	for(i=0;i<num_cols;i++)
	{
		if( types[i] == 'B' || types[i] == 'I')
		{
			integer_indices[total_num_integer] = i;
			total_num_integer++;
		}
	}
}

void add_split_pt	(split_pt *before,
			split_pt *after,
			double objval1,
			double objval2,
			double *x)
{
	if(before == NULL)
	{
		first_split_pt = (struct split_pt*) malloc( sizeof( struct split_pt ) ); 
		first_split_pt->f1 = objval1;
		first_split_pt->f2 = objval2;
		first_split_pt->x = x;
		first_split_pt->prev = NULL;
		first_split_pt->next = NULL;
	}
	else if(after == NULL)
	{
		last_split_pt = (struct split_pt*) malloc( sizeof( struct split_pt ) ); 
		last_split_pt->f1 = objval1;
		last_split_pt->f2 = objval2;
		last_split_pt->x = x;
		last_split_pt->prev = before;
		before->next = last_split_pt;
		last_split_pt->next = NULL;
	}
	else
	{
		split_pt *this_split_pt = (struct split_pt*) malloc( sizeof( struct split_pt ) );
		this_split_pt->f1 = objval1;
		this_split_pt->f2 = objval2;
		this_split_pt->x = x;
		this_split_pt->prev = before;
		before->next = this_split_pt;
		this_split_pt->next = after;
		after->prev = this_split_pt;
	}
} /* End of add_split_pt */

void delete_split_pts(split_pt *pt)
{
	if(pt->next) delete_split_pts(pt->next);
	free(pt);
}

double prob_width = 0.;

int print_in_cfe = 0;

int compute_feas_extremes	(CPXENVptr  env,
		    		CPXLPptr   lp1,
		    		CPXLPptr   lp2,
		    		double *obj_coef1,
				double *obj_coef2)
{
	weighted_coefs = (double*)calloc(cur_numcols+2, sizeof(double));
	
	CPXLPptr lpclone1 = NULL;
	
	/****************** Copy the Problem *****************************/

  	lpclone1 = CPXcloneprob (env, lp1, &status);
  	if ( status ) 
  	{
    		fprintf (stderr, "Failed to clone problem 1.\n");
    		exit(0);
  	}
  	
  	/******************************************************************/
  	
  	/****************** Set parameters ****************************/

/*	status = CPXsetintparam (env, CPX_PARAM_REPEATPRESOLVE, 3);*/
/*	status = CPXsetintparam (env, CPXPARAM_Emphasis_MIP,CPX_MIPEMPHASIS_BALANCED);*/
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
/*	status = CPXsetdblparam (env, CPX_PARAM_EPGAP, 1.0e-5);*/
/*	status = CPXsetdblparam (env, CPX_PARAM_TILIM, 180.);*/
/*	if ( status ) {*/
/*		printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);*/
/*		exit(0);*/
/*	}*/
	
	status = CPXsetintparam ((CPXENVptr) env, CPX_PARAM_SCRIND, CPX_OFF);
	
	/******************************************************************/
  	
  	/****************** Compute the extremes ****************************/
	
	x = (double *) malloc (cur_numcols*sizeof(double));
	int i=0;
	int j=0;
	
/*	chg_coefs_(env, lpclone1, indices, -10000000.);*/
	
/*	status = CPXwriteprob (env, lpclone1, "myprob.lp", "LP");*/
	
  	CPXmipopt (env, lpclone1);
  	num_mips_solved++;
/*	CPXpopulate(env,lpclone1);*/
    	status2 = CPXgetstat(env, lpclone1);
    	if(print_in_cfe) printf("status2: %d\n",status2);
    	
    	if(status2 == 107 || status2 == 108)
    	{
    		status = CPXsetintparam (env, CPXPARAM_Emphasis_MIP,CPX_MIPEMPHASIS_OPTIMALITY);
    		status = CPXsetdblparam (env, CPX_PARAM_TILIM, pow(10.,75.));
		if ( status ) {
			printf ("Failed to set solution time limit to 60 seconds, error %d.\n",status);
			exit(0);
		}
		CPXmipopt (env, lpclone1);
    	}
    	
/*    	status2 = CPXgetstat(env, lpclone1);*/
    	
    	
    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
    	{
		if(status2 == 118)
    		{
    			printf("The first test problem is unbounded. Exiting!\n");
    			exit(0);
    			exit(0);
    		}
    		else if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
    	}
    	
    	global_mip = CPXcloneprob (env, lpclone1, &status);
  	if ( status ) 
  	{
    		fprintf (stderr, "Failed to clone problem 1.\n");
    		exit(0);
  	}
    	
    	int num_starts = CPXgetnummipstarts (env, lpclone1);
/*    	printf("num starts: %d\n",num_starts);*/
  	global_num_starts = num_starts;
  	global_startspace = cur_numcols*num_starts;
/*  	printf("number of mip starts: %d\n",num_starts);*/
  	
/*  	printf("first allocation\n");*/
  	global_beg = (int *) malloc ((global_num_starts)*sizeof(int));
	global_varindices = (int *) malloc ((global_startspace)*sizeof(int));
	global_values = (double *) malloc ((global_startspace)*sizeof(double));
	global_effortlevel = (int *) malloc ((global_startspace)*sizeof(int));
	
	int nzcnt = 0, surplus = 0;
	
	status = CPXgetmipstarts (env, lpclone1, &nzcnt, global_beg, global_varindices, 
                           global_values, global_effortlevel, global_startspace,
                           &surplus, 0, global_num_starts-1);
                           
        status = CPXaddmipstarts (env, global_mip, global_num_starts, nzcnt, global_beg, global_varindices,
                           global_values, global_effortlevel, NULL);

	double obvals[2] = {0.,0.};
  	
  	status = CPXgetx (env, lpclone1, obvals, obj1_index, obj2_index);
  	if(status)
  	{
		printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
  		exit(0);
  	} 
  	
  	int add_check = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);
  	if(add_check == -1) 
	{
		status = CPXgetx (env, lpclone1, x, 0, cur_numcols-1);
		insert_level++;
		leftmost_val = SE_extreme_x + 10.;
		PSA_full(env,NULL,x,NULL,NULL);
	} 	
  		
  	SE_extreme_x = obvals[0];
  	SE_extreme_y = obvals[1];
  	
  	if(print_in_cfe) printf("SE_x: %lf \t SE_y: %lf\n",SE_extreme_x,SE_extreme_y );
  	
/*  	status = CPXgetx (env, lpclone1, x, 0, cur_numcols-1);*/
/*  	for(i=0;i<cur_numcols;i++)*/
/*  	{*/
/*  		printf("x%d: %lf\n",i,x[i]);*/
/*  	}*/
/*  	printf("\n\n");*/
  	
  	 
	chg_coefs(env, lpclone1, indices, 0.);
	
/*	num_starts = CPXgetnummipstarts (env, lpclone1);*/
/*    	printf("num starts: %d\n",num_starts);*/
/*    	exit(0);*/
	
/*	status = CPXwriteprob (env, lpclone1, "myprob1.lp", "LP");*/
	
	CPXmipopt (env, lpclone1);
	num_mips_solved++;
    	status2 = CPXgetstat(env, lpclone1);
    	if(print_in_cfe) printf("status2: %d\n",status2);
    	
/*    	global_mip = CPXcloneprob (env, lpclone1, &status);*/
/*  	if ( status ) */
/*  	{*/
/*    		fprintf (stderr, "Failed to clone problem 1.\n");*/
/*    		exit(0);*/
/*  	}*/

	num_starts = CPXgetnummipstarts (env, lpclone1);
/*	printf("num starts: %d\n",num_starts);*/
/*    	exit(0);*/
  	global_num_starts = num_starts;
  	global_startspace = cur_numcols*num_starts;
/*  	printf("number of mip starts: %d\n",num_starts);*/
  	
/*  	printf("first allocation\n");*/
  	global_beg = (int *) realloc (global_beg,(global_num_starts)*sizeof(int));
	global_varindices = (int *) realloc (global_varindices,(global_startspace)*sizeof(int));
	global_values = (double *) realloc (global_values,(global_startspace)*sizeof(double));
	global_effortlevel = (int *) realloc (global_effortlevel,(global_startspace)*sizeof(int));
	
	status = CPXgetmipstarts (env, lpclone1, &nzcnt, global_beg, global_varindices, 
                           global_values, global_effortlevel, global_startspace,
                           &surplus, 0, global_num_starts-1);
                           
        status = CPXaddmipstarts (env, global_mip, global_num_starts, nzcnt, global_beg, global_varindices,
                           global_values, global_effortlevel, NULL);
                           
        num_starts = CPXgetnummipstarts (env, global_mip);
    	
    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
    	{
		if(status2 == 118)
    		{
    			printf("The second test problem is unbounded. Exiting!\n");
    			exit(0);
    			exit(0);
    		}
    		else if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
    	}

	status = CPXgetx (env, lpclone1, obvals, obj1_index, obj2_index);
  	if(status)
  	{
		printf("Failed to get x-values from CPLEX. Line: %d\n",__LINE__);
		printf("status of second mip: %d\n",status2);
  		exit(0);
  	}
  	
  	add_check = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);
  	if(add_check == -1) 
	{
		status = CPXgetx (env, lpclone1, x, 0, cur_numcols-1);
		insert_level++;
		leftmost_val = SE_extreme_x + 10.;
		PSA_full(env,NULL,x,NULL,NULL);
	}
  	
  	NW_extreme_x = obvals[0];
  	NW_extreme_y = obvals[1];
  	
  	if(print_in_cfe) printf("NW_x: %lf \t NW_y: %lf\n",NW_extreme_x,NW_extreme_y );
  	
/*  	status = CPXgetx (env, lpclone1, x, 0, cur_numcols-1);*/
/*  	for(i=0;i<cur_numcols;i++)*/
/*  	{*/
/*  		printf("x%d: %lf\n",i,x[i]);*/
/*  	}*/
/*  	printf("\n\n");*/
  	    	
      	prob_width = SE_extreme_x - NW_extreme_x;
      	if(prob_width < .0001) 
      	{
      		printf("Non-conflicting objectives. Only Pareto point is ideal point: (%lf,%lf)\n",NW_extreme_x,NW_extreme_y);
      		finish_time = clock();
	   	double duration = (double)(finish_time - start_time) / CLOCKS_PER_SEC;
	   	printf("Total time: %lf\n",duration);
      		exit(0);
      	}

  	/******************************************************************/
  	
  	/****************** Code for TERMINATE *****************************/

/*  	TERMINATE:*/

  	free_and_null ((char **) &x);
  	free_and_null ((char **) &global_beg);
	free_and_null ((char **) &global_varindices);
	free_and_null ((char **) &global_values);
	free_and_null ((char **) &global_effortlevel);
  	CPXfreeprob(env, &lpclone1);
  	
  	return (status);
  	/******************************************************************/
  	
} /* End of compute_feas_extremes */

/*******************************************************************/
/* 	In this function we perform the main body of the
	algorithm, within the box in the objective space
	designated by the points (se_x,se_y) and (nw_x,nw_y).      */
/*******************************************************************/

int recursion_count = 0, max_recursion = 50000;
CPXLPptr lp_clone = NULL;

int print_stuff_in_mss = 0;
int these_two_needs_freed = 0;

int same_box_counter = 0, same_closest_counter = 0, same_WS_MIP_solution_counter = 0, same_CH_MIP_solution_counter = 0;

double prev_se_x, prev_se_y, prev_nw_x, prev_nw_y, prev_closest_se_x, prev_closest_se_y, prev_closest_nw_x, prev_closest_nw_y, prev_WS_1, prev_WS_2, prev_CH_1, prev_CH_2;

void mip_solve_sequential (CPXENVptr env, box *b1, box *b2) //(CPXENVptr env, double se_x, double se_y, double nw_x, double nw_y)
{

/*	if(recursion_count > 7450) print_stuff_in_mss = 1;*/
	
	double se_x = b1->se_x;
	double nw_x = b1->nw_x;
	double se_y = b1->se_y;
	double nw_y = b1->nw_y;
	
	if(prev_se_x == se_x && prev_se_y == se_y && prev_nw_x == nw_x && prev_nw_y == nw_y)
	{
		same_box_counter++;
	}
	else same_box_counter = 0;
	
	if(same_box_counter > 4)
	{
		printf("Encountered the same box 5 times. Exitting!\n");
		printf("Current iteration count was: %d\n",recursion_count);
		exit(0);
	}
	
	prev_se_x = se_x;
	prev_nw_x = nw_x;
	prev_se_y = se_y;
	prev_nw_y = nw_y;
	
/*	if(b2) changed = 0;*/

	these_two_needs_freed = 0;
	recursion_count++;
	
	if(print_stuff_in_mss) printf("\n\n Number of times through recursion: %d\n\n",recursion_count);
	
	if(recursion_count > max_recursion)
	{
		printf("reached max number of iterations!!\n");
		exit(0);
	}
	
	if(print_stuff_in_mss)
	{
		printf("___________________________________________\n");
		print_inorder(tree,2);
		printf("___________________________________________\n");
	}
	
	int status = 0;
	if(recursion_count == 1)
	{
		int nzcnt = 0, surplus = 0;
		int num_starts = CPXgetnummipstarts (env, global_mip);
/*		printf("num starts: %d\n",num_starts);*/
	    	
	  	global_num_starts = num_starts;
	  	global_startspace = cur_numcols*num_starts;
	/*  	printf("number of mip starts: %d\n",num_starts);*/
	  	
	/*  	printf("first allocation\n");*/
	  	global_beg = (int *) realloc (global_beg,(global_num_starts)*sizeof(int));
		global_varindices = (int *) realloc (global_varindices,(global_startspace)*sizeof(int));
		global_values = (double *) realloc (global_values,(global_startspace)*sizeof(double));
		global_effortlevel = (int *) realloc (global_effortlevel,(global_startspace)*sizeof(int));
	
		status = CPXgetmipstarts (env, global_mip, &nzcnt, global_beg, global_varindices, 
		                   global_values, global_effortlevel, global_startspace,
		                   &surplus, 0, global_num_starts-1);
	
		lp_clone = CPXcloneprob (env, lp1, &status);
	  	if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to clone problem 1.\n");
	    		exit(0);
	  	}
	  	
	  	status = CPXaddmipstarts (env, lp_clone, global_num_starts, nzcnt, global_beg, global_varindices,
		                   global_values, global_effortlevel, NULL);
  	}
  	
  	if(print_stuff_in_mss) 
  	{
  		printf("se_x: %lf \t se_y: %lf \t nw_x: %lf \t nw_y: %lf\n",se_x,se_y,nw_x,nw_y);
  		printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",nw_x,se_x,nw_y,nw_y);
  		printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",nw_x,se_x,se_y,se_y);
  		printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",nw_x,nw_x,se_y,nw_y);
  		printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",se_x,se_x,se_y,nw_y);			
  	}
  	
/*  	status = CPXwriteprob (env, lp_clone, "myprob.lp", "LP");*/
  	
  	if( nw_x >= se_x || se_y >= nw_y)
  	{
  		if(print_stuff_in_mss) printf("here\n");
  		 return;
  	}
  	
  	if(recursion_count == 1 || se_x != SE_extreme_x || se_y != SE_extreme_y || nw_x != NW_extreme_x || nw_y != NW_extreme_y)
  	{
  		int ind[4] = {obj1_index, obj1_index, obj2_index, obj2_index};
  		char lu[4] = {'L','U','L','U'};
  		double bd[4] = {nw_x, se_x, se_y, nw_y};
  		status = CPXchgbds (env, lp_clone, 4, ind, lu, bd);
  		if ( status ) 
	  	{
	    		fprintf (stderr, "Failed to change bounds.\n");
	    		exit(0);
	  	}
  	}
  	
/*  	status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/

  	
/*  	print_preorder(tree, NULL);*/

	RECALC_NODES:
	;
  	
  	closest_nodes *these_two = find_two_nodes_right_of_val(nw_x, nw_y, tree);
  	these_two_needs_freed = 1;
  	
  	if(prev_closest_se_x == these_two->closest->se_x && prev_closest_se_y == these_two->closest->se_y && prev_closest_nw_x == these_two->closest->nw_x && prev_closest_nw_y == these_two->closest->nw_y)
	{
		same_closest_counter++;
	}
	else same_closest_counter = 0;
	
	if(same_closest_counter > 4)
	{
		printf("Encountered the same closest solution with a box 5 times. Exitting!\n");
		printf("Current iteration count was: %d\n",recursion_count);
		exit(0);
	}
	
	prev_closest_se_x = these_two->closest->se_x;
	prev_closest_nw_x = these_two->closest->nw_x;
	prev_closest_se_y = these_two->closest->se_y;
	prev_closest_nw_y = these_two->closest->nw_y;
  	
  	REDO_THIS:
  	
  	;
  	
  	if(print_stuff_in_mss) printf("at recalc\n");
  	
  	if(!these_two) 
  	{
  		these_two_needs_freed = 0;
  		exit(0);
  	}
  	
  	double cl_x = these_two->closest->nw_x;
	double cl_y = these_two->closest->nw_y;
	double ne_x = these_two->next->se_x;
	double ne_y = these_two->next->se_y;
  	
/*  	printf("segment info: %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n\n",these_two->closest->type,these_two->next->type,fabs(these_two->closest->se_x - these_two->next->se_x),fabs(these_two->closest->se_y - these_two->next->se_y),fabs(these_two->closest->se_x - these_two->next->nw_x),fabs(these_two->closest->se_y - these_two->next->nw_y),fabs(these_two->next->se_x - these_two->closest->se_x),fabs(these_two->next->se_y - these_two->closest->se_y),fabs(these_two->next->se_x - these_two->closest->nw_x),fabs(these_two->next->se_y - these_two->closest->nw_y));*/
  	
  	if(these_two->closest && these_two->next)
  	{
  		if(print_stuff_in_mss)
  		{
  			printf("closest type: %d\tvals: %lf, %lf\n",these_two->closest->type,fabs(these_two->closest->se_x - these_two->closest->nw_x),fabs(these_two->closest->se_y - these_two->closest->nw_y));
  		}
  		if(fabs(these_two->closest->se_x - these_two->closest->nw_x) < .000001 && fabs(these_two->closest->se_y - these_two->closest->nw_y) < .000001)
  		{
  			these_two->closest->type = 1;
  			these_two->closest->se_x = these_two->closest->nw_x;
  			these_two->closest->se_y = these_two->closest->nw_y;
  		}
  		if(fabs(these_two->next->se_x - these_two->next->nw_x) < .000001 && fabs(these_two->next->se_y - these_two->next->nw_y) < .000001)
  		{
  			these_two->next->type = 1;
  			these_two->next->se_x = these_two->next->nw_x;
  			these_two->next->se_y = these_two->next->nw_y;
  		}
  		if(these_two->closest->type == 1)
  		{
  			if(print_stuff_in_mss)
	  		{
	  			printf("closest type is 1. Vals: %lf, %lf, %lf, %lf\n",these_two->closest->se_x - these_two->next->se_x,fabs(these_two->closest->se_y - these_two->next->se_y),these_two->closest->se_x - these_two->next->nw_x, fabs(these_two->closest->se_y - these_two->next->nw_y));
	  		}
  			if((these_two->closest->se_x < these_two->next->se_x + .000001 && fabs(these_two->closest->se_y - these_two->next->se_y) < .000001) || (these_two->closest->se_x < these_two->next->nw_x + .000001 && fabs(these_two->closest->se_y - these_two->next->nw_y) < .000001))
  			{	// There is a weakly-dominated point stored. Eliminate it.
  				delete_node(these_two->closest);
  				these_two = find_two_nodes_right_of_val(nw_x, nw_y, tree);
  				goto REDO_THIS;
  			}
  		}
  		if(these_two->next->type == 1)
  		{
  			if(print_stuff_in_mss)
	  		{
	  			printf("next type is 1. Vals: %lf, %lf, %lf, %lf\n",these_two->next->se_x - these_two->closest->se_x,fabs(these_two->next->se_y - these_two->closest->se_y),these_two->next->se_x - these_two->closest->nw_x, fabs(these_two->next->se_y - these_two->closest->nw_y));
	  		}
  			if((these_two->next->se_x < these_two->closest->se_x + .000001 && fabs(these_two->next->se_y - these_two->closest->se_y) < .000001) || (these_two->next->se_x < these_two->closest->nw_x + .000001 && fabs(these_two->next->se_y - these_two->closest->nw_y) < .000001))
  			{	// There is a weakly-dominated point stored. Eliminate it.
  				delete_node(these_two->next);
  				these_two = find_two_nodes_right_of_val(nw_x, nw_y, tree);
  				goto REDO_THIS;
  			}
  		}
  	} 
  	
  	if(print_stuff_in_mss)
  	{
	  	printf("plot([%lf,%lf],[%lf,%lf],'-go');\n",x_ideal-these_two->closest->nw_x,x_ideal-these_two->closest->se_x,y_ideal-these_two->closest->nw_y,y_ideal-these_two->closest->se_y);
	  	printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",x_ideal-these_two->next->nw_x,x_ideal-these_two->next->se_x,y_ideal-these_two->next->nw_y,y_ideal-these_two->next->se_y);
	}
	
	if(0)//nw_y > y_ideal - these_two->closest->se_y + .000001) //use epsilon constraint to reduce y-val
	{
		chg_coefs(env, lp_clone, indices, 0.);
		
		double target_obj_val = fmax(y_ideal-these_two->closest->se_y,se_y);

		if(print_stuff_in_mss) printf("Target obj val: %lf\n",target_obj_val );
	  	
	  	CPXmipopt (env, lp_clone);
	  	num_mips_solved++;
	    	status2 = CPXgetstat(env, lp_clone);
	    	if(print_stuff_in_mss) printf("status2: %d\n",status2);    	
	    	
	    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
	    	{
			if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
	 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
	    	}

		double obvals[2] = {0.,0.};
	  	
	  	status = CPXgetx (env, lp_clone, obvals, obj1_index, obj2_index);
	  	if(status)
	  	{
/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
			printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
			exit(0);
	  		exit(0);
	  	}  	
	  	
	  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);
	  	
	  	double objval = 0.;
	  	
	  	status = CPXgetobjval (env, lp_clone, &objval);
	  	
	  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);
	  	
	  	if(recursion_count > max_recursion) exit(0);
	  	
	  	if( objval < nw_y ) // the box can be reduced
	  	{
	  		if( fabs(objval - target_obj_val) < .000001 || rightmost_val == NW_extreme_x - 10.)
	  		{
/*	  			mip_solve_sequential (env, se_x, se_y, nw_x, objval);*/
	  			exit(0);
	  		}
	  	
/*	  		node *temp = these_two->closest;*/
	  		free(these_two);
	  		these_two_needs_freed = 0;
	  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);
	  		
	  		double temp_rmv = fmin(rightmost_val,se_x);
			double temp_rmvy = fmax(rightmost_val_y,se_y);
			
/*			mip_solve_sequential (env, se_x, se_y, temp_rmv, temp_rmvy - prob_width/1000.);*/
/*			mip_solve_sequential (env, se_x, temp_rmvy - prob_width/1000., nw_x, nw_y);*/
	  		
	  		exit(0);
	  	}
	}
	
	if(these_two->closest == these_two->next || x_ideal-these_two->closest->se_x - (x_ideal-these_two->next->nw_x) > -.000001 || x_ideal-these_two->closest->se_y - (x_ideal-these_two->next->nw_y) < .00001)
	{
		if(print_stuff_in_mss) printf("The segments are the same, or next is left of closest, we are nearing the end!\n");
		goto BOTH_SAME;
	}
  	
  	if(these_two->next && ( these_two->next == these_two->closest || (fabs(these_two->closest->nw_x - these_two->next->se_x) < .000001 && 
  		fabs(these_two->closest->nw_y - these_two->next->se_y) < .000001))) //segments touch
  	{
  		if(print_stuff_in_mss) printf("The segments touch!\n");
  		
  		if(these_two->next != these_two->closest && fabs(these_two->closest->slope - these_two->next->slope) < .00001)
  		{
  			if(print_stuff_in_mss) printf("The segments also have same slope, combine them!\n");
  			these_two->closest->nw_x = these_two->next->nw_x;
  			these_two->closest->nw_y = these_two->next->nw_y;
  			delete_node(these_two->next);
  			goto RECALC_NODES;
  		}
  		
/*  		printf("vals: %lf,%lf,%lf,%lf\n",x_ideal - these_two->next->se_x , se_x , y_ideal - these_two->next->se_y , se_y);*/
  		
  		if( (x_ideal - these_two->next->se_x) - se_x > -.000001 || y_ideal - these_two->next->se_y - se_y < .000001) //next segment is outside of current box
  		{
  			if(print_stuff_in_mss) printf("next segment is outside of current box!\n");
  			
  			if(y_ideal - these_two->closest->se_y - nw_y < -.000001 && x_ideal - these_two->closest->se_x - nw_x > .000001) //there is separation beween edge of box and current segment, reduce using epsilon constraint
  		
	  		{	
	  			if(print_stuff_in_mss) printf("there is separation beween top of box and current segment\n");
	  			printf("vals: %lf, %lf, %lf, %lf\n", x_ideal - these_two->closest->se_y,x_ideal - these_two->closest->nw_y, nw_y, se_y);
	  			
	  			chg_coefs(env, lp_clone, indices, 0.);
		  	
				double target_obj_val = fmax(y_ideal-these_two->closest->se_y,se_y);
	
				if(print_stuff_in_mss) printf("Target obj val: %lf\n",target_obj_val );
			  	
			  	sol_added = 0;
			  	CPXmipopt (env, lp_clone);
			  	if(print_stuff_in_mss) printf("sol added: %d\n",sol_added);
			  	num_mips_solved++;
			    	status2 = CPXgetstat(env, lp_clone);
			    	if(print_stuff_in_mss) printf("status2: %d\n",status2);    	
			    	
			    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
			    	{
					if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
			 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
			    	}

				double obvals[2] = {0.,0.};
			  	
			  	status = CPXgetx (env, lp_clone, obvals, obj1_index, obj2_index);
			  	if(status)
			  	{
	/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
					printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
					exit(0);
			  		exit(0);
			  	}  	
			  	
			  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);
			  	
			  	double objval = 0.;
			  	
			  	status = CPXgetobjval (env, lp_clone, &objval);
			  	
			  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);
			  	
			  	if(recursion_count > max_recursion) exit(0);
			  	
			  	free(these_two);
		  		these_two_needs_freed = 0;
		  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);
		  		
		  		printf("reached line %d\n, fix!",__LINE__);
		  		exit(0);
		  		
/*		  		mip_solve_sequential (env, se_x, se_y, nw_x, objval);*/
		  		exit(0);
		  	}
  			
  			goto BOTH_SAME; //BOTH_SAME is also location where we use weighted sum method
  		}
  		
  		if(these_two->closest == these_two->next || x_ideal-these_two->closest->se_x - (x_ideal-these_two->next->nw_x) > -.000001 || x_ideal-these_two->closest->se_y - (x_ideal-these_two->next->nw_y) < .000001) 
  		{
  			if(print_stuff_in_mss) printf("The segments are the same, or next is left of closest, we are nearing the end!\n");
  			goto BOTH_SAME;
  		}
  		else //split 
  		{
  			free(these_two);
  			these_two_needs_freed = 0;
  			box *new_box1 = (struct box*) malloc( sizeof( struct box ) ); 
  			box *new_box2 = (struct box*) malloc( sizeof( struct box ) );
  			new_box1->se_x = se_x;
  			new_box1->se_y = y_ideal-cl_y;
  			new_box1->nw_x = nw_x;
  			new_box1->nw_y = nw_y;
  			new_box2->se_x = se_x;
  			new_box2->se_y = se_y;
  			new_box2->nw_x = x_ideal-cl_x;
  			new_box2->nw_y = y_ideal-cl_y;
  			
  			if(print_stuff_in_mss) 
		  	{
		  		printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",new_box2->nw_x,new_box2->se_x,new_box2->nw_y,new_box2->nw_y);
		  		printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",new_box2->nw_x,new_box2->se_x,new_box2->se_y,new_box2->se_y);
		  		printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",new_box2->nw_x,new_box2->nw_x,new_box2->se_y,new_box2->nw_y);
		  		printf("plot([%lf,%lf],[%lf,%lf],'-ro');\n",new_box2->se_x,new_box2->se_x,new_box2->se_y,new_box2->nw_y);	
		  	}
		  	
		  	if(new_box1->se_y > new_box1->nw_y)
  			{
  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
  				exit(0);
  			}
  			
  			mip_solve_sequential (env, new_box1, new_box2);		//Upper box extends up and right
  			
  			if(new_box2->se_y > new_box2->nw_y)
  			{
  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
  				exit(0);
  			}
			mip_solve_sequential (env, new_box2, NULL);
			
			free(new_box1);
			free(new_box2);
			
/*			mip_solve_sequential (env, x_ideal-cl_x, y_ideal-cl_y, nw_x, nw_y);	//Upper box only extends up*/
/*			mip_solve_sequential (env, se_x, se_y, x_ideal-cl_x, nw_y);*/
  		}
  		
/*	  	status = CPXwriteprob (env, lp_clone, "myprob.lp", "LP");*/
/*	  	exit(0);*/
/*	  	*/
/*	  	chg_coefs(env, lp_clone, indices, these_two->closest->slope);*/

/*		double target_obj_val = (x_ideal-these_two->closest->nw_x) - (1./these_two->closest->slope)*(y_ideal-these_two->closest->nw_y);*/
/*	*/
/*		if(print_stuff_in_mss) printf("Target obj val: %lf\n",target_obj_val );*/
/*	  	*/
/*	  	CPXmipopt (env, lp_clone);*/
/*	    	status2 = CPXgetstat(env, lp_clone);*/
/*	    	if(print_stuff_in_mss) printf("status2: %d\n",status2);    	*/
/*	    	*/
/*	    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)*/
/*	    	{*/
/*			if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)*/
/*	 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");*/
/*	    	}*/

/*		double obvals[2] = {0.,0.};*/
/*	  	*/
/*	  	status = CPXgetx (env, lp_clone, obvals, obj1_index, obj2_index);*/
/*	  	if(status)*/
/*	  	{*/
/*			printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);*/
/*			exit(0);*/
/*	  		exit(0);*/
/*	  	}  	*/
/*	  	*/
/*	  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);*/
/*	  	*/
/*	  	double objval = 0.;*/
/*	  	*/
/*	  	status = CPXgetobjval (env, lp_clone, &objval);*/
/*	  	*/
/*	  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);*/
/*	  	*/
/*	  	if(recursion_count > max_recursion) exit(0);*/
/*	  	*/
/*	  	if( fabs(target_obj_val - objval) < .000001 || (fabs(target_obj_val - objval) < 5 && leftmost_val == SE_extreme_x + 10.))//leftmost_val == SE_extreme_x + 10.)//(fabs(target_obj_val - objval) < 5 && leftmost_val == SE_extreme_x + 10.)) // the current line segment is Pareto, move on*/
/*	  	{*/
/*	  		free(these_two);*/
/*	  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);*/
/*	  		mip_solve_sequential (env, se_x, se_y, x_ideal-cl_x, y_ideal-cl_y);*/
/*	  	}*/
/*	  	else // the current line segment may not be Pareto, split the frame*/
/*	  	{*/
/*	  	*/
/*		*/
/*			if(print_stuff_in_mss) printf("leftmost_val: %lf\n",leftmost_val);*/
/*			*/
/*			if(leftmost_val == SE_extreme_x + 10.)*/
/*			{*/
/*			*/
/*				last_insert_num = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);			*/
/*				find_leftmost_with_insert_num(last_insert_num, tree);*/
/*				*/
/*				if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-oo');\n",leftmost_val,leftmost_val_y);*/

/*				double temp_lmv = fmax(leftmost_val,nw_x);*/
/*				double temp_lmvy = fmin(leftmost_val_y,nw_y);*/
/*				*/
/*				if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*				mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*				mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*			}*/
/*			else //new solutions were found, split based on these.*/
/*			{*/
/*			*/
/*				if(leftmost_val < x_ideal - these_two->closest->se_x)*/
/*				{*/
/*					free(these_two);*/
/*					goto RECALC_NODES;*/
/*				}*/
/*			*/
/*				double temp_lmv = fmax(leftmost_val,nw_x);*/
/*				double temp_lmvy = fmin(leftmost_val_y,nw_y);*/
/*				*/
/*				if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*				mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*				mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*			}*/
/*	  	*/
/*	  	}*/
  	}
  	else //segments don't touch
  	{
  		if(print_stuff_in_mss) printf("The segments don't touch!\n");
/*	  	status = CPXwriteprob (env, lp_clone, "myprob.lp", "LP");*/
	  	
	  	if(x_ideal-these_two->closest->se_x >= x_ideal-these_two->next->nw_x || x_ideal-these_two->closest->se_y <= x_ideal-these_two->next->nw_y)
	  	{
	  		printf("Next is left of closest, we are nearing the end!\n");
	  		goto BOTH_SAME;
	  	}
	  	
/*	  	exit(0);*/
	  	
/*	  	printf("nadir point of segments should be:\n");*/
/*	  	*/
/*	  	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x_ideal-these_two->closest->nw_x,x_ideal-these_two->closest->nw_x,y_ideal-these_two->next->se_y,y_ideal-these_two->next->se_y);*/
/*	  	*/
/*	  	//create the chebychev problem*/
/*	  	*/
/*	  	status = CPXsetintparam ((CPXENVptr) env, CPX_PARAM_SCRIND, CPX_ON);*/
/*	  	*/
/*	  	CPXLPptr cheby_prob = CPXcloneprob (env, lp_clone, &status);*/
/*	  	double beta = (NW_extreme_y - (y_ideal-these_two->next->se_y))/(NW_extreme_y - (y_ideal-these_two->next->se_y) + SE_extreme_x - (x_ideal-these_two->closest->nw_x));*/
/*	  	*/
/*	  	double neg_one = -1.;*/
/*	  	status = CPXaddcols (env, cheby_prob, 1, 0, &neg_one, NULL, NULL, NULL, NULL, NULL, NULL);*/
/*	  	if(status)*/
/*	  	{*/
/*			printf("Failed to add column to cheby prob. Line: %d. Error Code: %d\n",__LINE__,status);*/
/*			exit(0);*/
/*	  		exit(0);*/
/*	  	}  */
/*	  	*/
/*	  	double *vals = (double*)calloc(cur_numcols, sizeof(double));*/
/*	  	*/
/*	  	status = CPXchgobj (env, cheby_prob, cur_numcols, indices, vals);*/
/*	  	if(status)*/
/*	  	{*/
/*			printf("Failed to change objective coefs of cheby prob. Line: %d. Error Code: %d\n",__LINE__,status);*/
/*			exit(0);*/
/*	  		exit(0);*/
/*	  	}*/
/*	  	*/
/*	  	double rhs[2] = {-beta*(SE_extreme_x),-(1.-beta)*(NW_extreme_y)};*/
/*	  	char sense[2] = {'L','L'};*/
/*	  	int rmatbeg[2] = {0,2};*/
/*	  	int rmatind[4] = {cur_numcols,obj1_index,cur_numcols,obj2_index};*/
/*	  	double rmatval[4] = {-1.,-beta,-1.,beta-1.};*/
/*	  	*/
/*	  	status = CPXaddrows (env, cheby_prob, 0, 2, 4, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);*/
/*	  	*/
/*	  	status = CPXwriteprob (env, cheby_prob, "cheby_prob.lp", "LP");*/
/*	  	*/
/*	  	free(vals);*/
/*	  	*/
/*	  	double target_obj_val = beta*(NW_extreme_y - (y_ideal-these_two->next->se_y));*/
/*	  	*/
/*	  	if(print_stuff_in_mss) printf("Target obj val: %lf\n",target_obj_val );*/
/*			  	*/
/*	  	CPXmipopt (env, cheby_prob);*/
/*	    	status2 = CPXgetstat(env, cheby_prob);*/
/*	    	if(print_stuff_in_mss) printf("status2: %d\n",status2);  */
/*	    	*/
/*	    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)*/
/*	    	{*/
/*			if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)*/
/*	 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");*/
/*	    	}*/

/*		double obvals[2] = {0.,0.};*/
/*	  	*/
/*	  	status = CPXgetx (env, cheby_prob, obvals, obj1_index, obj2_index);*/
/*	  	if(status)*/
/*	  	{*/
/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
/*			printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);*/
/*			exit(0);*/
/*	  		exit(0);*/
/*	  	}  	*/
/*	  	*/
/*	  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);*/
/*	  	*/
/*	  	double objval = 0.;*/
/*	  	*/
/*	  	status = CPXgetobjval (env, cheby_prob, &objval);*/
/*	  	*/
/*	  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);*/
/*	  	*/
/*	  	if(recursion_count > max_recursion) exit(0);*/
/*	  	*/
/*	  	printf("___________________________________________\n");*/
/*		print_inorder(tree,2);*/
/*		printf("___________________________________________\n");*/
/*	  	*/
/*	  	if( fabs(target_obj_val - objval) < .000001)// || leftmost_val == SE_extreme_x + 10.) // chebychev method has shown the gap to be empty*/
/*	  											      // solve two new subproblems*/
/*	  	{*/
/*	  		node *temp = these_two->closest;*/
/*	  		free(these_two);*/
/*	  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);*/
/*	  		mip_solve_sequential (env, x_ideal-temp->nw_x, y_ideal-temp->nw_y, nw_x, nw_y);*/
/*	  	}*/
/*	  	else // new solutions have been found, rerun*/
/*	  	{*/
/*	  		printf("%d write code for rerun\n", __LINE__);*/
/*	  		exit(0);*/
/*	  	}*/
/*	  	*/
/*	  	exit(0);*/
		
/*		printf("vals: %lf, %lf, %lf, %lf\n",x_ideal - these_two->next->se_x, se_x,y_ideal - these_two->next->se_y , se_y);*/
	  	
  		if(x_ideal - these_two->next->se_x - se_x > .0000001 || y_ideal - these_two->next->se_y - se_y < -.0000001) //next segment is outside of current box
  		{
  			if(print_stuff_in_mss) printf("next segment is outside of current box\n");
  		
  			if(0)//x_ideal - these_two->closest->nw_x < se_x) //there is separation beween edge of box and current segment, reduce using epsilon constraint
  		
	  		{	
	  			if(print_stuff_in_mss) printf("there is separation beween edge of box and current segment\n");
	  			
	  			chg_coefs(env, lp_clone, indices, -100000000.);
		  	
			/*  	status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
			/*  	exit(0);*/

				double target_obj_val = fmin(x_ideal-these_two->closest->nw_x,se_x);
	
				if(print_stuff_in_mss) printf("Target obj val: %lf\n",target_obj_val );
			  	
			  	sol_added = 0;
			  	CPXmipopt (env, lp_clone);
			  	num_mips_solved++;
			    	status2 = CPXgetstat(env, lp_clone);
			    	if(print_stuff_in_mss) printf("status2: %d\n",status2);    	
			    	
			    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
			    	{
					if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
			 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
			    	}

				double obvals[2] = {0.,0.};
			  	
			  	status = CPXgetx (env, lp_clone, obvals, obj1_index, obj2_index);
			  	if(status)
			  	{
	/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
					printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
					exit(0);
			  		exit(0);
			  	}  	
			  	
			  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);
			  	
			  	double objval = 0.;
			  	
			  	status = CPXgetobjval (env, lp_clone, &objval);
			  	
			  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);
			  	
			  	if(recursion_count > max_recursion) exit(0);
			  	
			  	int add_check = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);
	/*		  	printf("add check val: %d\n",add_check);*/
			  	if(add_check == -1) 
				{
					x = (double *) malloc (cur_numcols*sizeof(double));
					status = CPXgetx (env, lp_clone, x, 0, cur_numcols-1);
					if(status)
				  	{
			/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
						printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
						exit(0);
				  		exit(0);
				  	}  
					insert_level++;
					leftmost_val = SE_extreme_x + 10.;
					PSA_full(env,NULL,x,NULL,NULL);
					free(x);
				}
			  	
			  	if(sol_added && fabs(target_obj_val - objval) < .000001)// || leftmost_val == SE_extreme_x + 10.) // the box can be reduced
			  	{
/*			  		node *temp = these_two->closest;*/
			  		free(these_two);
			  		these_two_needs_freed = 0;
			  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);
			  		box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
		  			new_box1->se_x = x_ideal-cl_x;
		  			new_box1->se_y = y_ideal-cl_y;
		  			new_box1->nw_x = nw_x;
		  			new_box1->nw_y = nw_y;
  					
  					if(new_box1->se_y > new_box1->nw_y)
		  			{
		  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
		  				exit(0);
		  			}
			  		mip_solve_sequential (env, new_box1, NULL);
			  		
			  		free(new_box1);
			  	}
			  	else // the current line segment may not be Pareto, split the frame
			  	{	
			  		free(these_two);
			  		these_two_needs_freed = 0;
			  		box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
		  			new_box1->se_x = se_x;
		  			new_box1->se_y = se_y;
		  			new_box1->nw_x = nw_x;
		  			new_box1->nw_y = nw_y;
		  			
		  			if(new_box1->se_y > new_box1->nw_y)
		  			{
		  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
		  				exit(0);
		  			}
			  		mip_solve_sequential (env, new_box1, NULL); //new solutions were found, rerun
			  		
			  		free(new_box1);
			  			
/*					if(leftmost_val == SE_extreme_x + 10.)*/
/*					{*/
/*			*/

/*						last_insert_num = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);*/
/*						find_leftmost_with_insert_num(last_insert_num, tree);*/
/*				*/

/*						double temp_lmv = fmax(leftmost_val,nw_x);*/
/*						double temp_lmvy = fmin(leftmost_val_y,nw_y);*/
/*				*/
/*						if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*						mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*						mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*					}*/
/*					else //new solutions were found, split based on these.*/
/*					{*/
/*						if(leftmost_val < x_ideal - these_two->closest->se_x)*/
/*						{*/
/*							free(these_two);*/
/*							goto RECALC_NODES;*/
/*						}*/
/*				*/
/*						double temp_lmv = fmax(leftmost_val,nw_x);*/
/*						double temp_lmvy = fmin(leftmost_val_y,nw_y);*/
/*					*/
/*						if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*						mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*						mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*					}*/
			  	}
		  	}
		  	else //edges of box coincide with segment boundaries, use slope to reduce
		  	{
		  		if(print_stuff_in_mss) printf("edges of box coincide with segment boundaries, use slope to reduce\n");
		  		
		  		BOTH_SAME:
		  		
		  		if(these_two->closest->type == 1) 
		  		{
		  			printf("we were going to try to use weighted sum for a point\n");
/*		  			exit(0);*/
		  			goto CHEBY;
		  		}
		  		
		  		chg_coefs(env, lp_clone, indices, these_two->closest->slope);
	  	
			/*  	status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
			/*  	exit(0);*/

				double target_obj_val = (x_ideal-these_two->closest->nw_x) - (1./these_two->closest->slope)*(y_ideal-these_two->closest->nw_y);

				if(print_stuff_in_mss) printf("(%d) Target obj val: %lf\n",__LINE__,target_obj_val );
			  	
			  	sol_added = 0;
			  	
/*			  	int nzcnt = 0, surplus = 0;*/
/*			  	global_num_starts = CPXgetnummipstarts (env, global_mip);*/
/*			  	printf("global num starts: %d\n",global_num_starts);*/
/*			  	*/
/*			  	status = CPXgetmipstarts (env, global_mip, &nzcnt, global_beg, global_varindices, */
/*						   global_values, global_effortlevel, global_startspace,*/
/*						   &surplus, 0, global_num_starts-1);*/
/*						   */
/*				status = CPXaddmipstarts (env, lp_clone, global_num_starts, nzcnt, global_beg, global_varindices,*/
/*						   global_values, global_effortlevel, NULL);*/
/*			  	*/
			  	int num_starts = CPXgetnummipstarts (env, lp_clone);
/*			  	printf("num starts: %d\n",num_starts);*/
/*			  	exit(0);*/
			  	CPXmipopt (env, lp_clone);
			  	if(print_stuff_in_mss) printf("sol added: %d\n",sol_added);
			  	num_mips_solved++;
			    	status2 = CPXgetstat(env, lp_clone);
			    	if(print_stuff_in_mss) printf("status2: %d\n",status2);    	
			    	
			    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
			    	{
					if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
			 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
			    	}

				double obvals[2] = {0.,0.};
			  	
			  	status = CPXgetx (env, lp_clone, obvals, obj1_index, obj2_index);
			  	if(status)
			  	{
		/*	  		status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
					printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
					exit(0);
			  		exit(0);
			  	}  
			  	
			  	if(prev_WS_1 == obvals[0] && prev_WS_2 == obvals[1])
				{
					same_WS_MIP_solution_counter++;
				}
				else same_WS_MIP_solution_counter = 0;
	
				if(same_WS_MIP_solution_counter > 4)
				{
					printf("Encountered the same WS MIP solution 5 times. Exitting!\n");
					printf("Current iteration count was: %d\n",recursion_count);
					exit(0);
				}
	
				prev_WS_1 = obvals[0];
				prev_WS_2 = obvals[1];	
			  	
			  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);
			  	
			  	double objval = 0.;
			  	
			  	status = CPXgetobjval (env, lp_clone, &objval);
			  	
			  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);
			  	
			  	if(recursion_count > max_recursion) exit(0);
			  	
/*			  	int add_check = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);*/
/*			  	if(print_stuff_in_mss)  printf("add check val: %d\n",add_check);*/
/*			  	if(add_check == -1) */
/*				{*/
/*					x = (double *) malloc (cur_numcols*sizeof(double));*/
/*					status = CPXgetx (env, lp_clone, x, 0, cur_numcols-1);*/
/*					if(status)*/
/*				  	{*/
/*						printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);*/
/*						exit(0);*/
/*				  		exit(0);*/
/*				  	}  */
/*					insert_level++;*/
/*					leftmost_val = SE_extreme_x + 10.;*/
/*					PSA_full(env,NULL,x,NULL,NULL);*/
/*					free(x);*/
/*				}*/
			  	
			  	if(!sol_added || fabs(target_obj_val - objval) < .0005 || (fabs(target_obj_val - objval) < 5 && leftmost_val == SE_extreme_x + 10.))//leftmost_val == SE_extreme_x + 10.)//(fabs(target_obj_val - objval) < 5 && leftmost_val == SE_extreme_x + 10.)) // the current line segment is Pareto, move on
			  	{
			  		free(these_two);
			  		these_two_needs_freed = 0;
/*			  		node *temp = these_two->closest;*/
/*			  		free(these_two);*/
/*			  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);*/
/*			  		mip_solve_sequential (env, se_x, se_y, x_ideal-temp->nw_x, y_ideal-temp->nw_y);*/
			  	}
			  	else // the current line segment may not be Pareto, split the frame
			  	{
			  		free(these_two);
			  		these_two_needs_freed = 0;
			  		box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
			  		
			  		if(b2)
			  		{
			  			if(print_stuff_in_mss) printf("(%d) modifying the next box so that fewer iterations will be needed.\n",__LINE__);
			  			
/*			  			printf("---------------------------\n");*/
/*			  			print_inorder(tree,2);*/
/*			  			printf("---------------------------\n");*/
			  			
/*			  			printf("original b2 vals: %lf,%lf\n",b2->nw_x,b2->nw_y);*/
/*			  			*/
/*			  			b2->nw_x = PSA_x_val;*/
/*			  			b2->nw_y = PSA_y_val;*/
/*			  			*/
/*			  			printf("after changing using PSA vals: %lf,%lf\n",b2->nw_x,b2->nw_y);*/
			  			
			  			change_box_corner(b2,obvals[0],obvals[1],tree);
			  			
			  			if(print_stuff_in_mss) printf("after changing using function: %lf,%lf\n",b2->nw_x,b2->nw_y);
			  		}
			  		
		  			new_box1->se_x = se_x;
		  			if(b2) new_box1->se_y = b2->nw_y;
		  			else new_box1-> se_y = obvals[1];
		  			new_box1->nw_x = nw_x;
		  			new_box1->nw_y = nw_y;
		  			
		  			if(new_box1->se_y > new_box1->nw_y)
		  			{
		  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
/*		  				exit(0);*/
		  			}
			  		mip_solve_sequential (env, new_box1, NULL); //new solutions were found, rerun
			  		
			  		free(new_box1);
			  		
/*			  		if(b2)*/
/*			  		{*/
/*			  			printf("(%d) modifying the next box so that fewer iterations will be needed.\n",__LINE__);*/
/*			  			*/
/*			  			printf("---------------------------\n");*/
/*			  			print_inorder(tree,2);*/
/*			  			printf("---------------------------\n");*/
/*			  			*/
/*			  			printf("original b2 vals: %lf,%lf\n",b2->nw_x,b2->nw_y);*/
/*			  			*/
/*			  			b2->nw_x = PSA_x_val;*/
/*			  			b2->nw_y = PSA_y_val;*/
/*			  			*/
/*			  			printf("after changing using PSA vals: %lf,%lf\n",b2->nw_x,b2->nw_y);*/
/*			  			*/
/*			  			change_box_corner(b2,obvals[0],obvals[1],tree);*/
/*			  			*/
/*			  			printf("after changing using function: %lf,%lf\n",b2->nw_x,b2->nw_y);*/
/*			  		}*/
	
/*					if(print_stuff_in_mss) printf("leftmost_val: %lf\n",leftmost_val);*/
/*		*/
/*					if(leftmost_val == SE_extreme_x + 10.)*/
/*					{*/
/*		*/
/*		*/
/*						last_insert_num = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);*/
/*						find_leftmost_with_insert_num(last_insert_num, tree);*/
/*			*/

/*						double temp_lmv = fmax(leftmost_val,nw_x);*/
/*						double temp_lmvy = fmin(leftmost_val_y,nw_y);*/
/*			*/
/*						if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*						mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*						mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*					}*/
/*					else //new solutions were found, split based on these.*/
/*					{*/
/*		*/
/*						if(leftmost_val < x_ideal - these_two->closest->se_x)*/
/*						{*/
/*							free(these_two);*/
/*							goto RECALC_NODES;*/
/*						}*/
/*		*/
/*						double temp_lmv = fmax(leftmost_val,nw_x);*/
/*						double temp_lmvy = fmin(leftmost_val_y,nw_y);*/
/*			*/
/*						if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y,fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*						mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*						mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*					}*/
			  	
			  	}	
		  	}
  		}
  		else // both segments are in the current box. Use chebychev to see if we can split.
  		{
  			if(print_stuff_in_mss) printf("both segments are in the current box. Use chebychev to see if we can split.\n");
/*  			printf("nadir point of segments should be:\n");*/
	  	
/*		  	printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",x_ideal-these_two->closest->nw_x,x_ideal-these_two->closest->nw_x,y_ideal-these_two->next->se_y,y_ideal-these_two->next->se_y);*/
		  	
		  	//create the chebychev problem
		  	
/*		  	status = CPXsetintparam ((CPXENVptr) env, CPX_PARAM_SCRIND, CPX_ON);*/

			CHEBY:
			;
		  	
/*		  	if( fabs(these_two->closest->nw_y - these_two->next->se_y) < .00001 ) //no need to solve cheby prob since y-vals are*/
/*		  									      //essentially equal*/
/*		  	{*/
/*		  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);*/
/*		  		box *new_box1 = (struct box*) malloc( sizeof( struct box ) );*/
/*	  			new_box1->se_x = se_x;*/
/*	  			new_box1->se_y = se_y;*/
/*	  			new_box1->nw_x = ne_x;*/
/*	  			new_box1->nw_y = ne_y;*/
/*	  			*/
/*	  			if(new_box1->se_y > new_box1->nw_y)*/
/*	  			{*/
/*	  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);*/
/*	  				exit(0);*/
/*	  			}*/
/*	  			*/
/*		  		mip_solve_sequential (env, new_box1, NULL);*/
/*		  		*/
/*		  		free(new_box1);*/
/*		  		*/
/*			  	exit(0); */
/*			}*/
		  	
		  	CPXLPptr cheby_prob = CPXcloneprob (env, lp_clone, &status);
		  	double beta = (NW_extreme_y - (y_ideal-these_two->next->se_y))/(NW_extreme_y - (y_ideal-these_two->next->se_y) + SE_extreme_x - (x_ideal-these_two->closest->nw_x));
		  	
		  	double neg_one = -1.;
		  	status = CPXaddcols (env, cheby_prob, 1, 0, &neg_one, NULL, NULL, NULL, NULL, NULL, NULL);
		  	if(status)
		  	{
				printf("Failed to add column to cheby prob. Line: %d. Error Code: %d\n",__LINE__,status);
				exit(0);
		  		exit(0);
		  	}  
		  	
		  	double *vals = (double*)calloc(cur_numcols, sizeof(double));
		  	
		  	status = CPXchgobj (env, cheby_prob, cur_numcols, indices, vals);
		  	if(status)
		  	{
				printf("Failed to change objective coefs of cheby prob. Line: %d. Error Code: %d\n",__LINE__,status);
				exit(0);
		  		exit(0);
		  	}
		  	
		  	double rhs[2] = {-beta*(SE_extreme_x),-(1.-beta)*(NW_extreme_y)};
		  	char sense[2] = {'L','L'};
		  	int rmatbeg[2] = {0,2};
		  	int rmatind[4] = {cur_numcols,obj1_index,cur_numcols,obj2_index};
		  	double rmatval[4] = {-1.,-beta,-1.,beta-1.};
		  	
		  	status = CPXaddrows (env, cheby_prob, 0, 2, 4, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		  	
		  	int ind[4] = {obj1_index, obj1_index, obj2_index, obj2_index};
	  		char lu[4] = {'L','U','L','U'};
	  		double bd[4] = {NW_extreme_x, SE_extreme_x, SE_extreme_y, NW_extreme_y};
	  		status = CPXchgbds (env, cheby_prob, 4, ind, lu, bd);
	  		if ( status ) 
		  	{
		    		fprintf (stderr, "Failed to change bounds.\n");
		    		exit(0);
		  	}
		  	
/*		  	status = CPXwriteprob (env, cheby_prob, "cheby_prob.lp", "LP");*/
		  	
		  	free(vals);
		  	
		  	double target_obj_val = (beta-1.)*(NW_extreme_y - (y_ideal-these_two->next->se_y));
		  	
		  	if(print_stuff_in_mss) printf("Target obj val: %lf\n",target_obj_val );
			
			sol_added = 0;	
			
			int num_starts = CPXgetnummipstarts (env, lp_clone);
/*			printf("num starts: %d\n",num_starts); */
			
			int nzcnt = 0, surplus = 0;
					    	
		  	global_num_starts = num_starts;
		  	global_startspace = cur_numcols*num_starts;
		/*  	printf("number of mip starts: %d\n",num_starts);*/
		  	
		/*  	printf("first allocation\n");*/
		  	global_beg = (int *) realloc (global_beg,(global_num_starts)*sizeof(int));
			global_varindices = (int *) realloc (global_varindices,(global_startspace)*sizeof(int));
			global_values = (double *) realloc (global_values,(global_startspace)*sizeof(double));
			global_effortlevel = (int *) realloc (global_effortlevel,(global_startspace)*sizeof(int));
	
			status = CPXgetmipstarts (env, lp_clone, &nzcnt, global_beg, global_varindices, 
				           global_values, global_effortlevel, global_startspace,
				           &surplus, 0, global_num_starts-1);
	
			 	
		  	status = CPXaddmipstarts (env, cheby_prob, global_num_starts, nzcnt, global_beg, global_varindices,
		                   global_values, global_effortlevel, NULL);
			 	
/*			printf("num starts: %d\n",num_starts); */
			 	
		  	CPXmipopt (env, cheby_prob);

			num_starts = CPXgetnummipstarts (env, cheby_prob);
/*			printf("num mip starts after running cheby: %d\n",num_starts);*/
/*			printf("num cheby starts: %d\n",num_starts);*/
			
			if(num_starts > 1)//num_starts > global_num_starts)
			{
				global_num_starts = 2*(global_num_starts+num_starts+1);
			  	global_startspace = cur_numcols*global_num_starts;
/*			  	printf("number of mip starts: %d\n",num_starts);*/
			  	
			/*  	printf("first allocation\n");*/
			  	global_beg = (int *) realloc (global_beg,(global_num_starts)*sizeof(int));
				global_varindices = (int *) realloc (global_varindices,(global_startspace)*sizeof(int));
				global_values = (double *) realloc (global_values,(global_startspace)*sizeof(double));
				global_effortlevel = (int *) realloc (global_effortlevel,(global_startspace)*sizeof(int));
	
				status = CPXgetmipstarts (env, cheby_prob, &nzcnt, global_beg, global_varindices, 
						   global_values, global_effortlevel, global_startspace,
						   &surplus, 0, num_starts-1);
						   
/*				printf("status: %d, nzcnt: %d, surplus: %d\n",status,nzcnt,surplus);*/
	
				 	
			  	status = CPXaddmipstarts (env, lp_clone, num_starts, nzcnt+1, global_beg, global_varindices,
		                   global_values, global_effortlevel, NULL);
		                   
/*		                printf("status: %d, nzcnt: %d, surplus: %d\n",status,nzcnt,surplus);*/
		                   
		                num_starts = CPXgetnummipstarts (env, lp_clone);
/*				printf("num mip starts after adding to lp_clone: %d\n",num_starts);*/
			}

		  	num_mips_solved++;
		    	status2 = CPXgetstat(env, cheby_prob);
		    	if(print_stuff_in_mss) printf("status2: %d\n",status2);  
		    	
		    	if(status2 != 101 && status2 != 102 && status2 != 129 && status2 != 130)
		    	{
				if(status2 == 104 || status2 == 105 || status2 == 107 || status2 == 109 || status2 == 111 || status2 == 113)
		 			printf("Warning: CPLEX failed to reach an optimal solution. A suboptimal integer-feasible solution is being used.\n");
		    	}

			double obvals[2] = {0.,0.};
		  	
		  	status = CPXgetx (env, cheby_prob, obvals, obj1_index, obj2_index);
		  	if(status)
		  	{
	/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
				printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
				exit(0);
		  		exit(0);
		  	}
		  	
		  	if(prev_CH_1 == obvals[0] && prev_CH_2 == obvals[1])
			{
				same_CH_MIP_solution_counter++;
			}
			else same_CH_MIP_solution_counter = 0;

			if(same_CH_MIP_solution_counter > 4)
			{
				printf("Encountered the same Chebychev MIP solution 5 times. Exitting!\n");
				printf("Current iteration count was: %d\n",recursion_count);
				exit(0);
			}

			prev_CH_1 = obvals[0];
			prev_CH_2 = obvals[1];	  	
		  	
		  	if(print_stuff_in_mss) printf("plot([%lf],[%lf],'-mo');\n",obvals[0],obvals[1]);
		  	
		  	double objval = 0.;
		  	
		  	status = CPXgetobjval (env, cheby_prob, &objval);
		  	
		  	int add_check = mock_insert(1,obvals[0],obvals[1],0,0,0,&tree);
/*		  	printf("add check val: %d\n",add_check);*/
		  	if(add_check == -1) 
			{
				x = (double *) malloc (cur_numcols*sizeof(double));
				status = CPXgetx (env, cheby_prob, x, 0, cur_numcols-1);
				if(status)
			  	{
		/*				status = CPXwriteprob (env, lp_clone, "myprob2.lp", "LP");*/
					printf("Failed to get x-values from CPLEX. Line: %d. Error Code: %d\n",__LINE__,status);
					exit(0);
			  		exit(0);
			  	}  
				insert_level++;
				leftmost_val = SE_extreme_x + 10.;
				PSA_full(env,NULL,x,NULL,NULL);
				free(x);
			}
			if(print_stuff_in_mss) printf("sol added: %d\n",sol_added);
		  	
		  	if(print_stuff_in_mss) printf("(%d) objval: %lf\n",__LINE__,objval);
		  	
		  	if(recursion_count > max_recursion) exit(0);
		  	
/*		  	printf("___________________________________________\n");*/
/*			print_inorder(tree,2);*/
/*			printf("___________________________________________\n");*/
		  	
		  	if(!sol_added || fabs(target_obj_val - objval) < .0005)// || leftmost_val == SE_extreme_x + 10.) // chebychev method has shown the gap to be empty
		  											      // solve two new subproblems
		  	{
		  		free(these_two);
		  		these_two_needs_freed = 0;
		  		if(print_stuff_in_mss) printf("reducing box size at %d\n",__LINE__);
		  		box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
	  			new_box1->se_x = x_ideal-cl_x;
	  			new_box1->se_y = y_ideal-cl_y;
	  			new_box1->nw_x = nw_x;
	  			new_box1->nw_y = nw_y;
	  			box *new_box2 = (struct box*) malloc( sizeof( struct box ) );
	  			new_box2->se_x = se_x;
	  			new_box2->se_y = se_y;
	  			new_box2->nw_x = y_ideal-ne_x;
	  			new_box2->nw_y = y_ideal-ne_y;
	  			
	  			if(print_stuff_in_mss) printf("splitting and creating two boxes\n");
	  			if(print_stuff_in_mss) 
			  	{
			  		printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",new_box1->nw_x,new_box1->se_x,new_box1->nw_y,new_box1->nw_y);
			  		printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",new_box1->nw_x,new_box1->se_x,new_box1->se_y,new_box1->se_y);
			  		printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",new_box1->nw_x,new_box1->nw_x,new_box1->se_y,new_box1->nw_y);
			  		printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",new_box1->se_x,new_box1->se_x,new_box1->se_y,new_box1->nw_y);	
			  	}
			  	if(print_stuff_in_mss) 
			  	{
			  		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",new_box2->nw_x,new_box2->se_x,new_box2->nw_y,new_box2->nw_y);
			  		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",new_box2->nw_x,new_box2->se_x,new_box2->se_y,new_box2->se_y);
			  		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",new_box2->nw_x,new_box2->nw_x,new_box2->se_y,new_box2->nw_y);
			  		printf("plot([%lf,%lf],[%lf,%lf],'-ko');\n",new_box2->se_x,new_box2->se_x,new_box2->se_y,new_box2->nw_y);	
			  	}
	  			
	  			if(new_box1->se_y > new_box1->nw_y)
	  			{
	  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
/*	  				exit(0);*/
	  			}
	  			
		  		mip_solve_sequential (env, new_box1, NULL);
		  		
		  		if(new_box2->se_y > new_box2->nw_y)
	  			{
	  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
/*	  				exit(0);*/
	  			}
		  		mip_solve_sequential (env, new_box2, NULL);
		  		
		  		free(new_box1);
		  		free(new_box2);
		  	}
		  	else // new solutions have been found, rerun
		  	{
		  		if(objval != prev_cheby_sol)
		  		{
		  			prev_cheby_sol = objval;
		  			free(these_two);
		  			these_two_needs_freed = 0;
		  			box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
		  			new_box1->se_x = se_x;
		  			new_box1->se_y = se_y;
		  			new_box1->nw_x = nw_x;
		  			new_box1->nw_y = nw_y;
	  			
	  				if(new_box1->se_y > new_box1->nw_y)
		  			{
		  				printf("at line %d we generated a new box that will not be oriented correctly\n",__LINE__);
		  				exit(0);
		  			}
		  			mip_solve_sequential (env, new_box1, NULL);
		  			
		  			free(new_box1);
		  		}
		  		else
		  		{
		  			printf("%d same cheby sol twice, figure out what to do!\n",__LINE__);
		  			exit(0);
		  		}
/*		  		printf("%d write code for rerun\n", __LINE__);*/
/*		  		exit(0);*/
		  	}
		  	
/*		  	printf("exitting %d\n",__LINE__);*/
/*		  	exit(0);*/
  		
/*  			double temp_lmv = x_ideal - these_two->next->se_x;*/
/*			double temp_lmvy = y_ideal - these_two->next->se_y;*/
/*			*/
/*			if(temp_lmv == x_ideal - these_two->closest->nw_x) //x-values match*/
/*			{*/
/*				if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, temp_lmv, nw_y, temp_lmv, temp_lmvy, nw_x, nw_y);*/
/*				mip_solve_sequential (env, se_x, se_y, temp_lmv, temp_lmvy + prob_width/1000.);*/
/*				mip_solve_sequential (env, se_x, temp_lmvy + prob_width/1000., nw_x, nw_y);*/
/*			}*/
/*			else*/
/*			{*/
/*			*/
/*				if(print_stuff_in_mss) printf("(%d) splitting into two sub boxes\nFirst: %lf,%lf,%lf,%lf\nSecond: %lf,%lf,%lf,%lf\n",__LINE__,se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*				mip_solve_sequential (env, se_x, se_y, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), nw_y);*/
/*				mip_solve_sequential (env, fmax(temp_lmv - prob_width/1000., (nw_x+temp_lmv)/2.), temp_lmvy, nw_x, nw_y);*/
/*			}*/
	 	}
  	}
  	
  	if(these_two_needs_freed) free(these_two);
  	
  	if(print_stuff_in_mss) printf("here 1\n");
/*  	exit(0);*/
  	
  	
  	;
/*  	TERMINATE:*/
  	if(print_stuff_in_mss) printf("here 2\n");
/*  	free(these_two);*/
/*  	exit(0);*/
  	
  	return;
} /* End of mip_solve_sequential */

int main(int argc, char **argv)
{
	/* initialize Cplex environment *********************************/

  	env = CPXopenCPLEX(&status);

  	if(env==NULL)
  	{
    		char errmsg[1024];
    		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
    		CPXgeterrorstring (env, status, errmsg);
    		printf ("%s", errmsg);
    		exit(0);
  	}

/*  	if ((output = fopen (outfilename, "w+"))==NULL)*/
/*    	{*/
/*      		fprintf (stderr, "could not open output file, exiting\n");*/
/*      		exit(1);*/
/*    	}*/
    	/******************************************************************/

  	/************* Set to 1 thread **********************************/
  	
	status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	if ( status ) {
		printf ("Failure to set threads to 1, error %d.\n",status);
		exit(0);
	}
  	
  	/******************************************************************/
  	
  	/************* Set any desired CPLEX parameters here **************/
  	
	status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_AbsGap, 0.0);
	if ( status ) {
		printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		exit(0);
	}
	status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_RelGap, 0.0);
	if ( status ) {
		printf ("Failed to set solution pool gap to 0, error %d.\n",status);
		exit(0);
	}
	
	status = CPXsetincumbentcallbackfunc(env, userincumbent, NULL);
  	
  	/******************************************************************/
  	
  	/** Create test Problems *******************/

  	lp1 = CPXcreateprob (env,&status,argv[1]);
  	if(lp1==NULL) 
  	{
    		printf("CPXcreateprob, Failed to create LP1, error code %d\n", status);
    		exit(0);
    	}

  	lp2 = CPXcreateprob (env,&status,argv[2]);
  	if(lp2==NULL) 
  	{
    		printf("CPXcreateprob, Failed to create LP2, error code %d\n", status);
    		exit(0);
    	}
	/******************************************************************/

  	/** read the testP1 and testP2 from the model file **********/

  	status = CPXreadcopyprob(env,lp1,argv[1],NULL);

  	if ( status ) 
  	{
    		printf ("Could not read input 1, exiting. Error code: %d\n", status);
    		exit(0);
    	}

  	status = CPXreadcopyprob(env,lp2,argv[2],NULL);

  	if ( status ) 
  	{
    		printf ("Could not read input 2, exiting. Error code: %d\n", status);
    		exit(0);
    	}
	/******************************************************************/
  
  	/******************************************************************/
  	// Make sure the problems are in maximization form with less than
  	// or equal to and/or equality constraints.
  	/******************************************************************/

  	start_time = clock();
  	
/*  	status = CPXwriteprob (env, lp1, "myprob.lp", "LP");*/
	
	cur_numcols = CPXgetnumcols (env, lp1);
	cur_numrows = CPXgetnumrows (env, lp1);
	
	indices = (int*)calloc(cur_numcols*2, sizeof(int));
	
	char *row_sense;
    	row_sense = (char *) malloc ((cur_numrows)*sizeof(char));
    	status = CPXgetsense (env, lp1, row_sense, 0, cur_numrows-1);
    	
    	int i = 0, j = 0, equalities = 0;
    	for(i = 0;i<cur_numrows;i++)
    	{
    		if(row_sense[i] == 'G')
    		{
    			double coef = 0.;
    			double rhs[1] = {0.};
    			char new_sense[1] = {'L'};
    			for(j=0;j<cur_numcols;j++)
    			{
    				status = CPXgetcoef (env, lp1, i, j, &coef);
    				if(coef != 0.)
    				{
    					status = CPXchgcoef (env, lp1, i, j, -coef);
    					status = CPXchgcoef (env, lp2, i, j, -coef);
    				}
    			}
    			status = CPXgetrhs (env, lp1, rhs, i, i);
    			rhs[0] = -rhs[0];
    			status = CPXchgrhs (env, lp1, 1, &i, rhs);
    			status = CPXchgsense (env, lp1, 1, &i, new_sense);
    			status = CPXchgrhs (env, lp2, 1, &i, rhs);
    			status = CPXchgsense (env, lp2, 1, &i, new_sense);
    		}
    		if(!equalities && row_sense[i] == 'E')
    		{
    			equalities = 1;
    		}
    	}
    	for(i = 0;i < cur_numcols*2; i++) indices[i] = i;

  	obj_coef1 = (double *) malloc (cur_numcols*sizeof(double));

 	status = CPXgetobj(env, lp1, obj_coef1, 0, cur_numcols -1);

  	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP1, error code %d\n", status);
    		exit(0);
    	}

  	obj_coef2 = (double *) malloc (cur_numcols*sizeof(double));

  	status = CPXgetobj(env,lp2, obj_coef2, 0, cur_numcols -1);

  	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP2, error code %d\n", status);
    		exit(0);
    	}
    	
    	status = CPXgetobjsen (env, lp1);
    	
    	if ( status == 0) 
  	{
    		printf ("When trying to get objective sense, found that no problem object exists, error code %d\n", status);
    		exit(0);
    	}
    	else if( status == 1)	// Problem is minimization type - switch to maximization
    	{	
    		for(i=0;i<cur_numcols;i++)
    		{
    			obj_coef1[i] = -obj_coef1[i];
    			obj_coef2[i] = -obj_coef2[i];
    		}
    		status = CPXchgobj (env, lp1, cur_numcols, indices, obj_coef1);
    		if ( status ) 
	  	{
	    		printf ("Error changing objective 1, error code %d\n", status);
	    		exit(0);
	    	}
    		status = CPXchgobj (env, lp2, cur_numcols, indices, obj_coef2);
    		if ( status ) 
	  	{
	    		printf ("Error changing objective 2, error code %d\n", status);
	    		exit(0);
	    	}
    		status = CPXchgobjsen (env, lp1, CPX_MAX);
    		if ( status ) 
	  	{
	    		printf ("Error changing sense of objective 1, error code %d\n", status);
	    		exit(0);
	    	}
    		status = CPXchgobjsen (env, lp2, CPX_MAX);
    		if ( status ) 
	  	{
	    		printf ("Error changing sense of objective 2, error code %d\n", status);
	    		exit(0);
	    	}
    	}
    	
/*    	status = CPXwriteprob (env, lp1, "myprob2.lp", "LP");*/
/*    	exit(0);*/
    	
    	
    	/******************************************************************/
  	
  	/********* Turn Parallel Mode On or Off ***************************/
  	
/*  	int parallel = atoi(argv[3]);*/
/*  	//printf("parallel: %d\n",parallel);*/
/*  	*/
/*  	if(parallel != 0 && parallel != 1) */
/*  	{*/
/*  		printf("Third command-line argument should be:\n" */
/*  		       "0 - indicating non-parallel mode\n or\n"*/
/*  		       "1 - indicating parallel mode\n\n Error!!\n");*/
/*  		exit(0);*/
/*  	}*/
  	
  	/******************************************************************/
  	
  	char *xctype = (char *) malloc ((cur_numcols+2)*sizeof(char));
  	
  	status = CPXgetctype (env, lp1, xctype, 0, cur_numcols-1);
    	if ( status ) 
  	{
    		printf ("(%d) CPXgetctype, Failed to get variable types, error code %d\n", __LINE__,status);
    		exit(0);
    	}
    	
    	provide_xctype(env, xctype,cur_numcols);
  	
  	/****************************************************
    		 Add two new variables to keep track of
    		 the objective function values.
    	*****************************************************/
	
	double lb[2] = {-CPX_INFBOUND,-CPX_INFBOUND};
	double ub[2] = {CPX_INFBOUND,CPX_INFBOUND};

	double x1_[2] = {0.,0.};
	double x2_[2] = {0.,0.};
	
	status = (CPXnewcols (env, lp1, 2, NULL, lb, ub, NULL, NULL) || CPXnewcols (env, lp2, 2, NULL, lb, ub, NULL, NULL)); 
	if ( status ) 
  	{
    		printf ("CPXnewcols, Failed to add additional variables to the model, error code %d\n", status);
    		exit(0);
    	}
	
	cur_numcols += 2;
	obj1_index = cur_numcols - 2;
	obj2_index = cur_numcols - 1;
	
	char sense[4] = {'E','E'};
	
	status = (CPXnewrows (env, lp1, 2, NULL, sense, NULL, NULL) || CPXnewrows (env, lp2, 2, NULL, sense, NULL, NULL) ); 
	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
    		exit(0);
    	}
    	
    	int obj1_row_index = cur_numrows;
    	int obj2_row_index = cur_numrows+1;
    	
    	cur_numrows += 2;
    	
    	int *rowlist;
    	int *collist;
    	double *vallist;
    	
    	rowlist = (int *) malloc (cur_numcols*sizeof(int));
    	collist = (int *) malloc (cur_numcols*sizeof(int));
    	vallist = (double *) malloc (cur_numcols*sizeof(double));
    	
    	/******************************* Set coefficients for first new row *****************/
    	for(i=0;i<cur_numcols;i++)
    	{
    		rowlist[i] = cur_numrows - 2;
    		collist[i] = i;
    		if(i < cur_numcols - 2) vallist[i] = obj_coef1[i];
    		else if(i == cur_numcols - 2) vallist[i] = -1;
    		else vallist[i] = 0;
    	}
    	
    	status = (CPXchgcoeflist (env, lp1, cur_numcols, rowlist, collist, vallist) || CPXchgcoeflist (env, lp2, cur_numcols, rowlist, collist, vallist));
    	if ( status ) 
  	{
    		printf ("CPXgchgcoeflist, Failed to set coefficients for first new row, error code %d\n", status);
    		exit(0);
    	}
    	
    	double rh1 = 0.;
    	double rh2 = 0.;
    	
    	status = (CPXchgrhs (env, lp1, 1, &rowlist[0], &rh1) || CPXchgrhs (env, lp2, 1, &rowlist[0], &rh1));
    	/*************************************************************************************/
    	/******************************* Set coefficients for second new row *****************/
    	for(i=0;i<cur_numcols;i++)
    	{
    		rowlist[i] = cur_numrows - 1;
    		collist[i] = i;
    		if(i < cur_numcols - 2) vallist[i] = obj_coef2[i];
    		else if(i == cur_numcols - 1) vallist[i] = -1;
    		else vallist[i] = 0;
    	}
    	
    	status = ( CPXchgcoeflist (env, lp1, cur_numcols, rowlist, collist, vallist) || CPXchgcoeflist (env, lp2, cur_numcols, rowlist, collist, vallist));
    	if ( status ) 
  	{
    		printf ("CPXgchgcoeflist, Failed to set coefficients for second new row, error code %d\n", status);
    		exit(0);
    	}
    	
    	status = (CPXchgrhs (env, lp1, 1, &rowlist[0], &rh2) || CPXchgrhs (env, lp2, 1, &rowlist[0], &rh2));
    	
    	/******************************************************************/
  	
  	/********* Initialization: Compute Extreme Pareto Points **********/
  	
  	compute_feas_extremes   (env,
  				lp1,
  				lp2,
  				obj_coef1,
  				obj_coef2);
  				
/*  	print_inorder(tree,2);*/
/*  	exit(0);*/
  				
  	/******************************************************************/
  	
  	/***************************************************************************************/
  	/* 
  	 _  _  __  ____      ____   __   ____  ____  ____    ____   __   __    _  _  ____  ____
  	( \/ )(  )(  _ \ ___(  _ \ / _\ / ___)(  __)(    \  / ___) /  \ (  )  / )( \(  __)(  _ \      
	/ \/ \ )(  ) __/(___)) _ (/    \\___ \ ) _)  ) D (  \___ \(  O )/ (_/\\ \/ / ) _)  )   /      
	\_)(_/(__)(__)      (____/\_/\_/(____/(____)(____/  (____/ \__/ \____/ \__/ (____)(__\_) 
	  
	*/
	/***************************************************************************************/
/*  	if(parallel == 0) */
	box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
	new_box1->se_x = SE_extreme_x;
	new_box1->se_y = SE_extreme_y;
	new_box1->nw_x = NW_extreme_x;
	new_box1->nw_y = NW_extreme_y;
	  			
  	mip_solve_sequential (env, new_box1, NULL);
  	
  	free(new_box1);
/*              		      nPoints,*/
/*              		      &nNadir);*/
/*        else */
/*        mip_solve_parallel   (env,*/
/*              		      lp1,*/
/*              		      lp2,*/
/*              		      obj_coef1,*/
/*              		      obj_coef2);*/
/*              		      nPoints,*/
/*              		      &nNadir);*/

	printf("IT FINISHED!!\n\n");
	finish_time = clock();
	double duration = (double)(finish_time - start_time) / CLOCKS_PER_SEC;
	   	printf("Total time: %lf\n",duration);
	printf("Iterations: %d\n",recursion_count);
	printf("Num MIPs solved: %d\n\n\n",num_mips_solved);
	int n = get_num_nodes(tree);
	printf("Number of nodes in tree: %d\n",n);
	print_inorder(tree,1);

        /**************************************************************************************/
  	
  	/**********************Code for TERMINATE*************************/
/*  	TERMINATE:*/

  	/* Free the problem as allocated by CPXcreateprob and
     	CPXreadcopyprob, if necessary *************************************/

  	if ( lp1 != NULL || lp2 != NULL ) 
  	{
    		status = CPXfreeprob (env, &lp1) || CPXfreeprob (env, &lp2);
    		if ( status ) 
    		{
      			printf ("CPXfreeprob failed, error code %d.\n",status);
    		}
  	}
  	/******************************************************************/

  	/* Free the CPLEX environment, if necessary ***********************/

  	if(env!=NULL) status=CPXcloseCPLEX(&env);
  	/******************************************************************/

	/******************* Print Solutions to File **********************/
	
/*  	print_preorder (tree, output);*/
  	/******************************************************************/

	/***************** Close Files ************************************/

/*  	if (output) fclose(output);*/
/*  	if (init_nadir) fclose(init_nadir);*/
  	/*  fclose(all_inserted);*/
  	/*  fclose(testing_struct);*/
/*  	if (init_sol) fclose(init_sol);*/
  	/*  close_files();*/
  	/******************************************************************/
  	
  	/******************* Free Variables *******************************/

/*  	free_and_null ((char **) &nadirX_p);*/
/*  	free_and_null ((char **) &ParetoX_p);*/
/*  	//free_and_null ((char **) &Theta_k_ind);*/
  	destroy_tree (tree);
/*  	free_and_null ((char **) &ctype);*/
  	free_and_null ((char **) &obj_coef1);
	free_and_null ((char **) &obj_coef2);
	free_and_null ((char **) &row_sense);
	free_and_null ((char **) &xctype);
	free_and_null ((char **) &rowlist);
    	free_and_null ((char **) &collist);
    	free_and_null ((char **) &vallist);
    	free_and_null ((char **) &integer_indices);
    	free_and_null ((char **) &indices);
    	free_and_null ((char **) &weighted_coefs);
/* 	free_and_null ((char **) &obj_coef3);*/
/* 	 free_and_null ((char **) &theta);*/
/*  	free_and_null ((char **) &indices_new);*/
/*  	free_and_null ((char **) &indices_still_to_check);*/
/*  	free_and_null ((char **) &indices_f_n_l);*/
/*	delete_split_pts(first_split_pt);*/
  	/******************************************************************/

  	return 0;
  	/************** End of TERMINATE **********************************/
}

/**************************************************************/
/** free and nulll                                           **/
/**************************************************************/

void free_and_null (char **ptr)
{
  if ( *ptr != NULL ) {
    free (*ptr);
    *ptr = NULL;
  }
} /* END free_and_null */
