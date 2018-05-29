#include<stdlib.h>
#include<stdio.h>
#include <math.h>
#include<time.h>

#include "cplex.h"
#include "callbacks.h"
#include "bb-bicriteria.h"
#include "max_tree.h"

CPXLPptr   global_mip=NULL;
int *global_beg;
int *global_varindices; 
double *global_values; 
int *global_effortlevel;
int global_num_starts;
int global_startspace;
int last_insert_num;

int CPXPUBLIC
userincumbent (CPXCENVptr env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle,
               double     objval,
               double     *x,
               int        *isfeas_p,
               int        *useraction_p)
{

	*useraction_p = CPX_CALLBACK_DEFAULT;
	int status = 0;
	
	CPXLPptr nodelp = NULL;
	status = CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
	if ( status ) {
		printf ("CPXgetcallbacknodelp, Failed to catch nodelp, error code %d\n", status);
	}
/*	printf("in userincumbent, wherefrom is: %d\n",wherefrom);*/
/*	printf("in userincumbent, checking: %lf, %lf\n",x[obj1_index],x[obj2_index]);*/
	
	leftmost_val = SE_extreme_x + 10.;
	rightmost_val = NW_extreme_x - 10.;
	int add_check = mock_insert(1,x[obj1_index],x[obj2_index],0,0,0,&tree);
	int i;
	
/*	printf("add_check value: %d\n",add_check);*/

	if(add_check == -1) 
	{
		insert_level++;
		leftmost_val = SE_extreme_x + 10.;
		PSA_full(env,NULL,x,NULL,NULL);
	}
	last_insert_num = add_check;
	
	return (status);

} /* END userincumbent */

void chg_coefs(CPXCENVptr env, CPXLPptr prob, int *indices, double slope)
{
	int i,status;
	if(slope == 0.) 
	{
/*		printf("doing this\n");*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs[i] = obj_coef2[i];
		}
		weighted_coefs[obj1_index] = 0.;
		weighted_coefs[obj2_index] = 0.;
	}
	else if(slope < -10000000.)
	{
/*		printf("doing this2\n");*/
		for(i=0;i<obj1_index;i++)
		{
			weighted_coefs[i] = obj_coef1[i];
		}
		weighted_coefs[obj1_index] = 0.;
		weighted_coefs[obj2_index] = 0.;
	}
	else
	{
		for(i=0;i<obj1_index;i++) weighted_coefs[i] = 0.;
		weighted_coefs[obj1_index] = 1.;
		weighted_coefs[obj2_index] = -1./slope;
	}

/*	if(prob == lp_1 || prob == lp_2) status = CPXchgobj (env, prob, obj1_index, indices, weighted_coefs);*/
/*	else status = CPXchgobj (env, prob, obj1_index+2, indices, weighted_coefs);*/
	status = CPXchgobj (env, prob, obj1_index+2, indices, weighted_coefs);
	if ( status ) {
		printf ("Failed to change obj coef. Error code %d\n", status);
	}
/*	status = CPXwriteprob (env, prob, "myprob2.lp", "LP");*/
}

/*********************************************************************************************************************** 

	This function is used for generating integer feasible line segments that are part of the Pareto
	set of a slice problem. (Works from right to left)
	
***********************************************************************************************************************/

void PSA_full_right(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	if(!changed) printf("we'll be changing a corner of a box!\n");*/

/*	print_on = 1;*/
/*	printf("PSA_full right\n");*/

	int iter_cnt = 0, status = 0;
	int changed_locally = 0;
	
	if( !prob ) 
	{
		prob = CPXcloneprob (env, lp1, &status);
		status = CPXchgprobtype(env, prob, 0);
	}
	
	double first_obj1_val, first_obj2_val;
	
	/*************** Getting started ***************/
	
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	
	/*************** Fix the integer variables ***************/
	
	int i;
	char *b;
	double *integer_x;
/*	int *indices;*/
	b = (char *) malloc (total_num_integer*sizeof(char));
	integer_x = (double *) malloc (total_num_integer*sizeof(double));
/*	indices = (int *) malloc (cur_numcols * sizeof (int));*/
			
	for(i=0;i<total_num_integer;i++)
	{
		b[i] = 'B';
		integer_x[i] = x_orig[integer_indices[i]];
	}
	
/*	printf("changing bounds 1\n");*/
	status = CPXchgbds (env, prob1, total_num_integer, integer_indices, b, integer_x);
	
/*	char up[1] = {'U'};*/
	
/*	status = CPXchgbds (env, prob1, 1, &obj1_index, up, &x_orig[obj1_index]);*/
/*	if ( status ){*/
/*		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);*/
/*		 exit(0);}*/
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
/*  	for(i=0;i<cur_numcols;i++) indices[i] = i;*/
  	
  	chg_coefs(env,prob1,indices,-10000000.);
  	chg_coefs(env,prob2,indices,-.000000001);
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		exit(0);
/*    		exit(0);*/
  	}
  	
  	int lpstat = CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
/*  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");*/
  			PSA_full_right(env,NULL,x_orig,NULL,NULL);
  			exit(0);
/*  			exit(0);*/
  		}
  	}
  	else if(lpstat == 3) exit(0);
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;
  	
/*  	prev_obj_vals[0] = x_orig[obj1_index];*/
/*  	prev_obj_vals[1] = x_orig[obj2_index];*/

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
/*	if(!changed)*/
/*	{*/
/*		PSA_x_val = prev_obj_vals[0];*/
/*		PSA_y_val = prev_obj_vals[1];*/
/*	}*/
	
/*	printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		exit(0);
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
/*  	printf("plot([%lf],[%lf],'ro');\n",obj_vals[0],obj_vals[1]);*/
/*  	printf("differences: %lf,%lf\n",fabs(obj_vals[0] - prev_obj_vals[0]),fabs(obj_vals[1] - prev_obj_vals[1]));*/
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) > .0001 && fabs(obj_vals[1] - prev_obj_vals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",obj_vals[0],obj_vals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
  		insert(1,prev_obj_vals[0],prev_obj_vals[1],0.,0.,0.,&tree,NULL);
  		exit(0);
  	}
	
	int same_cnt = 0;
	double prev_coef_lb = -100000000000000000., add_val = .001;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .0001 )
	{
		iter_cnt++;
/*		printf("iteration: %d\n",iter_cnt);*/
		if(iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prev_obj_vals[0], obj_vals[0]);
			printf("same count: %d\n",same_cnt);
/*			printf("%d\n",__LINE__);*/
			exit(0);
		}
/*		printf("prev_obj_vals[0]: %lf, obj_vals[0]: %lf\n",prev_obj_vals[0],obj_vals[0]);*/
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(coef_ub > 10000000000.) 
		{	
/*			printf("coef_ub too big\n");*/
			break;
		}
		if(fabs(prev_coef_lb - coef_lb) < .00000001)
		{
			coef_ub += add_val;
			add_val = add_val*2;
/*			retval = 0;*/
/*			exit(0);*/
		}
		prev_coef_lb = coef_lb;
		chg_coefs(env,prob1,indices, -1./(coef_ub + .0000001*pow(100,(same_cnt+1))));
		RESOLVE:
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			exit(0);
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
	  	if(lpstat == 12)
	  	{
	  		double upper_limit,lower_limit,what_it_should_be;
	  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
	  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
	  		what_it_should_be = pow(10.,75.);
	  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
	  		{
	  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
	  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
/*	  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");*/
	  			goto RESOLVE;
	  		}
	  	}
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/

/*		printf("plot([%lf],[%lf],'ko');\n",new_obj_vals[0],new_obj_vals[1]);*/

/*		if(!changed)*/
/*		{*/
/*			PSA_x_val = new_obj_vals[0];*/
/*			PSA_y_val = new_obj_vals[1];*/
/*		}*/
		
		if(new_obj_vals[0] < prev_obj_vals[0])
		{
/*			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);*/
/*			if(check_for_stopping_PSA_full)*/
/*			{*/
/*				int insert_check = mock_insert(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],*/
/*				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree);*/
/*				if(!insert_check) */
/*				{*/
/*					printf("Stopping PSA full right early\n");	*/
/*					exit(0);*/
/*				}*/
/*			}*/
			insert(2,prev_obj_vals[0],prev_obj_vals[1],new_obj_vals[0],new_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree,NULL);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else 
		{
/*			printf("same cnt: %d\n",same_cnt);*/
			same_cnt++;
		}
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) 
		{	
			insert(2,prev_obj_vals[0],prev_obj_vals[1],obj_vals[0],obj_vals[1],
				(prev_obj_vals[1]-obj_vals[1])/(prev_obj_vals[0]-obj_vals[0]),&tree,NULL);
			break;
		}
	}
 	
/* 	TERMINATE:*/
 	
 	changed = 1;
 	
 	free_and_null ((char **) &b);
	free_and_null ((char **) &integer_x);
/*	free_and_null ((char **) &indices);*/
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

/*********************************************************************************************************************** 

	This function is used for generating integer feasible line segments that are part of the Pareto
	set of a slice problem. (Works from left to right)
	
***********************************************************************************************************************/

void PSA_full_left(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	print_on = 1;*/
/*	printf("PSA_full left\n");*/

	int iter_cnt = 0, status = 0;
	
	if( !prob )
	{
		prob = CPXcloneprob (env, lp1, &status);
		status = CPXchgprobtype(env, prob, 0);
	}
	
	double first_obj1_val, first_obj2_val;
	
	/*************** Getting started ***************/
	
	CPXLPptr prob1 = CPXcloneprob (env, prob, &status);
	
	/*************** Fix the integer variables ***************/
	
	int i;
	char *b;
	double *integer_x;
/*	int *indices;*/
	b = (char *) malloc (total_num_integer*sizeof(char));
	integer_x = (double *) malloc (total_num_integer*sizeof(double));
/*	indices = (int *) malloc (cur_numcols * sizeof (int));*/
	
/*	printf("_________----------___________\n");*/
/*	for(i=0;i<cur_numcols;i++) printf("x_%d: %lf\n",i,x_orig[i]);*/
/*	printf("_________----------___________\n");*/
			
	for(i=0;i<total_num_integer;i++)
	{
		b[i] = 'B';
		integer_x[i] = x_orig[integer_indices[i]];
/*		printf("x%d: %lf\n",integer_indices[i],x_orig[integer_indices[i]]);*/
	}
	
/*	printf("changing bounds 1\n");*/
	status = CPXchgbds (env, prob1, total_num_integer, integer_indices, b, integer_x);
	
	char low[1] = {'L'};
	
	status = CPXchgbds (env, prob1, 1, &obj1_index, low, &x_orig[obj1_index]);
	if ( status ){
		 printf("CPXchgbds, Failed to change bounds, error code %d\n", status);
		 exit(0);}
	
	/*************** Build prob 2 ***************/
	
	CPXLPptr prob2 = CPXcloneprob (env, prob1, &status);
/*  	for(i=0;i<cur_numcols;i++) indices[i] = i;*/
  	
  	chg_coefs(env,prob1,indices,-.0000001);
  	chg_coefs(env,prob2,indices,-1000000000.);
	
	/*************** Solve LP1 ***************/
		
	status = CPXlpopt (env, prob1);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		exit(0);
  	}
  	
  	int lpstat = 0;//CPXgetstat (env, prob1);
  	
  	if(lpstat == 12)
  	{
  		double upper_limit,lower_limit,what_it_should_be;
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_UpperObj, &upper_limit);
  		status = CPXgetdblparam (env, CPXPARAM_Simplex_Limits_LowerObj, &lower_limit);
  		what_it_should_be = pow(10.,75.);
  		if(upper_limit < what_it_should_be || lower_limit > -what_it_should_be)
  		{
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_UpperObj, what_it_should_be);
  			status = CPXsetdblparam ( (CPXENVptr) env, CPXPARAM_Simplex_Limits_LowerObj, -what_it_should_be);
/*  			printf("Warning: Within PSA_full encountered an issue where objective value limits were set to strange values. This was NOT a user error. Resetting these values to default and restarting PSA_full\n");*/
  			PSA_full_left(env,NULL,x_orig,NULL,NULL);
  			exit(0);
  		}
  	}
  	else if(lpstat == 3) exit(0);
  	
  	double prev_obj_vals[2] = {0.,0.};
  	double new_obj_vals[2] = {0.,0.};
  	double obj_vals[2] = {0.,0.};
  	double coef_lb = 0., coef_ub = 0.;
/*  	printf("the status of the solve: %d\n",lpstat);*/

	status = CPXgetx (env, prob1, prev_obj_vals, obj1_index, obj2_index);
	
/*	printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
	
	status = CPXlpopt (env, prob2);
 	if ( status ) {
    		printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
    		exit(0);
  	}
  	
  	status = CPXgetx (env, prob2, obj_vals, obj1_index, obj2_index);
  	
  	if( fabs(obj_vals[0] - prev_obj_vals[0]) > .0001 && fabs(obj_vals[1] - prev_obj_vals[1]) > .0001 )
  	{
/*  		printf("obj1 and obj2 solution are different. Keep going\n");*/
/*  		printf("plot([%lf],[%lf],'go');\n",obj_vals[0],obj_vals[1]);*/
  	}
  	else
  	{
/*  		printf("obj1 and obj2 solutions are the same. Stop\n");*/
/*		printf("plot([%lf],[%lf],'go');\n",prev_obj_vals[0],prev_obj_vals[1]);*/
  		insert(1,prev_obj_vals[0],prev_obj_vals[1],0.,0.,0.,&tree,NULL);
  		exit(0);
  	}
	
	int same_cnt = 0;
	double prev_coef_ub = 100000000000000000., add_val = .001;
	while( fabs(prev_obj_vals[0] - obj_vals[0]) > .01 )
	{
		iter_cnt++;
		if(iter_cnt > 5000)
		{
			printf("%lf, %lf\n",prev_obj_vals[0], obj_vals[0]);
			printf("same count: %d\n",same_cnt);
/*			printf("it_cnt: %d, (%d)\n",iter_cnt,__LINE__);*/
			exit(0);
		}
/*		printf("prev_obj_vals[0]: %lf, obj_vals[0]: %lf\n",prev_obj_vals[0],obj_vals[0]);*/
		status = CPXobjsa(env,prob1,obj2_index,obj2_index,&coef_lb,&coef_ub);
/*		printf("bounds on coefficient: %lf to %lf\n",coef_lb,coef_ub);*/
		if(coef_lb < -10000000000.) break;
		if(fabs(prev_coef_ub - coef_ub) < .00000001)
		{
			coef_lb -= add_val;
			add_val = add_val*2;
/*			retval = 0;*/
/*			exit(0);*/
		}
		prev_coef_ub = coef_ub;
		chg_coefs(env,prob1,indices, -1./(coef_lb - .0000001));
		status = CPXlpopt (env, prob1);
		if ( status ) {
			printf ("%s(%d): CPXprimopt, Failed to solve relaxation,  error code %d\n", __FILE__, __LINE__, status);
			exit(0);
		}
		lpstat = CPXgetstat (env, prob1);
/*	  	printf("the status of the solve: %d\n",lpstat);*/
		status = CPXgetx (env, prob1, new_obj_vals, obj1_index, obj2_index);
/*		printf("the next point might be: (%lf,%lf)\n",new_obj_vals[0],new_obj_vals[1]);*/
		
		if(new_obj_vals[0] > prev_obj_vals[0])
		{
			printf("plot([%lf,%lf],[%lf,%lf],'-mo');\n",prev_obj_vals[0],new_obj_vals[0],prev_obj_vals[1],new_obj_vals[1]);
/*			if(check_for_stopping_PSA_full)*/
/*			{*/
/*				int insert_check = mock_insert(2,new_obj_vals[0],new_obj_vals[1],prev_obj_vals[0],prev_obj_vals[1],*/
/*					(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree);*/
/*				if(!insert_check) */
/*				{*/
/*					printf("Stopping PSA full left early\n");*/
/*					exit(0);*/
/*				}*/
/*			}*/
			insert(2,new_obj_vals[0],new_obj_vals[1],prev_obj_vals[0],prev_obj_vals[1],
				(prev_obj_vals[1]-new_obj_vals[1])/(prev_obj_vals[0]-new_obj_vals[0]),&tree,NULL);
			prev_obj_vals[0] = new_obj_vals[0];
			prev_obj_vals[1] = new_obj_vals[1];
			same_cnt = 0;
		}
		else same_cnt++;
		if( fabs(prev_obj_vals[0] - obj_vals[0]) < .1 && same_cnt >=5 ) break;
	}
 	
/* 	TERMINATE:*/
 	
 	free_and_null ((char **) &b);
	free_and_null ((char **) &integer_x);
/*	free_and_null ((char **) &indices);*/
	CPXfreeprob(env,&prob1);
	CPXfreeprob(env,&prob2);
}

void PSA_full(CPXCENVptr env, CPXLPptr prob, double *x_orig, int *basis_col_info_orig, int *basis_row_info_orig)
{
/*	printf("plot([%lf],[%lf],'ko');\n",x_orig[obj1_index],x_orig[obj2_index]);*/
/*	printf("calling PSA full right\n");*/
	PSA_full_right(env, prob, x_orig, basis_col_info_orig, basis_row_info_orig);
/*	printf("calling PSA full left\n");*/
/*	PSA_full_left(env, prob, x_orig, basis_col_info_orig, basis_row_info_orig);*/
}
