
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <vector>
#include <string>

#include "cplex.h"

using namespace std;

/*******************************************************************/
/*              CPLEX VARIABLES                                    */
/*******************************************************************/

CPXENVptr  env=NULL;
CPXLPptr   lp = NULL;
CPXLPptr   lp1 = NULL;
vector<CPXLPptr>  lps;
int numObj = 0;
int status = 0;
int cur_numcols;
int cur_numrows;

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
  	
  	/******************************************************************/
  	
  	numObj = atoi(argv[1]);
  	
  	/** Create test Problems *******************/

    int k = 0;
/*  	for(int k = 0; k < numObj; k++)*/
/*	{*/
	  	lp = CPXcreateprob (env,&status,argv[k+2]);
	  	if(lp==NULL) 
	  	{
	    		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", k+1, status);
	    		exit(0);
	    }
	    	
	    status = CPXreadcopyprob(env,lp,argv[k+2],NULL);

	  	if ( status ) 
	  	{
	    		printf ("Could not read input %d, exiting. Error code: %d\n", k+1, status);
	    		exit(0);
	    }
	    
	    status = CPXwriteprob (env, lp, "myprob.lp", "LP");
	    
	    exit(0);
	    
/*	    lps.push_back(lp);*/
/*    }*/
	/******************************************************************/
  
  	/******************************************************************/
  	// Make sure the problems are in maximization form with less than
  	// or equal to and/or equality constraints.
  	/******************************************************************/

    string s = "lp0.mps";
    for(int k = 0; k < numObj; k++)
    {   
        s[2]++; 
    	lp1 = lps[0];
    	lps.erase (lps.begin());
	    cur_numcols = CPXgetnumcols (env, lp1);
	    cur_numrows = CPXgetnumrows (env, lp1);
	
	    int *indices = (int*)calloc(cur_numcols*2, sizeof(int));
	
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
    				}
    			}
    			status = CPXgetrhs (env, lp1, rhs, i, i);
    			rhs[0] = -rhs[0];
    			status = CPXchgrhs (env, lp1, 1, &i, rhs);
    			status = CPXchgsense (env, lp1, 1, &i, new_sense);
    		}
    		if(!equalities && row_sense[i] == 'E')
    		{
    			equalities = 1;
    		}
    	}
    	for(i = 0;i < cur_numcols*2; i++) indices[i] = i;

      	double *obj_coef1 = (double *) malloc (cur_numcols*sizeof(double));

     	status = CPXgetobj(env, lp1, obj_coef1, 0, cur_numcols -1);

      	if ( status ) 
      	{
        		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP1, error code %d\n", status);
        		exit(0);
        	}
    	
    	status = CPXgetobjsen (env, lp1);
    	
    	if ( status == 0) 
  	    {
    		printf ("When trying to get objective sense, found that no problem object exists, error code %d\n", status);
    		exit(0);
    	}
    	else if( status == -1)	// Problem is maximization type - switch to minimization
    	{	
    		for(i=0;i<cur_numcols;i++)
    		{
    			obj_coef1[i] = -obj_coef1[i];
    		}
    		status = CPXchgobj (env, lp1, cur_numcols, indices, obj_coef1);
    		if ( status ) 
	  	{
	    		printf ("Error changing objective 1, error code %d\n", status);
	    		exit(0);
	    	}
    		status = CPXchgobjsen (env, lp1, CPX_MIN);
    		if ( status ) 
	  	{
	    		printf ("Error changing sense of objective 1, error code %d\n", status);
	    		exit(0);
	    	}
    	}
    	
    	lps.push_back(lp1);
    	free(indices);
    	free(row_sense);
    	free(obj_coef1);
/*    	status = CPXwriteprob (env, lp1, &s[0], "LP");*/
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
    	
  	
  	/****************************************************
    		 Add two new variables to keep track of
    		 the objective function values.
    	*****************************************************/
	
	double lb[2] = {-CPX_INFBOUND,-CPX_INFBOUND};
	double ub[2] = {CPX_INFBOUND,CPX_INFBOUND};

	double x1_[2] = {0.,0.};
	double x2_[2] = {0.,0.};
	
	char sense[4] = {'G'};
	
	lp1 = lps[0];
	vector< char* > v;
	char* t;
	
	s = "OBJ2";
	for(int k = 0; k < numObj - 1; k++)
	{
	     t = (char *)s.c_str();
	     v.push_back(t);
	     
	     status = CPXnewrows (env, lp1, 1, NULL, sense, NULL, &v[0]); 
	    if ( status ) 
      	{
        		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
        		exit(0);
        	}
	     
	     s[3]++;
	}
	
    	
    	int *rowlist;
    	int *collist;
    	double *vallist;
    	
    	rowlist = (int *) malloc (cur_numcols*sizeof(int));
    	collist = (int *) malloc (cur_numcols*sizeof(int));
    	vallist = (double *) malloc (cur_numcols*sizeof(double));
    	
    	/******************************* Set coefficients for new rows *****************/
    	
    	for(int k = 1; k < numObj; k++)
    	{
    	    double *obj_coef1 = (double *) malloc (cur_numcols*sizeof(double));

         	status = CPXgetobj(env, lps[k], obj_coef1, 0, cur_numcols -1);

          	if ( status ) 
          	{
            		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP1, error code %d\n", status);
            		exit(0);
            	}
            	
	        for(int i=0;i<cur_numcols;i++)
	        {
		        rowlist[i] = cur_numrows + k - 1;
		        collist[i] = i;
		        vallist[i] = obj_coef1[i];
	        }
	
	        status = CPXchgcoeflist (env, lp1, cur_numcols, rowlist, collist, vallist);
	        if ( status ) 
        {
		        printf ("CPXgchgcoeflist, Failed to set coefficients for first new row, error code %d\n", status);
		        exit(0);
	        }
    	}
    	
    	status = CPXwriteprob (env, lp1, "multiobj_prob_MIP.mps", "MPS");
    	
    	status = CPXchgprobtype(env, lp1, CPXPROB_LP);
      	if ( status ) 
      	{
        		printf ("Failed to change the problem from a MIP to an LP.\n");
        		exit(0);
      	}
      	
      	status = CPXwriteprob (env, lp1, "multiobj_prob_LP.mps", "MPS");
    	
    	/*************************************************************************************/
    	/******************************* Set coefficients for second new row *****************/
    	
/*    	for(i=0;i<cur_numcols;i++)*/
/*    	{*/
/*    		rowlist[i] = cur_numrows - 1;*/
/*    		collist[i] = i;*/
/*    		if(i < cur_numcols - 2) vallist[i] = obj_coef2[i];*/
/*    		else if(i == cur_numcols - 1) vallist[i] = -1;*/
/*    		else vallist[i] = 0;*/
/*    	}*/
/*    	*/
/*    	status = ( CPXchgcoeflist (env, lp1, cur_numcols, rowlist, collist, vallist) || CPXchgcoeflist (env, lp2, cur_numcols, rowlist, collist, vallist));*/
/*    	if ( status ) */
/*  	{*/
/*    		printf ("CPXgchgcoeflist, Failed to set coefficients for second new row, error code %d\n", status);*/
/*    		exit(0);*/
/*    	}*/
/*    	*/
/*    	status = (CPXchgrhs (env, lp1, 1, &rowlist[0], &rh2) || CPXchgrhs (env, lp2, 1, &rowlist[0], &rh2));*/
    	
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
