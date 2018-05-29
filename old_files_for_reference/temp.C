/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#include "cplex.h"
#include "problem_class.h"
#include "multiobjective_solver.h"

int main(int argc, char **argv)
{
	/* initialize Cplex environment *********************************/
	bool debug = true;
	int status = 0;
	string phrase = "";
	MultiobjectiveProblem myProb;
	
  	env = CPXopenCPLEX(&status);

  	if(env==NULL)
  	{
    		char errmsg[1024];
    		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
    		CPXgeterrorstring (env, status, errmsg);
    		printf ("%s", errmsg);
    		exit(0);
  	}
  	
  	myProb.SetEnv(env);
  	
    	/******************************************************************/

  	/************* Set to 1 thread **********************************/
  	
	status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	if ( status ) {
		printf ("Failure to set threads to 1, error %d.\n",status);
		exit(0);
	}
  	
  	/******************************************************************/
  	
  	/************* Get Number of Objectives ***************************/
  	
	myProb.SetNumObj(atoi(argv[1]));
  	
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
	
//	status = CPXsetincumbentcallbackfunc(env, userincumbent, NULL);
  	
  	/******************************************************************/
  	
  	/** Create test Problems and read them in from the model files ****/

	for(int i = 0; i < myProb.GetNumObj(); i++)
	{
	  	lp = CPXcreateprob (env,&status,argv[i+2]);
	  	if(lp==NULL) 
	  	{
	    		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", i+1, status);
	    		exit(0);
	    	}
	    	
	    	status = CPXreadcopyprob(env,lp,argv[i+2],NULL);

	  	if ( status ) 
	  	{
	    		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	    		exit(0);
	    	}

		myProb.AddLP(lp);
    	}
    	
    	if(debug)
    	{
    		phrase = "myprob1.lp";
    		for(int i = 0; i < myProb.GetNumObj(); i++) 
    		{
    			phrase[6] = i+1;
    			status = CPXwriteprob (env, myProb.GetLP(i), phrase.c_str(), "LP");
    		}
    	}
	/******************************************************************/
  
  	/******************************************************************/
  	// Make sure the problems are in minimization form with less than
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
    	for(i = 0; i < cur_numcols*2; i++) indices[i] = i;

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
//    	
//    	
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
//  	
//  	char *xctype = (char *) malloc ((cur_numcols+2)*sizeof(char));
//  	
//  	status = CPXgetctype (env, lp1, xctype, 0, cur_numcols-1);
//    	if ( status ) 
//  	{
//    		printf ("(%d) CPXgetctype, Failed to get variable types, error code %d\n", __LINE__,status);
//    		exit(0);
//    	}
//    	
//    	provide_xctype(env, xctype,cur_numcols);
//  	
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
//  	
//  	/***************************************************************************************/
//  	/* 
//  	 _  _  __  ____      ____   __   ____  ____  ____    ____   __   __    _  _  ____  ____
//  	( \/ )(  )(  _ \ ___(  _ \ / _\ / ___)(  __)(    \  / ___) /  \ (  )  / )( \(  __)(  _ \      
//	/ \/ \ )(  ) __/(___)) _ (/    \\___ \ ) _)  ) D (  \___ \(  O )/ (_/\\ \/ / ) _)  )   /      
//	\_)(_/(__)(__)      (____/\_/\_/(____/(____)(____/  (____/ \__/ \____/ \__/ (____)(__\_) 
//	  
//	*/
//	/***************************************************************************************/
///*  	if(parallel == 0) */
//	box *new_box1 = (struct box*) malloc( sizeof( struct box ) );
//	new_box1->se_x = SE_extreme_x;
//	new_box1->se_y = SE_extreme_y;
//	new_box1->nw_x = NW_extreme_x;
//	new_box1->nw_y = NW_extreme_y;
//	  			
//  	mip_solve_sequential (env, new_box1, NULL);
//  	
//  	free(new_box1);
///*              		      nPoints,*/
///*              		      &nNadir);*/
///*        else */
///*        mip_solve_parallel   (env,*/
///*              		      lp1,*/
///*              		      lp2,*/
///*              		      obj_coef1,*/
///*              		      obj_coef2);*/
///*              		      nPoints,*/
///*              		      &nNadir);*/

//	printf("IT FINISHED!!\n\n");
//	finish_time = clock();
//	double duration = (double)(finish_time - start_time) / CLOCKS_PER_SEC;
//	   	printf("Total time: %lf\n",duration);
//	printf("Iterations: %d\n",recursion_count);
//	printf("Num MIPs solved: %d\n\n\n",num_mips_solved);
//	int n = get_num_nodes(tree);
//	printf("Number of nodes in tree: %d\n",n);
//	print_inorder(tree,1);

//        /**************************************************************************************/
//  	
//  	/**********************Code for TERMINATE*************************/
///*  	TERMINATE:*/

//  	/* Free the problem as allocated by CPXcreateprob and
//     	CPXreadcopyprob, if necessary *************************************/

//  	if ( lp1 != NULL || lp2 != NULL ) 
//  	{
//    		status = CPXfreeprob (env, &lp1) || CPXfreeprob (env, &lp2);
//    		if ( status ) 
//    		{
//      			printf ("CPXfreeprob failed, error code %d.\n",status);
//    		}
//  	}
//  	/******************************************************************/

//  	/* Free the CPLEX environment, if necessary ***********************/

//  	if(env!=NULL) status=CPXcloseCPLEX(&env);
//  	/******************************************************************/

//	/******************* Print Solutions to File **********************/
//	
///*  	print_preorder (tree, output);*/
//  	/******************************************************************/

//	/***************** Close Files ************************************/

///*  	if (output) fclose(output);*/
///*  	if (init_nadir) fclose(init_nadir);*/
//  	/*  fclose(all_inserted);*/
//  	/*  fclose(testing_struct);*/
///*  	if (init_sol) fclose(init_sol);*/
//  	/*  close_files();*/
//  	/******************************************************************/
//  	
//  	/******************* Free Variables *******************************/

///*  	free_and_null ((char **) &nadirX_p);*/
///*  	free_and_null ((char **) &ParetoX_p);*/
///*  	//free_and_null ((char **) &Theta_k_ind);*/
//  	destroy_tree (tree);
///*  	free_and_null ((char **) &ctype);*/
//  	free_and_null ((char **) &obj_coef1);
//	free_and_null ((char **) &obj_coef2);
//	free_and_null ((char **) &row_sense);
//	free_and_null ((char **) &xctype);
//	free_and_null ((char **) &rowlist);
//    	free_and_null ((char **) &collist);
//    	free_and_null ((char **) &vallist);
//    	free_and_null ((char **) &integer_indices);
//    	free_and_null ((char **) &indices);
//    	free_and_null ((char **) &weighted_coefs);
///* 	free_and_null ((char **) &obj_coef3);*/
///* 	 free_and_null ((char **) &theta);*/
///*  	free_and_null ((char **) &indices_new);*/
///*  	free_and_null ((char **) &indices_still_to_check);*/
///*  	free_and_null ((char **) &indices_f_n_l);*/
///*	delete_split_pts(first_split_pt);*/
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
