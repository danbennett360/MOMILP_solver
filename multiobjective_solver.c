/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#include "cplex.h"
#include "point_class.h"
#include "sydneys_class.h"
#include "simplex_class.h"
#include "problem_class.h"
#include "multiobjective_solver.h"

bool DEBUG = false;
bool SCAN_FOR_REPEATS = false;
bool SAVE_POINTS = true;

int main(int argc, char **argv)
{
	/* initialize Cplex environment *********************************/
	int status = 0;
	string phrase = "";
	MultiobjectiveProblem myProb;
	CPXENVptr  env = NULL;
    CPXLPptr   lp = NULL;
	
  	env = CPXopenCPLEX(&status);

  	if(env==NULL)
  	{
    		char errmsg[1024];
    		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
    		CPXgeterrorstring (env, status, errmsg);
    		printf ("%s", errmsg);
    		exit(0);
  	}
  	
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
	
	myProb.SetEnv(env);
	
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
    	
	/******************************************************************/
	
	/******************************************************************/
	// Make sure the problems are in minimization form with less than
  	// or equal to and/or equality constraints.
  	/******************************************************************/

    myProb.SetNumRowsAndCols();
    myProb.ConvertLPs();
 
    if(DEBUG)
	{
		phrase = "myprob0.lp";
		for(int i = 0; i < myProb.GetNumObj(); i++) 
		{
			phrase[6]++;
			status = CPXwriteprob (env, myProb.GetLP(i), phrase.c_str(), "LP");
		}
	}

     /******************************************************************/
  
  	/******************************************************************/
  	// Take in any other command line flags and set appropriate
  	// parameter values
  	/******************************************************************/

    if(DEBUG) cout << "setting param vals" << endl;
    myProb.SetParamVals(argc, argv);

    /******************************************************************/

     /****************************************************
    		 Add new variables to keep track of
    		 the objective function values.
    	*****************************************************/
	
	
	if(myProb.StoreObjectivesInMainProb()) 
	{
	    if(DEBUG) cout << "adding rows for objectives" << endl;
	    myProb.AddRowsForObjectives();
	}
  	
  	/********* Initialization *****************************************/
  	
    //Implement some sort of a heuristic here that generates starting solutions.
    
/*    Simplex mySimplex(myProb.GetNumObj());*/
  				
  	/******************************************************************/

    /********* TEMPORARY CODE *****************************************/
    
    if(DEBUG) status = CPXwriteprob (env, myProb.GetMainLP(), "overall_prob.lp", "LP");
  	
  	if(DEBUG) cout << "finding initial simplices" << endl;
    vector<Simplex> t = myProb.DichotomicSearch();
  				
  	/******************************************************************/
    

  	return 0;
}

