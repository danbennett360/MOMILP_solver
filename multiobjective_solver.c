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
	string phrase2 = "";
	string filename1 = "", filename2 = "";
	string temp = "";
	MultiobjectiveProblem myProb;
	CPXENVptr  env = NULL;
    CPXLPptr   lp = NULL;
    ifstream fin;
    ofstream fout;
    bool done = false;
	
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
    
    phrase = argv[2];
    phrase = phrase.substr(phrase.find("."));
    filename1 = argv[2];
    filename2 = "temp0.mop";
    
    if(phrase[1] == 'm') // assume we are working with a single *.mops file
    {
        for(int i = 0; i < myProb.GetNumObj(); i++)
	    {
	        if(i != 0) 
	        {
	            done = false;
	            fin.open(filename1);
	            fout.open(filename2);
	            while (getline(fin, phrase))
                {
                    phrase2 = regex_replace(phrase, regex("^ +"), "");
                    if(!done && ((phrase2[0] == 'N' || phrase2[0] == 'n') && (phrase2[1] == ' ' || phrase2[1] == '\t'))) 
                    {
/*                        cout << phrase2 << endl;*/
                        temp = phrase;
                        done = true;
                    }
                    else
                    {
                        if(temp.length() && (phrase2[0] != 'N' && phrase2[0] != 'n' && phrase2[0] != 'L' && phrase2[0] != 'l' && phrase2[0] != 'G' && phrase2[0] != 'g'))
                        {
                            fout << temp << endl;
                            temp = "";
                        }
                        fout << phrase << endl;
                    }
                }
                fin.close();
                fout.close();
                filename1 = filename2;
                filename2[4]++;
	        }
	      	lp = CPXcreateprob (env,&status,filename1.c_str());
	      	if(lp==NULL) 
	      	{
	        		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", i+1, status);
	        		exit(0);
	        }
	        	
	        status = CPXreadcopyprob(env,lp,filename1.c_str(),NULL);
	      	if ( status ) 
	      	{
	        		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	        		exit(0);
	        }

		    myProb.AddLP(lp);
        }
        filename2[4] = '0';
        for(int i = 0; i < myProb.GetNumObj() - 1; i++)
        {
            remove(filename2.c_str());
            filename2[4]++;
        }
    }
    else // assume we are working with multiple *.lp files
    {
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

