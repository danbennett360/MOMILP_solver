/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#ifdef CPLEX
    #include "cplex.h"
#else
    #include <glpk.h>
#endif

#include "point_class.h"
#include "sydneys_class.h"
#include "simplex_class.h"
#include "problem_class.h"
#include "multiobjective_solver.h"

bool DEBUG = false;
bool SCAN_FOR_REPEATS = false;
bool SAVE_POINTS = true;
bool SCAN_FOR_NEGATIVE_NORMAL = false;
double EPSILON = .000000001;

int main(int argc, char **argv)
{
	/* initialize Cplex environment *********************************/
	int status = 0;
	string phrase = "";
	string phrase2 = "";
	string phraseCopy = "";
	string filename1 = "", filename2 = "";
	string temp = "";
	MultiobjectiveProblem myProb;
	bool minProb = true;
	#ifdef CPLEX
	    CPXENVptr  env = NULL;
        CPXLPptr   lp = NULL;
    #else 
        glp_prob *lp;
    #endif
    ifstream fin;
    ofstream fout;
    bool done = false;
	
	#ifdef CPLEX
      	env = CPXopenCPLEX(&status);

      	if(env==NULL)
      	{
        		char errmsg[1024];
        		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
        		CPXgeterrorstring (env, status, errmsg);
        		printf ("%s", errmsg);
        		exit(0);
      	}
  	#endif
  	
	/******************************************************************/

  	/************* Set to 1 thread **********************************/
  	
  	#ifdef CPLEX
	    status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	    if ( status ) {
		    printf ("Failure to set threads to 1, error %d.\n",status);
		    exit(0);
	    }
	#endif
  	
  	/******************************************************************/
  	
  	/************* Get Number of Objectives ***************************/
  	
	myProb.SetNumObj(atoi(argv[1]));
  	
  	/******************************************************************/
  	
  	
  	/************* Set any desired CPLEX parameters here **************/
  	
  	#ifdef CPLEX
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
	#else
	    glp_term_out(GLP_OFF); // terminal output is on, but I will want it off (on is default)
	#endif
	
//	status = CPXsetincumbentcallbackfunc(env, userincumbent, NULL);
  	
  	/******************************************************************/
  	
  	/** Create test Problems and read them in from the model files ****/
    
    phrase = argv[2];
    phrase = phrase.substr(phrase.find("."));
    filename1 = argv[2];
    filename2 = "temp0.mop";
    
    if(phrase[1] == 'm') // assume we are working with a single *.mops file
    {
        #ifndef CPLEX
/*        cout << "modifying file to not contain OBJSENSE" << endl;*/
            fin.open(filename1);
	        fout.open(filename2);
            while (getline(fin, phrase))
            {
                phrase2 = regex_replace(phrase, regex("^ +"), "");
                if((phrase2[0] == 'O' || phrase2[0] == 'o') && (phrase2[1] == 'B' || phrase2[1] == 'b') && (phrase2[2] == 'J' || phrase2[2] == 'j'))
                {
                    getline(fin, phrase);
                    phrase2 = regex_replace(phrase, regex("^ +"), "");
                    if((phrase2[0] == 'M' || phrase2[0] == 'm') && (phrase2[1] == 'A' || phrase2[1] == 'a') && (phrase2[2] == 'X' || phrase2[2] == 'x')) minProb = false;
                }
                else
                {
                    fout << phrase << endl;
                }
            }
            fin.close();
            fout.close();
            filename1 = filename2;
            filename2[4]++;
        #endif
        
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
/*                    cout << phrase2 << endl;*/
                    if(!done && ((phrase2[0] == 'N' || phrase2[0] == 'n') && (phrase2[1] == ' ' || phrase2[1] == '\t'))) 
                    {
/*                        cout << phrase2 << endl;*/
                        temp = phrase;
                        done = true;
                    }
                    else
                    {
                        if(temp.length() && (phrase2[0] != 'N' && phrase2[0] != 'n' && phrase2[0] != 'L' && phrase2[0] != 'l' && phrase2[0] != 'G' && phrase2[0] != 'g' && phrase2[0] != 'E' && phrase2[0] != 'e'))
                        {
/*                            cout << "here" << endl;*/
                            fout << temp << endl;
                            temp = "";
                        }
/*                        cout << phrase << endl;*/
                        fout << phrase << endl;
                    }
                }
                fin.close();
                fout.close();
                filename1 = filename2;
                filename2[4]++;
	        }
	        #ifdef CPLEX
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
	        #else
	            lp = glp_create_prob();
	            status = glp_read_mps(lp, GLP_MPS_FILE, NULL, filename1.c_str());
	            if ( status ) 
	          	{
	            		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	            		exit(0);
	            }
	            if(minProb) glp_set_obj_dir(lp, GLP_MIN);
	            else glp_set_obj_dir(lp, GLP_MAX);
	        #endif
	        
	        myProb.AddLP(lp);

        }
        filename2[4] = '0';
        
        #ifdef CPLEX
            for(int i = 0; i < myProb.GetNumObj() - 1; i++)
            {
                remove(filename2.c_str());
                filename2[4]++;
            }
        #else
            for(int i = 0; i < myProb.GetNumObj(); i++)
            {
                remove(filename2.c_str());
                filename2[4]++;
            }
        #endif
    }
    else // assume we are working with multiple *.lp files
    {
	    for(int i = 0; i < myProb.GetNumObj(); i++)
	    {
	        #ifdef CPLEX
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
	        #else
	            lp = glp_create_prob();
	        
	            status = glp_read_lp(lp, NULL, argv[i+2]);
	            if ( status ) 
	          	{
	            		printf ("Could not read input %d, exiting. Error code: %d\n", i+1, status);
	            		exit(1);
	            }
	        #endif
	        
	        myProb.AddLP(lp);
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
			#ifdef CPLEX
			    status = CPXwriteprob (env, myProb.GetLP(i), phrase.c_str(), "LP");
			#else 
			    status = glp_write_lp(myProb.GetLP(i), NULL, phrase.c_str());
			    cout << status << endl;
			#endif
		}
	}

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
    
    if(DEBUG) 
    {
        #ifdef CPLEX
            status = CPXwriteprob (env, myProb.GetMainLP(), "overall_prob.lp", "LP");
        #endif
    }
  	
  	if(DEBUG) cout << "finding initial simplices" << endl;
    vector<Simplex> t = myProb.DichotomicSearch();
  				
  	/******************************************************************/
    

  	return 0;
}

