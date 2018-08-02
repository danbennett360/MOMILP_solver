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

#include "problem_class.h"
#include "point_class.h"

void TellAboutCurrent(vector<Simplex *> active) {

   return;
   size_t i;
   Simplex * adj;

   cout << endl;
   cout << "There are " << active.size() << " active simplices" << endl;
   for(auto x: active) {
      cout <<"Simplex " << x->ID() << endl;

      // print the pointsd
      i = 0;
      cout << "\tPoints:" << endl;
      for (auto & point: x->GetExtremePoints()) {
         cout << "\t\t" << i <<":  (";
	 for(auto cord: point) {
	    cout << cord << " ";
	 }
	 cout <<")" << endl;
	 i++;
      }

      cout << "\tNormal:" << endl;
      cout << "\t\t Has " << x->GetNormal().size() << " components" << endl;
      cout << "\t\t<";
      for(auto point: x->GetNormal()) {
          cout << point << " ";
      }
      cout << ")" << endl;
      cout << endl;

      // print the adjacent
      cout << "\tAdjacent: ";
      for(i = 0; i < size_t(x->GetDimension()); i++) {
         adj = x->Adjacent(i); 
	 if (adj != nullptr) {
	    cout << adj->ID() << " ";
	 } else {
	    cout << " X ";
	 }
      }
      cout << endl;
      cout << endl;
   }
   cout << endl;
   return;
}


#ifdef CPLEX
    void MultiobjectiveProblem::SetEnv(CPXENVptr e)
    {
	    env = e;
	    return;
    }

    CPXENVptr MultiobjectiveProblem::GetEnv()
    {
	    return env;
    }
#endif

void MultiobjectiveProblem::SetNumObj(int a)
{
	numObjectives = a;
	return;
}

#ifdef CPLEX
    void MultiobjectiveProblem::AddLP(CPXLPptr lp)
    {
	    lps.push_back(lp);
	    return;
    }
#else
    void MultiobjectiveProblem::AddLP(glp_prob *lp)
    {
	    lps.push_back(lp);
	    return;
    }
#endif

int MultiobjectiveProblem::GetNumObj()
{
	return numObjectives;
}

#ifdef CPLEX
    CPXLPptr MultiobjectiveProblem::GetLP(int i)
    {
	    return lps[i];
    }
#else
    glp_prob *MultiobjectiveProblem::GetLP(int i)
    {
	    return lps[i];
    }
#endif

void MultiobjectiveProblem::SetNumRowsAndCols()
{
    #ifdef CPLEX
         numCols = CPXgetnumcols (env, lps[0]);
         numRows = CPXgetnumrows (env, lps[0]);
    #else
         numCols = glp_get_num_cols(lps[0]);
         numRows = glp_get_num_rows(lps[0]);
    #endif
}

void MultiobjectiveProblem::ConvertLPs()
{
	char *row_sense;
	row_sense = new char[numRows];
	int status = 0;
	double coef = 0.;
	double *objCoefs;
	objCoefs = new double[numCols];
	int *colIndices;
	colIndices = new int[numCols];
	double *rowCoefs;
	rowCoefs = new double[numCols];
	double rhs[1] = {0.};
	char new_sense[1] = {'L'};
	int k = 0;
	vector<double> temp(numCols);
	
	#ifdef CPLEX
	    status = CPXgetsense (env, lps[k], row_sense, 0, numRows-1);
	#else
	    for(int i = 0; i < numRows; i++)
	    {
	        status = glp_get_row_type(lps[k], i+1);
	        if(status == GLP_LO) row_sense[i] = 'G';
	        else if(status == GLP_UP) row_sense[i] = 'L';
	        else if(status == GLP_FX) row_sense[i] = 'E';
	        else
	        {
	            cout << "There is a row of type " << status << ". Write code to deal with this!\n";
	            exit(1);
	        }
	    }
	#endif
	
	for(int i = 0; i < numRows; i++)
	{
	    #ifndef CPLEX
	        status = glp_get_mat_row(lps[k], i+1, colIndices, rowCoefs);
	    #endif
		if(row_sense[i] == 'G')
		{
		    #ifdef CPLEX
			    for(int j=0; j < numCols; j++)
			    {
				    status = CPXgetcoef (env, lps[k], i, j, &coef);
				    if(coef != 0.)
				    {
					        status = CPXchgcoef (env, lps[k], i, j, -coef);
				    }
			    }
			    status = CPXgetrhs (env, lps[k], rhs, i, i);
			    rhs[0] = -rhs[0];
			    status = CPXchgrhs (env, lps[k], 1, &i, rhs);
			    status = CPXchgsense (env, lps[k], 1, &i, new_sense);
			#else
			    for(int j = 0; j < status; j++)
			    {
			        rowCoefs[j] = -rowCoefs[j];
			    }
			    rhs[0] = glp_get_row_lb(lps[k], i+1);
			    rhs[0] = -rhs[0];
			    glp_set_mat_row(lps[k], i+1, status, colIndices, rowCoefs);
			    //glp_set_row_bnds(lps[k], i+1, GLP_UP, NULL, rhs[0]);
			    glp_set_row_bnds(lps[k], i+1, GLP_UP, 0, rhs[0]);
			#endif
		}
		if(!containsEqualities && row_sense[i] == 'E')
		{
			containsEqualities = true;
		}
	}
	for(int i = 0; i < numCols*2; i++) indices.push_back(i);

    for(int i = 0; i < numObjectives; i++)
    {
        #ifdef CPLEX
           	status = CPXgetobj(env, lps[i], objCoefs, 0, numCols -1);
           	if ( status ) 
           	{
             		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP1, error code %d\n", status);
             		exit(0);
         	}
        #else
/*            glp_write_lp(lps[i], NULL, "prob.lp");*/
/*            exit(0);*/
/*            cout << glp_version() << endl;*/
/*            exit(0);*/
/*            cout << "------------------------\n" << glp_get_num_cols(lps[i]) << glp_get_obj_coef(lps[i], 1) << endl;*/
            for(int j = 0; j < numCols; j++)
            {
                objCoefs[j] = glp_get_obj_coef(lps[i], j+1);
            }
     	#endif
     	
     	#ifdef CPLEX
     	    status = CPXgetobjsen (env, lps[i]);
     	#else
     	    status = glp_get_obj_dir(lps[i]);
     	    if(status == GLP_MIN) status = 1;
     	    else if(status == GLP_MAX) status = -1;
     	    else status = 0;
     	#endif
     	if ( status == 0 ) 
       	{
		printf ("When trying to get objective sense, found that no problem object exists, error code %d\n", status);
       		exit(0);
     	}
     	else if( status == -1 )	// Problem is maximization type - switch to minimization
     	{	
     		for(int j=0; j < numCols; j++)
     		{
     			objCoefs[j] = -objCoefs[j];
     		}
     		
     		#ifdef CPLEX
         		status = CPXchgobj (env, lps[i], numCols, &indices[0], objCoefs);
         		if ( status ) 
           	    {
             		printf ("Error changing objective, error code %d\n", status);
             		exit(0);
             	}
         		status = CPXchgobjsen (env, lps[i], CPX_MIN);
         		if ( status ) 
           	    {
             		printf ("Error changing sense of objective 1, error code %d\n", status);
             		exit(0);
             	}
            #else
                for(int j = 0; j < numCols; j++)
                {
                    glp_set_obj_coef(lps[i], j+1, objCoefs[j]);
                }
                glp_set_obj_dir(lps[i], GLP_MIN);
         	#endif
         	
         	memcpy(&temp[0], &objCoefs[0], numCols*sizeof(double));
         	objectiveCoefs.push_back(temp);
     	}
     	else
     	{
     	    memcpy(&temp[0], &objCoefs[0], numCols*sizeof(double));
     	    objectiveCoefs.push_back(temp); 
     	}
	}
	
	#ifdef CPLEX
	    mainProb = CPXcloneprob (env, lps[k], &status);
      	if ( status ) 
      	{
        		printf ("Failed to clone problem 1.\n");
        		exit(0);
      	}
    #else
        mainProb = glp_create_prob();
        glp_copy_prob(mainProb, lps[k], GLP_ON);
  	#endif
  	
/*  	status = CPXwriteprob (env, mainProb, "mainprob.lp", "LP");*/
  	
  	delete[] row_sense; 
  	delete[] objCoefs;
  	delete[] rowCoefs;
  	delete[] colIndices;
    	
	return;
}

void MultiobjectiveProblem::AddRowsForObjectives()
{
     vector<double> lb(numObjectives, -infinity);
     vector<double> ub(numObjectives, infinity);
     vector<char> sense(numObjectives, 'E');
     vector<int> rowlist(numCols + numObjectives + 1, numRows - 2);
	 vector<int> collist(numCols + numObjectives + 1, 0);
	 vector<double> vallist(numCols + numObjectives + 1, 0.);  	
     int status = 0;
     double dZero = 0.;
     
     for(int i = 0; i < numObjectives; i++)
     {
          objectiveColIndices.push_back(numCols + i);
          objectiveRowIndices.push_back(numRows + i);
/*          cout << "i: " << i << endl;*/
     }

    #ifdef CPLEX
	    status = CPXnewcols (env, mainProb, numObjectives, NULL, &lb[0], &ub[0], NULL, NULL); 
	    if ( status ) 
      	{
        		printf ("CPXnewcols, Failed to add additional variables to the model, error code %d\n", status);
        		exit(0);
        }
    #else
        status = glp_add_cols(mainProb, numObjectives);
        for(int i = status; i < status + numObjectives; i++) 
        {
            //glp_set_col_bnds(mainProb, i, GLP_FR, NULL, NULL);
            glp_set_col_bnds(mainProb, i, GLP_FR, 0, 0);
        }
    #endif
	
	numCols += numObjectives;
	
	#ifdef CPLEX
	    status = CPXnewrows (env, mainProb, numObjectives, NULL, &sense[0], NULL, NULL) ; 
	    if ( status ) 
      	{
        		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
        		exit(0);
        }
    #else
        status = glp_add_rows(mainProb, numObjectives);
        //for(int i = status; i < status + numObjectives; i++) glp_set_row_bnds(mainProb, i, GLP_FX, 0., NULL);
        for(int i = status; i < status + numObjectives; i++) glp_set_row_bnds(mainProb, i, GLP_FX, 0., 0);
    #endif
    	
	numRows += numObjectives;
	
	/******************************* Set coefficients for new rows *****************/
	
	if(DEBUG) cout << "Setting coefficients for new rows" << endl;
	
	for(int j = 0; j < numObjectives; j++)
	{
/*	    cout << "j: " << j << endl;*/
     	for(int i = 0; i < numCols; i++)
     	{
/*     	    cout << "i: " << i << endl;*/
     	    rowlist[i] = objectiveRowIndices[j];
     	    
/*     	    cout << objectiveColIndices[j] << endl;*/
            #ifdef CPLEX
     		    collist[i] = i;
     		#else 
     		    collist[i+1] = i+1;
     		#endif
     		
     		if(i < objectiveColIndices[j] && i < numCols - numObjectives)
     		{
/*     		     cout << vallist[i] << endl;*/
/*     		     cout << objectiveCoefs[j][i] << endl;*/
                 #ifdef CPLEX
     		        vallist[i] = objectiveCoefs[j][i];
     		     #else
     		        vallist[i+1] = objectiveCoefs[j][i];
     		     #endif
     		}
     		else if(i == objectiveColIndices[j]) 
     		{
     		    #ifdef CPLEX
     		        vallist[i] = -1;
     		    #else
     		        vallist[i+1] = -1;
     		    #endif
     		}
     		else 
     		{
     		    #ifdef CPLEX
     		        vallist[i] = 0;
     		    #else
     		        vallist[i+1] = 0;
     		    #endif
     		}
/*     		#ifndef CPLEX*/
/*     		    vallist[objectiveColIndices[j]+1] = objectiveCoefs[j][[objectiveColIndices[j]];*/
/*     		    vallist[objectiveColIndices[j]+2] = -1;*/
     	}
     	
     	#ifdef CPLEX
         	status = CPXchgcoeflist (env, mainProb, numCols, &rowlist[0], &collist[0], &vallist[0]);
         	if ( status ) 
       	    {
         		printf ("CPXgchgcoeflist, Failed to set coefficients for first new row, error code %d\n", status);
         		exit(0);
         	}
         	
         	status = CPXchgrhs (env, mainProb, 1, &rowlist[0], &dZero);
         	if ( status ) 
       	    {
         		printf ("CPXchgrhs, Failed to change rhs, error code %d\n", status);
         		exit(0);
         	}
        #else
/*            cout << "------------------------\n" << glp_get_num_cols(mainProb) << "\t" << glp_get_num_rows(mainProb) << endl;*/
            glp_set_mat_row(mainProb, rowlist[0]+1, collist.size() - 1, &collist[0], &vallist[0]);
     	#endif
	}
/*	status = glp_write_lp(mainProb, NULL, "prob.lp");*/
/*	exit(0);*/
/*	status = CPXwriteprob (env, mainProb, "mainprob.lp", "LP");*/
	
	/*************************************************************************************/

    return;
}

void MultiobjectiveProblem::SetParamVals(int argc, char **argv)
{
     // bennett 7/4/18
     // Change needed because the number of input files is no longer fixed.
     // We always need to have the flag start with a - and consume all arguments.
     //int i = numObjectives + 2;
     int i = 2;
     int limit = min (argc, numObjectives+2);
     int debugLevel = 0;
     while (i < limit and argv[i][0] != '-') {
         i++;
     }

     while(i < argc)
     {
/*        cout << argv[i] << endl;*/
        if(argv[i][0] != '-')
        {
            cout << "Invalid Command Line Argument (missing '-'). Exiting." << endl;
	        cout << "Got " << argv[i] << endl;
            exit(1);
        }
        else
        {
            if(!strcmp(argv[i],"-store"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    storeObjectivesInMainProb = true;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    storeObjectivesInMainProb = false;
                    objectiveColIndices.push_back(numCols-1);
                }
                else
                {
                    cout << "Invalid value for flag '-store', ignoring. Valid values are 'T' and 'F'." << endl;
                }
                i++;
            }
	    // bennett 7/18
            else if(!strcmp(argv[i],"-epsilon")) {
                i++;
		        if (argv[i] != nullptr) {
		            try {
		                epsilon = stod(argv[i]);
		                EPSILON = epsilon;
		                i++;
		            } catch(const invalid_argument & ia) {
		           
                        cout << "Invalid value for flag '-epsilon'" <<  argv[i] 
			            << " , ignoring. Valid values are positive doubles." << endl;
		            }
		        } else {
                            cout << "Invalid value for flag '-epsilon'" 
		                 << ", ignoring. Valid values are positive doubles." << endl;
		        }
            }
/*            else if(!strcmp(argv[i],"-normalize"))*/
/*            {*/
/*                i++;*/
/*                if(toupper(argv[i][0]) == 'T')*/
/*                {*/
/*                    normalizeObjectiveMultipliers = true;*/
/*                }*/
/*                else if(toupper(argv[i][0]) == 'F')*/
/*                {*/
/*                    normalizeObjectiveMultipliers = false;*/
/*                }*/
/*                else*/
/*                {*/
/*                    cout << "Invalid value for flag '-normalize', ignoring. Valid values are 'T' and 'F'." << endl;*/
/*                }*/
/*                i++;*/
/*            }*/
            else if(!strcmp(argv[i],"-lexopt"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    useLexicographicOptimization = true;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    useLexicographicOptimization = false;
                }
                else
                {
                    cout << "Invalid value for flag '-lexopt', ignoring. Valid values are 'T' and 'F'." << endl;
                }
                i++;
            }
            else if(!strcmp(argv[i],"-interior"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    interiorPoint = true;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    interiorPoint = false;
                }
                else
                {
                    cout << "Invalid value for flag '-interior', ignoring. Valid values are 'T' and 'F'." << endl;
                }
                i++;
            }
            else if(!strcmp(argv[i],"-showprogress"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    showProgress = true;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    showProgress = false;
                }
                else
                {
                    cout << "Invalid value for flag '-showprogress', ignoring. Valid values are 'T' and 'F'." << endl;
                }
                i++;
            }
            else if(!strcmp(argv[i],"-progressvalue"))
            {
                i++;
                showProgressIterations = atoi(argv[i]);
		        i++;
            }
            else if(!strcmp(argv[i],"-maxiter"))
            {
                i++;
                maxIterations = atoi(argv[i]);
                i++;
            }
            else if(!strcmp(argv[i],"-debug"))
            {
                i++;
                debugLevel = atoi(argv[i]);
                if(debugLevel < 0 || debugLevel > 2)
                {
                    cout << "Invalid value for flag '-debug', ignoring. Valid values are '0', '1' and '2'." << endl;
                }
                else
                {
                    if(debugLevel >= 1) 
                    {
                        SCAN_FOR_REPEATS = true;
                        SCAN_FOR_NEGATIVE_NORMAL = true;
                    }
                    if(debugLevel >= 2) DEBUG = true;
                }
                i++;
/*                cout << "Value of debug flag: " << debugLevel << "\tValue of DEBUG: " << DEBUG << endl;*/
            }
            else if(!strcmp(argv[i],"-reldist"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    relativeDistance = true;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    relativeDistance = false;
                }
                else
                {
                    cout << "Invalid value for flag '-reldist', ignoring. Valid values are 'T' and 'F'." << endl;
                }
                i++;
            }
            else
            {
                cout << "Invalid command line flag. Exiting!\nFor a list of valid flags, see 'HELP.txt'" << endl;
                exit(0);
            }
        }
     }
}

bool MultiobjectiveProblem::StoreObjectivesInMainProb()
{
    return storeObjectivesInMainProb;
}

void MultiobjectiveProblem::DichotomicSearch(SimplexStore & store)
{
    int status = 0;
    SimplexStore retVal;
    
    #ifdef CPLEX
        tempProb = CPXcloneprob (env, mainProb, &status);
        /*    cout << status << endl;*/
        if ( status ) 
      	{
       		printf ("Failed to clone main problem.\n");
       		exit(0);
      	}
      	
      	status = CPXchgprobtype(env, tempProb, CPXPROB_LP);
      	if ( status ) 
      	{
       		printf ("Failed to change the problem from a MIP to an LP.\n");
       		exit(0);
      	}
    #else
        tempProb = glp_create_prob();
        glp_copy_prob(tempProb, mainProb, GLP_ON);
    #endif
  	
    /*  status = CPXwriteprob (env, tempProb, "tempprob.lp", "LP");*/
  	
    MeatOfDichotomicSearch(store);
    return;
}

// compute primal points and save in the point stack.
// This responsiblitiy should probably be moved to the point class.
Point MultiobjectiveProblem::MakePoint(vector<double> point) {

    Point rv;
    double * x;
    int status = 0;
    x = new double[numCols];

    #ifdef CPLEX
        status = CPXgetx (env, tempProb, x, 0, numCols-1);
        if ( status ) {
                printf ("%s(%d): CPXgetx, Failed to get x values,  error code %d\n", __FILE__, __LINE__, status);
                exit(1);
        }
    #else
        if(!interiorPoint) for(int j = 0; j < numCols; j++) x[j] = glp_get_col_prim(tempProb, j+1);
        else for(int j = 0; j < numCols; j++) x[j] = glp_ipt_col_prim(tempProb, j+1);
    #endif
    if(storeObjectivesInMainProb) {
         rv = Point(point, x,  numCols-numObjectives);
    }
    else {
        rv = Point(point, x, numCols);
    }

    delete[] x;
    return rv;
}

// Get initial simplex using lexicographic minimization
void MultiobjectiveProblem::AddFirstSimplex(SimplexStore & store, vector<vector<double> > & extremes, vector<Simplex *> & simplexStack) {

    vector<double> mins(numObjectives, infinity);
    vector<double> maxs(numObjectives, -infinity);
    vector<int> simplexPoints;
    vector<double> point(numObjectives, 0.);
    Simplex * temp=nullptr;

    for(int i = 0; i < numObjectives; i++)
    {
        extremes.push_back(LexicographicMinimization(i));
        if(i != 0) CheckForDomination(extremes, epsilon);
    }
    
    for(int i = 0; i < numObjectives; i++)
    {
        for(int j = 0; j < numObjectives; j++)
        {
            mins[i] = min(mins[i], extremes[j][i]);
            maxs[i] = max(maxs[i], extremes[j][i]);
        }
    }
   
    int location;
    for(int i = 0; i < numObjectives; i++)
    {
        for(int j = 0; j < numObjectives; j++)
        {
            if(i == j) point[j] = infinity; 
            else point[j] = mins[j];
        }
	Point p = MakePoint(point);
	location = store.AddPoint(p);
	simplexPoints.push_back(location);
    }

    temp = store.AddSimplex(simplexPoints, NormalizeObjectiveMultipliers());
   
    simplexStack.push_back(temp);

    // why scan for repeats, we only have one simplex present?
    if(SCAN_FOR_REPEATS)
    {
        scanForRepeats(simplexStack);
    }
    
    if(DEBUG)
    {
        cout << "**************************\nInitial Simplex: " << endl;
        for(unsigned int i = 0; i < simplexStack.size(); i++) simplexStack[i]->WriteOctaveCodeToPlotSimplex(true);
        cout << "**************************" << endl;
    }
    return;
}


void MultiobjectiveProblem::DoInitialSplits(SimplexStore & simplexStore, vector<vector<double> > & extremes, vector<Simplex *> & simplexStack) {

    Simplex * currentSimplex;
    vector<Simplex * > simplicesToSplit;
    size_t k;

    for(int i = 0; i < numObjectives; i++) {
        if(DEBUG)
        {
            cout << "Adding point: ";
            for(unsigned int j = 0; j < extremes[i].size(); j++) cout << extremes[i][j] << "\t";
            cout << endl;
        }
    
        if(i == 0)
        {
            currentSimplex = simplexStack.back();
            simplexStack.pop_back();
        }
        else
        {
            k = 0;
            while( (!relativeDistance && simplexStack[k]->MultiplyPointsByNormal(extremes[i]) >= simplexStack[k]->PlaneVal() - epsilon) ||
                   (relativeDistance && 
                        (simplexStack[k]->MultiplyPointsByNormal(extremes[i]) - simplexStack[k]->PlaneVal())/abs(simplexStack[k]->PlaneVal()) >= - epsilon) )
            {
	        // changed the order of k++ and test so the loop doesn't bomb when k = size.
                k++;
                if(k >= simplexStack.size())
                {
                    cout << "Error, accessing elements past the end of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                    exit(0);
                }
            }
            currentSimplex = simplexStack[k];
            simplexStack.erase(simplexStack.begin() + k);
        }

	// at this point, currentSimplex is the simplex we want to split and 
	//                extremes[i] is the point we want to split it on.
        
        simplicesToSplit.push_back(currentSimplex);
	currentSimplex->MarkForSplit();

        CheckIfAdjacentsAreShadowed(simplexStack, simplicesToSplit, currentSimplex, extremes[i], numObjectives, relativeDistance, 0.000000001, simplexStore);
	//epsilon*epsilon);
        
        if(DEBUG) 
        {
            cout << "Size of simplicesToSplit: " << simplicesToSplit.size() << endl;
            cout << "The simplices that will be split are: " << endl;
        }
       
        DoSplit(simplexStack, simplexStore, simplicesToSplit, extremes[i]);

	TellAboutCurrent(simplexStack);
       
        size_t l = simplexStack.size();

        // should not need this any longer.
        if(simplicesToSplit.size() > 1) deleteRepeats(simplexStack, l);
       
        simplicesToSplit.resize(0);
        
        if(DEBUG)
        {
            cout << "**************************\nStored simplices after adding " << i+1 << "th initial extreme point: " << endl;
            for(unsigned int j = 0; j < simplexStack.size(); j++) simplexStack[j]->WriteOctaveCodeToPlotSimplex(true);
            cout << "**************************" << endl;
        }
        
        if(SCAN_FOR_NEGATIVE_NORMAL)
        {
            scanForNegativeNormal(simplexStack);
        }
    }
    return;
}

void MultiobjectiveProblem::MeatOfDichotomicSearch(SimplexStore & simplexStore)
{
    int status = 0;
    vector<Simplex * > simplexStack ;
    vector< vector<double> > extremes;
    vector<Simplex *> simplicesToSplit;
    Simplex * currentSimplex;

    vector<double> point(numObjectives, 0.);
    vector<double> point2(numObjectives, 0.);
    int iterator = 0;
    vector<int> removeTheseIndicesFromSimplexStack;
    long int k = 0;
    long int  l = 0;
    int numSaved = 0;
    long maxIter = maxIterations;
    
/*    bool DEBUG = true;*/
    
    // Get initial simplex using lexicographic minimization
    AddFirstSimplex(simplexStore, extremes, simplexStack);
    TellAboutCurrent(simplexStack);

    DoInitialSplits(simplexStore, extremes, simplexStack);
   
    cout << "Finished adding initial simplices, entering body of the simplex generation process." << endl;
    
    k = max(1, int(simplexStack.size())-1);

    while(k > 0 && iterator < maxIter)
    {
/*        if(iterator == 5080) DEBUG = true;*/
        iterator++;
        
        if(showProgress && iterator % showProgressIterations == 0)
        {
            for(unsigned int j = 0; j < simplexStack.size(); j++) 
            {
                if(simplexStack[j]->GetSaveForSolution()) numSaved++;
/*                simplexStack[j]->WriteOctaveCodeToPlotSimplex(true);*/
            }
            cout << "Summary after iteration " << iterator << ":" << endl;
            cout << "Number of Simplices in List: " << simplexStack.size() << "\tNumber Saved: " << numSaved << "\tNumber left to search: " << simplexStack.size() - numSaved << endl;
            cout << "------------------------------------------------------\n";
            numSaved = 0;
        }
        
        if(DEBUG)
        {
            cout << "-----------------------------------------------------------------------------------" << endl;
            cout << "Iteration: " << iterator << "\tSize of stack: " << simplexStack.size() <<endl;
            cout << "**************************\nSimplices in stack: " << endl;
            for(unsigned int i = 0; i < simplexStack.size(); i++) simplexStack[i]->WriteOctaveCodeToPlotSimplex(true);
            cout << "\%**************************\n\%Simplices stored in return vector:" << endl;
            //for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex(true);
            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
        }
        
        k = simplexStack.size()-1;
        while(k >= 0 and simplexStack[k]->GetSaveForSolution()) k--;

	//cout << "k = " << k <<  " and the stack size is " << simplexStack.size() << endl;

	//DEBUG = true;
        if(k > 0)
        {
            currentSimplex = simplexStack[k];
            simplexStack.erase(simplexStack.begin() + k);
    /*        simplexStack.pop_back();*/
            
    /*        currentSimplex->PrintData();*/
            
    /*        status = CPXwriteprob (env, tempProb, "prob_before.lp", "LP");*/
    /*        exit(0);*/
        
            ChangeTempObjCoefs(currentSimplex->GetNormal());
	    //cout << "Here i am  1" << endl;
            
    /*        status = CPXwriteprob (env, tempProb, "prob.lp", "LP");*/
    /*        exit(0);*/
            
            #ifdef CPLEX
                status = CPXlpopt (env, tempProb);
             	if ( status ) {
                		printf ("%s(%d): CPXlpopt, Failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
                		exit(0);
              	}
            #else
                if(interiorPoint) status = glp_interior(tempProb, NULL);
                else status = glp_simplex(tempProb, NULL);
                if ( status ) {
                		printf ("%s(%d): GLPK failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
                		exit(1);
              	}
            #endif
          	
            if(DEBUG)
            {
              	point = currentSimplex->GetNormal();
              	cout << "new current normal is: <";
              	for(int i = 0; i < numObjectives; i++) cout << point[i] << "\t";
              	cout << ">" << endl;
            }
          	
          	point = GetObjectiveValues(tempProb);
          	
          	if(DEBUG)
          	{
              	    cout << "new point is: (";
              	    for(int i = 0; i < numObjectives; i++) cout << point[i] << "\t";
              	    cout << ")" << endl;
              	    cout << currentSimplex->MultiplyPointsByNormal(point) << "\t" << currentSimplex->PlaneVal() << endl;
                    currentSimplex->WriteOctaveCodeToPlotSimplex(true);
              	    cout << "Distance to simplex: " << currentSimplex->MultiplyPointsByNormal(point) - (currentSimplex->PlaneVal() - epsilon) << "\nIn Front of An Adjacent: " << PointIsInFrontOfAnAdjacent(simplexStack, currentSimplex, point, epsilon) << endl;
          	}
          	
          	if( (!relativeDistance && currentSimplex->MultiplyPointsByNormal(point) < currentSimplex->PlaneVal() - epsilon) ||
         	   (relativeDistance && (currentSimplex->MultiplyPointsByNormal(point) - currentSimplex->PlaneVal())/abs(currentSimplex->PlaneVal()) < - epsilon) )
          	{
                
          	if(DEBUG) {
          	        /*cout << "Distance to simplex: " << currentSimplex->MultiplyPointsByNormal(point) - (currentSimplex->PlaneVal() - epsilon) 
                              << "\nIn Front of An Adjacent: " << PointIsInFrontOfAnAdjacent(simplexStack, currentSimplex, point, epsilon) << endl;*/
          	        cout << "Adding three new simplices to stack" << endl;
          	    }
          	    
          	simplicesToSplit.push_back(currentSimplex);
		currentSimplex->MarkForSplit();
                CheckIfAdjacentsAreShadowed(simplexStack, simplicesToSplit, currentSimplex, point, numObjectives, relativeDistance, 0.000000001, simplexStore);
		//epsilon*epsilon);
               
                if(DEBUG) cout << "Size of simplicesToSplit: " << simplicesToSplit.size() << endl;
                if(DEBUG) cout << "Size of simplexStack: " << simplexStack.size() << endl;
                if(DEBUG) 
                {
                    cout << "-----------------------\nThe simplices to be split are: " << endl;
                    for(unsigned int j = 0; j < simplicesToSplit.size(); j++) 
                    {
                        simplicesToSplit[j]->WriteOctaveCodeToPlotSimplex(false);
                    }
                    cout << "-----------------------\n";
                }
                
                DoSplit(simplexStack, simplexStore, simplicesToSplit, point);
	        TellAboutCurrent(simplexStack);
                
                if(DEBUG) cout << "Size of simplexStack after adding: " << simplexStack.size() << endl;
                
                if(simplicesToSplit.size() > 1) deleteRepeats(simplexStack, l);
                
                if(DEBUG) cout << "Size of simplexStack after deleting: " << simplexStack.size() << endl;
                
                simplicesToSplit.resize(0);
          	    
          	}
          	else
          	{
    /*      	    retVec.push_back(currentSimplex);*/
                    currentSimplex->SaveForSolution();
                    simplexStack.push_back(currentSimplex);
          	    
          	    if(DEBUG) cout << "No new point found" << endl;
          	}
      	}
	//DEBUG = false;
      	
      	if(SCAN_FOR_NEGATIVE_NORMAL)
        {
            scanForNegativeNormal(simplexStack);
        }
    }

    if(iterator == maxIter)
    {
        cout << "Max iterations (" << maxIter << ") reached in meat of dichotomic search. Consider increasing value. Exiting!" << endl;
        
        cout << "**************************\nSimplices still on stack: " << endl;
        for(unsigned int j = 0; j < simplexStack.size(); j++) 
        {
            if(simplexStack[j]->GetSaveForSolution()) numSaved++;
            simplexStack[j]->WriteOctaveCodeToPlotSimplex(true);
        }
        cout << "\%**************************" << endl;
        
        cout << "Number of Simplices in List: " << simplexStack.size() << "\tNumber Saved: " << numSaved << "\tNumber left to search: " << simplexStack.size() - numSaved << endl;
        
/*        cout << "\%**************************\n\%Simplices discovered so far: " << endl;*/
/*        for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex(true);*/
/*        cout << "**************************" << endl;*/
        
        exit(0);
    }
    
/*    if(DEBUG) */
/*    {*/
        cout << "**************************\nSimplices after dichotomic search (which took " << iterator << " iterations): " << endl;
        for(unsigned int i = 0; i < simplexStack.size(); i++) simplexStack[i]->WriteOctaveCodeToPlotSimplex(true);
        cout << "**************************\nTotal number in list is: " << simplexStack.size() << endl;
        
        cout << "*************************************\nThe individual extreme points discovered are:\n";
        WritePoints(simplexStack); 
/*    }*/
    
/*    PostProcessDichotomicSearch(retVec);*/
/*    */
/*    if(DEBUG) */
/*    {*/
/*     */
/*        for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex(true);*/
/*        cout << "**************************" << endl;*/
/*    }*/

    // bennett 6/18
    if(SAVE_POINTS) {
       DumpPoints(simplexStore);
    }
    
    return;
}

// this will eventually go away.
void MultiobjectiveProblem::DumpPoints(SimplexStore & store) {
    vector<string> varNames;

    #ifdef CPLEX
       varNames = GetVarNames(env, tempProb, numCols);
    #else
       varNames = GetVarNames(tempProb, numCols);
    #endif

    int i;
    bool extreme = false;
    Point tmpPoint;
    bool append = false;

    for(i=0;i<store.PointCount(); i++) {
        tmpPoint = store.GetPoint(i);
	extreme = false;
	for(auto x: tmpPoint.Data()) {
            if (x >= infinity) {
	       extreme = true;
	    }
	}

	if (not extreme) {
            tmpPoint.WritePointToFile("points.txt", varNames, append);
	    append = true;
	}
    }
    return;
}

void MultiobjectiveProblem::ChangeTempObjCoefs(int i)
{
    #ifdef CPLEX
        int status = CPXchgobj (env, tempProb, objectiveColIndices[0], &indices[0], &objectiveCoefs[i][0]);
	    if ( status ) {
		    printf ("Failed to change obj coef. Error code %d\n", status);
	    }
	#else
	    if(storeObjectivesInMainProb) for(int j = 0; j < numCols - numObjectives; j++) glp_set_obj_coef(tempProb, j+1, objectiveCoefs[i][j]);
	    else for(int j = 0; j < numCols; j++) glp_set_obj_coef(tempProb, j+1, objectiveCoefs[i][j]);
	#endif
	
/*	#ifdef CPLEX*/
/*        status = CPXwriteprob (env, tempProb, "temp_prob_cplex.lp", "LP");*/
/*    #else*/
/*        int status = glp_write_lp(tempProb, NULL, "temp_prob_glpk.lp");*/
/*    #endif*/
/*    */
/*    cout << i << endl;*/
/*    exit(0);*/
	
	return;
}

void MultiobjectiveProblem::ChangeTempObjCoefs(const vector<double> & v)
{
    int status = 0;
    
/*    cout << weightedCoefs.size() << "\t" << objectiveColIndices[0] << endl;*/
    if(int(weightedCoefs.size()) < objectiveColIndices[0]) weightedCoefs.resize(objectiveColIndices[0] + 1);
    
    for(int i = 0; i < objectiveColIndices[0] + 1; i++)
    {
        weightedCoefs[i] = 0.;
        for(int j = 0; j < numObjectives; j++)
        {
            weightedCoefs[i] += objectiveCoefs[j][i]*v[j];
        }
/*        cout << weightedCoefs[i] << endl;*/
    }
    
/*    status = CPXwriteprob (env, tempProb, "prob_before.lp", "LP");*/
   
    #ifdef CPLEX
        status = CPXchgobj (env, tempProb, objectiveColIndices[0], &indices[0], &weightedCoefs[0]);
	    if ( status ) {
		    printf ("Failed to change obj coef. Error code %d\n", status);
	    }
	#else
	    for(int i = 0; i < objectiveColIndices[0] + 1; i++) glp_set_obj_coef(tempProb, indices[i] + 1, weightedCoefs[i]);
	#endif
	
/*	status = CPXwriteprob (env, tempProb, "prob_after.lp", "LP");*/
/*	 exit(0);*/
	return;
}

#ifdef CPLEX
    vector<double> MultiobjectiveProblem::GetObjectiveValues(const CPXLPptr & lp)
    {
        vector<double> retVec(numObjectives);
        double *temp;
        temp = new double[numObjectives];
        double *x;
        
        int status = 0;
        
        if(StoreObjectivesInMainProb())
        {
            status = CPXgetx (env, lp, temp, objectiveColIndices[0], objectiveColIndices[numObjectives-1]);
            if ( status ) {
	            printf ("Failed to get obj values. Error code %d\n", status);
            }
        }
        else
        {
            x = new double[numCols];
            status = CPXgetx (env, lp, x, 0, numCols-1);
            if ( status ) {
	            printf ("Failed to get x. Error code %d\n", status);
            }
	        
	        for(int i = 0; i < numObjectives; i++)
	        {
	            temp[i] = 0.;
	            for(int j = 0; j < numCols; j++)
	            {
	                temp[i] += x[j]*objectiveCoefs[i][j];
	            }
	        }
	        delete[] x;
        }
        
        memcpy(&retVec[0], &temp[0], numObjectives*sizeof(double));
        
        delete[] temp;
        return retVec;
    }
#else
    vector<double> MultiobjectiveProblem::GetObjectiveValues(glp_prob *lp)
    {
        vector<double> retVec(numObjectives);
        double *temp;
        temp = new double[numObjectives];
        double *x;
        
        int status = 0;
        
        if(StoreObjectivesInMainProb())
        {
            if(interiorPoint) for(int i = 0; i < numObjectives; i++) temp[i] = glp_ipt_col_prim(lp, objectiveColIndices[0] + i + 1);
            else for(int i = 0; i < numObjectives; i++) temp[i] = glp_get_col_prim(lp, objectiveColIndices[0] + i + 1);
        }
        else
        {
            x = new double[numCols];
            if(interiorPoint) for(int i = 0; i < numCols; i++) x[i] = glp_ipt_col_prim(lp, i + 1);
            else for(int i = 0; i < numCols; i++) x[i] = glp_get_col_prim(lp, i + 1);
	        
	        for(int i = 0; i < numObjectives; i++)
	        {
	            temp[i] = 0.;
	            for(int j = 0; j < numCols; j++)
	            {
	                temp[i] += x[j]*objectiveCoefs[i][j];
	            }
	        }
	        delete[] x;
        }
        
        memcpy(&retVec[0], &temp[0], numObjectives*sizeof(double));
        
        delete[] temp;
        return retVec;
    }
#endif

bool MultiobjectiveProblem::NormalizeObjectiveMultipliers()
{
    return normalizeObjectiveMultipliers;
}

vector<double> MultiobjectiveProblem::LexicographicMinimization(int i)
{
    string s = "lp1.lp";
    int status = 0;
    vector<double> point(numObjectives, 0.);
    vector<double> v(numObjectives, 1.);
    char u = 'U';
    char sense = 'L';
    double val = 0.;
    
    ChangeTempObjCoefs(i);
    
    #ifdef CPLEX
        status = CPXlpopt (env, tempProb);
     	if ( status ) {
        	printf ("%s(%d): CPXlpopt, Failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
        	exit(1);
      	}
    #else  	
       if(interiorPoint) status = glp_interior(tempProb, NULL);
       else status = glp_simplex(tempProb, NULL);
       if ( status ) {
        	printf ("%s(%d): GLPK failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
        	exit(1);
      	}
  	#endif
  	
  	point = GetObjectiveValues(tempProb);
  	
/*  	for(int k = 0; k < numObjectives; k++) cout << point[k] << "\t" << endl;*/
  	
  	//add bounds to set max value of objective i
  	if(UseLexicographicOptimization())
  	{
      	if(StoreObjectivesInMainProb())
      	{
      	    #ifdef CPLEX
      	        status = CPXchgbds (env, tempProb, 1, &objectiveColIndices[i], &u, &point[i]);
      	    #else
      	        val = glp_get_col_lb(tempProb, objectiveColIndices[i] + 1);
      	        glp_set_col_bnds(tempProb, objectiveColIndices[i] + 1, GLP_DB, val, point[i]);
      	    #endif
      	}
      	else
      	{
      	    vector<int> rowIndices(objectiveColIndices[0],numRows);
      	    
      	    #ifdef CPLEX
              	status = CPXnewrows (env, tempProb, 1, &point[i], &sense, NULL, NULL) ; 
                if ( status ) 
              	{
                		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
                		exit(0);
                }
                
                status = CPXchgcoeflist (env, tempProb, objectiveColIndices[0], &rowIndices[0], &indices[0], &objectiveCoefs[i][0]);
             	if ( status ) 
                {
             		printf ("CPXgchgcoeflist, Failed to set coefficients for first new row, error code %d\n", status);
             		exit(0);
             	}
            #else
                status = glp_add_rows(tempProb, 1);
                //glp_set_row_bnds(tempProb, status, GLP_UP, NULL, point[i]);
                glp_set_row_bnds(tempProb, status, GLP_UP, 0, point[i]);
                indices.insert(indices.begin(), 0);
                glp_set_mat_row(tempProb, status, rowIndices.size(), &indices[0], &objectiveCoefs[i][0]);
                indices.erase(indices.begin());
         	#endif
        }
      	
      	for(int j = 0; j < numObjectives; j++)
      	{
            //optimize over additional objectives
            if(j != i)
            {
                ChangeTempObjCoefs(j);
                   
                #ifdef CPLEX     
                    status = CPXlpopt (env, tempProb);
                 	if ( status ) {
                    		printf ("%s(%d): CPXlpopt, Failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
                    		exit(0);
                  	}
                #else
                    if(interiorPoint) status = glp_interior(tempProb, NULL);
                   else status = glp_simplex(tempProb, NULL);
                   if ( status ) {
                    		printf ("%s(%d): GLPK failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
                    		exit(1);
                  	}
              	#endif
              	
              	point = GetObjectiveValues(tempProb);
              	
              	if(StoreObjectivesInMainProb())
              	{
              	    #ifdef CPLEX
              	        status = CPXchgbds (env, tempProb, 1, &objectiveColIndices[j], &u, &point[j]);
              	    #else
              	        val = glp_get_col_lb(tempProb, objectiveColIndices[j] + 1);
              	        glp_set_col_bnds(tempProb, objectiveColIndices[j] + 1, GLP_DB, val, point[j]);
              	    #endif
              	}
              	else
              	{
              	    vector<int> rowIndices(objectiveColIndices[0],numRows+j+1);
              	    
              	    #ifdef CPLEX
                      	status = CPXnewrows (env, tempProb, 1, &point[j], &sense, NULL, NULL) ; 
                        if ( status ) 
                      	{
                        		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
                        		exit(0);
                        }
                        
                        status = CPXchgcoeflist (env, tempProb, objectiveColIndices[0], &rowIndices[0], &indices[0], &objectiveCoefs[j][0]);
                     	if ( status ) 
                        {
                     		printf ("CPXgchgcoeflist, Failed to set coefficients for first new row, error code %d\n", status);
                     		exit(0);
                     	}
                    #else
                        status = glp_add_rows(tempProb, 1);
                        //glp_set_row_bnds(tempProb, status, GLP_UP, NULL, point[j]);
                        glp_set_row_bnds(tempProb, status, GLP_UP, 0, point[j]);
                        indices.insert(indices.begin(), 0);
                        glp_set_mat_row(tempProb, status, rowIndices.size(), &indices[0], &objectiveCoefs[j][0]);
                        indices.erase(indices.begin());
                 	#endif
                }
            }

      	}
      	
      	if(StoreObjectivesInMainProb())
      	{
      	    vector<char> up(numObjectives, 'U');
      	    vector<double> bd(numObjectives, infinity);
      	    
      	    #ifdef CPLEX
      	        status = CPXchgbds (env, tempProb, numObjectives, &objectiveColIndices[0], &up[0], &bd[0]);
      	    #else
      	        for(int j = 0; j < numObjectives; j++)
      	        {
          	        val = glp_get_col_lb(tempProb, objectiveColIndices[j] + 1);
                  	//glp_set_col_bnds(tempProb, objectiveColIndices[j] + 1, GLP_LO, val, NULL);
                  	glp_set_col_bnds(tempProb, objectiveColIndices[j] + 1, GLP_LO, val, 0);
              	}
      	    #endif
      	}
      	else
      	{
      	    #ifdef CPLEX
      	        status = CPXdelrows (env, tempProb, numRows, numRows+numObjectives);
      	    #else
      	        vector<int> num(numObjectives+1, 0);
      	        for(int j = 0; j < numObjectives; j++) num[j+1] = numRows + j + 1;
      	        glp_del_cols(tempProb, numObjectives, &num[0]);
      	    #endif
        }
    }
  	
  	return point;
}

void  MultiobjectiveProblem::DoSplit(vector<Simplex *> & simplexStack, SimplexStore & simplexStore, vector<Simplex *> & simplicesToSplit, vector<double> & point) {
    
    int pos;
    vector<Simplex *> newList;
    Simplex * tmp;

    // save point 
    Point p =MakePoint(point);
    pos = simplexStore.AddPoint(p);

    /*
     for(auto x: simplicesToSplit) {
        cout << "Splitting " << x->ID() << endl;
     }
     cout << endl;
     */
    // call split
    newList = simplexStore.SplitSimplex(simplicesToSplit, pos, NormalizeObjectiveMultipliers());

    // put results back into simplexStack
    while(newList.size() > 0) {
       tmp = newList.back();
       newList.pop_back();
       simplexStack.push_back(tmp);
    }

    simplicesToSplit.clear();
    return;
}

bool MultiobjectiveProblem::UseLexicographicOptimization()
{
    return useLexicographicOptimization;
}

double GetAngleBetween(const Simplex * s1, const Simplex * s2, bool normalize)
{
    vector<double> v1 = s1->GetNormal();
    double dotProduct = s2->MultiplyPointsByNormal(v1);
    cout << "dot product: " << dotProduct << endl;
    double cosOfAngle = 0.;
    if(!normalize) cosOfAngle = dotProduct/(GetVectorMagnitude(v1)*GetVectorMagnitude(s2->GetNormal()));
    else cosOfAngle = dotProduct;
    cout << "cosine of angle: " << cosOfAngle << endl;
    return acos(cosOfAngle);
}

double GetVectorMagnitude(const vector<double> & v)
{
    double retVal = 0.;
    for(unsigned int i = 0; i < v.size(); i++)
    {
        retVal += v[i]*v[i];
    } 
    cout << "vector magnitude: " << sqrt(retVal) << endl;
    return sqrt(retVal);
}


/*void MultiobjectiveProblem::PostProcessDichotomicSearch(vector<Simplex> & simplices)*/
/*{*/
/*    vector<Simplex> potentials;*/
/*    vector<bool> areDummies;*/
/*    Simplex temp(numObjectives);*/
/*    Simplex temp2(numObjectives);*/
/*    int num = 0, index = 0, simplexIndex = 0, newPointIndex = 0, dim = 0;*/
/*    unsigned int i = 0;*/
/*    bool keepGoing = true;*/
/*    */
/*    while(i < simplices.size())*/
/*    {*/
/*        num = simplices[i].GetNumDummyPoints();*/
/*        if(num >= 1)*/
/*        {*/
/*            if(num == 1) potentials.push_back(simplices[i]);*/
/*            simplices.erase(simplices.begin() + i);*/
/*            i--;*/
/*        }*/
/*        */
/*        i++;*/
/*    }*/
/*    */
/*    cout << "**************************\npotentials: " << endl;*/
/*    for(unsigned int k = 0; k < potentials.size(); k++) potentials[k].WriteOctaveCodeToPlotSimplex(true);*/
/*    cout << "**************************" << endl;*/
/*    */
/*    i = 0;*/
/*    while(i < potentials.size())*/
/*    {*/
/*        cout << "i: " << i << endl;*/
/*        areDummies = potentials[i].ExtremesAreDummies();*/
/*        index = 0;*/
/*        while(!areDummies[index]) */
/*        {*/
/*            cout << "areDummies[" << index << "]: " << areDummies[index] << endl;*/
/*            index++;*/
/*        }*/
/*        */
/*        for(int j = 0; j < numObjectives; j++)*/
/*        {*/
/*            cout << "j: " << j << endl;*/
/*            if(j != index) */
/*            {*/
/*                cout << "j was not index" << endl;*/
/*                dim = 1;*/
/*                keepGoing = true;*/
/*                while(dim > 0 && keepGoing)*/
/*                {*/
/*                    cout << "dim > 0 and keepGoing" << endl;*/
/*                    simplexIndex = i;*/
/*                    newPointIndex = index;*/
/*                    cout << "*************************************************************************" << endl;*/
/*                    temp = potentials[i].FindAdjacentContainingOriginalPoints(potentials, simplexIndex, newPointIndex, j, true);*/
/*                    dim = temp.GetDimension();*/
/*                    cout << "dim: " << dim << endl;*/
/*                    if(temp.GetDimension() > 0)*/
/*                    {*/
/*                        temp2.Reset();*/
/*                        for(int k = 0; k < numObjectives; k++)*/
/*                        {*/
/*                            cout << "k: " << k << "\tindex: " << index << endl;*/
/*                            if(k != index) temp2.AddExtreme(potentials[i].GetExtremePoint(k), NormalizeObjectiveMultipliers());*/
/*                        }*/
/*                        temp2.AddExtreme(temp.GetExtremePoint(newPointIndex), NormalizeObjectiveMultipliers());*/
/*                        temp2.WriteOctaveCodeToPlotSimplex(true);*/
/*                        temp2.PrintData();*/
/*                        exit(0);*/
/*                    }*/
/*                }*/
/*            }*/
/*        }*/
/*        */
/*        i++;*/
/*    }*/
/*    */
/*    return;*/
/*}*/


void CheckIfAdjacentsAreShadowed( vector<Simplex *> & simplexStack, vector<Simplex *> & simplicesToSplit, const Simplex * currentSimplex,
                                  const vector<double> & newPoint, const int & numObjectives, bool relativeDistance, double epsilon, SimplexStore & simplexStore)
{
    Simplex * temp, *temp2;
    vector<Simplex *>::iterator pos;
    
/*    bool DEBUG = true;*/

    if(DEBUG)
    {
         cout << "Checking the following simplex for shadowed adjacents: " << endl;
         currentSimplex->WriteOctaveCodeToPlotSimplex(true);
    }
    
    for(int i = 0; i < numObjectives; i++)
    {
        // get the adjacent at facet i
	temp2 = currentSimplex->Adjacent(i);

	// if it is really there (not null) and the point is in front
        if( temp2 != nullptr &&  not temp2->IsMarkedForSplit() and 
          	    ((!relativeDistance && temp2->MultiplyPointsByNormal(newPoint) <= temp2->PlaneVal() - epsilon) || 
          	     (relativeDistance && (temp2->MultiplyPointsByNormal(newPoint) - temp2->PlaneVal())/abs(temp2->PlaneVal()) <= - epsilon) ) )
        {
            if(DEBUG)
            {
                cout << "The following simplex needs split  ID=" <<  temp2->ID() << endl;
                temp2->WriteOctaveCodeToPlotSimplex(true);
            }

            
	    // add it to the list to split
            simplicesToSplit.push_back(temp2);
	    temp2->MarkForSplit();

	    pos = find(simplexStack.begin(), simplexStack.end(), temp2);
	    if (pos != simplexStack.end()) {
	        simplexStack.erase(pos);
	    }

	    // check it's adjacents.
            CheckIfAdjacentsAreShadowed( simplexStack, simplicesToSplit, temp2, newPoint, numObjectives, relativeDistance, epsilon, simplexStore);
            if(DEBUG)
            {
                 cout << "After returning from a recursion, the current simplex is: " << endl;
                 currentSimplex->WriteOctaveCodeToPlotSimplex(true);
            }
        }
    }
            
    return;
}

#ifdef CPLEX
    CPXLPptr MultiobjectiveProblem::GetMainLP()
    {
        return mainProb;
    }
#else
    glp_prob *MultiobjectiveProblem::GetMainLP()
    {
        return mainProb;
    }
#endif

void scanForRepeats(const vector<Simplex *> & simplexStack)
{
    vector< vector<double> > v;
    int count = 0;
    
    for(unsigned int i = 0; i < simplexStack.size() - 1; i++)
    {
        v.resize(0);
        for(int k = 0; k < simplexStack[0]->GetDimension(); k++)
        {
            v.push_back(simplexStack[simplexStack.size()-1]->GetExtremePoint(k));
        }
    
        for(int j = 0; j < simplexStack[0]->GetDimension(); j++)
        {
            v.push_back(simplexStack[i]->GetExtremePoint(j));
        }
/*        for(unsigned int j = 0; j < v.size(); j++)*/
/*        {*/
/*            for(int k = 0; k < simplexStack[0].GetDimension(); k++)*/
/*            {*/
/*                cout << v[j][k] << "\t";*/
/*            }*/
/*            cout << endl;*/
/*        }*/
/*        cout << endl;*/
        sort(v.begin(), v.end());
/*        for(unsigned int j = 0; j < v.size(); j++)*/
/*        {*/
/*            for(int k = 0; k < simplexStack[0].GetDimension(); k++)*/
/*            {*/
/*                cout << v[j][k] << "\t";*/
/*            }*/
/*            cout << endl;*/
/*        }*/
/*        cout << endl;*/
        count = unique(v.begin(), v.end()) - v.begin();
/*        for(unsigned int j = 0; j < v.size(); j++)*/
/*        {*/
/*            for(int k = 0; k < simplexStack[0].GetDimension(); k++)*/
/*            {*/
/*                cout << v[j][k] << "\t";*/
/*            }*/
/*            cout << endl;*/
/*        }*/
/*        cout << endl;*/
        if(count <= simplexStack[0]->GetDimension())
        {
            cout << "Error, a repeated simplex was just added to the stack! Exiting!\n";
            if(DEBUG)
            {
                cout << "**************************\nSimplices: " << endl;
                for(unsigned int j = 0; j < simplexStack.size(); j++) simplexStack[j]->WriteOctaveCodeToPlotSimplex(true);
                cout << "**************************" << endl;
            }
            exit(0);
        }
    }
}

void scanForNegativeNormal(const vector<Simplex *> & simplexStack)
{
    for(size_t i = 0; i < simplexStack.size(); i++)
    {
        if(!simplexStack[i]->IsPositive())
        {
            cout << "There is a simplex on the stack that is not oriented correctly! Exiting!\n";
            cout << "The culprit is: \n";
            simplexStack[i]->WriteOctaveCodeToPlotSimplex(true);
            exit(1);
        }
    }
}


// bennett 7/18
// changed completely due to adjacency list
bool PointIsInFrontOfAnAdjacent(const vector<Simplex *> & simplexStack, const Simplex * simp, const vector<double> & point, const double & epsilon)
{
    bool retBool = false;
    int dim = simp->GetDimension();
    Simplex * temp;
   
    int i;
    for(i=0;i<dim;i++) {
       temp = simp->Adjacent(i); 
       if (temp != nullptr) {
           if(temp->GetDimension() > 0 && temp->MultiplyPointsByNormal(point) < temp->PlaneVal() - epsilon)
	      return true;
       }
    }
    return false;
}

/*
bool PointIsInFrontOfAnAdjacent(const vector<Simplex> & simplexStack, const Simplex & simp, const vector<double> & point, const double & epsilon)
{
    bool retBool = false;
    int dim = simp.GetDimension();
    Simplex temp(dim);
    int simplexIndex = 0, newPointIndex = 0, i = 0;
    
    while(i < dim && !retBool)
    {
        temp = simp.FindAdjacentContainingOriginalPoints(simplexStack, simplexIndex, newPointIndex, i, true);
        if(temp.GetDimension() > 0 && temp.MultiplyPointsByNormal(point) < temp.PlaneVal() - epsilon) retBool = true;
        i++;
    }
    
    return retBool;
}
*/

void WritePoints(const vector<Simplex *> & simplexStack)
{
    vector< vector<double> > v;
    
    for(unsigned int i = 0; i < simplexStack.size(); i++)
    {
        for(int j = 0; j < simplexStack[0]->GetDimension(); j++) 
        {
            if(!simplexStack[i]->ExtremesAreDummies()[j]) v.push_back(simplexStack[i]->GetExtremePoint(j));
        }
    }
    sort(v.begin(), v.end());
    v.resize(distance(v.begin(), unique(v.begin(), v.end())));
    
    for(unsigned int i = 0; i < v.size(); i++)
    {
        cout << "[ ";
        for(int j = 0; j < simplexStack[0]->GetDimension(); j++) cout << v[i][j] << " ";
        cout << "]" << endl;
    }
    cout << "***********************************\nTotal number of extreme points: " << v.size() << endl;
}

#ifdef CPLEX
    vector<string> GetVarNames(const CPXENVptr & env, const CPXLPptr & lp, int numCols)
    {
        vector<string> retVec;
        int           status = 0;
        char          **cur_colname = NULL;
        char          *cur_colnamestore = NULL;
        int           cur_colnamespace;
        int           surplus;
        
        status = CPXgetcolname (env, lp, NULL, NULL, 0, &surplus, 0, numCols-1);

        if (( status != CPXERR_NEGATIVE_SURPLUS ) && ( status != 0 ))  
        {
          printf ("Could not determine amount of space for column names.\n");
          exit(0);
        }

        cur_colnamespace = - surplus;
        cur_colname      = (char **) malloc (sizeof(char *)*numCols);
        cur_colnamestore = (char *)  malloc (cur_colnamespace);
        if ( cur_colname == NULL || cur_colnamestore == NULL ) 
        {
           printf ("Failed to get memory for column names.\n");
           exit(0);
        }
        status = CPXgetcolname (env, lp, cur_colname, cur_colnamestore, cur_colnamespace, &surplus, 0, numCols-1);
        if ( status ) 
        {
           printf ("CPXgetcolname failed.\n");
           exit(0);
        }
        
        for (int j = 0; j < numCols; j++) 
        {
          retVec.push_back(cur_colname[j]);
        }
        
        free(cur_colname);
        free(cur_colnamestore);
        
        return retVec;  
    }
#else
    vector<string> GetVarNames(glp_prob *lp, int numCols)
    {
        vector<string> retVec;
        const char     *colname = NULL;
        
        for(int i = 0; i < numCols; i++) 
        {
            colname = glp_get_col_name(lp, i+1);
/*            cout << colname << endl;*/
            if(colname) 
            {
                retVec.push_back(colname);
/*                cout << retVec[i] << endl;*/
            }
        }
        
        return retVec;  
    }
#endif

// bennett 7/18
void MultiobjectiveProblem::Epsilon(double e) {
    epsilon= e;
    return;
}

double MultiobjectiveProblem::Epsilon(void) const{
    return epsilon;
}

void deleteRepeats(vector<Simplex *> & simplexStack, int startingScanIndex)
{
    unsigned int i = startingScanIndex;
    
    while(i < simplexStack.size())
    {
        if(simplexStack[i]->deleteRepeats(simplexStack, i + 1)){
	    cerr << " ***********************************************" << endl;
	    cerr << "Deleting a repeat " << endl;
	    cerr << " ***********************************************" << endl;
	    simplexStack.erase(simplexStack.begin() + i);
	} else {
	   i++;
	}
    }
    
    return;
}

void CheckForDomination(const vector< vector<double> > & points, const double & epsilon)
{
    unsigned int k = points.size() - 1;
    bool oldDominatesNew = true, newDominatesOld = true;
    
    for(unsigned int i = 0; i < points.size() - 1; i++)
    {
        oldDominatesNew = true;
        newDominatesOld = true;
        if(DEBUG) cout << "Comparing extreme points " << i + 1 << " and " << k + 1 << ".\n";
        for(unsigned int j = 0; j < points[i].size(); j++)
        {
            if(DEBUG) cout << "\t " << j+1 << "th elements: " << points[i][j] << "\t" << points[k][j] << endl;
            if(points[k][j] > points[i][j] + epsilon) newDominatesOld = false;
            if(points[i][j] > points[k][j] + epsilon) oldDominatesNew = false;
            if(DEBUG) cout << "\t\t newDominatesOld: " << newDominatesOld << "\toldDominatesNew: " << oldDominatesNew << endl;
            if(!newDominatesOld && !oldDominatesNew) break;
            if(j == points.size() - 1)
            {
                cout << "There is no conflict between objectives " << i+1 << " and " << k+1 << ". It does not make sense to continue. ";
                cout << "You may opt to remove one of these objectives and resolve. Exiting!\n";
                exit(0);
            }
        }
    }
    
    return;
}
