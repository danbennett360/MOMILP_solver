/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#include "cplex.h"
#include "problem_class.h"
#include "point_class.h"

void MultiobjectiveProblem::SetEnv(CPXENVptr e)
{
	env = e;
	return;
}

CPXENVptr MultiobjectiveProblem::GetEnv()
{
	return env;
}

void MultiobjectiveProblem::SetNumObj(int a)
{
	numObjectives = a;
	return;
}

void MultiobjectiveProblem::AddLP(CPXLPptr lp)
{
	lps.push_back(lp);
	return;
}

int MultiobjectiveProblem::GetNumObj()
{
	return numObjectives;
}

CPXLPptr MultiobjectiveProblem::GetLP(int i)
{
	return lps[i];
}

void MultiobjectiveProblem::SetNumRowsAndCols()
{
     numCols = CPXgetnumcols (env, lps[0]);
     numRows = CPXgetnumrows (env, lps[0]);
}

void MultiobjectiveProblem::ConvertLPs()
{
	char *row_sense;
	row_sense = new char[numRows];
	int status = 0;
	double coef = 0.;
	double *objCoefs;
	objCoefs = new double[numCols];
	double rhs[1] = {0.};
	char new_sense[1] = {'L'};
	int k = 0;
	vector<double> temp(numCols);
	
	status = CPXgetsense (env, lps[k], row_sense, 0, numRows-1);
	
	for(int i = 0; i < numRows; i++)
	{
		if(row_sense[i] == 'G')
		{
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
		}
		if(!containsEqualities && row_sense[i] == 'E')
		{
			containsEqualities = true;
		}
	}
	for(int i = 0; i < numCols*2; i++) indices.push_back(i);

    for(int i = 0; i < numObjectives; i++)
    {
       	status = CPXgetobj(env, lps[i], objCoefs, 0, numCols -1);
       	if ( status ) 
       	{
         		printf ("CPXgetobj, Failed to read the objective coef. of reduced myBOMIP1, error code %d\n", status);
         		exit(0);
     	}
     	
     	status = CPXgetobjsen (env, lps[i]);
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
         	memcpy(&temp[0], &objCoefs[0], numCols*sizeof(double));
         	objectiveCoefs.push_back(temp);
     	}
     	else
     	{
     	    memcpy(&temp[0], &objCoefs[0], numCols*sizeof(double));
     	    objectiveCoefs.push_back(temp); 
     	}
	}
	
	mainProb = CPXcloneprob (env, lps[k], &status);
  	if ( status ) 
  	{
    		printf ("Failed to clone problem 1.\n");
    		exit(0);
  	}
  	
/*  	status = CPXwriteprob (env, mainProb, "mainprob.lp", "LP");*/
  	
  	delete[] row_sense; 
  	delete[] objCoefs;
    	
	return;
}

void MultiobjectiveProblem::AddRowsForObjectives()
{
     vector<double> lb(numObjectives, -infinity);
     vector<double> ub(numObjectives, infinity);
     vector<char> sense(numObjectives, 'E');
     vector<int> rowlist(numCols + numObjectives, numRows - 2);
	 vector<int> collist(numCols + numObjectives, 0);
	 vector<double> vallist(numCols + numObjectives, 0.);  	
     int status = 0;
     double dZero = 0.;
     
     for(int i = 0; i < numObjectives; i++)
     {
          objectiveColIndices.push_back(numCols + i);
          objectiveRowIndices.push_back(numRows + i);
/*          cout << "i: " << i << endl;*/
     }

	status = CPXnewcols (env, mainProb, numObjectives, NULL, &lb[0], &ub[0], NULL, NULL); 
	if ( status ) 
  	{
    		printf ("CPXnewcols, Failed to add additional variables to the model, error code %d\n", status);
    		exit(0);
    	}
	
	numCols += numObjectives;
	
	status = CPXnewrows (env, mainProb, numObjectives, NULL, &sense[0], NULL, NULL) ; 
	if ( status ) 
  	{
    		printf ("CPXgetobj, Failed to add additional constraints to the model, error code %d\n", status);
    		exit(0);
    }
    	
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
     		collist[i] = i;
     		
     		if(i < objectiveColIndices[j] && i < numCols - numObjectives)
     		{
/*     		     cout << vallist[i] << endl;*/
/*     		     cout << objectiveCoefs[j][i] << endl;*/
     		     vallist[i] = objectiveCoefs[j][i];
     		}
     		else if(i == objectiveColIndices[j]) vallist[i] = -1;
     		else vallist[i] = 0;
     	}
     	
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
	}
	
/*	status = CPXwriteprob (env, mainProb, "mainprob.lp", "LP");*/
	
	/*************************************************************************************/

    return;
}

void MultiobjectiveProblem::SetParamVals(int argc, char **argv)
{
     int i = numObjectives + 2;
     while(i < argc)
     {
        if(argv[i][0] != '-')
        {
            cout << "Invalid Command Line Argument (missing '-'). Exiting." << endl;
            exit(0);
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
		       i++;
		    } catch(const invalid_argument & ia) {
		       
                        cout << "Invalid value for flag '-epsilon'" <<  argv[i] 
			     << " , ignoring. Valid values are positive dobules." << endl;
		    }
		} else {
                    cout << "Invalid value for flag '-epsilon'" 
		         << ", ignoring. Valid values are positive dobules." << endl;
		}
            }
            else if(!strcmp(argv[i],"-normalize"))
            {
                i++;
                if(toupper(argv[i][0]) == 'T')
                {
                    normalizeObjectiveMultipliers = true;
                }
                else if(toupper(argv[i][0]) == 'F')
                {
                    normalizeObjectiveMultipliers = false;
                }
                else
                {
                    cout << "Invalid value for flag '-normalize', ignoring. Valid values are 'T' and 'F'." << endl;
                }
                i++;
            }
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

vector<Simplex> MultiobjectiveProblem::DichotomicSearch()
{
    int status = 0;
    vector<Simplex> retVec;
    
    tempProb = CPXcloneprob (env, mainProb, &status);
    cout << status << endl;
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
  	
/*  	status = CPXwriteprob (env, tempProb, "tempprob.lp", "LP");*/
  	
  	retVec = MeatOfDichotomicSearch();
  	return retVec;
}

vector<Simplex> MultiobjectiveProblem::MeatOfDichotomicSearch()
{
    int status = 0;
    vector<Simplex> retVec;
    vector<Simplex> simplexStack;
    vector<Point> pointStack;
    Simplex temp(numObjectives);
    Simplex temp2(numObjectives);
    Simplex currentSimplex(numObjectives);
    vector<double> point(numObjectives, 0.);
    vector<double> point2(numObjectives, 0.);
    vector<double> mins(numObjectives, infinity);
    vector<double> maxs(numObjectives, -infinity);
    vector< vector<double> > extremes;
    vector<string> varNames;
    int iterator = 0;
    int simplexIndex = 0, newPointIndex = numObjectives-1;
    vector<int> removeTheseIndicesFromSimplexStack;
    unsigned int k = 0;
    int numSaved = 0;
    long maxIter = 100000;
    double *x;
    
/*    bool DEBUG = true;*/
    
    // Get initial simplex using lexicographic minimization
    for(int i = 0; i < numObjectives; i++)
    {
        extremes.push_back(LexicographicMinimization(i));
        if(SAVE_POINTS) 
        {
            x = new double[numCols];
            status = CPXgetx (env, tempProb, x, 0, numCols-1);
            if ( status ) {
        		printf ("%s(%d): CPXgetx, Failed to get x values,  error code %d\n", __FILE__, __LINE__, status);
        		exit(0);
      	    }
            pointStack.push_back(Point(extremes[i], x, numCols));
            delete[] x;
        }
    }
    
/*    cout << extremes.size() << endl;*/
    
    for(int i = 0; i < numObjectives; i++)
    {
        for(int j = 0; j < numObjectives; j++)
        {
            mins[i] = min(mins[i], extremes[j][i]);
            maxs[i] = max(maxs[i], extremes[j][i]);
            cout << extremes[i][j] << "\t";
        }
        cout << endl;
        cout << mins[i] << "\t" << maxs[i] << endl;
    }
    
    for(int i = 0; i < numObjectives; i++)
    {
        for(int j = 0; j < numObjectives; j++)
        {
            if(i == j) point[j] = infinity; 
            else point[j] = mins[j];
            cout << point[j] << "\t";
        }
        cout << endl;
        temp.AddExtreme(point, NormalizeObjectiveMultipliers());
    }
    
    simplexStack.push_back(temp);
    if(SCAN_FOR_REPEATS)
    {
        scanForRepeats(simplexStack);
    }
    
    if(DEBUG)
    {
        cout << "**************************\nInitial Simplex: " << endl;
        for(unsigned int i = 0; i < simplexStack.size(); i++) simplexStack[i].WriteOctaveCodeToPlotSimplex();
        cout << "**************************" << endl;
    }
    
    for(int i = 0; i < numObjectives; i++)
    {
        if(i == 0)
        {
            currentSimplex = simplexStack[simplexStack.size()-1];
            simplexStack.pop_back();
        }
        else
        {
            k = 0;
            while(simplexStack[k].MultiplyPointsByNormal(extremes[i]) >= simplexStack[k].PlaneVal() - epsilon) 
            {
                if(k >= simplexStack.size())
                {
                    cout << "Error, accessing elements past the end of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                    exit(0);
                }
                k++;
            }
            currentSimplex = simplexStack[k];
            simplexStack.erase(simplexStack.begin() + k);
            
        }
        AddNewSimplices(simplexStack, currentSimplex, extremes[i], NormalizeObjectiveMultipliers(), false, epsilon);
        CheckForSimplicesThatNeedReplaced(simplexStack, simplexIndex, newPointIndex, numObjectives, extremes[i], NormalizeObjectiveMultipliers(), epsilon);
        if(DEBUG)
        {
            cout << "**************************\nStored simplices after adding " << i+1 << "th initial extreme point: " << endl;
            for(unsigned int j = 0; j < simplexStack.size(); j++) simplexStack[j].WriteOctaveCodeToPlotSimplex();
            cout << "**************************" << endl;
        }
    }
    
    k = simplexStack.size()-1;
    while(k > 0 && iterator < maxIter)
    {
        iterator++;
        if(DEBUG)
        {
            cout << "-----------------------------------------------------------------------------------" << endl;
            cout << "Iteration: " << iterator << "\tSize of stack: " << simplexStack.size() <<endl;
            cout << "**************************\nSimplices in stack: " << endl;
            for(unsigned int i = 0; i < simplexStack.size(); i++) simplexStack[i].WriteOctaveCodeToPlotSimplex();
            cout << "\%**************************\n\%Simplices stored in return vector:" << endl;
            for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex();
            cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
        }
        
        k = simplexStack.size()-1;
        while(simplexStack[k].GetSaveForSolution() && k > 0) k--;
        currentSimplex = simplexStack[k];
        simplexStack.erase(simplexStack.begin() + k);
/*        simplexStack.pop_back();*/
        
/*        currentSimplex.PrintData();*/
        
/*        status = CPXwriteprob (env, tempProb, "prob_before.lp", "LP");*/
/*        exit(0);*/
        
        ChangeTempObjCoefs(currentSimplex.GetNormal());
        
/*        status = CPXwriteprob (env, tempProb, "prob.lp", "LP");*/
/*        exit(0);*/
        
        status = CPXlpopt (env, tempProb);
     	if ( status ) {
        		printf ("%s(%d): CPXlpopt, Failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
        		exit(0);
      	}
      	
      	if(DEBUG)
      	{
          	point = currentSimplex.GetNormal();
          	cout << "new current normal is: <";
          	for(int i = 0; i < numObjectives; i++) cout << point[i] << "\t";
          	cout << ">" << endl;
      	}
      	
      	point = GetObjectiveValues(tempProb);
      	if(SAVE_POINTS) 
        {
            x = new double[numCols];
            status = CPXgetx (env, tempProb, x, 0, numCols-1);
            if ( status ) {
        		printf ("%s(%d): CPXgetx, Failed to get x values,  error code %d\n", __FILE__, __LINE__, status);
        		exit(0);
      	    }
            pointStack.push_back(Point(point, x, numCols));
            delete[] x;
        }
      	
      	if(DEBUG)
      	{
          	cout << "new point is: (";
          	for(int i = 0; i < numObjectives; i++) cout << point[i] << "\t";
          	cout << ")" << endl;
          	cout << currentSimplex.MultiplyPointsByNormal(point) << "\t" << currentSimplex.PlaneVal() << endl;
          	currentSimplex.WriteOctaveCodeToPlotSimplex();
          	cout << "Distance to simplex: " << currentSimplex.MultiplyPointsByNormal(point) - (currentSimplex.PlaneVal() - epsilon) << "\nIn Front of An Adjacent: " << PointIsInFrontOfAnAdjacent(simplexStack, currentSimplex, point, epsilon) << endl;
      	}
      	
      	if(currentSimplex.MultiplyPointsByNormal(point) < currentSimplex.PlaneVal() - epsilon)// || PointIsInFrontOfAnAdjacent(simplexStack, currentSimplex, point, epsilon))
      	{
      	    if(DEBUG)
      	    {
/*      	        cout << "Distance to simplex: " << currentSimplex.MultiplyPointsByNormal(point) - (currentSimplex.PlaneVal() - epsilon) << "\nIn Front of An Adjacent: " << PointIsInFrontOfAnAdjacent(simplexStack, currentSimplex, point, epsilon) << endl;*/
      	        cout << "Adding three new simplices to stack" << endl;
      	    }
      	    
      	    AddNewSimplices(simplexStack, currentSimplex, point, NormalizeObjectiveMultipliers(), false, epsilon);

            CheckForSimplicesThatNeedReplaced(simplexStack, simplexIndex, newPointIndex, numObjectives, point, NormalizeObjectiveMultipliers(), epsilon);
      	}
      	else
      	{
/*      	    retVec.push_back(currentSimplex);*/
            currentSimplex.SaveForSolution();
            simplexStack.push_back(currentSimplex);
      	    
      	    if(DEBUG) cout << "No new point found" << endl;
      	}
    }
    if(iterator == maxIter)
    {
        cout << "Max iterations (" << maxIter << ") reached in meat of dichotomic search. Consider increasing value. Exiting!" << endl;
        
        cout << "**************************\nSimplices still on stack: " << endl;
        for(unsigned int j = 0; j < simplexStack.size(); j++) 
        {
            if(simplexStack[j].GetSaveForSolution()) numSaved++;
            simplexStack[j].WriteOctaveCodeToPlotSimplex();
        }
        cout << "\%**************************" << endl;
        
        cout << "Number of Simplices in List: " << simplexStack.size() << "\tNumber Saved: " << numSaved << "\tNumber left to search: " << simplexStack.size() - numSaved << endl;
        
/*        cout << "\%**************************\n\%Simplices discovered so far: " << endl;*/
/*        for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex();*/
/*        cout << "**************************" << endl;*/
        
        exit(0);
    }
    
/*    if(DEBUG) */
/*    {*/
  	    cout << "**************************\nSimplices after dichotomic search (which took " << iterator << " iterations): " << endl;
        for(unsigned int i = 0; i < simplexStack.size(); i++) simplexStack[i].WriteOctaveCodeToPlotSimplex();
        cout << "**************************\nTotal number in list is: " << simplexStack.size() << endl;
        
        cout << "*************************************\nThe individual extreme points discovered are:\n";
        WritePoints(simplexStack); 
/*    }*/
    
/*    PostProcessDichotomicSearch(retVec);*/
/*    */
/*    if(DEBUG) */
/*    {*/
/*  	    cout << "**************************\nSimplices after post processing: " << endl;*/
/*        for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex();*/
/*        cout << "**************************" << endl;*/
/*    }*/

    // bennett 6/18
    if(SAVE_POINTS) 
    {
       varNames = GetVarNames(env, tempProb, numCols);

       vector<Point>tmp;
       bool append = false;

       for (auto & tmpPt : pointStack) {
          if (find(tmp.begin(),tmp.end(),tmpPt) == tmp.end())  {
	      tmp.push_back(tmpPt);
	      tmpPt.WritePointToFile("points.txt", varNames, append);
	      append = true;
	  }

       }
    }
    
    return simplexStack;
}

void MultiobjectiveProblem::ChangeTempObjCoefs(int i)
{
    int status = CPXchgobj (env, tempProb, objectiveColIndices[0], &indices[0], &objectiveCoefs[i][0]);
	if ( status ) {
		printf ("Failed to change obj coef. Error code %d\n", status);
	}
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
   
    
    status = CPXchgobj (env, tempProb, objectiveColIndices[0], &indices[0], &weightedCoefs[0]);
	if ( status ) {
		printf ("Failed to change obj coef. Error code %d\n", status);
	}
	
/*	status = CPXwriteprob (env, tempProb, "prob_after.lp", "LP");*/
/*	 exit(0);*/
	return;
}

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
    
    ChangeTempObjCoefs(i);
    
    status = CPXlpopt (env, tempProb);
 	if ( status ) {
    		printf ("%s(%d): CPXlpopt, Failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
    		exit(0);
  	}
  	
  	point = GetObjectiveValues(tempProb);
  	
/*  	for(int k = 0; k < numObjectives; k++) cout << point[k] << "\t" << endl;*/
  	
  	//add bounds to set max value of objective i
  	if(UseLexicographicOptimization())
  	{
      	if(StoreObjectivesInMainProb())
      	{
      	    status = CPXchgbds (env, tempProb, 1, &objectiveColIndices[i], &u, &point[i]);
      	}
      	else
      	{
      	    vector<int> rowIndices(objectiveColIndices[0],numRows);
      	    
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
        }
      	
      	for(int j = 0; j < numObjectives; j++)
      	{
            //optimize over additional objectives
            if(j != i)
            {
                ChangeTempObjCoefs(j);
                        
                status = CPXlpopt (env, tempProb);
             	if ( status ) {
                		printf ("%s(%d): CPXlpopt, Failed to solve tempProb,  error code %d\n", __FILE__, __LINE__, status);
                		exit(0);
              	}
              	
              	point = GetObjectiveValues(tempProb);
              	
/*              	for(int k = 0; k < numObjectives; k++) cout << point[k] << "\t" << endl;*/
              	
              	if(StoreObjectivesInMainProb())
              	{
              	    status = CPXchgbds (env, tempProb, 1, &objectiveColIndices[j], &u, &point[j]);
              	}
              	else
              	{
              	    vector<int> rowIndices(objectiveColIndices[0],numRows+j+1);
              	    
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
                }
            }

      	}
      	
      	if(StoreObjectivesInMainProb())
      	{
      	    vector<char> up(numObjectives, 'U');
      	    vector<double> bd(numObjectives, infinity);
      	    status = CPXchgbds (env, tempProb, numObjectives, &objectiveColIndices[0], &up[0], &bd[0]);
      	}
      	else
      	{
      	    status = CPXdelrows (env, tempProb, numRows, numRows+numObjectives);
        }
    }
  	
  	return point;
}

void AddNewSimplices(   vector<Simplex> & simplexStack, const Simplex & currentSimplex, const vector<double> & point, bool normalize, bool useAdjacent, 
                        double epsilon)
{
    vector< vector<double> > extremePoints = currentSimplex.GetExtremePoints();
    int dim = currentSimplex.GetDimension();
    Simplex temp(dim);
    Simplex temp2(dim);
    int simplexIndex = 0;
    int newPointIndex = 0;
    
    epsilon *= epsilon;
    
/*    bool DEBUG = true;*/
    
    for(int i = 0; i < dim; i++)
    {
        temp.Reset();
        
        if(useAdjacent)
        {
          	temp2 = currentSimplex.FindAdjacentContainingOriginalPoints(simplexStack, simplexIndex, newPointIndex, i, false);
          	if(temp2.GetDimension() > 0 && temp2.MultiplyPointsByNormal(point) < temp2.PlaneVal() - epsilon)
          	{
          	    if(DEBUG) 
          	    {
              	    cout << "\%There is a 'dominated' simplex which needs replaced. Its index is " << simplexIndex << ". The new point index is: " << newPointIndex << ". The simplex is shown below: " << endl;
              	    temp2.WriteOctaveCodeToPlotSimplex();
          	    }
                for(int j = 0; j < dim; j++)
                {
                    if(j != newPointIndex) 
                    {
                        temp.Reset();
      	                for(int k = 0; k < dim; k++)
      	                {
      	                    if(k != j) temp.AddExtreme(temp2.GetExtremePoint(k), normalize);
      	                }
      	                temp.AddExtreme(point, normalize);
      	                if(DEBUG) temp.WriteOctaveCodeToPlotSimplex();
                        simplexStack.push_back(temp);
                        if(SCAN_FOR_REPEATS)
                        {
                            scanForRepeats(simplexStack);
                        }
                    }
                }
                if(simplexIndex >= int(simplexStack.size()))
                {
                    cout << "Error, accessing elements past the end of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                    exit(0);
                }
                else if( simplexIndex < 0 )
                {
                    cout << "Error, accessing elements before the start of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                    exit(0);
                }
                simplexStack.erase(simplexStack.begin() + simplexIndex);
            }
            else
            {
                for(int j = 0; j < dim; j++)
                {
                    if(j != i) temp.AddExtreme(extremePoints[j], normalize);
                }
                temp.AddExtreme(point, normalize);
                if(DEBUG) temp.WriteOctaveCodeToPlotSimplex();
                simplexStack.push_back(temp);
                if(SCAN_FOR_REPEATS)
                {
                    scanForRepeats(simplexStack);
                }
            }
        }
        else
        {
            for(int j = 0; j < dim; j++)
            {
                if(j != i) temp.AddExtreme(extremePoints[j], normalize);
            }
            temp.AddExtreme(point, normalize);
            if(DEBUG) temp.WriteOctaveCodeToPlotSimplex();
            simplexStack.push_back(temp);
            if(SCAN_FOR_REPEATS)
            {
                scanForRepeats(simplexStack);
            }
        }
    }    
    
/*    cout << "accessing last element of simplexStack, which has size: " << simplexStack.size() << endl;*/
/*    temp = simplexStack[simplexStack.size()-1];*/
    
    return;
}

bool MultiobjectiveProblem::UseLexicographicOptimization()
{
    return useLexicographicOptimization;
}

double GetAngleBetween(const Simplex & s1, const Simplex & s2, bool normalize)
{
    vector<double> v1 = s1.GetNormal();
    double dotProduct = s2.MultiplyPointsByNormal(v1);
    cout << "dot product: " << dotProduct << endl;
    double cosOfAngle = 0.;
    if(!normalize) cosOfAngle = dotProduct/(GetVectorMagnitude(v1)*GetVectorMagnitude(s2.GetNormal()));
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

void SplitSimplexInTwoUsingPoint(const Simplex & s, const vector<double> & point, vector<Simplex> & simplexStack, int newPointIndex, bool normalize)
{
    int dim = s.GetDimension();
    bool save = s.GetSaveForSolution();
    Simplex temp(dim);
    
/*    bool DEBUG = true;*/
    
    if(DEBUG) 
    {
        cout << "\%There is a 'dominated' simplex which needs replaced. The new point index is: " << newPointIndex << ". The simplex is shown below: " << endl;
        s.WriteOctaveCodeToPlotSimplex();
    }
    for(int j = 0; j < dim; j++)
    {
        if(j != newPointIndex) 
        {
            temp.Reset();
            if(save) temp.SaveForSolution();
            for(int k = 0; k < dim; k++)
            {
                if(k != j) temp.AddExtreme(s.GetExtremePoint(k), normalize);
            }
            temp.AddExtreme(point, normalize);
            if(DEBUG) 
            {
                temp.WriteOctaveCodeToPlotSimplex();
                cout << "size: " << simplexStack.size() << "\n";
                cout << "capacity: " << simplexStack.capacity() << "\n";
                cout << "max_size: " << simplexStack.max_size() << "\n";
            }
            simplexStack.push_back(temp);
            if(SCAN_FOR_REPEATS)
            {
                scanForRepeats(simplexStack);
            }
        }
    }
    return;
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
/*    for(unsigned int k = 0; k < potentials.size(); k++) potentials[k].WriteOctaveCodeToPlotSimplex();*/
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
/*                        temp2.WriteOctaveCodeToPlotSimplex();*/
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

void CheckForSimplicesThatNeedReplaced( vector<Simplex> & simplexStack, int & simplexIndex, int & newPointIndex, const int & numObjectives, 
                                        const vector<double> & newPoint, bool normalize, double epsilon)
{
    int k = simplexStack.size() - numObjectives;
    Simplex temp(numObjectives);
    Simplex temp2(numObjectives);
    
    epsilon *= epsilon;
    
/*    bool DEBUG = true;*/
    
    while(k < int(simplexStack.size()))
    {
        if(k >= int(simplexStack.size()))
        {
            cout << "Error, accessing elements past the end of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
            exit(0);
        }
        else if( k < 0 )
        {
            cout << "Error, accessing elements before the start of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
            exit(0);
        }
        temp = simplexStack[k];
        if(DEBUG)
        {
            cout << "$$$$$$$$$$$$$$$$$$$$$$$\n scanning this simplex: "<< endl;
            temp.WriteOctaveCodeToPlotSimplex();
            cout << "$$$$$$$$$$$$$$$$$$$$$$$"<< endl;
        }
        simplexIndex = k;
        newPointIndex = numObjectives - 1;
        temp2 = temp.FindAdjacentContainingOriginalPoints(simplexStack, simplexIndex, newPointIndex, numObjectives - 1, true);
        if(!temp.IsPositive())
        {
            if(DEBUG) cout << "normal vector is not oriented correctly" << endl;
           
            SplitSimplexInTwoUsingPoint(temp2, newPoint, simplexStack, newPointIndex, normalize);
            if(DEBUG)
            {
                cout << "\%erasing simplex: " << endl;
                simplexStack[k].WriteOctaveCodeToPlotSimplex();
                cout << "\%erasing simplex: " << endl;
                simplexStack[simplexIndex].WriteOctaveCodeToPlotSimplex();
            }
            simplexStack.erase(simplexStack.begin() + k);
            if(simplexIndex >= int(simplexStack.size()))
            {
                cout << "Error, accessing elements past the end of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                exit(0);
            }
            else if( simplexIndex < 0 )
            {
                cout << "Error, accessing elements before the start of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                exit(0);
            }
            simplexStack.erase(simplexStack.begin() + simplexIndex);
            k -= 2;
        }
        else
        {
/*            if(temp2.GetDimension() > 0) cout << temp2.MultiplyPointsByNormal(newPoint) << "\t" << temp2.PlaneVal() - epsilon << "\t" << temp2.MultiplyPointsByNormal(newPoint) - (temp2.PlaneVal() - epsilon) << endl;*/
            if(temp2.GetDimension() > 0 && temp2.MultiplyPointsByNormal(newPoint) < temp2.PlaneVal() - epsilon)
            {
                if(DEBUG)
                {
/*                            temp.PrintData();*/
/*                            temp2.PrintData();*/
/*                            temp.WriteOctaveCodeToPlotSimplex();*/
/*                            temp2.WriteOctaveCodeToPlotSimplex();*/
                    cout << "found two simplices which cause nonconvexity" << endl;
                    cout << "\%erasing simplex: " << endl;
                    simplexStack[k].WriteOctaveCodeToPlotSimplex();
/*                    temp.PrintData();*/
                    cout << "\%erasing simplex: " << endl;
                    simplexStack[simplexIndex].WriteOctaveCodeToPlotSimplex();
/*                    temp2.PrintData();*/
                }
                
                SplitSimplexInTwoUsingPoint(temp2, newPoint, simplexStack, newPointIndex, normalize);
                simplexStack.erase(simplexStack.begin() + k);
                if(simplexIndex >= int(simplexStack.size()))
                {
                    cout << "Error, accessing elements past the end of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                    exit(0);
                }
                else if( simplexIndex < 0 )
                {
                    cout << "Error, accessing elements before the start of the stack! File: " << __FILE__ << "  Line: " << __LINE__ << "  Exiting!\n";
                    exit(0);
                }
                simplexStack.erase(simplexStack.begin() + simplexIndex);
                k -= 2;
            }
        }
        k++;
    }
            
    return;
}

CPXLPptr MultiobjectiveProblem::GetMainLP()
{
    return mainProb;
}

void scanForRepeats(const vector<Simplex> & simplexStack)
{
    vector< vector<double> > v;
    int count = 0;
    
    for(unsigned int i = 0; i < simplexStack.size() - 1; i++)
    {
        v.resize(0);
        for(int k = 0; k < simplexStack[0].GetDimension(); k++)
        {
            v.push_back(simplexStack[simplexStack.size()-1].GetExtremePoint(k));
        }
    
        for(int j = 0; j < simplexStack[0].GetDimension(); j++)
        {
            v.push_back(simplexStack[i].GetExtremePoint(j));
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
        if(count <= simplexStack[0].GetDimension())
        {
            cout << "Error, a repeated simplex was just added to the stack! Exiting!\n";
            if(DEBUG)
            {
                cout << "**************************\nSimplices: " << endl;
                for(unsigned int j = 0; j < simplexStack.size(); j++) simplexStack[j].WriteOctaveCodeToPlotSimplex();
                cout << "**************************" << endl;
            }
            exit(0);
        }
    }
}

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

void WritePoints(const vector<Simplex> & simplexStack)
{
    vector< vector<double> > v;
    
    for(unsigned int i = 0; i < simplexStack.size(); i++)
    {
        for(int j = 0; j < simplexStack[0].GetDimension(); j++) 
        {
            if(!simplexStack[i].ExtremesAreDummies()[j]) v.push_back(simplexStack[i].GetExtremePoint(j));
        }
    }
    sort(v.begin(), v.end());
    v.resize(distance(v.begin(), unique(v.begin(), v.end())));
    
    for(unsigned int i = 0; i < v.size(); i++)
    {
        cout << "[ ";
        for(int j = 0; j < simplexStack[0].GetDimension(); j++) cout << v[i][j] << " ";
        cout << "]" << endl;
    }
    cout << "***********************************\nTotal number of extereme points: " << v.size() << endl;
}

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

// bennett 7/18
void MultiobjectiveProblem::Epsilon(double e) {
    epsilon= e;
    return;
}

double MultiobjectiveProblem::Epsilon(void) const{
    return epsilon;
}
