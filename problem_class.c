/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#include "cplex.h"
#include "problem_class.h"

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
	}
	
	mainProb = CPXcloneprob (env, lps[k], &status);
  	if ( status ) 
  	{
    		printf ("Failed to clone problem 1.\n");
    		exit(0);
  	}
  	
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
	
	for(int j = 0; j < numObjectives; j++)
	{
     	for(int i = 0; i < numCols; i++)
     	{
     	     rowlist[i] = objectiveRowIndices[j];
     		collist[i] = i;
     		if(i < objectiveColIndices[j] && i < numCols - numObjectives) vallist[i] = objectiveCoefs[j][i];
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
  	
  	retVec = MeatOfDichotomicSearch();
  	return retVec;
}

vector<Simplex> MultiobjectiveProblem::MeatOfDichotomicSearch()
{
    int status = 0;
    vector<Simplex> retVec;
    vector<Simplex> simplexStack;
    Simplex temp(numObjectives);
    Simplex temp2(numObjectives);
    Simplex currentSimplex(numObjectives);
    vector<double> point(numObjectives, 0.);
    vector<double> point2(numObjectives, 0.);
    vector<double> mins(numObjectives, infinity);
    vector<double> maxs(numObjectives, -infinity);
    vector< vector<double> > extremes;
    int iterator = 0;
    int simplexIndex = 0, newPointIndex = numObjectives-1;
    vector<int> removeTheseIndicesFromSimplexStack;
    unsigned int k = 0;
    long maxIter = 50;
    
/*    bool DEBUG = true;*/
    
    // Get initial simplex using lexicographic minimization
    for(int i = 0; i < numObjectives; i++)
    {
        extremes.push_back(LexicographicMinimization(i));
    }
    
    for(int i = 0; i < numObjectives; i++)
    {
        for(int j = 0; j < numObjectives; j++)
        {
            mins[i] = min(mins[i], extremes[j][i]);
            maxs[i] = max(maxs[i], extremes[j][i]);
        }
    }
    
    for(int i = 0; i < numObjectives; i++)
    {
        for(int j = 0; j < numObjectives; j++)
        {
            if(i == j) point[j] = infinity; 
            else point[j] = mins[j];
        }
        temp.AddExtreme(point, NormalizeObjectiveMultipliers());
    }
    
    simplexStack.push_back(temp);
    
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
            while(simplexStack[k].MultiplyPointsByNormal(extremes[i]) >= simplexStack[k].PlaneVal() - epsilon) k++;
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
    
    while(simplexStack.size() > 0 && iterator < maxIter)
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
        
        currentSimplex = simplexStack[simplexStack.size()-1];
        simplexStack.pop_back();
        
        ChangeTempObjCoefs(currentSimplex.GetNormal());
        
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
      	
      	if(DEBUG)
      	{
          	cout << "new point is: (";
          	for(int i = 0; i < numObjectives; i++) cout << point[i] << "\t";
          	cout << ")" << endl;
          	cout << currentSimplex.MultiplyPointsByNormal(point) << "\t" << currentSimplex.PlaneVal() << endl;
          	currentSimplex.WriteOctaveCodeToPlotSimplex();
      	}
      	
      	if(currentSimplex.MultiplyPointsByNormal(point) < currentSimplex.PlaneVal() - epsilon)
      	{
      	    if(DEBUG) cout << "Adding three new simplices to stack" << endl;
      	    
      	    AddNewSimplices(simplexStack, currentSimplex, point, NormalizeObjectiveMultipliers(), false, epsilon);

            CheckForSimplicesThatNeedReplaced(simplexStack, simplexIndex, newPointIndex, numObjectives, point, NormalizeObjectiveMultipliers(), epsilon);
      	}
      	else
      	{
      	    retVec.push_back(currentSimplex);
      	    
      	    if(DEBUG) cout << "No new point found" << endl;
      	}
    }
    if(iterator == maxIter)
    {
        cout << "Max iterations (" << maxIter << ") reached in meat of dichotomic search. Consider increasing value. Exiting!" << endl;
        exit(0);
    }
    
    if(DEBUG) 
    {
  	    cout << "**************************\nSimplices after dichotomic search: " << endl;
        for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex();
        cout << "**************************" << endl;
    }
    
/*    PostProcessDichotomicSearch(retVec);*/
/*    */
/*    if(DEBUG) */
/*    {*/
/*  	    cout << "**************************\nSimplices after post processing: " << endl;*/
/*        for(unsigned int i = 0; i < retVec.size(); i++) retVec[i].WriteOctaveCodeToPlotSimplex();*/
/*        cout << "**************************" << endl;*/
/*    }*/
    
    return retVec;
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
    if(int(weightedCoefs.size()) < objectiveColIndices[0]) weightedCoefs.resize(objectiveColIndices[0]);
    
    for(int i = 0; i < objectiveColIndices[0] + 1; i++)
    {
        weightedCoefs[i] = 0.;
        for(int j = 0; j < numObjectives; j++)
        {
            weightedCoefs[i] += objectiveCoefs[j][i]*v[j];
        }
    }
    
    int status = CPXchgobj (env, tempProb, objectiveColIndices[0], &indices[0], &weightedCoefs[0]);
	if ( status ) {
		printf ("Failed to change obj coef. Error code %d\n", status);
	}
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
                    }
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
        }
    }    
    
/*    cout << "accessing last element of simplexStack, which has size: " << simplexStack.size() << endl;*/
    temp = simplexStack[simplexStack.size()-1];
    
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
            for(int k = 0; k < dim; k++)
            {
                if(k != j) temp.AddExtreme(s.GetExtremePoint(k), normalize);
            }
            temp.AddExtreme(point, normalize);
            if(DEBUG) temp.WriteOctaveCodeToPlotSimplex();
            simplexStack.push_back(temp);
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
    
/*    bool DEBUG = true;*/
    
    while(k < int(simplexStack.size()))
    {
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
            simplexStack.erase(simplexStack.begin() + simplexIndex);
            k -= 2;
        }
        else
        {
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
                    cout << "\%erasing simplex: " << endl;
                    simplexStack[simplexIndex].WriteOctaveCodeToPlotSimplex();
                }
                
                SplitSimplexInTwoUsingPoint(temp2, newPoint, simplexStack, newPointIndex, normalize);
                simplexStack.erase(simplexStack.begin() + k);
                simplexStack.erase(simplexStack.begin() + simplexIndex);
                k -= 2;
            }
        }
        k++;
    }
            
    return;
}
