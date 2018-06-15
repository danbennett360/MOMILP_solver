#ifndef PROBLEMCLASS
#define PROBLEMCLASS

/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <vector>
#include <string>
#include <iterator>
#include <iostream>
#include <algorithm>
#include "simplex_class.h"
#include "multiobjective_solver.h"

using namespace std;

class MultiobjectiveProblem 
{
    	CPXENVptr  env;
    	CPXLPptr mainProb;
    	CPXLPptr tempProb;
    	vector<CPXLPptr> lps;
    	vector< vector<double> > objectiveCoefs;
    	int numCols = 0;
    	int numRows = 0;
    	bool containsEqualities = false;
    	vector<int> indices;
    	vector<double> weightedCoefs;
    	vector<int> objectiveColIndices;
    	vector<int> objectiveRowIndices;
    	vector<Simplex> MeatOfDichotomicSearch();
//    	void PostProcessDichotomicSearch(vector<Simplex> & simplices);
    	void ChangeTempObjCoefs(int i);
    	void ChangeTempObjCoefs(const vector<double> & v);
    	double epsilon = .01;
//    	double epsilon2 = .0000001;
    	double infinity = CPX_INFBOUND;
/************************************************/
/*      Parameter values that are 
        user-modifyable                         */
/************************************************/
        int numObjectives = 2;
        bool storeObjectivesInMainProb = true;
        bool normalizeObjectiveMultipliers = true;
        bool useLexicographicOptimization = false;
  public:
    	void 	    SetEnv(CPXENVptr e);
    	CPXENVptr	GetEnv();
	    void 	    SetNumObj(int a);
	    void 	    AddLP(CPXLPptr lp);
	    int 		GetNumObj();
	    CPXLPptr 	GetLP(int i);
	    CPXLPptr 	GetMainLP();
	    void 	    ConvertLPs();
	    void 	    SetNumRowsAndCols();
	    void        AddRowsForObjectives();
	    void        SetParamVals(int argc, char **argv);
	    bool        StoreObjectivesInMainProb();
	    bool        NormalizeObjectiveMultipliers();
	    bool        UseLexicographicOptimization();
	    vector<Simplex> DichotomicSearch();
	    vector<Simplex> DichotomicSearch(const vector<int> & indices, const vector<int> & vals);
	    vector<double> GetObjectiveValues(const CPXLPptr & lp);
	    vector<double> LexicographicMinimization(int i);
};

void AddNewSimplices(   vector<Simplex> & simplexStack, const Simplex & currentSimplex, const vector<double> & point, bool normalize, bool useAdjacent, 
                        double epsilon);
                        
double GetAngleBetween(const Simplex & s1, const Simplex & s2, bool normalize);

double GetVectorMagnitude(const vector<double> & v);

void SplitSimplexInTwoUsingPoint(const Simplex & s, const vector<double> & point, vector<Simplex> & simplexStack, int newPointIndex, bool normalize);

void CheckForSimplicesThatNeedReplaced( vector<Simplex> & simplexStack, int & simplexIndex, int & newPointIndex, const int & numObjectives, 
                                        const vector<double> & newPoint, bool normalize, double epsilon);
                                        
void scanForRepeats(const vector<Simplex> & simplexStack);

bool PointIsInFrontOfAnAdjacent(const vector<Simplex> & simplexStack, const Simplex & simp, const vector<double> & point, const double & epsilon);

void WritePoints(const vector<Simplex> & simplexStack); 

vector<string> GetVarNames(const CPXENVptr & env, const CPXLPptr & lp, int numCols);

#endif
