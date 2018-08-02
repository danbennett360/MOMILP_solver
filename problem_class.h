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

#include "SimplexStore.h"

using namespace std;

class MultiobjectiveProblem 
{
        #ifdef CPLEX
        	CPXENVptr  env;
        	CPXLPptr mainProb;
        	CPXLPptr tempProb;
        	vector<CPXLPptr> lps;
        #else
            glp_prob *mainProb;
        	glp_prob *tempProb;
            vector<glp_prob *> lps;
    	#endif
    	
    	vector< vector<double> > objectiveCoefs;
    	int numCols = 0;
    	int numRows = 0;
    	bool containsEqualities = false;
    	vector<int> indices;
    	vector<double> weightedCoefs;
    	vector<int> objectiveColIndices;
    	vector<int> objectiveRowIndices;
    	//vector<Simplex *> MeatOfDichotomicSearch();
    	void MeatOfDichotomicSearch(SimplexStore &);
//    	void PostProcessDichotomicSearch(vector<Simplex> & simplices);
    	void ChangeTempObjCoefs(int i);
    	void ChangeTempObjCoefs(const vector<double> & v);
    	double epsilon = .05;
//    	double epsilon2 = .0000001;
        
        #ifdef CPLEX
    	    double infinity = CPX_INFBOUND;
    	#else
    	    double infinity = 100000000000000000000.;
    	#endif
/************************************************/
/*      Parameter values that are 
        user-modifyable                         */
/************************************************/
        int numObjectives = 2;
        bool storeObjectivesInMainProb = true;
        bool normalizeObjectiveMultipliers = true;
        bool useLexicographicOptimization = false;
        bool showProgress = false;
        int showProgressIterations = 500;
        int maxIterations = 100000;
        bool relativeDistance = false;
        bool interiorPoint = false;

	// bennett private member functions
	Point MakePoint(vector<double> point);
	void AddFirstSimplex(SimplexStore & store, vector<vector<double> > & extremes, vector<Simplex *> & simplexStack);
	void DoInitialSplits(SimplexStore & store, vector<vector<double> > & extremes, vector<Simplex *> & simplexStack);
	void DumpPoints(SimplexStore & store);
        void DoSplit(vector<Simplex *> & simplexStack, SimplexStore & simplexStore, vector<Simplex *> & simplicesToSplit, vector<double> & point);

  public:
        #ifdef CPLEX
        	void 	    SetEnv(CPXENVptr e);
        	CPXENVptr	GetEnv();
        	void 	    AddLP(CPXLPptr lp);
        	CPXLPptr 	GetLP(int i);
	        CPXLPptr 	GetMainLP();
	        vector<double> GetObjectiveValues(const CPXLPptr & lp);
	    #else
	        void 	    AddLP(glp_prob *lp);
	        glp_prob    *GetLP(int i);
	        glp_prob    *GetMainLP();
	        vector<double> GetObjectiveValues(glp_prob *lp);
	    #endif
	    
	    void 	    SetNumObj(int a);
	    int 		GetNumObj();
	    void 	    ConvertLPs();
	    void 	    SetNumRowsAndCols();
	    void        AddRowsForObjectives();
	    void        SetParamVals(int argc, char **argv);
	    bool        StoreObjectivesInMainProb();
	    bool        NormalizeObjectiveMultipliers();
	    bool        UseLexicographicOptimization();
	    void DichotomicSearch(SimplexStore &);
	    // does not seem to be implemented.  Will need to change anyway, I suspect it will need
	    // to have a SimplexStore passed in.
	    //SimplexStore DichotomicSearch(const vector<int> & indices, const vector<int> & vals);
	    vector<double> LexicographicMinimization(int i);

	    // bennett 7/18
	    void 	Epsilon(double e);
	    double	Epsilon(void) const;
};

/*
void AddNewSimplices(   vector<Simplex *> & simplexStack, const Simplex * currentSimplex, const vector<double> & point, bool normalize, bool useAdjacent, 
                        bool relativeDistance, double epsilon, SimplexStore & simplexStore);
*/
/*
bool SplitSimplexInTwoUsingPoint(const Simplex * s, const vector<double> & point, vector<Simplex *> & simplexStack, int newPointIndex, bool normalize, int startingScanIndex);
void CheckForSimplicesThatNeedReplaced( vector<Simplex *> & simplexStack, int & simplexIndex, int & newPointIndex, const int & numObjectives, 
                                        const vector<double> & newPoint, bool normalize, bool relativeDistance, double epsilon);
*/

bool PointIsInFrontOfAnAdjacent(const vector<Simplex *> & simplexStack, const Simplex * simp, const vector<double> & point, const double & epsilon);

double GetAngleBetween(const Simplex * s1, const Simplex * s2, bool normalize);

double GetVectorMagnitude(const vector<double> & v);
                                        
void CheckIfAdjacentsAreShadowed( vector<Simplex *> & simplexStack, vector<Simplex *> & simplicesToSplit, const Simplex * currentSimplex, 
                                  const vector<double> & newPoint, const int & numObjectives, bool relativeDistance, double epsilon, SimplexStore & simplexStore);
                                        
void scanForRepeats(const vector<Simplex *> & simplexStack);

void scanForNegativeNormal(const vector<Simplex *> & simplexStack);


void WritePoints(const vector<Simplex *> & simplexStack); 

#ifdef CPLEX
    vector<string> GetVarNames(const CPXENVptr & env, const CPXLPptr & lp, int numCols);
#else
    vector<string> GetVarNames(glp_prob *lp, int numCols);
#endif

void deleteRepeats(vector<Simplex *> & simplexStack, int startingScanIndex);

void CheckForDomination(const vector< vector<double> > & points, const double & epsilon);

#endif
