#ifndef SIMPLEXCLASS
#define SIMPLEXCLASS

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
//#include <unordered_map>
//#include <unordered_set>
#include <map>
#include <set>
#include <cmath>
#include "multiobjective_solver.h"

using namespace std;

bool CompareDoubleVectors(const vector<double> & v1, const vector<double> & v2);

struct DblVecCompare 
{
    bool operator() (const vector<double> & v1, const vector<double> & v2) const 
    {
        return CompareDoubleVectors(v1, v2);
    }
};

class Simplex
{
        int dimension = -1;
    	vector< vector<double> > extremePoints;
    	vector<double> normal;
    	vector<bool> extremesAreDummies;
    	void GenerateNormal();
    	double planeVal;
    	bool isPositive = true;
        double epsilon = .0000001;
        double infinity = CPX_INFBOUND;
        int numDummies = 0;
        bool showNormalsInPlots = false;
        bool saveForSolution = false;
  public:
    	Simplex(int i);
    	Simplex() : Simplex(2) {};

	// bennett 7/18
	double Epsilon(void) const;
	void Epsilon(double e);

    	void AddExtreme(const vector<double> & v, bool normalize);
    	void PrintData() const;
    	void NormalizeNormal();
    	void Reset();
    	vector<double> GetNormal() const;
    	double MultiplyPointsByNormal(const vector<double> & v) const;
    	double PlaneVal();
    	vector< vector<double> > GetExtremePoints() const;
    	vector<double> GetExtremePoint(int i) const;
    	void WriteOctaveCodeToPlotSimplex(bool a) const;
    	bool IsPositive();
    	Simplex FindAdjacentContainingOriginalPoints(const vector<Simplex> & simplices, int & simplexIndex, int & newPointIndex, int ignoreIndex, bool isIn) const;
    	int GetDimension() const;
    	int GetNumDummyPoints();
    	vector<bool> ExtremesAreDummies() const;
    	int GetPointIndex(const vector<double> & point);
    	void SaveForSolution();
    	bool GetSaveForSolution() const;
    	void PrintNormal();
    	bool isInSubsetOfStack(const vector<Simplex> & simplexStack, int startingScanIndex);
    	bool deleteRepeats(vector<Simplex> & simplexStack, int startingScanIndex);
};

double Determinant(const vector< vector<double> > & matrix);

#endif
