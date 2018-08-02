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

#include "point_class.h"

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
    	vector<double> normal;
    	vector<bool> extremesAreDummies;
    	double planeVal;
    	bool isPositive = true;

        int numDummies = 0;
        bool showNormalsInPlots = false;
        bool saveForSolution = false;
       
        // changed to static so they are shared across all instances
        // static, make sure you don't copy these in a copy constructor
        // or asignment operator
	static double epsilon;
	static double infinity;

    	void GenerateNormal();

	// bennett new
        vector<Simplex * > adj;
        int id;
        vector<int> points;
        bool split = false;

        // static, make sure you don't copy these in a copy constructor
        // or asignment operator
        static vector<Point> * pointVector;
        static int nextID;

        // bennett kill
        //vector< vector<double> > extremePoints;

        // new private member functions
        Simplex(int i);
        Simplex() : Simplex(2) {};
        Simplex(const Simplex * other);

	friend class SimplexStore;
  public:
	// bennett 7/18
	double Epsilon(void) const;
	void Epsilon(double e);

    	//void AddExtreme(const vector<double> & v, bool normalize);
	void AddExtreme(int , bool normalize);
    	void PrintData() const;
    	void NormalizeNormal();
    	void Reset();
    	vector<double> GetNormal() const;
    	double MultiplyPointsByNormal(const vector<double> & v) const;
    	double PlaneVal();
    	vector< vector<double> > GetExtremePoints() const;
    	vector<double> GetExtremePoint(int i) const;
	void WriteOctaveCodeToPlotSimplex(bool a) const;
	bool IsPositive() const;
    	Simplex * FindAdjacentContainingOriginalPoints(const vector<Simplex *> & simplices, int & simplexIndex, int & newPointIndex, int ignoreIndex, bool isIn) const;
    	int GetDimension() const;
    	int GetNumDummyPoints();
    	vector<bool> ExtremesAreDummies() const;
    	int GetPointIndex(const vector<double> & point);
    	void SaveForSolution();
    	bool GetSaveForSolution() const;

    	void PrintNormal();
    	bool isInSubsetOfStack(const vector<Simplex * > & simplexStack, int startingScanIndex);
    	bool deleteRepeats(vector<Simplex*> & simplexStack, int startingScanIndex);

	// bennett new
	void Adjacent(int i, Simplex * a);
	void Adjacent(Simplex * b);
	Simplex * Adjacent(int i) const;

	void RemoveAdjacent(Simplex * s);
	void Points(vector<int> p, bool normalize);
	void MarkForSplit(void);
	bool IsMarkedForSplit(void) const;

	int ID(void) const;

	vector<Point> * GetPointVector(void) const;
	void SetPointVector(vector<Point> * pv);
};

double Determinant(const vector< vector<double> > & matrix);

#endif
