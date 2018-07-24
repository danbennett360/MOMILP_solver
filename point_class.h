#ifndef POINTCLASS
#define POINTCLASS


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

using namespace std;

class Point
{
        vector<double> point;
        vector<int> nonzeroPrimalIndices;
        vector<double> nonzeroPrimalValues;
	    double pointCloseEpsilon;   // bennett for == operator

  public:
    	Point(const vector<double> & p, double *vals, int valsSize,
	                                      double epsilon=0.0000001);
    	void WritePointToFile(const string & filename, 
	                 const vector<string> & varNames, bool append);
        bool operator == (const Point & other);
	    double Epsilon(void) const;
	    void Epsilon(double e);
};

#endif
