
#include <iostream>
#include <iomanip>
#include <fstream>
#include "point_class.h"

Point::Point(const vector<double> & p, double *vals, int valsSize)
{
    point = p;
    
    for(int i = 0; i < valsSize; i++)
    {
        if(vals[i] != 0.) 
        {
            nonzeroPrimalIndices.push_back(i);
            nonzeroPrimalValues.push_back(vals[i]);
        }
    }
}

void Point::WritePointToFile(const string & filename, const vector<string> & varNames, bool append)
{
    ofstream fout;	
	if(append) fout.open(filename.c_str(), ofstream::app);
	else fout.open(filename.c_str(), ofstream::trunc);
	
	fout << "[ ";
	for(unsigned int i = 0; i < point.size(); i++) fout << point[i] << " ";
	fout << "] ";
	for(unsigned int i = 0; i < nonzeroPrimalIndices.size(); i++) fout << varNames[nonzeroPrimalIndices[i]] << " = " << nonzeroPrimalValues[i] << " ";
	fout << endl;
	
	fout.close();

    return;
}
