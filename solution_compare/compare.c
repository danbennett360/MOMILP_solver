
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <time.h>
#include <pthread.h>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "cplex.h"

using namespace std;

/*******************************************************************/
/*              CPLEX VARIABLES                                    */
/*******************************************************************/

CPXENVptr  env=NULL;
CPXLPptr   lp = NULL;
CPXLPptr   lp1 = NULL;
CPXLPptr   lp2 = NULL;
vector<CPXLPptr>  lps;
int numObj = 0;
int status = 0;
int cur_numcols;
int cur_numrows;

int main(int argc, char **argv)
{

    vector< vector <double> > points1;
    vector< vector <double> > points2;
    vector< vector <double> > diff1;
    vector< vector <double> > diff2;
    vector<double> temp;
    vector<double> lb;
    ifstream fin;	
    double val = 0.;
    string s;
    int ind = 0;
    int num = 0;
    double va = 0., min = CPX_INFBOUND, max = -CPX_INFBOUND;
    char l = 'L';
    vector<int> cmatbeg;
    vector<int> cmatind;
    vector<double> cmatval;

	/* initialize Cplex environment *********************************/

  	env = CPXopenCPLEX(&status);

  	if(env==NULL)
  	{
    		char errmsg[1024];
    		printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
    		CPXgeterrorstring (env, status, errmsg);
    		printf ("%s", errmsg);
    		exit(0);
  	}
  	
  	/************* Set any desired CPLEX parameters here **************/
  	
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
  	
  	/******************************************************************/
  	
  	fin.open(argv[1]);
  	
  	while(!fin.eof())	
	{
		fin >> s;
/*		cout << s << endl;*/
		temp.resize(0);
		while(fin.good())
		{
		    fin >> val;
/*		    cout << val << endl;*/
		    if(fin.good()) temp.push_back(val);
		}
		if(temp.size()) points1.push_back(temp);
		fin.clear();
		fin.ignore(100000,'\n');
	}
	
	fin.close();
	
	fin.open(argv[2]);
  	
  	while(!fin.eof())	
	{
		fin >> s;
/*		cout << s << endl;*/
		temp.resize(0);
		while(fin.good())
		{
		    fin >> val;
/*		    cout << val << endl;*/
		    if(fin.good()) temp.push_back(val);
		}
		if(temp.size()) points2.push_back(temp);
		fin.clear();
		fin.ignore(100000,'\n');
	}
	
	fin.close();
	
	cout << "Comparing two input files ...\n\n\n";
	cout << "The first input file contains " << points1.size() << " points." << endl;
	cout << "The second input file contains " << points2.size() << " points." << endl;
	
	sort(points1.begin(), points1.end());
	sort(points2.begin(), points2.end());
	
	memcpy ( &diff1, &points1, sizeof(points2) );
	
	diff1.resize(set_difference(points1.begin(), points1.end(), points2.begin(), points2.end(), diff1.begin()) - diff1.begin());
	
	memcpy ( &diff2, &points2, sizeof(points2) );
	
	diff2.resize(set_difference(points2.begin(), points2.end(), points1.begin(), points1.end(), diff2.begin()) - diff2.begin());
	
	cout << "There are " << diff1.size() << " points in the first file which are not found in the second file." << endl;
	cout << "There are " << diff2.size() << " points in the second file which are not found in the first file." << endl;
	
  	
  	/** Create Problem 1 *******************/

  	lp = CPXcreateprob (env,&status,"lp1.lp");
  	if(lp == NULL) 
  	{
    		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", 1, status);
    		exit(0);
    }
    
    for(unsigned int i = 0; i < points1.size(); i++) lb.push_back(-CPX_INFBOUND);
    
    status = CPXaddcols (env, lp, points1.size(), 0, NULL, NULL, NULL, NULL, &lb[0], NULL, NULL);
    status = CPXaddrows (env, lp, 0, points1[0].size() + 1 + points1.size(), 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
    for(unsigned int i = 0; i < points1.size(); i++)
    {
        for(unsigned int j = 0; j < points1[0].size() + 1; j++)
        {
            if(j != points1[0].size()) status = CPXchgcoef (env, lp, j, i, points1[i][j]);
            else status = CPXchgcoef (env, lp, j, i, 1.);
        }
    }
    ind = points1[0].size();
    va = 1.;
    status = CPXchgrhs (env, lp, 1, &ind, &va);
    
    num = points1[0].size() + 1;
    for(unsigned int i = points1[0].size() + 1; i < points1[0].size() + 1 + points1.size(); i++)
    {
        ind = i;
        status = CPXchgsense (env, lp, 1, &ind, &l);
        status = CPXchgcoef (env, lp, i, i - num, -1.);
/*        cout << i << "\t" << i - num << "\t" << status << endl;*/
    }
    
    for(unsigned int i = 0; i < points1.size(); i++)
    {
        cmatbeg.push_back(i);
        cmatind.push_back(i + points1[0].size() + 1);
        cmatval.push_back(1.);
    }
    status = CPXaddcols (env, lp, 1, points1.size(), &va, &cmatbeg[0], &cmatind[0], &cmatval[0], &lb[0], NULL, NULL);
    CPXchgobjsen (env, lp, CPX_MAX);
    
    status = CPXaddcols (env, lp, points1[0].size(), points1[0].size(), NULL, &cmatbeg[0], &cmatbeg[0], &cmatval[0], NULL, NULL, NULL);
 
    /* Solve problem 1 for each point in list 2 which is not in list 1 *******/
    
    cout << "\nDetermining distance from each point in file 2 to the polyhedron specified by the extreme points present in file 1 (extreme directions included).\n";

    for(unsigned int i = 0; i < diff2.size(); i++)
    {
        lp1 = CPXcloneprob (env, lp, &status);
/*        cout << "i: " << i << endl;*/
        for(unsigned int j = 0; j < diff2[0].size(); j++)
        {
            status = CPXchgrhs (env, lp1, diff2[0].size(), &cmatbeg[0], &diff2[i][0]);
            if(status) cout << status << endl;
        }
        
        status = CPXlpopt (env, lp1);
        if(status) cout << status << endl;
        status = CPXgetobjval (env, lp1, &va);
        if(status) cout << status << endl;
        if(va < min) min = va;
        if(va > max) max = va;
/*        if(va < 0) */
/*        {*/
/*            cout << va << endl;*/
/*            for(unsigned int j = 0; j < diff2[0].size(); j++)*/
/*            {*/
/*               cout << diff2[i][j] << " ";*/
/*            }*/
/*            cout << endl;*/
/*        }*/
    }
    
    cout << "Maximum distance: " << abs(min) << endl;
    
    /** Create Problem 2 *******************/

  	lp2 = CPXcreateprob (env,&status,"lp2.lp");
  	if(lp2 == NULL) 
  	{
    		printf("CPXcreateprob, Failed to create LP%d, error code %d\n", 1, status);
    		exit(0);
    }
    
    for(unsigned int i = 0; i < points2.size(); i++) lb.push_back(-CPX_INFBOUND);
    
    status = CPXaddcols (env, lp2, points2.size(), 0, NULL, NULL, NULL, NULL, &lb[0], NULL, NULL);
    status = CPXaddrows (env, lp2, 0, points2[0].size() + 1 + points2.size(), 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
    for(unsigned int i = 0; i < points2.size(); i++)
    {
        for(unsigned int j = 0; j < points2[0].size() + 1; j++)
        {
            if(j != points2[0].size()) status = CPXchgcoef (env, lp2, j, i, points2[i][j]);
            else status = CPXchgcoef (env, lp2, j, i, 1.);
        }
    }
    ind = points2[0].size();
    va = 1.;
    status = CPXchgrhs (env, lp2, 1, &ind, &va);
    
    num = points2[0].size() + 1;
    for(unsigned int i = points2[0].size() + 1; i < points2[0].size() + 1 + points2.size(); i++)
    {
        ind = i;
        status = CPXchgsense (env, lp2, 1, &ind, &l);
        status = CPXchgcoef (env, lp2, i, i - num, -1.);
/*        cout << i << "\t" << i - num << "\t" << status << endl;*/
    }
    
    for(unsigned int i = 0; i < points2.size(); i++)
    {
        cmatbeg.push_back(i);
        cmatind.push_back(i + points2[0].size() + 1);
        cmatval.push_back(1.);
    }
    status = CPXaddcols (env, lp2, 1, points2.size(), &va, &cmatbeg[0], &cmatind[0], &cmatval[0], &lb[0], NULL, NULL);
    CPXchgobjsen (env, lp2, CPX_MAX);
    
    status = CPXaddcols (env, lp2, points2[0].size(), points2[0].size(), NULL, &cmatbeg[0], &cmatbeg[0], &cmatval[0], NULL, NULL, NULL);
    
    cout << "\nDetermining distance from each point in file 1 to the polyhedron specified by the extreme points present in file 2 (extreme directions included).\n";
    min = CPX_INFBOUND;
    max = -CPX_INFBOUND;
    /* Solve problem 2 for each point in list 1 which is not in list 2 *******/

    for(unsigned int i = 0; i < diff1.size(); i++)
    {
        lp1 = CPXcloneprob (env, lp2, &status);
/*        cout << "i2: " << i << endl;*/
        for(unsigned int j = 0; j < diff1[0].size(); j++)
        {
            status = CPXchgrhs (env, lp1, diff1[0].size(), &cmatbeg[0], &diff1[i][0]);
            if(status) cout << status << endl;
        }
        
        status = CPXlpopt (env, lp1);
        if(status) cout << status << endl;
        status = CPXgetobjval (env, lp1, &va);
        if(status) cout << status << endl;
        if(va < min) min = va;
        if(va > max) max = va;
/*        if(va < 0) */
/*        {*/
/*            cout << va << endl;*/
/*            for(unsigned int j = 0; j < diff1[0].size(); j++)*/
/*            {*/
/*               cout << diff1[i][j] << " ";*/
/*            }*/
/*            cout << endl;*/
/*        }*/
    }
    
    cout << "Maximum distance: " << abs(min) << endl;

    exit(0);
  	return 0;
}

