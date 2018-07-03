
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
#include <sstream>
#include <map>

#include "cplex.h"

using namespace std;

const double EPSILON = 0.000001;

ostream & operator << (ostream & s, const vector<double> & pts) {
 for (auto x: pts) {
    s << x << " ";
 }
 return s;
}

// bad dan, make this  a class.
struct Point{
   vector<double> point;
   map<string, double> extra;

   bool operator < (const Point & other) {
      return point < other.point; 
   }

   bool operator == (const Point & other) {
       if (point.size() != other.point.size()) {
          return false;
       }
       for(size_t i =0;i<point.size();i++) {
          if (EPSILON < fabs(point[i]-other.point[i])) {
	    return false;
	  }
       }
       return true;
   }

};

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


Point ReadPoint(string line, char * fileName) {
    Point p;
    stringstream first, second;
    string junk, name;
    double val = 0.;

    size_t pos1 = line.find('[');
    size_t pos2 = line.find(']');

    // a tad bit of validation
    if (pos1 == string::npos or pos2 == string::npos or pos1 >= pos2) {
       cerr << "Bad line format in " << fileName << "  should include [ ] but is " << endl;
       cerr << "But is \"" << line << "\"" << endl;
       exit(1);
    }

    // extract the point
    first.str(line.substr(pos1+1,pos2-pos1-1));
    first >> val;
    while (first ) {
        p.point.push_back(val);
	first >> val;
    }

    // Oh PolySCIP, why can't you just put spaces around your = signs?
    for(auto & x: line) {
       if (x=='=') {
          x = ' ';
       }

    }

    // extract the other point values
    second.str(line.substr(pos2+1,string::npos));
    second >> name >> val;
    while (second) {
       p.extra[name] = val;
       second >> name >>  val;
    }
    return p;
}

void ReadPoints(char *  fileName, vector<Point> &  points){
    ifstream fin;	
    string line;
    Point p;

    fin.open(fileName);
    if (!fin) {
       cerr << "Failed to open " << fileName  << endl;
       exit(1);
    }
  	
    getline(fin, line);
    while(fin){
        p = ReadPoint(line, fileName);
	if (p.point.size() > 0) {
	    points.push_back(p);
	}
        getline(fin, line);
    }
	
    fin.close();
    return;
}

// returns the size of the difference between two sets of points.
int ShowDifference(vector<Point> a, vector<Point> b, char * fn1, char * fn2, bool showDiff){
    vector<Point> diff;
    set_difference(a.begin(),a.end(), b.begin(), b.end(), inserter(diff, diff.begin()));

    if (diff.size() > 0 and showDiff) {
       cout << "\tItems in " << fn1 << " but not in " << fn2 << ":" << endl;
       for(auto x: diff) {
          cout << "\t\t" << x.point << endl;
       }
       cout << endl;
    }

    return diff.size();
}

void CompareExtras(vector<double> point, map<string, double> a, map<string, double> b,
                 char * fn1, char * fn2) {

   map<string, double>::iterator i,j;

   bool header = false;
   for(i=a.begin(); i!= a.end(); i++) {
      j = b.find(i->first);
      if (not header) {
          cout << "\t\t" << point << " in " << fn1 << " has " << endl;
	  header = true;
      }

      if (j == b.end()) {
	  cout <<"\t\t\t" << i->first << " = " << i->second << " not in " << fn2 << endl;
      } else if (fabs(i->second-j->second) > EPSILON) {
	  cout <<"\t\t\t" << i->first << " = " << i->second << " but is " << j->second 
	       << " in "  << fn2  << endl;
      }
   }

   for(i=b.begin(); i!= b.end(); i++) {
      j = a.find(i->first);
      if (not header) {
          cout << "\t\t" << point << " in " << fn2 << " has " << endl;
	  header = true;
      }
      if (j == a.end()) {
	  cout <<"\t\t\t" << i->first << " = " << i->second << " not in " << fn1 << endl;
      }
   }

   return;
}

int ShowExtraDifference(vector<Point> a, vector<Point> b ,char * fn1, char * fn2){
   vector<Point> common;
   vector<Point>::iterator first, second;

   set_intersection(a.begin(), a.end(), b.begin(), b.end(),inserter(common, common.begin()));

   cout <<"\tThere are " << common.size() << " common points between the two"
       << endl << endl;
   for(auto x: common) {
      first = find(a.begin(), a.end(), x); 
      second = find(b.begin(), b.end(), x); 
      CompareExtras(x.point, first->extra, second->extra, fn1, fn2);
   }

   return 0;

}

// returns a bool if the two 
bool Compare(vector<Point> set1, vector<Point> set2, char * fn1, char * fn2, bool showDiff){

    int size1= 0, size2= 0;

    cout << "Comparing two input files ..." << endl << endl;

    if (showDiff) {
        cout << endl;
	cout << "\tPoints unique to each file " << endl;
    }

    sort(set1.begin(), set1.end());
    sort(set2.begin(), set2.end());
    size1 = ShowDifference(set1, set2, fn1, fn2, showDiff);
    size2 = ShowDifference(set2, set1, fn2, fn1, showDiff);

    cout << endl;
    cout << "\t" << fn1 << " contains " << set1.size() << " points." << endl;
    cout << "\t" << fn2 << " contains " << set2.size() << " points." << endl;
    cout << "\tThere are " << set1.size()-size1 << " common points." << endl;
    cout << endl;
   
    return (size1 == 0 and size2 == 0) ;
}

void CplexInit(void) {
    /* initialize Cplex environment *********************************/

    env = CPXopenCPLEX(&status);

    if(env==NULL) {
	char errmsg[1024];
	printf("CPXopenCPLEX, Couldn't open the CPLEX environment, error code %d\n", status);
    	CPXgeterrorstring (env, status, errmsg);
    	printf ("%s", errmsg);
    	exit(1);
    }
  	
    /************* Set any desired CPLEX parameters here **************/
  	
    status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_AbsGap, 0.0);
    if ( status ) {
	printf ("Failed to set solution pool gap to 0, error %d.\n",status);
	exit(1);
    }
    status = CPXsetdblparam (env, CPXPARAM_MIP_Pool_RelGap, 0.0);
    if ( status ) {
 	printf ("Failed to set solution pool gap to 0, error %d.\n",status);
	exit(1);
    }
  	
    /******************************************************************/
}

double FindDistance(vector<Point> a, vector<Point> b) {
    vector<double> lb;
    int ind = 0;
    int num = 0;
    double va = 0., min = CPX_INFBOUND, max = -CPX_INFBOUND;
    char l = 'L';
    vector<int> cmatbeg;
    vector<int> cmatind;
    vector<double> cmatval;

    size_t points = a.size();
    size_t dim = a[0].point.size();

    // need to find the diff between points in a and points in b.
    // store this in diff
    vector<Point> diff;
    set_difference(a.begin(),a.end(), b.begin(), b.end(), inserter(diff, diff.begin()));
    
    lp = CPXcreateprob (env,&status,"lp1.lp");
  	if(lp == NULL) {
     	    printf("CPXcreateprob, Failed to create LP%d, error code %d\n", 1, status);
      	    exit(1);
    }
    
    for(unsigned int i = 0; i < points; i++) lb.push_back(-CPX_INFBOUND);
    
    status = CPXaddcols (env, lp, a.size(), 0, NULL, NULL, NULL, NULL, &lb[0], NULL, NULL);
    status = CPXaddrows (env, lp, 0, dim + 1 + points, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
    for(unsigned int i = 0; i < points; i++)
    {
        for(unsigned int j = 0; j < dim + 1; j++)
        {
            if(j != dim) status = CPXchgcoef (env, lp, j, i, a[i].point[j]);
            else status = CPXchgcoef (env, lp, j, i, 1.);
        }
    }
    ind = dim; 
    va = 1.;
    status = CPXchgrhs (env, lp, 1, &ind, &va);
    
    num = dim + 1;
    for(unsigned int i = dim + 1; i < dim + 1 + points; i++)
    {
        ind = i;
        status = CPXchgsense (env, lp, 1, &ind, &l);
        status = CPXchgcoef (env, lp, i, i - num, -1.);
        //  cout << i << "\t" << i - num << "\t" << status << endl;
    }
    
    for(unsigned int i = 0; i < points; i++)
    {
        cmatbeg.push_back(i);
        cmatind.push_back(i + dim + 1);
        cmatval.push_back(1.);
    }
    status = CPXaddcols (env, lp, 1, points, &va, &cmatbeg[0], &cmatind[0], &cmatval[0], &lb[0], NULL, NULL);
    CPXchgobjsen (env, lp, CPX_MAX);
    
    status = CPXaddcols (env, lp, dim, dim, NULL, &cmatbeg[0], &cmatbeg[0], &cmatval[0], NULL, NULL, NULL);
 
    // Solve problem 1 for each point in list 2 which is not in list 1 *******
    

    for(unsigned int i = 0; i < diff.size(); i++)
    {
        lp1 = CPXcloneprob (env, lp, &status);
        //  cout << "i: " << i << endl;
        for(unsigned int j = 0; j < dim; j++)
        {
            status = CPXchgrhs (env, lp1, dim, &cmatbeg[0], &diff[i].point[0]);
            if(status) cout << status << endl;
        }
        
        status = CPXlpopt (env, lp1);
        if(status) cout << status << endl;
        status = CPXgetobjval (env, lp1, &va);
        if(status) cout << status << endl;
        if(va < min) min = va;
        if(va > max) max = va;

// nate had this commented out
/*        if(va < 0) 
        {
            cout << va << endl;
            for(unsigned int j = 0; j < dim; j++)
            {
               cout << diff.point[i][j] << " ";
            }
            cout << endl;
        }
*/
    }
    return min;
}

int main(int argc, char **argv)
{
    bool findDistances = true;
    bool compOtherPoints = true;
    bool showDiff = true;
    bool same = false;
    int i;

    if (argc < 3) {
        cout << argv[0] << " requires two file names" << endl;
	return 0;
    }

    i = 3; 
    while ( i < argc) {
       if (strcmp(argv[i], "-noDistance") == 0) {
           findDistances = false;
	   i++;
       } else if (strcmp(argv[i], "-noFull") == 0) {
          compOtherPoints = false;
	  i++;
       } else if (strcmp(argv[i], "-noDiff") == 0) {
          showDiff = false;
	  i++;
       } else {
          cout << "Unknown Argument " << argv[i] << endl;
	  i++;
       }
    }

    double min;
    vector<Point> set1, set2;

    ReadPoints(argv[1], set1);
    ReadPoints(argv[2], set2);

    same = Compare(set1, set2, argv[1], argv[2],showDiff);
    if ( compOtherPoints) {
       ShowExtraDifference(set1, set2, argv[1], argv[2]);
    }

    if (same and findDistances) {
       cout << "The points are the same, not checking distance between surfaces " << endl;
    }

    if (not same and findDistances) {
        CplexInit();
	// nate, check this to see if I have the verbage right.
	// I had changed everything before I realized what the output said.
        cout << "Determining distance from each point to " 
             << "the polyhedron specified by the extreme points in the other file"  << endl
	     << " (extreme directions included)." << endl;;

        // Create Problem 1 
        min =  FindDistance(set1, set2);
        cout << "\tFrom " << argv[1] << " to " << argv[2] << " : "
             << abs(min) << endl;

        //  Create Problem 2 *******************
        min = FindDistance(set2, set1);
        cout << "\tFrom " << argv[2] << " to " << argv[1] << " : "
             << abs(min) << endl;
    }

    cout << endl;
   
    return 0;
}

