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
bool TERSE = false;

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
      // bennett 7/5/18
      // set_intersection is based on < for =, man page says  !(a<b) and !(b<a)
      if (*this==other) {  // if they are within epsilon, they are the same.
        return false;
      }
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


// read a point from the INNER/BenSolve output file format
Point ReadPointIB(string line, char validC) {
    Point p;
    stringstream data;
    char type;
    double num;

    data.str(line); 

    data >> type;

    // make sure that it is a valid point.
    if (type == validC) {
       data >> num;
       while (data) {
           p.point.push_back(num);
           data >> num;
       }
    }

    return p;
}

// read a point in the Nate and PolyScip Output File Format
Point ReadPointN(string line, char * fileName) {
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

enum  FileType {NATE, INNER, BENSOLVE,  UNKNOWN};

void ReadPoints(char *  fileName, vector<Point> &  points){
    ifstream fin;	
    string line;
    FileType type = UNKNOWN;
    Point p;

    fin.open(fileName);
    if (!fin) {
       cerr << "Failed to open " << fileName  << endl;
       exit(1);
    }
  
    line = "";
    getline(fin, line);
    while (line.size() == 0) {
        getline(fin, line);
    }

    if (line[0] == '[') { 
       type = NATE; 
    } else if (line[0] == 'C' or line[0] == 'v' or line[0] == 'V') {
       type = INNER ;
    } else if (line[0] == '0' or line[0] == '1') {
       type = BENSOLVE ;
    }  else {
       cerr << "Unknown format for file " << fileName << endl;
       exit (1);
    }

    while(fin){
        switch (type) {
	   case NATE:
               p = ReadPointN(line, fileName);
	       break;
	   case INNER:
               p = ReadPointIB(line,'V');
	       break;
	   case BENSOLVE:
               p = ReadPointIB(line,'1');
	       break;
	   case UNKNOWN:
	   default:
	       // shouldn't make it here but ...
	       cerr << "Error reading input file." << endl;
	       exit (0);
	       break;
	}
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

int ShowExtraDifference(vector<Point> a, vector<Point> b, char * fn1, char * fn2){
   vector<Point> common;
   vector<Point>::iterator first, second;

   set_intersection(a.begin(), a.end(), b.begin(), b.end(),inserter(common, common.begin()));

   if (TERSE) {
       cout << "," << common.size();
   } else {
       cout <<"\tThere are " << common.size() << " common points between the two"
           << endl << endl;
   }
   for(auto x: common) {
      first = find(a.begin(), a.end(), x); 
      second = find(b.begin(), b.end(), x); 
      if (not TERSE) {
          CompareExtras(x.point, first->extra, second->extra, fn1, fn2);
      }
   }

   return 0;

}

// returns a bool if the two 
bool Compare(vector<Point> set1, vector<Point> set2, char * fn1, char * fn2, bool showDiff){

    int size1= 0, size2= 0;

    if (not TERSE) {
        cout << "Comparing two input files ..." << endl << endl;
    }

    if (showDiff) {
        cout << endl;
	cout << "\tPoints unique to each file " << endl;
    }

    size1 = ShowDifference(set1, set2, fn1, fn2, showDiff);
    size2 = ShowDifference(set2, set1, fn2, fn1, showDiff);

    if (not TERSE) {
        cout << endl;
        cout << "\t" << fn1 << " contains " << set1.size() << " points." << endl;
        cout << "\t" << fn2 << " contains " << set2.size() << " points." << endl;
        cout << "\tThere are " << set1.size()-size1 << " common points." << endl;
        cout << endl;
    } else {
       cout << set1.size() <<","<< set2.size() << "," << size1 << "," << size2 << ",";
    }
   
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
  	
    status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
	if ( status ) {
		printf ("Failed to turn screen printing on, error %d.\n",status);
		exit(1);
	}
  	
    /******************************************************************/
}

double FindDistance(vector<Point> a, vector<Point> b, bool relDist, double & minMagnitude) {
    vector<double> lb;
    int ind = 0;
    int num = 0;
    double va = 0., min = CPX_INFBOUND, max = -CPX_INFBOUND;
/*    char l = 'L';*/
/*    double d[1] = {0.};*/
    vector<int> cmatbeg;
    vector<int> rmatbeg;
    vector<int> cmatind;
    vector<double> cmatval;
    vector<char> sense;

    size_t points = a.size();
    size_t dim = a[0].point.size();

    // need to find the diff between points in a and points in b.
    // store this in diff
    vector<Point> diff;
    set_difference(b.begin(),b.end(), a.begin(), a.end(), inserter(diff, diff.begin()));
    if (diff.size() == 0 ) {
       if (not TERSE) {
           cout << "\tThe points are the same (or a subset)" << endl;
       }
       return 0;
    }
    
/*    for(unsigned int i = 0; i < diff.size(); i++)*/
/*    {*/
/*        for(unsigned int j = 0; j < diff[i].point.size(); j++)*/
/*        {*/
/*            cout << diff[i].point[j] << "\t";*/
/*        }*/
/*        cout << endl;*/
/*    } */
    
    lp = CPXcreateprob (env,&status,"lp1.lp");
    if(lp == NULL) {
     	    printf("CPXcreateprob, Failed to create LP%d, error code %d\n", 1, status);
      	    exit(1);
    }
    
    for(unsigned int i = 0; i < points; i++) lb.push_back(-CPX_INFBOUND);
    
    status = CPXaddcols (env, lp, a.size(), 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    status = CPXaddrows (env, lp, 0, dim + 1, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    
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
    
/*    num = dim + 1;*/
/*    for(unsigned int i = dim + 1; i < dim + 1 + points; i++)*/
/*    {*/
/*        ind = i;*/
/*        status = CPXchgsense (env, lp, 1, &ind, &l);*/
/*        status = CPXchgcoef (env, lp, i, i - num, -1.);*/
/*        //  cout << i << "\t" << i - num << "\t" << status << endl;*/
/*    }*/
    
    for(unsigned int i = 0; i < points; i++)
    {
        cmatbeg.push_back(i);
        cmatind.push_back(i + dim + 1);
        cmatval.push_back(1.);
    }
/*    status = CPXaddcols (env, lp, 1, points, &va, &cmatbeg[0], &cmatind[0], &cmatval[0], NULL, NULL, NULL);*/
    CPXchgobjsen (env, lp, CPX_MIN);
    
    status = CPXaddcols (env, lp, dim, dim, NULL, &cmatbeg[0], &cmatbeg[0], &cmatval[0], NULL, NULL, NULL);
    
    for(unsigned int i = 0; i < cmatval.size(); i++) cmatval[i] = -1.;
    
    status = CPXaddcols (env, lp, dim, dim, NULL, &cmatbeg[0], &cmatbeg[0], &cmatval[0], &lb[0], NULL, NULL);
    
    status = CPXchgprobtype (env, lp, CPXPROB_QP);
    
    num = CPXgetnumcols (env, lp);
    
    for(unsigned int i = 0; i < dim; i++) status = CPXchgqpcoef (env, lp, num - i - 1, num - i - 1, 2.);
    
 
    // Solve problem 1 for each point in list 2 which is not in list 1 *******
    
    cmatind.resize(dim);
    cmatval.resize(dim);

    for(unsigned int i = 0; i < diff.size(); i++)
    {
        
        if(relDist)
        {
            va = 0.;
            for(unsigned int j = 0; j < diff[i].point.size(); j++) va += diff[i].point[j]*diff[i].point[j];
            va = sqrt(va);
            if(va < minMagnitude) minMagnitude = va;
        }
    
        lp1 = CPXcloneprob (env, lp, &status);
/*        cout << "i: " << i << endl;*/
/*        for(unsigned int j = 0; j < dim; j++)*/
/*        {*/
/*            status = CPXchgrhs (env, lp1, dim, &cmatbeg[0], &diff[i].point[0]);*/
/*            if(status) cout << status << endl;*/
/*        }*/
        
/*        if(i == 0)*/
/*        {*/
/*            status = CPXwriteprob (env, lp1, "lp1.lp", "LP");*/
/*            cout << status << endl;*/
/*            exit(0);*/
/*        }*/

        for(unsigned int j = 0; j < dim; j++) 
        {
            cmatind[j] = num - dim + j;
            cmatval[j] = -2*diff[i].point[j];
        }

        status = CPXchgobj (env, lp1, dim, cmatind.data(), cmatval.data());
        
        status = CPXqpopt (env, lp1);
        if(status) 
        {
            cout << "Unable to solve problem. CPLEX error code: " << status << ". Exiting." << endl;
            exit(1);
        }
        status = CPXgetobjval (env, lp1, &va);
        if(status) 
        {
            cout << "Unable to get objective value. CPLEX error code: " << status << ". Exiting." << endl;
            exit(1);
        }
        
/*        for(unsigned int j = 0; j < num; j++) */
/*        {*/
/*            status = CPXgetx (env, lp1, d, j, j);*/
/*            cout << "x" << j+1 << ": " << d[0] << endl;*/
/*        }*/
        
        for(unsigned int j = 0; j < dim; j++) va += diff[i].point[j]*diff[i].point[j];
        
        va = sqrt(va);
        
        if(!status)
        {
            if(va < min) min = va;
            if(va > max) max = va;
/*            cout << va << endl;*/
/*            cout << "min: " << min << "\tmax: " << max << endl;*/
        }
        
/*        status = CPXwriteprob (env, lp1, "lp1.lp", "LP");*/
/*        exit(0);*/

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
    return max;
}

int main(int argc, char **argv)
{
    bool findDistances = true;
    bool compOtherPoints = true;
    bool showDiff = true;
    bool same = false;
    bool relDist = false;
    int i;
    double minMagnitude = CPX_INFBOUND;

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
       } else if (strcmp(argv[i], "-Terse") == 0) {
          TERSE = true;
	  i++;
       } else if (strcmp(argv[i], "-RelDist") == 0) {
          relDist = true;
	  i++;
       } else {
          cout << "Unknown Argument " << argv[i] << endl;
	  i++;
       }
    }

    if (TERSE) {
       showDiff = false;
    }

    double min;
    vector<Point> set1, set2;

    ReadPoints(argv[1], set1);
    ReadPoints(argv[2], set2);

    sort(set1.begin(), set1.end());
    sort(set2.begin(), set2.end());


    same = Compare(set1, set2, argv[1], argv[2],showDiff);
    if ( compOtherPoints) {
       ShowExtraDifference(set1, set2, argv[1], argv[2]);
    }

    if (same and findDistances) {
       if (not TERSE) {
           cout << "The points are the same, not checking distance between surfaces " << endl;
       } else {
           cout << "0,0" << endl;
	   return  0;
       }
    }

    if (not same and findDistances) {
        CplexInit();
	// nate, check this to see if I have the verbage right.
	// I had changed everything before I realized what the output said.
	if (not TERSE) {
        cout << "Determining distance from each point to " 
             << "the polyhedron specified by the extreme points in the other file"  << endl
	     << " (extreme directions included)." << endl;;
        }
        // Create Problem 1 
        min =  FindDistance(set1, set2, relDist, minMagnitude);
	if (not TERSE) {
            cout << "\tFrom " << argv[1] << " to " << argv[2] << " : "
                 << abs(min) << endl;
            if(relDist) cout << "\t\tAs a relative distance: " << abs(min)/minMagnitude << endl;
	} else {
	    cout << abs(min) << ",";
	}

        //  Create Problem 2 *******************
        min = FindDistance(set2, set1, false, minMagnitude);
	if (not TERSE) {
           cout << "\tFrom " << argv[2] << " to " << argv[1] << " : "
                 << abs(min) << endl;
           if(relDist) cout << "\t\tAs a relative distance: " << abs(min)/minMagnitude << endl;
        } else {
	    cout << abs(min) << endl;
	}
    }

    cout << endl;
   
    return 0;
}

