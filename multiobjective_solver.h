#ifndef MULTIOBJECTIVE
#define MULTIOBJECTIVE

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
#include <fstream>
#include <regex>

using namespace std;

//CPXENVptr  env=NULL;
//CPXLPptr   lp=NULL;

extern bool DEBUG;
extern bool SCAN_FOR_REPEATS;
extern bool SAVE_POINTS;

#endif
