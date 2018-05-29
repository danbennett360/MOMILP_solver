#ifndef SYDNEYSCLASS
#define SYDNEYSCLASS


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
#include "simplex_class.h"

using namespace std;

class SydneysClass
{
        vector<Simplex> data;
  public:
    	void Insert(Simplex s);
};

#endif
