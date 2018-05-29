
#include "cplex.h"
#include "sydneys_class.h"
#include "simplex_class.h"

void SydneysClass::Insert(Simplex s)
{
    data.push_back(s);
    
    return;
}
