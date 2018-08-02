#ifndef SIMPLEX_STORE
#define SIMPLEX_STORE

#include "point_class.h"

#ifdef CPLEX
    #include "cplex.h"
#else
    #include <glpk.h>
#endif

#include "simplex_class.h"

#include <vector>
#include <list>

class SimplexStore{
    public:
       SimplexStore();
       
       SimplexStore(SimplexStore & s) = delete;
       SimplexStore operator = (SimplexStore & o) = delete;

       ~SimplexStore();

       int AddPoint(Point p);
       Point GetPoint(int i)const;
       int PointCount(void) const;

       Simplex* AddSimplex(vector<int> Points, bool normalize);
       int SimplexCount(void) const;
       vector<Simplex *> SplitSimplex(vector<Simplex *> s, int point, bool normalize);

       void RemoveSimplex(Simplex * s);

    private:

       void DoSplit(Simplex * simplex, vector<Simplex * > & keepList,
                   vector<Simplex * > & mergeList, int point, bool normalize);
       void MakeNewSimplicies(Simplex * simplex, vector<Simplex * > & newList,
                 int point, bool normalize);
       void MakeAdjacent(vector<Simplex *> newList);
       bool HasSplitAdj(Simplex * simplex);

       std::vector<Point> points;
       std::list<Simplex *> simplicies;
};

#endif
