#include <iostream>

#include "SimplexStore.h"
#include "point_class.h"
#include "simplex_class.h"
#include <algorithm>

using namespace std;

SimplexStore::SimplexStore(){
    return;
}

SimplexStore::~SimplexStore(){
    Simplex * s;
    while (simplicies.size() > 0) {
       s = simplicies.back();
       simplicies.pop_back();
       delete s;
    }
    return;
}

int SimplexStore::AddPoint(Point p){
    points.push_back(p);
    return points.size()-1;
}

Point SimplexStore::GetPoint(int i)const{
    if (size_t(i) < points.size()) {
       return points[i];
    } else {
       cerr << "Error, point list has " << points.size() << " elements"
            << " but an attempt to access element " << i << " occured" << endl;
       exit(1);
    }
    return points[0];
}

int SimplexStore::PointCount(void) const{
    return points.size();
}


// simply add a new simplex to the store.
// Assume the points are already in the store, but we could check this if needed
// Assume that this is an isolated simplex, so there are no adjacents.

Simplex* SimplexStore::AddSimplex(vector<int> simplexPoints, bool normalize){
    Simplex * s;

    s = new Simplex(simplexPoints.size());

    // only need to do this once  but do it before anything else.
    if (nullptr == s->GetPointVector()){
        cerr << endl;
        cerr << "TODO:, make sure normalize is right SimplexStore::AddSymplex"
	     << endl;
        cerr << endl;
        s->SetPointVector(&points);
    }

    size_t i;
    //size_t p;

    for (i =0; i<simplexPoints.size()-1;i++) {
        //p = simplexPoints[i];
        s->AddExtreme(simplexPoints[i], normalize);
	s->Adjacent(i, nullptr);
    }

    //s->AddExtreme(simplexPoints[i], true);
    s->AddExtreme(simplexPoints[i], false);
    s->Adjacent(i, nullptr);

    simplicies.push_back(s);


    return s;
}

int SimplexStore::SimplexCount(void) const{
    return simplicies.size();
}

// for each point in the simplex
//      create a new simplex by replacing that point with the new point
void SimplexStore::MakeNewSimplicies(Simplex * simplex, 
                   vector<Simplex * > & newList, int point, bool normalize){
    Simplex * s;
    int pos;

    // I'm a friend so I can mess with the internal structure.
    for(int i=0;i<simplex->dimension;i++) {
       s = new Simplex(simplex);
       s->points[i] = point;
       s->GenerateNormal();
       if(normalize) {
          s->NormalizeNormal();
       }
       s->planeVal = s->MultiplyPointsByNormal(s->pointVector->at(s->points[0]).Data());

       sort(s->points.begin(),s->points.end());

       // damn it, I sorted it so it is possibly no longer at position i, in 
       // fact is is probably at position dimension-1, but who knows.
       pos = simplex->dimension-1;
       while (pos > 0 and s->points[pos] != point) {
	     pos--;
       }
       

       s->split = false;
       for(int j=0;j<int(s->adj.size());j++) {
          if (j != pos) {
              s->adj[j] = nullptr;
	  } else {
              // we are still adjacent to our parent's adjacent on this edge
	      // but our parent is no longer adjacet to this one.
	      Simplex * neighbor = simplex->adj[i];
	      if (neighbor != nullptr) {
		  simplex->adj[i] = nullptr;
		  // mark child as adjacent
		  neighbor->Adjacent(s);
	      }
	  }
       }
       newList.push_back(s);

    }

    return;
}

void SimplexStore::MakeAdjacent(vector<Simplex *> newList){
    int dim = newList[0]->dimension;

    int i,j;
    Simplex * s;

    // for each new simplex
    for(i=0;i<dim;i++) { 
       s = newList[i];
       for (j=i+1;j<dim;j++) {
           s->Adjacent(newList[j]);
       }
    }

    return;
}

bool SimplexStore::HasSplitAdj(Simplex * simplex) {
    for (auto x: simplex->adj) {
       if (x != nullptr and x->IsMarkedForSplit()) {
          return true;
       }
    }
    return false;
}

void SimplexStore::DoSplit(Simplex * simplex, vector<Simplex * > & keepList,
                   vector<Simplex * > & mergeList, int point, bool normalize) {

     vector<Simplex *> newList;

     // generate the new siplicies  from this simplex and the new point.
     MakeNewSimplicies(simplex, newList, point, normalize);

     // fix the adjacencies for each new simplex
     MakeAdjacent(newList);
     
     // put in proper list
     for(size_t i=0;i<newList.size();i++) {
        if (HasSplitAdj(newList[i])) {
	    mergeList.push_back(newList[i]);
	} else {
	    keepList.push_back(newList[i]);
	}
     }
     return;
}

vector<Simplex *> SimplexStore::SplitSimplex(vector<Simplex *> sList, int point,
                       bool normalize){
    vector<Simplex * > keepList, mergeList;
    Simplex * a, *b=nullptr, *aPrime,* bPrime;
    size_t i;
    bool found;

    // generate all new simplicies
    for (auto  simplex: sList) {
        DoSplit(simplex, keepList, mergeList, point,normalize);
    }

    // fix items in mergeList
    while (mergeList.size() > 0) {
        a = mergeList.back();
	mergeList.pop_back();

        // find the mirror  of A in the merge list B 
        i = 0;
	found = false;
	while (i < mergeList.size() and not found) {
	    if (a->points == mergeList[i] ->points) {
	       b = mergeList[i];
	       mergeList[i] = mergeList.back();
	       mergeList.pop_back();
	       found = true;
	    }
	    i++;
	}

	if (found) {
	    for (i=0;i<a->adj.size();i++) {
	        if (a->adj[i] != nullptr) {
		  aPrime = a->adj[i];
		  if (not aPrime->IsMarkedForSplit()) {
		      bPrime = b->adj[i];
		      aPrime->Adjacent(bPrime);
		  }
		}
	    }
	} else {
	    // probably should have an error mesage here.
	}

	RemoveSimplex(a);
	RemoveSimplex(b);

    }
    // now go through and delete all simplies in the sList
    for (auto simplex: sList) {
        RemoveSimplex(simplex);
    }

    return keepList;
}

// remove the simplex from the dlist
// Remove any adjacent pointers to it 
//      (assumes adjacents are pointing to each other)
// delete the simplex.
void SimplexStore::RemoveSimplex(Simplex * s) {
    list<Simplex *>::iterator i;
    Simplex * a;

    i = find(simplicies.begin(),simplicies.end(), s);
    if (i != simplicies.end()) {
       simplicies.erase(i);

       // deleting any entries in an adjacent simplex list
       for(int j=0;j<s->GetDimension();j++) {
          a = s->Adjacent(j);
	  if (a != nullptr) {
	     a->RemoveAdjacent(s);
	  }
       }
       delete s; 
    }

    return;
}
