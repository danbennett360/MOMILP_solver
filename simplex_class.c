/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#ifdef CPLEX
    #include "cplex.h"
#else
    #include <glpk.h>
#endif

#include "simplex_class.h"

int Simplex::nextID = 0;
double Simplex::epsilon = .0000001;
#ifdef CPLEX
    double Simplex::infinity = CPX_INFBOUND;
#else
    double Simplex::infinity = 100000000000000000000.;
#endif

// need a reference to the point vector  so we can
// gain access to our points.
vector<Point> * Simplex::pointVector = nullptr;

vector<Point> * Simplex::GetPointVector(void) const {
   return pointVector;
}

void Simplex::SetPointVector(vector<Point> * pv) {
   pointVector = pv;
   return;
}

Simplex::Simplex(int i)
{
    dimension = i;
    epsilon = EPSILON;
    id = nextID;
    nextID++;
}

int Simplex::ID(void) const {
   return id;
}

Simplex::Simplex(const Simplex * other) {
    id = nextID;
    nextID++;

    dimension = other->dimension;
    normal = other->normal;
    extremesAreDummies = other->extremesAreDummies;
    planeVal = other->planeVal;
    isPositive = other->isPositive;
    numDummies = other->numDummies;
    showNormalsInPlots = other->showNormalsInPlots;
    saveForSolution = other->saveForSolution;
    adj = other->adj;
    points = other->points;
    split = other->split;

    return;
}

void Simplex::AddExtreme(int pointPos, bool normalize) {
   bool isDummy = false;
   int i = 0;

   if(int(points.size()) < dimension)
    {
        points.push_back(pointPos);
        while(i < dimension && !isDummy)
        {
            if(pointVector->at(pointPos).Data()[i] >= infinity) isDummy = true;
            i++;
        }
        extremesAreDummies.push_back(isDummy);
        if(isDummy) numDummies++;

        if(points.size() == size_t(dimension))
        {
            GenerateNormal();
            if(normalize)  NormalizeNormal();
            planeVal = MultiplyPointsByNormal(pointVector->at(points[0]).Data());
        }
    }
    else
    {
        cout << "Simplex already contains maximum number of extremes, ignoring point!" << endl;
    }

}

/*
void Simplex::AddExtreme(const vector<double> & v, bool normalize)
{
   bool isDummy = false;
   int i = 0;
    
   if(int(extremePoints.size()) < dimension)
    {
        extremePoints.push_back(v);
        while(i < dimension && !isDummy)
        {
            if(v[i] >= infinity) isDummy = true;
            i++;
        }
        extremesAreDummies.push_back(isDummy);
        if(isDummy) numDummies++;
        
        if(int(extremePoints.size()) == dimension)
        {
            GenerateNormal();
            if(normalize)  NormalizeNormal();
            planeVal = MultiplyPointsByNormal(extremePoints[0]);
        }
    }
    else
    {
        cout << "Simplex already contains maximum number of extremes, ignoring point!" << endl;
    }

    return;
}
*/

void Simplex::GenerateNormal()
{
    vector< vector<double> > vectorsInHyperplane;
    vector< vector<double> > submat;
    vector<double> temp;
/*    bool allSameSign = true;*/
/*    bool isZero = true;*/
    double val = 0.;
    double t1, t2;
   

    for(int i = 0; i < dimension-1; i++)
    {
        temp.resize(0);
        for(int j = 0; j < dimension; j++)
        {
            t1 = pointVector->at(points[i]).Data()[j];
            t2 = pointVector->at(points[dimension-1]).Data()[j];
            temp.push_back(t1-t2);
        }
        vectorsInHyperplane.push_back(temp);
    }

    normal.resize(0);
    
    for(int k = 0; k < dimension; k++)
    {
        submat.resize(0);
        for(int i = 0; i < dimension - 1; i++)
        {
            temp.resize(0);
            for(int j = 0; j < dimension; j++)
            {
                if(k != j)
                {
                    temp.push_back(vectorsInHyperplane[i][j]);
                }
            }
            submat.push_back(temp);
        }
        if(k % 2 == 0)
        {
            val = Determinant(submat);
        }
        else
        {
            val = -1.*Determinant(submat);
        }
/*        cout << "Val: " << val << "\tNormal Size: " << normal.size() << "\tisZero: " << isZero << "\tisPositive: " << isPositive << "\tallSameSign: " << allSameSign << endl;*/
/*        if(normal.size() > 0 && !isZero && ((val < -epsilon && isPositive) || (val > epsilon && !isPositive))) allSameSign = false;*/
/*        if(val < -epsilon) isPositive = false;*/
/*        if(abs(val) > epsilon) isZero = false;*/
        normal.push_back(val); 
    }
    
/*    if(allSameSign && !isPositive)*/
/*    {*/
/*        for(int i = 0; i < dimension; i++) normal[i] = -1.*normal[i];*/
/*        isPositive = true;*/
/*    }*/
/*    else if(!allSameSign) cout << "There is an error! A normal vector did not have all positive weights!" << endl;*/

    if(DEBUG) 
    {
        cout << "\%Calculated normal: ";
        
        for(int i = 0; i < dimension; i++) cout << normal[i] << "\t";
        cout << endl << endl;
    }
    
    return;
}

double Determinant(const vector< vector<double> > & matrix)
{
    double retVal = 0.;
    int n = matrix.size();
    int subi, subj;
    vector< vector<double> > submat;
    vector<double> temp(n-1);
    
    for(int k = 0; k < n-1; k++)
    {
        submat.push_back(temp);
    }
    
    if (n == 2) 
    {
        retVal = (matrix[0][0] * matrix[1][1]) - (matrix[1][0] * matrix[0][1]);
    }
    else
    {  
        for(int c = 0; c < n; c++)
        {  
            subi = 0;  
            for(int i = 1; i < n; i++)
            {  
                subj = 0;
                for(int j = 0; j < n; j++)
                {    
                    if (j != c)
                    {
                        submat[subi][subj] = matrix[i][j];
                        subj++;
                    }
                }
                subi++;
            }
            retVal = retVal + (pow(-1 ,c) * matrix[0][c] * Determinant(submat));
        }
    }
    return retVal;
}

void Simplex::PrintData() const
{
    for(int i = 0; i < dimension; i++)
    {
        cout << "Extreme Point " << i+1 << ": " << endl;
        for(int j = 0; j < dimension; j++)
        {
		cout << pointVector->at(points[i]).Data()[j] << "\t";
        }
        cout << endl;
    }
    
    cout << endl << "Normal Vector: " << endl;
    for(int i = 0; i < dimension; i++)
    {
        cout << normal[i] << "\t";
    }
    cout << endl;
}

void Simplex::NormalizeNormal()
{
    double magnitude = 0.;
    int sign = 0;
    bool allSameSign = true;
    
    for(int i = 0; i < dimension; i++)
    {
        magnitude += normal[i]*normal[i];
    }
    magnitude = sqrt(magnitude);
    
    for(int i = 0; i < dimension; i++)
    {
        normal[i] = normal[i]/magnitude;
        if(isnan(normal[i]))
        {
            cout << "Error. There is a normal component that doesn't make sense. Exiting!\n";
            exit(0);
        }
        if(!sign)
        {
            if(normal[i] > epsilon) sign = 1;
            else if(normal[i] < -epsilon) sign = -1;
        }
        else if((sign > 0 && normal[i] < -epsilon) || (sign < 0 && normal[i] > epsilon)) allSameSign = false;
    }
    
    if(!allSameSign) isPositive = false;
    else if(sign < 0)
    {
        for(int i = 0; i < dimension; i++)
        {
            normal[i] = -normal[i];
        }
    }
    
    return;
}

void Simplex::Reset()
{
    points.resize(0);
    normal.resize(0);
    extremesAreDummies.resize(0);
    numDummies = 0;
    isPositive = true;
    saveForSolution = false;
    
    return;
}

vector<double> Simplex::GetNormal() const
{ 
    return normal;
}

double Simplex::MultiplyPointsByNormal(const vector<double> & v) const
{
    double retVal = 0.;
    
    for(int i = 0; i < dimension; i++)
    {
        retVal += normal[i]*v[i];
    }
    
    return retVal;
}

double Simplex::PlaneVal()
{ 
    return planeVal;
}

vector<double> Simplex::GetExtremePoint(int i) const
{
    return pointVector->at(points[i]).Data();
}

void Simplex::WriteOctaveCodeToPlotSimplex(bool a) const
{
    vector<double> midpoint(dimension, 0.);
    vector<double> norm = GetNormal(); 
    cout << "patch([";
    for(int i = 0; i < dimension; i++)
    {
        for(int j = 0; j < dimension; j++)
        {
		cout << pointVector->at(j).Data()[i] << " ";
        }
        cout << "], [";
    }
    for(int i = 0; i < 3; i++) 
    {
/*        if(a) cout << ((double) rand() / (RAND_MAX)) << " ";*/
        if(a) cout << ".75" << " ";
        else cout << "0 ";
    }
    
    cout << "]);" <<  endl;
    
    for(int i = 0; i < dimension; i++)
    {
        double t;
        for(int j = 0; j < dimension; j++) {
             t = pointVector->at(j).Data()[i];
             midpoint[i] += (1./3.)*t;
        }
    }
    
    if(showNormalsInPlots)
    {
        cout << "plot3(";
        for(int i = 0; i < dimension; i++) 
        {
            if(i!=dimension-1) cout << "[" << midpoint[i] << " " << midpoint[i] - 1000*norm[i] << "],";
            else cout << "[" << midpoint[i] << " " << midpoint[i] - 1000*norm[i] << "]";
        }
    /*    for(int i = 0; i < dimension; i++) */
    /*    {*/
    /*        if(i != dimension-1) cout << midpoint[i] - 10*norm[i] << ",";*/
    /*        else cout << midpoint[i] - 10*norm[i];*/
    /*    }  */
        cout << ")" << endl;
    }
    
    if(DEBUG)
    {
        cout << "\%";
        for(int i = 0; i < dimension; i++) 
        {
            cout << norm[i] << "\t";
        }
        cout << endl;
    }
    
    return;
}

bool Simplex::IsPositive() const
{
    return isPositive;
}


// bennett: I bet this needs rewritten
// I think we only need to look at the adjacent
Simplex *  Simplex::FindAdjacentContainingOriginalPoints(const vector<Simplex * > & simplices, int & simplexIndex, int & newPointIndex, int ignoreIndex, bool isIn) const {

    cerr << "*****************************************************" << endl;
    cerr << " In FindAdjacentContainingOriginalPoints, check to see if this routine needs rewritten " << endl;
    cerr << "*****************************************************" << endl;
    Simplex * retSimplex= nullptr;
    Simplex * temp;
    vector<double> p;
    int k = 0, l = 0, m = 0;
    bool foundIt = false;
    bool couldBeIt = true;
    bool tempBool = true;
    vector<bool> matches(dimension,false);
    
/*    bool DEBUG = true;*/
/*    bool DEBUG = false;*/
    
    newPointIndex = -1;
    
    if(DEBUG) 
    {
        cout << "\%Looking for an adjacent simplex containing the points: " << endl;
        
        for(int i = 0; i < dimension; i++)
        {
            if(i != ignoreIndex)
            {
                cout << "\%\t" ;
                for(int j = 0; j < dimension; j++) cout << GetExtremePoint(i)[j] << "\t";
                cout << endl;
            }
        }
        
        cout << "ignored index: " << ignoreIndex << endl;
    }
    
    while(k < int(simplices.size()) && !foundIt)
    {
        couldBeIt = true;
        if(isIn && k == simplexIndex) k++;
        if(k < int(simplices.size()))
        {
            temp = simplices[k];
/*            if(DEBUG) */
/*            {*/
/*                cout << "^^^^^^^^^^^^" << endl;*/
/*                PrintData();*/
/*                cout << "^^^^^^^^^^^^" << endl;*/
/*                temp->PrintData();*/
/*                cout << "^^^^^^^^^^^^" << endl;*/
/*            }*/
            for(int i = 0; i < dimension; i++) matches[i] = false;
            k++;
            l = 0;
            while(l < dimension && couldBeIt)
            {
                if(l == ignoreIndex) l++;
                if(l < dimension)
                {
                    couldBeIt = false;
                    m = 0;
                    while(m < dimension)
                    {
                        tempBool = (GetExtremePoint(l) == temp->GetExtremePoint(m));
                        couldBeIt = max(couldBeIt, tempBool);
                        if(tempBool && !matches[m]) matches[m] = true;
/*                        if(DEBUG) cout << "tempBool: " << tempBool << "\tl:" << l << "\tm:" << m << "\tmatches[m]: " << matches[m] << endl; */
                        m++;
                    }
                    l++;
                }
            }
            if(couldBeIt) 
            {
                foundIt = true;
                retSimplex = temp;
                simplexIndex = k-1;
            }
        }
    }
    
    k = 0;
    while(newPointIndex < 0 && k < dimension)
    {
        if(!matches[k]) newPointIndex = k;
        k++;
    }
    
    if(DEBUG) 
    {
        if(retSimplex!= nullptr)
        {
            cout << "\%The adjacent simplex contains: " << endl;
            
            for(int i = 0; i < dimension; i++)
            {
                cout << "\%\t" ;
                for(int j = 0; j < dimension; j++) cout << retSimplex->GetExtremePoint(i)[j] << "\t";
                cout << endl;
            }
        }
        else
        {
            cout << "\%No adjacent simplex containing these points" << endl;
        }
    }
    
    return retSimplex;
}

int Simplex::GetDimension() const
{
   return dimension; 
}

bool CompareDoubleVectors(const vector<double> & v1, const vector<double> & v2)
{
    bool retBool = false;
    int i = 0;
    
    if(v1.size() == v2.size())
    {
        while(i < int(v1.size()) && v1[i] == v2[i]) i++;
        if(i == int(v1.size())) retBool = true;
    }
    
    return retBool;
}

vector< vector<double> >  Simplex::GetExtremePoints() const
{
    vector<vector<double> > rv;
    for(int i=0;i<dimension;i++) {
       rv.push_back(pointVector->at(points[i]).Data());
    }
    return rv;

}

int Simplex::GetNumDummyPoints()
{
    return numDummies;
}

vector<bool> Simplex::ExtremesAreDummies() const
{
    return extremesAreDummies;
}

int Simplex::GetPointIndex(const vector<double> & point)
{
    size_t i = 0;
    while(i < points.size() and pointVector->at(points[i]).Data() != point) i++;
    return i;
}

void Simplex::SaveForSolution()
{
    saveForSolution = true;
    return;
}

bool Simplex::GetSaveForSolution() const
{
    return saveForSolution;
}

// bennett 7/18
void Simplex::Epsilon(double e) {
    epsilon = e;
}

double Simplex::Epsilon(void) const{
    return epsilon;
}

void Simplex::PrintNormal()
{
    for(unsigned int i = 0; i < normal.size(); i++) cout << normal[i] << "\t";
/*    cout << endl;*/
    return;
}

// new routines added by bennett for adjacency

// bennett 7/24/18


// This is mostly for taking care of myself.  I know which
// facet I am adding an ajdacent to.
//
// There will be another one which takes a list of point indexes
// and a simplex pointer and adds that as an adjacent.  This will take much
// more work.
//
void Simplex::Adjacent(int i, Simplex * a) {
    while ( adj.size() <= size_t(i)) {
       adj.push_back(nullptr);
    }
    adj[i] = a;
    return;
}

// just give me the ith adjacent simplex.
Simplex * Simplex::Adjacent(int i) const {
    return adj[i];
}

// Remove a simplex as adjacent. 
// Call this routine when the adjacent simplex is to be deleted.
void Simplex::RemoveAdjacent(Simplex * s){
    size_t i;
    for(i=0;i<adj.size();i++) {
       if (adj[i] == s) {
          adj[i] = nullptr;
       }
    }
    return;
}

// add the point indicies to the simplex.  
// This should eventually replace AddExtreme();
void Simplex::Points(vector<int> p, bool normalize){
    points = p;
    GenerateNormal();
    if(normalize)  NormalizeNormal();
    planeVal = MultiplyPointsByNormal(pointVector->at(points[0]).Data());

    // for set difference later.
    sort(points.begin(), points.end());
    return;
}

void Simplex::MarkForSplit(void){
    split = true;
    return;
}

bool Simplex::IsMarkedForSplit(void) const{
    return split;
}

vector<int> DoSetDiff(vector<int> & a, vector<int> & b) {
    vector<int> rv;
    set_difference(a.begin(), a.end(), b.begin(), b.end(), inserter(rv, rv.begin()));
    return rv;
}

int FindPos(int value, vector<int> & points) {
    size_t i = 0;
    while (i < points.size() and points[i] != value) {
       i++;
    }
    return i;

}

void Simplex::Adjacent( Simplex * b){
    vector<int> da, db;
    int pos;

    if (b == this) {
       cout <<"Refuse to mark myself as adjacent to myself" << endl;
       return;

    }

    da = DoSetDiff(points, b->points);
    if (da.size() == 1) {
       db = DoSetDiff(b->points, points);
       pos = FindPos(da[0],points);
       if (not split) {
          adj[pos] = b;
       }
       pos = FindPos(db[0],b->points);
       if (not b->split) {
           b->adj[pos] = this;
       }
    }

    return;
}


// below here needs to be fixed, extremePoints[] needs to go away


bool Simplex::isInSubsetOfStack(const vector<Simplex *> & simplexStack, int startingScanIndex)
{
    vector< vector<double> > v;
    int count = 0;
    bool retBool = false;
    
    for(unsigned int i = startingScanIndex; i < simplexStack.size(); i++)
    {
        v.resize(0);
        for(int k = 0; k < simplexStack[0]->GetDimension(); k++)
        {
            v.push_back(pointVector->at(k).Data()) ;
            //v.push_back(extremePoints[k]);
        }
    
        for(int j = 0; j < simplexStack[0]->GetDimension(); j++)
        {
            v.push_back(simplexStack[i]->GetExtremePoint(j));
        }
        sort(v.begin(), v.end());
        count = unique(v.begin(), v.end()) - v.begin();
        if(count <= simplexStack[0]->GetDimension())
        {
            cout << "The new simplex is already here! Don't add it!\n";
            retBool = true;
            break;
        }
    }
    
    return retBool;
}

bool Simplex::deleteRepeats(vector<Simplex*> & simplexStack, int startingScanIndex)
{
    vector< vector<double> > v;
    int count = 0;
    bool retBool = false;
    unsigned int i = startingScanIndex;
    
    while(i < simplexStack.size())
    {
        v.resize(0);
        for(int k = 0; k < simplexStack[0]->GetDimension(); k++)
        {
            v.push_back(pointVector->at(k).Data()) ;
            //v.push_back(extremePoints[k]);
        }
    
        for(int j = 0; j < simplexStack[0]->GetDimension(); j++)
        {
            v.push_back(simplexStack[i]->GetExtremePoint(j));
        }
        sort(v.begin(), v.end());
        count = unique(v.begin(), v.end()) - v.begin();
        if(count <= simplexStack[0]->GetDimension())
        {
            if(DEBUG) cout << "Deleting a repeated simplex.\n";
            simplexStack.erase(simplexStack.begin() + i);
            retBool = true;
        }
        else i++;
    }
    
    return retBool;
}
