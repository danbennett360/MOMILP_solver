/* File created by Dr. Nathan Adelgren, Assistant Professor at Edinboro University of PA.
Collaborators are currently Dr. Dan Bennett and Sydney Lesseski.
Started: 5/12/2018 
Finished: N/A
This work is a start toward the solution of multiobjective mixed-integer linear programs. 
Initially we will just build a data structure for storing the (minimally excessive) set 
of nondominated solutions.*/

#include "cplex.h"
#include "simplex_class.h"

Simplex::Simplex(int i)
{
    dimension = i;
}

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

void Simplex::GenerateNormal()
{
    vector< vector<double> > vectorsInHyperplane;
    vector< vector<double> > submat;
    vector<double> temp;
    bool allSameSign = true;
    bool isZero = true;
    double val = 0.;
    
    for(int i = 0; i < dimension-1; i++)
    {
        temp.resize(0);
        for(int j = 0; j < dimension; j++)
        {
            temp.push_back(extremePoints[i][j] - extremePoints[dimension-1][j]);
        }
        vectorsInHyperplane.push_back(temp);
    }
    
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
        if(normal.size() > 0 && !isZero && ((val < 0 && isPositive) || (val > 0 && !isPositive))) allSameSign = false;
        if(val < -epsilon) isPositive = false;
        if(abs(val) > epsilon) isZero = false;
        normal.push_back(val); 
    }
    
    if(allSameSign && !isPositive)
    {
        for(int i = 0; i < dimension; i++) normal[i] = -1.*normal[i];
        isPositive = true;
    }
/*    else if(!allSameSign) cout << "There is an error! A normal vector did not have all positive weights!" << endl;*/

    if(DEBUG) 
    {
        cout << "\%Calculated normal: ";
        
        for(int i = 0; i < dimension; i++) cout << normal[i] << "\t";
        cout << endl;
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
            cout << extremePoints[i][j] << "\t";
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
    }
    
    return;
}

void Simplex::Reset()
{
    extremePoints.resize(0);
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
    return extremePoints[i];
}

void Simplex::WriteOctaveCodeToPlotSimplex() const
{
    vector<double> midpoint(dimension, 0.);
    vector<double> norm = GetNormal(); 
    cout << "patch([";
    for(int i = 0; i < dimension; i++)
    {
        for(int j = 0; j < dimension; j++)
        {
            cout << extremePoints[j][i] << " ";
        }
        cout << "], [";
    }
    for(int i = 0; i < 3; i++) cout << ((double) rand() / (RAND_MAX)) << " ";
    
    cout << "]);" <<  endl;
    
    for(int i = 0; i < dimension; i++)
    {
        for(int j = 0; j < dimension; j++) midpoint[i] += (1./3.)*extremePoints[j][i];
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

bool Simplex::IsPositive()
{
    return isPositive;
}

Simplex Simplex::FindAdjacentContainingOriginalPoints(const vector<Simplex> & simplices, int & simplexIndex, int & newPointIndex, int ignoreIndex, bool isIn) const
{
    Simplex retSimplex(-1);
    Simplex temp;
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
/*                temp.PrintData();*/
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
                        tempBool = (GetExtremePoint(l) == temp.GetExtremePoint(m));
                        couldBeIt = max(couldBeIt, tempBool);
                        if(tempBool && !matches[m]) matches[m] = true;
                        if(DEBUG) cout << "tempBool: " << tempBool << "\tl:" << l << "\tm:" << m << "\tmatches[m]: " << matches[m] << endl; 
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
        if(retSimplex.GetDimension() > 0)
        {
            cout << "\%The adjacent simplex contains: " << endl;
            
            for(int i = 0; i < dimension; i++)
            {
                cout << "\%\t" ;
                for(int j = 0; j < dimension; j++) cout << retSimplex.GetExtremePoint(i)[j] << "\t";
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
    return extremePoints;
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
    int i = 0;
    while(extremePoints[i] != point) i++;
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

bool Simplex::isInSubsetOfStack(const vector<Simplex> & simplexStack, int startingScanIndex)
{
    vector< vector<double> > v;
    int count = 0;
    bool retBool = false;
    
    for(unsigned int i = startingScanIndex; i < simplexStack.size(); i++)
    {
        v.resize(0);
        for(int k = 0; k < simplexStack[0].GetDimension(); k++)
        {
            v.push_back(extremePoints[k]);
        }
    
        for(int j = 0; j < simplexStack[0].GetDimension(); j++)
        {
            v.push_back(simplexStack[i].GetExtremePoint(j));
        }
        sort(v.begin(), v.end());
        count = unique(v.begin(), v.end()) - v.begin();
        if(count <= simplexStack[0].GetDimension())
        {
            cout << "The new simplex is already here! Don't add it!\n";
            retBool = true;
            break;
        }
    }
    
    return retBool;
}

bool Simplex::deleteRepeats(vector<Simplex> & simplexStack, int startingScanIndex)
{
    vector< vector<double> > v;
    int count = 0;
    bool retBool = false;
    unsigned int i = startingScanIndex;
    
    while(i < simplexStack.size())
    {
        v.resize(0);
        for(int k = 0; k < simplexStack[0].GetDimension(); k++)
        {
            v.push_back(extremePoints[k]);
        }
    
        for(int j = 0; j < simplexStack[0].GetDimension(); j++)
        {
            v.push_back(simplexStack[i].GetExtremePoint(j));
        }
        sort(v.begin(), v.end());
        count = unique(v.begin(), v.end()) - v.begin();
        if(count <= simplexStack[0].GetDimension())
        {
            if(DEBUG) cout << "Deleting a repeated simplex.\n";
            simplexStack.erase(simplexStack.begin() + i);
            retBool = true;
        }
        else i++;
    }
    
    return retBool;
}
