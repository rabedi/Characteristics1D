#ifndef COMMON_FUNCTIONS__H
#define COMMON_FUNCTIONS__H

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include <vector>
#include <map>
#include <cstdlib>
using namespace std;

#include "Dims.h"
#include <float.h>

void gChangeSerialNum(string& filename);

char* ITOA( int value, char* result, int base);
void SplitFileName(const string& filename, string& path, string& name, string& ext, string& folderSpliter);
void MergeFileName(string& filename, const string& path, const string& name, const string& ext, string folderSpliter = "/");
// "" values means ext and root are not changed
void ChangeFileNameExtensionPath(string& newfilename, const string& filename, string newExt = "", string newPath = "", string newName = "", bool newPathAddedAtTheEndOfCurrenctPath = false, 
								string subtractEndOfFileName = "", string addEndOfFileName = "", string subtractBeginningOfFileName = "", string addBeginningFileName = "");

// existenceMode = 0, does not check existence of file
//				 = 1, checks if file exists, if not, exits the code
//				 = 2, checks if file exists, if not, uses original file name
void AddChangeSerialNumber2FileName(const string& filename, string& filenameNew, int serialNum, int existenceMode);
void AddChangeShortVSerialNumber2FileName(string& filename, int serialNum = -1, int existenceMode = 2);
// returns false if checkIfFileCanBeOpened && fileCannot be opened
// serialNo < 0 is not added
bool AddSerialNumber2FileName(const string& fileNameWOSerial, string& fileNameOut, int serialNo = -1, string txtBeforeSerial = "", bool checkIfFileCanBeOpened = true);

void print_state(const std::ios& stream);

// this is a global serial number that is used to change the serial number of certain files =>
// this allows running many problems with for example different random meshes ... useful for Monte Carlo method
// less than zero means it's inactive
extern int g_i_serial;
// 0 does not check existence of file, 1 checks it and exit if file does not exist, 2 does not change file, if serial number changed file does not exist - see AddChangeShortVSerialNumber2FileName
extern int g_existenceMode;

void gChangeSerialNum(string& filename);

template <class T>
void Integrate(const vector<double>& x, vector<T>& vals, vector<T>& integrals, int num_x = -1, bool uniform_xGrid = false)
{
	if (num_x < 0)
		num_x = MIN(x.size(), vals.size());
	if (num_x == 0)
		return;
	integrals.resize(num_x);
	integrals[0] = 0.0 * vals[0];
	if (num_x == 1)
		return;

	if (!uniform_xGrid || (num_x == 2))
	{
		vector<double> half_delx(num_x - 1);
		for (int i = 0; i < num_x - 1; ++i)
			half_delx[i] = 0.5 * (x[i + 1] - x[i]);

		for (int i = 1; i < num_x; ++i)
		{
			//T delIntegral = half_delx[i - 1] * (vals[i - 1] + vals[i]);
			T delIntegral = half_delx[i - 1] * vals[i - 1] ;
			delIntegral +=  half_delx[i - 1] * vals[i];
			integrals[i] = integrals[i - 1] + delIntegral;
		}
		return;
	}
	double h = x[1] - x[0];
	//integrals[1] = (0.5 * h) * (vals[0] + vals[1]);
	integrals[1] = (0.5 * h) * vals[0] ;
	integrals[1] += (0.5 * h) *vals[1];
	double sa = h / 3.0, sb = 4.0 * sa;
	//integrals[2] = sa * (vals[0] + vals[2]) + sb * vals[1];
	integrals[2] = sa * vals[0];
	integrals[2] += sa * vals[2];
	integrals[2] +=  sb * vals[1];

	if (num_x == 3)
		return;
	//integrals[3] = (3.0 / 8.0 * h) * ((vals[0] + vals[3]) + 3.0 * (vals[1] + vals[2]));
	double tmp1 = (3.0 / 8.0 * h), tmp2 = 3.0 * tmp1;
	integrals[3] =  tmp1 * vals[0];
	integrals[3] += tmp1 * vals[3];
	integrals[3] += tmp2 * vals[1] ;
	integrals[3] += tmp2 * vals[2];

	if (num_x == 4)
		return;
	for (int i = 4; i < num_x; ++i)
	{
		//T delIntegral = sa * (vals[i - 2] + vals[i]) + sb * vals[i - 1];
		T delIntegral = sa * vals[i - 2];
		delIntegral += sa *  vals[i];
		delIntegral +=  sb * vals[i - 1];

		//integrals[i] = delIntegral + integrals[i - 2];
		integrals[i] = delIntegral;
		integrals[i] += integrals[i - 2];
	}
}

/*
// Obsolete version of the function above
template <class T>
inline void ReadUnitReachSpecificTypeValue(istream& in, T& val)
{
	in >> val;
	while (in.fail() == true)
	{
//		print_state(in);
		in.clear();
		in.ignore(std::numeric_limits<streamsize>::max(),'\n');
		in >> val;
	}
}
*/

// Ian McNamara moved these two to commonClasses.h
/*template <class T>
ostream& operator<<(ostream& output, vector<T>& vec)
{
	int sz = vec.size();
	output<<"\nsz\t"<<sz<<endl;
	for (int i = 0; i < sz; ++i)
		output<<vec[i]<<endl;
	return output;
}


template <class T>
ostream& operator<<(ostream& output, const vector<T>& vec)
{
	int sz = vec.size();
	output<<"\nsz\t"<<sz<<endl;
	for (int i = 0; i < sz; ++i)
		output<<vec[i]<<endl;
	return output;
}*/

template <class T> 
inline int Add(vector<T>& vec, const T &value)
{
	int index = Find(vec, value);
	if (index < 0)
	{
		vec.push_back(value);
		index = vec.size() - 1;
	}
	return index;
}

template <class T> 
inline int Add(vector<T>& vec, const T &value, bool &found)
{
	int index = Find(vec, value);
	if (index < 0)
	{
		found = false;
		vec.push_back(value);
		index = vec.size() - 1;
	}
	else
		found = true;
	return index;
}


template <class T> 
inline int Union(vector<T>& vecA, vector<T>& vecB, vector<T>& vecUnion)
{
	vecUnion.clear();
	int szA = vecA.size();
	int szB = vecB.size();

	if (szA >= szB)
	{
		vecUnion = vecA;
		for (int i = 0; i < szB; ++i)
			Add(vecUnion, vecB[i]);
		return vecUnion.size();
	}
	vecUnion = vecB;
	for (int i = 0; i < szA; ++i)
		Add(vecUnion, vecA[i]);
	return vecUnion.size();
}

template <class T> 
inline int Union(vector< vector<T> >& vecs, vector<T>& vecUnion)
{
	int sz = vecs.size();
	if (sz == 0)
		return 0;
	if (sz == 1)
	{
		vecUnion = vecs[0];
		return 1;
	}
	Union(vecs[0], vecs[1], vecUnion);
	for (int i = 2; i < sz; ++i)
	{
		vector<T> temp;
		temp = vecUnion;
//		vecUnion.clear();
		return Union(vecs[i], temp, vecUnion);
	}
}

template <class T> 
inline bool HaveIntersection(vector<T>& vecA, vector<T>& vecB)
{
	int szA = vecA.size();
	int szB = vecB.size();

	if (szA <= szB)
	{
		for (int i = 0; i < szA; ++i)
			if (Find(vecB, vecA[i]) >= 0)
				return true;
		return false;
	}
	for (int i = 0; i < szB; ++i)
		if (Find(vecA, vecB[i]) >= 0)
			return true;
	return false;
}


template <class T> 
inline int Intersection(vector<T>& vecA, vector<T>& vecB, vector<T>& vecIntersection)
{
	vecIntersection.clear();
	int szA = vecA.size();
	int szB = vecB.size();

	if (szA <= szB)
	{
		for (int i = 0; i < szA; ++i)
			if (Find(vecB, vecA[i]) >= 0)
				vecIntersection.push_back(vecA[i]);
		return vecIntersection.size();
	}
	for (int i = 0; i < szB; ++i)
		if (Find(vecA, vecB[i]) >= 0)
			vecIntersection.push_back(vecB[i]);
	return vecIntersection.size();
}

template <class T> 
inline int Intersection(vector< vector<T> >& vecs, vector<T>& vecIntersection)
{
	int sz = vecs.size();
	if (sz == 0)
		return 0;
	if (sz == 1)
	{
		vecIntersection = vecs[0];
		return 1;
	}
	Intersection(vecs[0], vecs[1], vecIntersection);
	for (int i = 2; i < sz; ++i)
	{
		vector<T> temp;
		temp = vecIntersection;
//		vecIntersection.clear();
		Intersection(vecs[i], temp, vecIntersection);
	}
	return vecIntersection.size();
}

// check whether setA is subset of setB
template <class T> 
inline bool isSubset(vector<T>& setA, vector<T>& setB)
{
	for (int i = 0 ; i < setA.size(); ++i)
		if (Find(setB, setA[i]) < 0)
			return false;
	return true;
}


#if 0
template<class KEY, class TARGET>
void removeFromMap(TARGET targetVal, std::map<KEY, TARGET>& mp)
{
	std::map<KEY, TARGET>::iterator iter;
	for (iter = mp.begin();
		iter != mp.end();  /* no increment! */)
	{
		if (iter->second == targetVal)
			mp.erase(iter++);
		else
			++iter;
	}
}
#endif

int getNewComponentPosition(vector<double>& vec, double newComp, double tol, bool& foundInVector);
// if pos == -1 sent to the function, pos is first calculated from the function above
void addNewComponentPosition(vector<double>& vec, double newComp, double tol, int& pos, bool& foundInVector);

template<class T>
void addNewComponentPosition(vector<T>& vec, T& newComp, int pos)
{
	int sz = vec.size();
	vec.push_back(newComp);
	if (pos < sz)
	{
		for (int i = sz; i > pos; --i)
			vec[i] = vec[i - 1];
		vec[pos] = newComp;
	}
}



/*PC 3/30/2017*/
/*AUX FUNC*/
//=============================================
bool AreSame(const double a, const double b, double tolerance = -1.0);
bool LessThan(const double a, const double b, double tolerance = -1.0);
bool LessThanOrEqual(const double a, const double b, double tolerance = -1.0);
bool GreaterThan(const double a, const double b, double tolerance = -1.0);
bool GreaterThanOrEqual(const double a, const double b, double tolerance = -1.0);
bool AreSame(const vector<double>& a, const vector<double>& b, double tolerance = -1.0);
//=============================================

template<class T>
inline void findIntervalIndices(vector<int>& pos, vector<T>& vals, T value, int left = -1, int right = -1)
{
	if (left == -1)
		left = 0;
	else if (left < 0)
		THROW("left boundary index must be non negative real integer.");

	if (right == -1)
		right = (int)vals.size() - 1;
	else if (right >= (int)vals.size())
		THROW("right boundary index out of bounds.");

	//fill(pos.begin(), pos.end(), -1);
	pos.clear();

	int sz = (int)vals.size();

	while (left <= right) {
		int middle = left + (right - left) / 2; // (left + right) / 2;
		if (AreSame(vals[middle], value))
		{
			pos.push_back(middle);
			return;
		}
		else
		{
			if (vals[middle] > value)
			{
				if (LessThanOrEqual(vals[middle - 1], value) && LessThanOrEqual(value, vals[middle]))
				{
					pos.push_back(middle - 1);
					pos.push_back(middle);
					return;
				}
				else
					right = middle - 1;
			}
			else
			{
				if (LessThanOrEqual(vals[middle], value) && LessThanOrEqual(value, vals[middle + 1]))
				{
					pos.push_back(middle);
					pos.push_back(middle + 1);
					return;
				}
				else
					left = middle + 1;
			}
		}
	}
	//return -1;
};

// General Math Functions
double nFactorial(int n);
void multivariatePolynomialTerms(vector<double>& crd, int level, vector<int>& I, vector<int>& Imax, vector<double>& terms, int& arbitrary_counter, bool computeDerivative = false, int wrtJ = 0, int nthDiff = 1);
double PolyValue(double x, vector<double>& coef, int p = -1);
double genPolynomialFunc(vector<double>& crd, vector<int> polyOrder, vector<double>& coefficients, bool computeDerivative = false, int wrtJ = 0, int nthDiff = 1);

double genPolynomialDerivative(vector<double>& crd, vector<int> polyOrder, vector<double>& coefficients, int j, int n);

// each term is in form a*exp(b*x+c) ...
double genUnivariateExponentialFunc(vector<double>& crd, int numTerms, vector<double>& coefficients);

double genUnivariateExponentialDerivative(vector<double>& crd, int numTerms, vector<double>& coefficients, int wrtJ, int nthDiff);

// each term is in form a*sin(b*x+c)
double genUnivariateHarmonicFunc(vector<double>& crd, int numTerms, vector<double>& coefficients);

double genUnivariateHarmonicDerivative(vector<double>& crd, int numTerms, vector<double>& coefficients, int wrtJ, int nthDiff);

// each term is in form a*ln(b*x+c) ...
double genUnivariateLogFunc(vector<double>& crd, int numTerms, vector<double>& coefficients);

double genUnivariateLogDerivative(vector<double>& crd, int numTerms, vector<double>& coefficients, int wrtJ, int nthDiff);

// -1, 0, 1 value for x < 0; x = 0, x > 0 
int Sign(double x);

// if Max(abs(valA), abs(valB)) < relTol, absolute error check is used, otherwise relative
bool IsRelativeErrorAccepable(double valA, double valB, double relTol, double& relError);


void ReadSetInteger(istream& in, std::set<int>& dat);
void ReadVectorInteger(istream& in, vector<int>& dat);
void ReadVectorString(istream& in, vector<string>& dat);
void ReadMapInteger2Integer(istream& in, std::map<int, int>& dat);
void ReadMapString2String(istream& in, std::map<string, string>& dat);
void ReadSetDouble(istream& in, std::set<double>& dat);
void ReadVectorDouble(istream& in, vector<double>& dat);
void ReadMapInteger2Double(istream& in, std::map<int, double>& dat);

void WriteSetInteger(ostream& out, const std::set<int>& dat);
void WriteVectorInteger(ostream& out, const vector<int>& dat);
void WriteVectorString(ostream& out, const vector<string>& dat);
void WriteMapInteger2Integer(ostream& out, const std::map<int, int>& dat);
void WriteSetDouble(ostream& out, const std::set<double>& dat);
void WriteVectorDouble(ostream& out, const vector<double>& dat);
void WriteMapInteger2Double(ostream& out, const std::map<int, double>& dat);


// http://en.wikipedia.org/wiki/Non-analytic_smooth_function
double getInfinitelySmoothTransitionFunctionX_0_2_Infinity_Y_0_2_1(double x, double& dydx, bool compute_dydx = false);
double getInfinitelySmoothRampFunctionX_0_2_1_Y_0_2_1(double x, double& dydx, bool compute_dydx = false);
double getInfinitelySmoothBumpFunctionX_0_2_1_Y_0_2_1_2_0(double x, double& dydx, bool compute_dydx = false);
// the function is regularized in [-diracPositivieWidth, diracPositivieWidth]
double getInfinitlySmoothDeltaDirac(double x, double& dydx, double diracPositivieWidth = 0.01, bool compute_dydx = false);

// approximate Dirac function: from Jing 1995: Jing_1995_Dynamic contact force model for contactable cracks with static and kinetic friction
// contact/papers/numericalExamples 
// g(t) = 16 * [g4(tau) - 4 * g4(tau - 1/4) + 6 g4(tau - 1/2) -4 g4(tau - 3/4) + g4(tau - 1)]
// tau = t/T
// g4(x) = x^3 H(x)
// T = duration of Delta Dirac approximation
// this goes from [0, T] to center it around zero send t - T/2 to the function or modify the function input parameter
//	normalizeOpt		0			-> maximum function value is 1
//	normalizeOpt		1			-> integral under the curve is 1
double Dirac(double t, double T, int normalizeOpt);


// x <= xLimMin -> y = functionValMin
// x >= xLimMax -> y = functionvalMax
// infinitely smooth transition in between
double getInfinitlySmoothRampFunction(double x, double& dydx, double xLimMax = 1.0, double functionValMax = 1.0, double xLimMin = 0.0, double functionValMin = 0.0, bool compute_dydx = false);

#endif

