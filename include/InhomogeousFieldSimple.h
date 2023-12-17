#ifndef INHOMOGENEOUS_FIELD_SIMPLE__H
#define INHOMOGENEOUS_FIELD_SIMPLE__H

#include <iostream> /* for input and output using >> and << operator */
#include <fstream>  /* file streams for input and output */
#include <cstring>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <assert.h>     /* assert */
#include <math.h>     /* math functions */
#include <array>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <complex>
#include <cmath>

#ifndef THROW
#if VCPP
#define THROW(msg){{char tmpstr[255];sprintf_s(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);	cerr<<tmpstr;getchar(); getchar(); throw(tmpstr); }}
#else
#define THROW(msg){{char tmpstr[255];sprintf(tmpstr, "In %s, line %d : %s \n", __FILE__, __LINE__, msg);	cerr<<tmpstr;getchar(); getchar(); throw(tmpstr); }}
#endif
#endif

using namespace std;

template<class T>
void toStringTmp(T& dat, string& name)
{
	ostringstream convert;
	convert << dat;
	name = convert.str();
}

class InhomogeneousSimple
{
public:
	InhomogeneousSimple();
	void Compute_InhomogeneousSimple(int serialNoIn = 0, double llcIn = -3.5, double deltaIn = 0.9, int shapeIn = 2, double resolutionFactorIn = 1, string ssofFSIn = "min", int baseResolution_pIn = 14);
	// inputs
	int serialNo;
	double llc;
	int shape;
	double delta;
	// 2, 4, ....
	double resolutionFactor;
	string ssoFS; // min, mean_arithmetic, mean_harmonic, valStart
	int baseResolution_p; // 2^14 base resolution meshes

	// computed
	vector<double> snValues;
	vector<double> mappedValues; // between (1 - delta, 1 + delta)
	vector<double> mappedValues_coarsened; // coarsened values from the list above
	unsigned int coarseningRat;
	unsigned int sz_mappedValues_coarsened;
	double min_mappedValues_coarsened, mean_mappedValues_coarsened;
	int index_min_mappedValues_coarsened;

	static string root;

private:
	bool isUniform;
	double inv_shape;
	double GenSymTriangle_getInverseCDF(double sn);
};

double getStatvalueSimple(const vector<double>& vals, string ssoFS, int& index);


void Test_InhomogeneousSimple();

#endif
