#ifndef MICRO_CROSSING_SOLVER_PERIODIC_D_FRAGMENT__H
#define MICRO_CROSSING_SOLVER_PERIODIC_D_FRAGMENT__H

#include "MinCrossingSolver1D.h"
#include "SLInterfaceCalculator.h"

// entry of log(a) for calculation
// la_noNormalization: actual values of log(a)		is read
// la_2Normalization:  tauC = 0.5 is used log(a'), a' = a * 0.5,		is read
// la_GNormalization: normalization based on G, a' = a * G, log(a') is read; G = 0.5 for Ortiz, exp(1) for Xu-Needleman
typedef enum {la_noNormalization, la_2Normalization, la_GNormalization, la_NormalizationType_SIZE} la_NormalizationType;
string getName(la_NormalizationType dat);
bool name2Type(string& name, la_NormalizationType& typeVal);
ostream& operator<<(ostream& out, la_NormalizationType dat);
istream& operator>>(istream& in, la_NormalizationType& dat);

class gFx2yPDF : public gFx2y
{
public:
	gFx2yPDF();
	virtual void InitializeValues(Solver1D* conf, const map<string, string>& str_mapIn, const vector<unsigned int>& indices4ParasIn, const vector<double>& parasIn, double& xMin, double& xMax, double& tol_x, vector<double>& primary_xs, vector<double>& secondary_xs, int& num_y, vector<double>& tol_ys);
	virtual bool ComputeValue(genIndexVal& giv);
	virtual void Print_YHeader(ostream& out, int num_y = -1) const;
	SLFF_TSRType	tsrModel;
	// entry of log(a) for calculation
	// la_noNormalization: actual values of log(a)		is read
	// la_2Normalization:  tauC = 0.5 is used log(a'), a' = a * 0.5,		is read
	// la_GNormalization: normalization based on G, a' = a * G, log(a') is read; G = 0.5 for Ortiz, exp(1) for Xu-Needleman
	la_NormalizationType laNormT;

	// writes log 10 files for a and l
	void Writel10_al_TSR_file(double l10l = 1.0, int index_primary_l = 0, int index_secondary_l = 0);
};

void MAIN_SolvePeriodicFragmentSize(string confSolve = "config/OneInterface/Periodic_segment/Periodic_segmen_confSolver1D.txt");

#endif

